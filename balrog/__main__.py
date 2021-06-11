import os
import tarfile
import time
import pickle
import numpy as np
from Bio.Seq import Seq
from scipy.special import expit
from scipy.special import logit

import torch
import torch.nn.functional as F

""" Get directories for model and seengenes """
module_dir = os.path.dirname(os.path.realpath(__file__))
model_dir = os.path.join(module_dir, "balrog_models")

""" Print what the program is doing."""
verbose = True

""" Use kmer prefilter to increase gene sensitivity. 
May not play nice with very high GC genomes."""
protein_kmer_filter = False

""" Nucleotide to amino acid translation table. 11 for most bacteria/archaea.
4 for Mycoplasma/Spiroplasma."""
translation_table = 11
# translation_table = 4

""" Batch size for the temporal convolutional network used to score genes.
Small batches and big batches slow down the model. Very big batches may crash the 
GPU. """
gene_batch_size = 200
TIS_batch_size = 1000

""" All following are internal parameters. Change at your own risk."""
weight_gene_prob = 0.9746869839852076
weight_TIS_prob = 0.25380288790532707
score_threshold = 0.47256101519707244
weight_ATG = 0.84249804151264
weight_GTG = 0.7083689705744909
weight_TTG = 0.7512400826652517
unidirectional_penalty_per_base = 3.895921717182765  # 3' 5' overlap
convergent_penalty_per_base = 4.603432608883688  # 3' 3' overlap
divergent_penalty_per_base = 3.3830814940689975  # 5' 5' overlap
k_seengene = 10
multimer_threshold = 2

nuc_encode = {"A": 0,
              "T": 1,
              "G": 2,
              "C": 3,
              "N": 0,
              "M": 0,
              "R": 0,
              "Y": 0,
              "W": 0,
              "K": 0}

start_enc = {"ATG": 0,
             "GTG": 1,
             "TTG": 2}

aa_table = {"L": 1,
            "V": 2,
            "I": 3,
            "M": 4,
            "C": 5,
            "A": 6,
            "G": 7,
            "S": 8,
            "T": 9,
            "P": 10,
            "F": 11,
             "Y": 12,
            "W": 13,
            "E": 14,
            "D": 15,
            "N": 16,
            "Q": 17,
            "K": 18,
            "R": 19,
            "H": 20,
            "*": 0,
            "X": 0}


# generate ORF sequences from coordinates
@profile
def generate_sequence(graph_vector, nodelist, node_coords, overlap):
    sequence = ""
    for i in range(0, len(nodelist)):
        id = nodelist[i]
        coords = node_coords[i]

        # calculate strand based on value of node (if negative, strand is false)
        strand = True if id >= 0 else False

        if strand:
            unitig_seq = graph_vector[abs(id) - 1].seq
        else:
            unitig_seq = str(Seq(graph_vector[abs(id) - 1].seq).reverse_complement())

        if len(sequence) == 0:
            substring = unitig_seq[coords[0]:(coords[1] + 1)]
        else:
            if coords[1] >= overlap:
                substring = unitig_seq[overlap:(coords[1] + 1)]

        sequence += substring

    return sequence


@profile
def tokenize_aa_seq(aa_seq):
    """ Convert amino acid letters to integers."""
    tokenized = torch.tensor([aa_table[aa] for aa in aa_seq])
    return tokenized


@profile
def get_ORF_info(ORF_vector, graph, overlap):
    ORF_seq_list = []
    TIS_seqs = []
    # iterate over list of ORFs
    for ORFNodeVector in ORF_vector:
        # need to determine ORF sequences from paths
        ORF_nodelist = ORFNodeVector[0]
        ORF_node_coords = ORFNodeVector[1]
        TIS_nodelist = ORFNodeVector[3]
        TIS_node_coords = ORFNodeVector[4]

        # generate ORF_seq, as well as upstream and downstream TIS seq
        ORF_seq = graph.generate_sequence(ORF_nodelist, ORF_node_coords, overlap)
        upstream_TIS_seq = graph.generate_sequence(TIS_nodelist, TIS_node_coords, overlap)
        downstream_TIS_seq = ORF_seq[0:19]

        # generate Seq class for translation
        seq = Seq(ORF_seq)

        # translate once per frame, then slice. Note, do not include start or stop codons
        aa = str(seq[3:-3].translate(table=translation_table, to_stop=False))

        ORF_seq_list.append(aa)
        TIS_seqs.append((upstream_TIS_seq, downstream_TIS_seq))

    # convert amino acids into integers
    ORF_seq_enc = [tokenize_aa_seq(x) for x in ORF_seq_list]

    return ORF_seq_enc, TIS_seqs


@profile
def predict(model, X):
    model.eval()
    with torch.no_grad():
        if torch.cuda.device_count() > 0:
            X_enc = F.one_hot(X, 21).permute(0, 2, 1).float().cuda()
            probs = expit(model(X_enc).cpu())
            del X_enc
            torch.cuda.empty_cache()
        else:
            X_enc = F.one_hot(X, 21).permute(0, 2, 1).float()
            probs = expit(model(X_enc).cpu())

    return probs


@profile
def predict_tis(model_tis, X):
    model_tis.eval()
    with torch.no_grad():
        if torch.cuda.device_count() > 0:
            X_enc = F.one_hot(X, 4).permute(0, 2, 1).float().cuda()
        else:
            X_enc = F.one_hot(X, 4).permute(0, 2, 1).float()
        probs = expit(model_tis(X_enc).cpu())
    return probs


@profile
def kmerize(seq, k):
    kmerset = set()
    for i in range(len(seq) - k + 1):
        kmer = tuple(seq[i: i + k].tolist())
        kmerset.add(kmer)
    return kmerset


def load_kmer_model():
    # check if directory exists. If not, unzip file
    if not os.path.exists(model_dir):
        tar = tarfile.open(model_dir + ".tar.gz", mode="r:gz")
        tar.extractall(module_dir)
        tar.close()

    """Load k-mer filters"""
    genexa_kmer_path = os.path.join(model_dir, "10mer_thresh2_minusARF_all.pkl")

    with open(genexa_kmer_path, "rb") as f:
        aa_kmer_set = pickle.load(f)

    return aa_kmer_set


def load_gene_models():
    # check if directory exists. If not, unzip file
    if not os.path.exists(model_dir):
        tar = tarfile.open(model_dir + ".tar.gz", mode="r:gz")
        tar.extractall(module_dir)
        tar.close()

    torch.hub.set_dir(model_dir)
    # print("Loading convolutional model...")
    if torch.cuda.device_count() > 0:
        # print("GPU detected...")
        model = torch.hub.load(model_dir, "geneTCN", source='local').cuda()
        model_tis = torch.hub.load(model_dir, "tisTCN", source='local').cuda()
        time.sleep(0.5)
    else:
        # print("No GPU detected, using CPU...")
        model = torch.hub.load(model_dir, "geneTCN", source='local')
        model_tis = torch.hub.load(model_dir, "tisTCN", source='local')
        time.sleep(0.5)

    return (model, model_tis)


@profile
def score_genes(ORF_vector, graph_vector, minimum_ORF_score, overlap, model, model_tis, aa_kmer_set):
    # get sequences and coordinates of ORFs
    # print("Finding and translating open reading frames...")

    ORF_seq_enc, TIS_seqs = get_ORF_info(ORF_vector, graph_vector, overlap)

    # seengene check
    if protein_kmer_filter:
        seengene = []
        for s in ORF_seq_enc:
            kmerset = kmerize(s, k_seengene)
            # s = [x in aa_kmer_set for x in kmerset]
            s = np.isin(list(kmerset), aa_kmer_set)
            seen = np.count_nonzero(s) >= multimer_threshold
            seengene.append(seen)

    # score
    # print("Scoring ORFs with temporal convolutional network...")

    # sort by length to minimize impact of batch padding
    ORF_lengths = np.asarray([len(x) for x in ORF_seq_enc])
    length_idx = np.argsort(ORF_lengths)
    ORF_seq_sorted = [ORF_seq_enc[i] for i in length_idx]

    # pad to allow creation of batch matrix
    prob_list = []
    for i in range(0, len(ORF_seq_sorted), gene_batch_size):
        batch = ORF_seq_sorted[i:i + gene_batch_size]
        seq_lengths = torch.LongTensor(list(map(len, batch)))
        seq_tensor = torch.zeros((len(batch), seq_lengths.max())).long()

        for idx, (seq, seqlen) in enumerate(zip(batch, seq_lengths)):
            seq_tensor[idx, :seqlen] = torch.LongTensor(seq)

        pred_all = predict(model, seq_tensor)

        pred = []
        for j, length in enumerate(seq_lengths):
            subseq = pred_all[j, 0, 0:int(length)]
            predprob = float(expit(torch.mean(logit(subseq))))
            pred.append(predprob)

        prob_list.extend(pred)
    prob_arr = np.asarray(prob_list, dtype=float)

    # unsort
    unsort_idx = np.argsort(length_idx)
    ORF_prob = prob_arr[unsort_idx]

    # recombine ORFs
    idx = 0
    ORF_gene_score = [None] * len(ORF_seq_enc)
    for k, coord in enumerate(ORF_gene_score):
        ORF_gene_score[k] = float(ORF_prob[idx])
        idx += 1

    # print("Scoring translation initiation sites...")

    # extract nucleotide sequence surrounding potential start codons
    ORF_TIS_seq_flat = []
    ORF_TIS_seq_idx = []
    ORF_TIS_prob = [None] * len(TIS_seqs)
    ORF_start_codon = [None] * len(ORF_seq_enc)

    for i, TIS in enumerate(TIS_seqs):
        # unpack tuple. Note, downsteam includes start codon, which needs to be removed
        upstream, downstream = TIS
        if len(upstream) == 16:
            TIS_seq = torch.tensor([nuc_encode[c] for c in (upstream + downstream[3:])[::-1]],
                                   dtype=int)  # model scores 3' to 5' direction
            ORF_TIS_seq_flat.append(TIS_seq)
            ORF_TIS_seq_idx.append(i)
        else:
            ORF_TIS_prob[i] = 0.5

        # encode start codon
        start_codon = start_enc[downstream[0:3]]
        ORF_start_codon[i] = start_codon

    # batch score TIS
    TIS_prob_list = []
    for i in range(0, len(ORF_TIS_seq_flat), TIS_batch_size):
        batch = ORF_TIS_seq_flat[i:i + TIS_batch_size]
        TIS_stacked = torch.stack(batch)
        pred = predict_tis(model_tis, TIS_stacked)

        TIS_prob_list.extend(pred)
    y_pred_TIS = np.asarray(TIS_prob_list, dtype=float)

    # reindex batched scores
    for i, prob in enumerate(y_pred_TIS):
        idx = ORF_TIS_seq_idx[i]
        ORF_TIS_prob[idx] = float(prob)

    # combine all info into single score for each ORF
    if protein_kmer_filter:
        ORF_score_flat = []
        for i, geneprob in enumerate(ORF_gene_score):
            if not geneprob:
                ORF_score_flat.append(None)
                continue
            seengene_idx = 0
            # calculate length by multiplying number of amino acids by 3, then adding 6 for start and stop
            length = (len(ORF_seq_enc[i]) * 3) + 6
            TIS_prob = ORF_TIS_prob[i]
            start_codon = ORF_start_codon[i]
            ATG = start_codon == 0
            GTG = start_codon == 1
            TTG = start_codon == 2

            combprob = geneprob * weight_gene_prob \
                       + TIS_prob * weight_TIS_prob \
                       + ATG * weight_ATG \
                       + GTG * weight_GTG \
                       + TTG * weight_TTG
            maxprob = weight_gene_prob + weight_TIS_prob + max(weight_ATG, weight_TTG, weight_GTG)
            probthresh = score_threshold * maxprob
            score = (combprob - probthresh) * length + 1e6 * seengene[seengene_idx]
            seengene_idx += 1

            ORF_score_flat.append(score)

    else:
        ORF_score_flat = []
        for i, geneprob in enumerate(ORF_gene_score):
            if not geneprob:
                ORF_score_flat.append(None)
                continue

            # calculate length by multiplying number of amino acids by 3, then adding 6 for start and stop
            length = len(ORF_seq_enc[i]) * 3
            TIS_prob = ORF_TIS_prob[i]
            start_codon = ORF_start_codon[i]
            ATG = start_codon == 0
            GTG = start_codon == 1
            TTG = start_codon == 2

            combprob = geneprob * weight_gene_prob \
                       + TIS_prob * weight_TIS_prob \
                       + ATG * weight_ATG \
                       + GTG * weight_GTG \
                       + TTG * weight_TTG
            maxprob = weight_gene_prob + weight_TIS_prob + max(weight_ATG, weight_TTG, weight_GTG)
            probthresh = score_threshold * maxprob
            score = (combprob - probthresh) * length

            ORF_score_flat.append(score)

    # update initial dictionary, removing low scoring ORFs and create score mapping score within a tuple
    ORF_score_dict = {}
    for i, score in enumerate(ORF_score_flat):
        # if score greater than minimum, add to the ORF_score_dict
        if score >= minimum_ORF_score:
            ORF_score_dict[i] = score

    return ORF_score_dict
