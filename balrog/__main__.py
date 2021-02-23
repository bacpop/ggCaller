import os
import tarfile
import time
import pickle
import numpy as np
from tqdm.auto import tqdm
from Bio.Seq import Seq
from scipy.special import expit
from scipy.special import logit

import torch
#import torch.nn as nn
import torch.nn.functional as F
#from torch.nn.utils import weight_norm

""" Get directories for model and seengenes """
module_dir = os.path.dirname(os.path.realpath(__file__))
model_dir = os.path.join(module_dir, "balrog_models")

""" Print what the program is doing."""
verbose = True

""" Use kmer prefilter to increase gene sensitivity. 
May not play nice with very high GC genomes."""
protein_kmer_filter = True

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
unidirectional_penalty_per_base = 3.895921717182765 # 3' 5' overlap
convergent_penalty_per_base = 4.603432608883688 # 3' 3' overlap
divergent_penalty_per_base = 3.3830814940689975 # 5' 5' overlap
k_seengene = 10
multimer_threshold = 2

#find genes
def tokenize_aa_seq(aa_seq):
    """ Convert amino acid letters to integers."""
    table = {"L": 1,
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
    tokenized = torch.tensor([table[aa] for aa in aa_seq])
    return tokenized


def get_start_codon(seq, orfcoords, strand):
    if strand == 1:
        # forward strand
        startcoord = orfcoords[0]
        return seq[startcoord - 3:startcoord]
    else:
        # reverse strand
        startcoord = orfcoords[1]
        return seq[startcoord:startcoord + 3].reverse_complement()

def get_ORF_info(ORF_colour_ID_map):
    ORF_colour = []
    ORF_IDs = []
    ORF_seq = []
    ORF_coord = []
    ORF_nucseq = []
    # iterate over list of ORFs
    for colour, seq_dict in ORF_colour_ID_map.items():
        for ORF_ID, ORF in seq_dict.items():

            seq = Seq(ORF)
            length = len(ORF)

            ### find the ORF 16bp upstream of start of sequence
            orf_coords = [16, length - 3]

            # standardize coords
            full_ORF_nuccord = (orf_coords[0] + 3, orf_coords[1])

            # translate once per frame, then slice
            aa = str(seq[full_ORF_nuccord[0]:full_ORF_nuccord[1]].translate(table=translation_table, to_stop=False))

            ORF_IDs.append(ORF_ID)
            ORF_colour.append(colour)
            ORF_seq.append(aa)
            ORF_coord.append(full_ORF_nuccord)
            ORF_nucseq.append(ORF)

    return ORF_IDs, ORF_colour, ORF_seq, ORF_nucseq, ORF_coord

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


def predict_tis(model_tis, X):
    model_tis.eval()
    with torch.no_grad():
        if torch.cuda.device_count() > 0:
            X_enc = F.one_hot(X, 4).permute(0, 2, 1).float().cuda()
        else:
            X_enc = F.one_hot(X, 4).permute(0, 2, 1).float()
        probs = expit(model_tis(X_enc).cpu())
    return probs


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


def tensor_to_seq(tensor):
    table = {0: "X",
             1: "L",
             2: "V",
             3: "I",
             4: "M",
             5: "C",
             6: "A",
             7: "G",
             8: "S",
             9: "T",
             10: "P",
             11: "F",
             12: "Y",
             13: "W",
             14: "E",
             15: "D",
             16: "N",
             17: "Q",
             18: "K",
             19: "R",
             20: "H"}
    return "".join([table[x] for x in tensor])

def kmerize(seq, k):
    kmerset = set()
    for i in range(len(seq) - k + 1):
        kmer = tuple(seq[i: i + k].tolist())
        kmerset.add(kmer)
    return kmerset

def score_genes(ORF_dict, ORF_colour_ID_map, minimum_ORF_score, num_threads):
    # set up load gene and TIS models
    """ extract and load balrog model """

    # check if directory exists. If not, unzip file
    if not os.path.exists(model_dir):
        tar = tarfile.open(model_dir + ".tar.gz", mode="r:gz")
        tar.extractall(module_dir)
        tar.close()

    torch.hub.set_dir(model_dir)
    print("Loading convolutional model...")
    if torch.cuda.device_count() > 0:
        print("GPU detected...")
        model = torch.hub.load(model_dir, "geneTCN", source='local').cuda()
        model_tis = torch.hub.load(model_dir, "tisTCN", source='local').cuda()
        time.sleep(0.5)
    else:
        print("No GPU detected, using CPU...")
        torch.set_num_threads(num_threads)
        model = torch.hub.load(model_dir, "geneTCN", source='local')
        model_tis = torch.hub.load(model_dir, "tisTCN", source='local')
        time.sleep(0.5)

    """Load k-mer filters"""
    genexa_kmer_path = os.path.join(model_dir, "10mer_thresh2_minusARF_all.pkl")

    # load kmer filter
    with open(genexa_kmer_path, "rb") as f:
        aa_kmer_set = pickle.load(f)

    # get sequences and coordinates of ORFs
    if verbose:
        print("Finding and translating open reading frames...")

    ORF_ID_list, ORF_colour_list, ORF_seq_list, ORF_nucseq_list, ORF_coord_list = get_ORF_info(ORF_colour_ID_map)

    # encode amino acids as integers
    if verbose:
        print("Encoding amino acids...")
    ORF_seq_enc = [tokenize_aa_seq(x) for x in ORF_seq_list]

    # seengene check
    if protein_kmer_filter:
        if verbose:
            print("Applying protein kmer filter...")
        seengene = []
        for s in ORF_seq_enc:
            kmerset = kmerize(s, k_seengene)
            s = [x in aa_kmer_set for x in kmerset]
            seen = np.sum(s) >= multimer_threshold

            seengene.append(seen)

    # score
    if verbose:
        print("Scoring ORFs with temporal convolutional network...")

    # sort by length to minimize impact of batch padding
    ORF_lengths = np.asarray([len(x) for x in ORF_seq_enc])
    length_idx = np.argsort(ORF_lengths)
    ORF_seq_sorted = [ORF_seq_enc[i] for i in length_idx]

    # pad to allow creation of batch matrix
    prob_list = []
    for i in tqdm(range(0, len(ORF_seq_sorted), gene_batch_size), unit=" batch"):
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
    ORF_gene_score = [None] * len(ORF_coord_list)
    for k, coord in enumerate(ORF_gene_score):
        ORF_gene_score[k] = float(ORF_prob[idx])
        idx += 1

    if verbose:
        print("Scoring translation initiation sites...")

    # extract nucleotide sequence surrounding potential start codons
    ORF_TIS_seq = ORF_coord_list[:]
    ORF_start_codon = [None] * len(ORF_coord_list)

    for i, contig in enumerate(ORF_TIS_seq):
        n = 0  # count to index into flat structure # TODO make sure this works as expected

        nucseq = ORF_nucseq_list[i]  # easier to use coords relative to single nucseq

        if any(contig):
            coords = ORF_coord_list[i]
            n += 1
            fiveprime = coords[0]
            if fiveprime >= 16 + 3:  # NOTE 16 HARD CODED HERE
                downstream = nucseq[fiveprime: fiveprime + 16]
                upstream = nucseq[fiveprime - 16 - 3: fiveprime - 3]
                start_codon = start_enc[nucseq[fiveprime - 3: fiveprime]]
                TIS_seq = torch.tensor([nuc_encode[c] for c in (upstream + downstream)[::-1]],
                                       dtype=int)  # model scores 3' to 5' direction
            else:
                TIS_seq = -1  # deal with gene fragments later
                start_codon = 2

            ORF_TIS_seq[i] = TIS_seq
            ORF_start_codon[i] = start_codon

    # flatten TIS for batching
    ORF_TIS_prob = [None] * len(ORF_TIS_seq)

    ORF_TIS_seq_flat = []
    ORF_TIS_seq_idx = []
    for i, contig in enumerate(ORF_TIS_seq):
        if type(contig) == int:  # fragment
            ORF_TIS_prob[i] = 0.5  # HOW BEST TO DEAL WITH FRAGMENT TIS?
        elif len(contig) != 32:
            ORF_TIS_prob[i] = 0.5
        else:
            ORF_TIS_seq_flat.append(contig)
            ORF_TIS_seq_idx.append(i)

    # batch score TIS
    TIS_prob_list = []
    for i in tqdm(range(0, len(ORF_TIS_seq_flat), TIS_batch_size), unit=" batch"):
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
            length = ORF_coord_list[i][1] - ORF_coord_list[i][0] + 1
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

            length = ORF_coord_list[i][1] - ORF_coord_list[i][0] + 1
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

    # update initial dictionary with strand and score within a tuple
    for i, score in enumerate(ORF_score_flat):
        # if score greater than minimum, add to the ORF_dict
        if score >= minimum_ORF_score:
            ORF_dict[ORF_colour_list[i]][ORF_ID_list[i]] = score
        # else, remove the ORF from the full ORF list
        else:
            del ORF_dict[ORF_colour_list[i]][ORF_ID_list[i]]

    return ORF_dict