import re
import sys
import subprocess
from functools import partial
from .generate_output import get_unannotated_nodes, output_aa_sequence, output_dna_sequence
from Bio import SeqIO
import pandas as pd
from Bio import SearchIO
from collections import defaultdict

def check_diamond_install():
    command = ["diamond", "help"]

    p = str(
        subprocess.run(command,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.DEVNULL))

    present = False

    find_ver = re.search(r'diamond v\d+\.\d+', p)
    if find_ver != None:
        present = True

    if present == False:
        sys.stderr.write("Need diamond to be installed and available in PATH " + "\n")
        sys.exit(1)

    return present


def check_HMMER_install():
    command = ["hmmscan", "-h"]

    p = str(
        subprocess.run(command,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.DEVNULL))

    present = False

    find_ver = re.search(r'HMMER \d+\.\d+', p)
    if find_ver != None:
        present = True

    if present == False:
        sys.stderr.write("Need HMMER3 to be installed and available in PATH  " + "\n")
        sys.exit(1)

    return present


def generate_HMMER_index(infile):
    command = ["hmmpress", infile]
    result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if result.returncode != 0:
        raise Exception("hmmpress failed to run on file: " + infile)
        sys.exit(1)

def generate_diamond_index(infile):
    outfile = infile.split(".")[0] + ".dmnd"
    command = ["diamond", "makedb", "--in", infile, "-d",
               infile.split(".")[0]]
    result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if result.returncode != 0:
        raise Exception("Diamond make-db failed to run on file: " + infile)
        sys.exit(1)

    return outfile


def run_diamond_search(G, shd_arr_tup, overlap, annotation_temp_dir, annotate, annotation_db, evalue, pool):
    # list of sequence records
    all_centroid_dna = []

    # get unannotated nodes
    unannotated_nodes = get_unannotated_nodes(G)

    # Multithread writing amino acid sequences to disk (temp directory)
    for centroid_dna in pool.map(partial(output_dna_sequence, shd_arr_tup=shd_arr_tup, overlap=overlap),
                                unannotated_nodes):
        all_centroid_dna = all_centroid_dna + centroid_dna

    # write all sequences to single file
    all_centroid_aa = (x for x in all_centroid_dna)
    SeqIO.write(all_centroid_aa, annotation_temp_dir + "dna_d.fasta", 'fasta')

    if annotate == "fast":
        command = ["diamond", "blastx", "--iterate", "--evalue", str(evalue), "-d",
                   annotation_db, "--outfmt", "6", "qseqid", "sseqid", "bitscore", "stitle", "-q",
                   annotation_temp_dir + "dna_d.fasta", "-o", annotation_temp_dir + "dna_d.tsv"]
    else:
        command = ["diamond", "blastx", "--sensitive", "--evalue", str(evalue), "-d",
                   annotation_db, "--outfmt", "6", "qseqid", "sseqid", "bitscore", "stitle", "-q",
                   annotation_temp_dir + "dna_d.fasta", "-o", annotation_temp_dir + "dna_d.tsv"]

    #result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise Exception("Diamond search failed!")
        sys.exit(1)

    # read in file, map highest scoring annotation and bitscore to query
    try:
        df = pd.read_csv(annotation_temp_dir + "dna_d.tsv", sep='\t', header=None)
    except pd.errors.EmptyDataError:
        return G
    df = df.set_axis(['query_id', 'id', 'bitscore', 'description'], axis=1, copy=False)

    # split node in which query is found
    df['node'] = df['query_id'].str.split(';').str[0]

    # only take highest value for each node
    df = df.sort_values('bitscore').drop_duplicates("node", keep='last')

    # pull information from dataframe
    node_info = [(int(w), x, y, z) for w, x, y, z in
                 zip(df['node'], df['id'], df['bitscore'], df['description'])]

    # add information to graph. If multiple centroids aligned, take the annotation with the highest bitscore
    for entry in node_info:
        G.nodes[entry[0]]['annotation'] = entry[1]
        G.nodes[entry[0]]['bitscore'] = entry[2]
        G.nodes[entry[0]]['description'] = entry[3]

    return G


def run_HMMERscan(G, shd_arr_tup, overlap, annotation_temp_dir, annotation_db, evalue, pool, n_cpu):
    # list of sequence records
    all_centroid_aa = []

    # get unannotated nodes
    unannotated_nodes = get_unannotated_nodes(G)

    # Multithread writing amino acid sequences to disk (temp directory)
    for centroid_aa in pool.map(partial(output_aa_sequence, shd_arr_tup=shd_arr_tup, overlap=overlap),
                                unannotated_nodes):
        all_centroid_aa = all_centroid_aa + centroid_aa

    # write all sequences to single file
    all_centroid_aa = (x for x in all_centroid_aa)
    SeqIO.write(all_centroid_aa, annotation_temp_dir + "aa_h.fasta", 'fasta')

    command = ["hmmscan", "-E", str(evalue), "--tblout",
               annotation_temp_dir + "aa_h.tsv",
               "--cpu", str(n_cpu), annotation_db, annotation_temp_dir + "aa_h.fasta"]

    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise Exception("HMMscan failed!")
        sys.exit(1)

    # read in file, map highest scoring annotation and bitscore to query
    attribs = ['id', 'query_id', 'bitscore', 'description']
    hits = defaultdict(list)
    with open(annotation_temp_dir + "aa_h.tsv") as handle:
        for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
            for hit in queryresult.hits:
                for attrib in attribs:
                    hits[attrib].append(getattr(hit, attrib))

    # generate pandas dataframe for determining most signficant hit per node
    df = pd.DataFrame.from_dict(hits)

    # check if empty
    if df.empty:
        return G

    # split node in which query is found
    df['node'] = df['query_id'].str.split(';').str[0]

    # only take highest value for each node
    df = df.sort_values('bitscore').drop_duplicates("node", keep='last')

    # pull information from dataframe
    node_info = [(int(w), x, y, z) for w, x, y, z in
                 zip(df['node'], df['id'], df['bitscore'], df['description'])]

    # add information to graph. If multiple centroids aligned, take the annotation with the highest bitscore
    for entry in node_info:
        G.nodes[entry[0]]['annotation'] = entry[1]
        G.nodes[entry[0]]['bitscore'] = entry[2]
        G.nodes[entry[0]]['description'] = entry[3]

    return G


def iterative_annotation_search(G, shd_arr_tup, overlap, annotation_temp_dir, annotation_db, hmm_db,
                                evalue, annotate, n_cpu, pool):
    # run initial iterative search
    G = run_diamond_search(G, shd_arr_tup, overlap, annotation_temp_dir,
                                              annotate, annotation_db, evalue, pool)

    # run ultra-sensitive search
    if annotate == "ultrasensitive":
        G = run_HMMERscan(G, shd_arr_tup, overlap, annotation_temp_dir,
                                             hmm_db, evalue, pool, n_cpu)

    return G
