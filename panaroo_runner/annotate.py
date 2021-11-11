import os
import re
import sys
import subprocess
from .generate_output import get_unannotated_nodes, output_aa_sequence
from Bio import SeqIO
import pandas as pd

def check_diamond_install():
    # remove
    command = ["/home/sth19/miniconda3/envs/ggCaller/bin/diamond", "help"]
    # command = ["diamond", "help"]

    p = str(
        subprocess.run(command,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.DEVNULL))

    present = False

    find_ver = re.search(r'diamond v\d+\.\d+\.\d+', p)
    if find_ver != None:
        present = True

    if present == False:
        sys.stderr.write("Need diamond to be installed and available in PATH " + "\n")
        sys.exit(1)

    return present


def check_HMMER_install():
    # remove
    command = ["/home/sth19/miniconda3/envs/ggCaller/bin/hmmscan", "-h"]
    # command = ["diamond", "help"]

    p = str(
        subprocess.run(command,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.DEVNULL))

    present = False

    find_ver = re.search(r'HMMER \d+\.\d+\.\d+', p)
    if find_ver != None:
        present = True

    if present == False:
        sys.stderr.write("Need HMMER3 to be installed and available in PATH  " + "\n")
        sys.exit(1)

    return present


def generate_HMMER_index(infile):
    # remove
    command = ["/home/sth19/miniconda3/envs/ggCaller/bin/hmmpress", infile]
    # command = ["hmmpress", infile]
    result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if result.returncode != 0:
        raise Exception("hmmpress failed to run on file: " + infile)
        sys.exit(1)

def generate_diamond_index(infile):
    outfile = infile.split["."][0] + ".dmnd"
    # remove
    command = ["/home/sth19/miniconda3/envs/ggCaller/bin/diamond", "makedb", "--in", "infile", "-d",
               infile.split["."][0]]
    # command = ["diamond", "makedb", "--in", "infile", "-d", infile.split["."][0]]
    result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if result.returncode != 0:
        raise Exception("Diamond make-db failed to run on file: " + infile)
        sys.exit(1)

    return outfile


def run_diamond_search(G, annotation_temp_dir, annotation_db, nb_colours, evalue, pool):
    # list of sequence records
    all_centroid_aa = []

    # list of dictionary of annotations for each sequence if available
    annotation_list = [{}] * nb_colours

    # get unannotated nodes
    unannotated_nodes = get_unannotated_nodes(G)

    # Multithread writing amino acid sequences to disk (temp directory)
    for centroid_aa in pool.map(output_aa_sequence, unannotated_nodes):
        all_centroid_aa = all_centroid_aa + centroid_aa

    # write all sequences to single file
    all_centroid_aa = (x for x in all_centroid_aa)
    SeqIO.write(all_centroid_aa, annotation_temp_dir + "aa_d.fasta", 'fasta')

    command = ["/home/sth19/miniconda3/envs/ggCaller/bin/diamond", "blastp", "--iterate", "--evalue", str(evalue), "-d",
               annotation_db, "-q",
               annotation_temp_dir + "aa_d.fasta", "-o", annotation_temp_dir + "aa_d.tsv"]

    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise Exception("Snakemake failed with diamond search!")
        sys.exit(1)

    # read in file, map highest scoring annotation and bitscore to query
    df = pd.read_csv(annotation_temp_dir + "aa_d.tsv", sep='\t', header=None)
    df = pd.concat([df.iloc[:, 0:2], df.iloc[:, 10:]], axis=1)
    df.set_axis(['query', 'target', 'evalue', 'bitscore'], axis=1, inplace=True)

    # split node in which query is found
    df['node'] = df['query'].str.split(';').str[0]

    # only take highest value for each node
    df = df.sort_values('bitscore').drop_duplicates("node", keep='last')

    # pull information from dataframe
    node_info = [(w, x, y, z) for w, x, y, z in
                 zip(df['query'], df['node'], df['target'], df['bitscore'])]

    # add information to graph. If multiple centroids aligned, take the annotation with the highest bitscore
    for entry in node_info:
        G.nodes[entry[1]]['annotation'] = entry[2]
        G.nodes[entry[1]]['bitscore'] = entry[3]

        # add entries for node to annotation_list
        gene_id = entry[0].split(";")[-1]
        genome = gene_id.split("_")[0]

        # add list for each seqid entry for annotation
        for entry in G.nodes[entry[1]]['seqIDs']:
            annotation_list[genome][entry] = ["diamond", entry[2], entry[3]]

    return G, annotation_list


def run_HMMERscan(G, annotation_list, annotation_temp_dir, annotation_db, evalue, n_cpu, pool):
    # list of sequence records
    all_centroid_aa = []

    # get unannotated nodes
    unannotated_nodes = get_unannotated_nodes(G)

    # Multithread writing amino acid sequences to disk (temp directory)
    for centroid_aa in pool.map(output_aa_sequence, unannotated_nodes):
        all_centroid_aa = all_centroid_aa + centroid_aa

    # write all sequences to single file
    all_centroid_aa = (x for x in all_centroid_aa)
    SeqIO.write(all_centroid_aa, annotation_temp_dir + "aa_h.fasta", 'fasta')

    command = ["/home/sth19/miniconda3/envs/ggCaller/bin/hmmscan", "-E", str(evalue), "--tblout",
               annotation_temp_dir + "aa_h.tsv",
               "--cpu", str(n_cpu), annotation_db, annotation_temp_dir + "aa_h.fasta"]

    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise Exception("Snakemake failed with HMMERscan!")
        sys.exit(1)

    # read in file, map highest scoring annotation and bitscore to query
    df = pd.read_csv(annotation_temp_dir + "aa_h.tsv", delim_whitespace=True, header=None,
                     comment='#')
    df = pd.concat([df.iloc[:, 0], df.iloc[:, 2], df.iloc[:, 4:6]], axis=1)
    df.set_axis(['target', 'query', 'evalue', 'bitscore'], axis=1, inplace=True)

    # split node in which query is found
    df['node'] = df['query'].str.split(';').str[0]

    # only take highest value for each node
    df = df.sort_values('bitscore').drop_duplicates("node", keep='last')

    # pull information from dataframe
    node_info = [(int(x.split(";")[0]), x.split(";")[1], y, z) for x, y, z in
                 zip(df['query'], df['target'], df['bitscore'])]

    # add information to graph. If multiple centroids aligned, take the annotation with the highest bitscore
    for entry in node_info:
        G.nodes[entry[0]]['annotation'] = entry[2]
        G.nodes[entry[0]]['bitscore'] = entry[3]

        # add entries for node to annotation_list
        gene_id = entry[0].split(";")[-1]
        genome = gene_id.split("_")[0]

        # add list for each seqid entry for annotation
        for entry in G.nodes[entry[1]]['seqIDs']:
            annotation_list[genome][entry] = ["diamond", entry[2], entry[3]]

    return G


def iterative_annotation_search(G, annotation_temp_dir, annotation_db, hmm_db, evalue, nb_colours, n_cpu, pool):
    # run initial iterative search
    G, annotation_list = run_diamond_search(G, annotation_temp_dir, annotation_db, evalue, nb_colours, pool)

    # run ultra-sensitive search
    G, annotation_list = run_HMMERscan(G, annotation_list, annotation_temp_dir, hmm_db, evalue, n_cpu, pool)

    return G, annotation_list
