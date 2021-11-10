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
        sys.stderr.write("Need diamond to be installed " + "\n")
        sys.exit(1)

    return present


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


def run_diamond_search(G, annotate, annotation_temp_dir, annotation_db, threshold, pool):
    # first iteration of annotation
    all_centroid_aa = []

    # get unannotated nodes
    unannotated_nodes = get_unannotated_nodes(G)

    # Multithread writing amino acid sequences to disk (temp directory)
    for centroid_aa in pool.map(output_aa_sequence, unannotated_nodes):
        all_centroid_aa = all_centroid_aa + centroid_aa

    # write all sequences to single file
    all_centroid_aa = (x for x in all_centroid_aa)
    SeqIO.write(all_centroid_aa, annotation_temp_dir + "aa_" + annotate[0] + ".fasta", 'fasta')

    # set working directory and reference for snakefile
    os.environ["ANNOWORKDIR"] = annotation_temp_dir
    os.environ["ANNODB"] = annotation_db
    os.environ['DIAMONDSETTING'] = annotate[0]
    os.environ['EVALUETHRESHOLD'] = str(threshold)

    # set off Snakemake pipeline
    # remove
    snakemake = "/home/sth19/miniconda3/envs/ggCaller/bin/snakemake"
    if annotate == "fast":
        command = [snakemake, "-c1", "-R", "diamond_search_fast"]
    elif annotate == "medium":
        command = [snakemake, "-c1", "-R", "diamond_search_med"]
    elif annotate == "sensitive":
        command = [snakemake, "-c1", "-R", "diamond_search_sens"]
    elif annotate == "iterative":
        command = [snakemake, "-c1", "-R", "diamond_search_iter"]

    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise Exception("Snakemake failed with diamond search!")
        sys.exit(1)

    # read in file, map highest scoring annotation and bitscore to query
    df = pd.read_csv(annotation_temp_dir + "aa_" + annotate[0] + ".tsv", sep='\t', header=None)
    df = pd.concat([df.iloc[:, 0:2], df.iloc[:, 10:]], axis=1)
    df.set_axis(['query', 'target', 'evalue', 'bitscore'], axis=1, inplace=True)

    df = df.sort_values('bitscore').drop_duplicates("query", keep='last')

    # pull information from dataframe
    node_info = [(int(x.split(";")[0]), x.split(";")[1], y, z) for x, y, z in
                 zip(df['query'], df['target'], df['bitscore'])]

    # add information to graph. If multiple centroids aligned, take the annotation with the highest bitscore
    for entry in node_info:
        if G.nodes[entry[0]]['bitscore'] < entry[3]:
            G.nodes[entry[0]]['annotation'] = entry[2]
            G.nodes[entry[0]]['bitscore'] = entry[3]

    return G


def iterative_diamond_search(G, annotation_temp_dir, annotation_db, threshold, pool):
    # run initial iterative search
    G = run_diamond_search(G, "iterative", annotation_temp_dir, annotation_db, threshold, pool)

    # run ultra-sensitive search
    G = run_diamond_search(G, "sensitive", annotation_temp_dir, annotation_db, threshold, pool)

    return G
