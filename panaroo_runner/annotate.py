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


def run_diamond_search(annotate, annotation_temp_dir, annotation_db):
    # set working directory and reference for snakefile
    os.environ["ANNOWORKDIR"] = annotation_temp_dir
    os.environ["ANNODB"] = annotation_db
    os.environ['DIAMONDSETTING'] = annotate[0]

    # set off Snakemake pipeline
    # remove
    snakemake = "/home/sth19/miniconda3/envs/ggCaller/bin/snakemake"
    if annotate == "fast":
        command = [snakemake, "-c1", "-R", "diamond_search_90perc"]
    elif annotate == "medium":
        command = [snakemake, "-c1", "-R", "diamond_search_70perc"]
    elif annotate == "sensitive":
        command = [snakemake, "-c1", "-R", "diamond_search_0perc"]

    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise Exception("Snakemake failed with diamond search!")
        sys.exit(1)


def iterative_diamond_search(G, annotation_temp_dir, annotation_db, threshold, pool):
    # first iteration of annotation
    all_centroid_aa = []

    # get unannotated nodes
    unannotated_nodes = get_unannotated_nodes(G, threshold)

    # Multithread writing amino acid sequences to disk (temp directory)
    for centroid_aa in pool.map(output_aa_sequence, unannotated_nodes):
        all_centroid_aa = all_centroid_aa + centroid_aa

    # write all sequences to single file
    all_centroid_aa = (x for x in all_centroid_aa)
    SeqIO.write(all_centroid_aa, annotation_temp_dir + "all_aa_f.fasta", 'fasta')

    run_diamond_search("fast", annotation_temp_dir, annotation_db)

    # read in file, map highest scoring annotation and escore
    df = pd.read_csv(annotation_temp_dir + "all_aa_f_f.tsv", sep='\t')
    df = pd.concat([df.iloc[:, 0:2], df.iloc[:, 10:]], axis=1)
    df.rename(index={0: "query", 1: "target", 2: "evalue"})

    df.sort_values('evalue').drop_duplicates("query", keep='last')

    return G
