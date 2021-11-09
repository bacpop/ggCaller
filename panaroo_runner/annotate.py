import os
import re
import sys
import subprocess


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


def run_diamond_search(annotate, annotation_temp_dir, annotation_db, n_cpu, verbose):
    if verbose:
        print("annotating gene families...")

    # set working directory and reference for snakefile
    os.environ["ANNOWORKDIR"] = annotation_temp_dir
    os.environ["ANNODB"] = annotation_db
    os.environ['DIAMONDSETTING'] = annotate[0]

    # set off Snakemake pipeline
    # remove
    snakemake = "/home/sth19/miniconda3/envs/ggCaller/bin/snakemake"
    if annotate == "fast":
        command = [snakemake, "--cores", str(n_cpu), "-R", "diamond_search_90perc"]
    elif annotate == "medium":
        command = [snakemake, "--cores", str(n_cpu), "-R", "diamond_search_70perc"]
    elif annotate == "sensitive":
        command = [snakemake, "--cores", str(n_cpu), "-R", "diamond_search_0perc"]

    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise Exception("Snakemake failed with diamond search!")
        sys.exit(1)

    return
