import os
import subprocess
import sys
import re

from joblib import Parallel, delayed
from functools import partial
from tqdm import tqdm
from Bio.Align.Applications import MafftCommandline
from multiprocessing import Pool


def run_snpsites_dir(annotation_dir, vcf_dir, no_vc_set, n_cpus):
    with Pool(processes=n_cpus, maxtasksperchild=1) as pool:
        pool.map(partial(run_snpsites, annotation_dir=annotation_dir, vcf_dir=vcf_dir),
                 [file for file in os.listdir(annotation_dir) if file not in no_vc_set])

    return True


def run_snpsites(file, annotation_dir, vcf_dir):
    outfile = os.path.join(vcf_dir, file.split(".")[0] + ".vcf")
    file = os.path.join(annotation_dir, file)

    command = ["snp-sites", "-v", "-o", outfile, file]

    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if result.returncode != 0:
        find_err = re.search(r'No SNPs were detected', str(result))
        if find_err is None:
            raise Exception("Snp-sites failed to run on file: " + file)

    return


def check_snpsites_install():
    command = ["snp-sites", "-V"]

    p = str(
        subprocess.run(command,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE))
    present = False
    find_ver = re.search(r'snp-sites 2\.\d+\.\d+', p)
    if find_ver != None:
        present = True

    if present == False:
        sys.stderr.write("Need snp-sites v2 or greater to be installed " + "\n")
        sys.exit(1)

    return present


def check_aligner_install():
    """Checks for the presence of the specified aligned in $PATH

    Args:
        check_aligner_install(str)
            str = specified aligner

    Returns:
        presence (bool)
            True/False aligner present
    """
    command = ["mafft", "--help"]

    p = str(
        subprocess.run(command,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE))
    present = False

    find_ver = re.search(r'MAFFT v7\.\d+', p)
    if find_ver != None:
        present = True

    if present == False:
        sys.stderr.write("Need Mafft version 7 or later to be installed " + "\n")
        sys.exit(1)

    return present


def get_alignment_commands(fastafile_name, outdir, aligner):
    if aligner == "def":
        command = MafftCommandline(input=fastafile_name,
                                   auto=True,
                                   amino=True)
    elif aligner == "ref":
        ref_file, seq_file = fastafile_name
        outfile = outdir + "aligned_gene_sequences/" + seq_file.split('/')[-1].split('.')[0] + '.aln.fas'
        command = ["mafft", "--amino", "--6merpair", "--addfragments",
                   seq_file, ref_file, outfile]

    return (command, fastafile_name)


def align_sequences(command, outdir, aligner):
    if aligner == "def":
        name = str(command[0]).split()[-1].split('/')[-1].split('.')[0]
        stdout, stderr = command[0]()
        with open(outdir + name + '.aln.fas', 'w+') as handle:
            handle.write(stdout)
    elif aligner == "ref":
        with open(command[0][-1], 'w+') as handle:
            result = subprocess.run(command[0][:-1], stdout=handle, stderr=subprocess.DEVNULL)
            if result.returncode != 0:
                raise Exception("Mafft-ref failed to run on file: " + command[0][-1].split("/")[-1])

    if isinstance(command[1], tuple):
        for file in command[1]:
            try:
                os.remove(file)
            except FileNotFoundError:
                None
    else:
        try:
            os.remove(command[1])
        except FileNotFoundError:
            None
    return True


def multi_align_sequences(commands, outdir, threads, aligner, quiet):
    alignment_results = Parallel(n_jobs=threads, prefer="threads")(
        delayed(align_sequences)(x, outdir, aligner) for x in tqdm(commands, disable=quiet))

    return True


def write_alignment_header(alignment_list, outdir):
    out_entries = []
    # Set the tracking variables for gene positions
    gene_start = 1
    gene_end = 0
    for gene in alignment_list:
        # Get length and name from one sequence in the alignment
        # Set variables that need to be set pre-output
        gene_end += gene[2]
        gene_name = gene[0]
        # Create the 3 line feature entry
        gene_entry1 = "FT   feature         " + str(gene_start) + ".." + str(
            gene_end) + '\n'
        gene_entry2 = "FT                   /label=" + gene_name + '\n'
        gene_entry3 = "FT                   /locus_tag=" + gene_name + '\n'
        gene_entry = gene_entry1 + gene_entry2 + gene_entry3
        # Add it to the output list
        out_entries.append(gene_entry)
        # Alter the post-output variables
        gene_start += gene[2]
    # Create the header and footer
    header = ("ID   Genome standard; DNA; PRO; 1234 BP.\nXX\nFH   Key" +
              "             Location/Qualifiers\nFH\n")
    footer = ("XX\nSQ   Sequence 1234 BP; 789 A; 1717 C; 1693 G; 691 T;" +
              " 0 other;\n//\n")
    # open file and output
    with open(outdir + "core_alignment_header.embl", "w+") as outhandle:
        outhandle.write(header)
        for entry in out_entries:
            outhandle.write(entry)
        outhandle.write(footer)
    return True
