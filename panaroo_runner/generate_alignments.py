import os
import subprocess
import sys
import re

from joblib import Parallel, delayed
from tqdm import tqdm

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Bio.Align.Applications import MafftCommandline


def check_aligner_install():
    """Checks for the presence of the specified aligned in $PATH

    Args:
        check_aligner_install(str)
            str = specified aligner

    Returns:
        presence (bool)
            True/False aligner present
    """
    command = "mafft --help"

    p = str(
        subprocess.run(command,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       shell=True))
    present = False

    find_ver = re.search(r'MAFFT v7\.\d+', p)
    if find_ver != None:
        present = True

    if present == False:
        sys.stderr.write("Need Mafft version 7 or later to be installed " + "\n")
        sys.exit(1)

    return present


def output_sequence(node, isolate_list, temp_directory, outdir):
    # Get the name of the sequences for the gene of interest
    sequence_ids = node["seqIDs"]
    output_sequences = []
    # Counter for the number of sequences to
    isolate_no = 0
    # Look for gene sequences among all genes (from disk)
    for seq in SeqIO.parse(outdir + "combined_DNA_CDS.fasta", 'fasta'):
        isolate_num = int(seq.id.split('_')[0])
        isolate_name = isolate_list[isolate_num].replace(";",
                                                         "") + ";" + seq.id
        if seq.id in sequence_ids:
            output_sequences.append(
                SeqRecord(seq.seq, id=isolate_name, description=""))
            isolate_no += 1
    # Put gene of interest sequences in a generator, with corrected isolate names
    output_sequences = (x for x in output_sequences)
    # set filename to gene name, if more than one sequence to be aliged
    if isolate_no > 1:
        outname = temp_directory + node["name"] + ".fasta"
    else:
        # If only one sequence, output it to aliged directory and break
        outname = outdir + "/aligned_gene_sequences/" + node["name"] + ".fasta"
        SeqIO.write(output_sequences, outname, 'fasta')
        return None
    # check to see if filename is too long
    if len(outname) >= 248:
        outname = outname[:248] + ".fasta"
    # Write them to disk
    SeqIO.write(output_sequences, outname, 'fasta')
    return outname


def get_alignment_commands(fastafile_name, outdir, aligner):
    if aligner == "def":
        command = MafftCommandline(input=fastafile_name,
                                   auto=True,
                                   nuc=True)
    elif aligner == "ref":
        ref_file, seq_file = fastafile_name
        outfile = outdir + "aligned_gene_sequences/" + seq_file.split('/')[-1].split('.')[0] + '.aln.fas'
        command = ["mafft", "--6merpair", "--addfragments",
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
