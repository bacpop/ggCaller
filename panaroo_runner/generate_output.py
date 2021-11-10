from functools import partial
from ggCaller.shared_memory import *
import networkx as nx
import numpy as np
from Bio import AlignIO
import itertools as iter

from .generate_alignments import *

def output_aa_sequence(node_pair):
    # unpack node_pair
    node = node_pair[1]

    ref_output_sequences = []

    # iterate over centroids to generate fasta files
    for i in range(0, len(node["centroid"])):
        name = str(node_pair[0]) + ";" + node["centroid"][i]
        ref_output_sequences.append(SeqRecord(Seq(node["protein"][i]), id=name, description=""))

    return ref_output_sequences


def output_dna_sequence(node_pair, isolate_list, temp_directory, outdir, shd_arr_tup, high_scoring_ORFs, overlap,
                        ref_aln):
    # load shared memory items
    existing_shm = shared_memory.SharedMemory(name=shd_arr_tup.name)
    shd_arr = np.ndarray(shd_arr_tup.shape, dtype=shd_arr_tup.dtype, buffer=existing_shm.buf)

    # unpack node_pair
    node = node_pair[1]

    # Counter for the number of sequences to
    seq_no = 0

    # get outname for reference alignment file
    ref_outname = None

    # Get the name of the sequences for the gene of interest
    sequence_ids = node["seqIDs"]
    centroid_sequence_ids = set(node["centroid"])
    ref_output_sequences = []
    output_sequences = []

    # if reference-guided alignment being done, separate centroids and other sequences
    if ref_aln:
        for i in range(0, len(node["centroid"])):
            member = int(node["centroid"][i].split("_")[0])
            isolate_name = isolate_list[member].replace(";",
                                                        "") + ";" + node["centroid"][i]
            ref_output_sequences.append(SeqRecord(Seq(node["dna"][i]), id=isolate_name, description=""))
        centroid_no = len(ref_output_sequences)
        # Put gene of interest sequences in a generator, with corrected isolate names
        ref_output_sequences_gen = (x for x in ref_output_sequences)
        if len(sequence_ids) == centroid_no and centroid_no == 1:
            # If only one sequence, output it to aligned directory and break
            # if no other sequences, then just output with no alignment
            ref_outname = outdir + "/aligned_gene_sequences/" + node["name"] + ".aln.fas"
            if len(ref_outname) >= 248:
                ref_outname = ref_outname[:248] + ".fasta"
            SeqIO.write(ref_output_sequences_gen, ref_outname, 'fasta')
            return None, None
        else:
            # if centroid is on it's own, give name aln for aligned, otherwise ref
            if centroid_no > 1:
                ref_outname = temp_directory + node["name"] + "_ref.fasta"
            else:
                ref_outname = temp_directory + node["name"] + "_ref.aln.fas"
        if len(ref_outname) >= 248:
            ref_outname = ref_outname[:248] + ".fasta"
        SeqIO.write(ref_output_sequences_gen, ref_outname, 'fasta')

    # Look for gene sequences among all genes (from memory)
    for seq in sequence_ids:
        # check if reference true and seqID is centroid. Is so, pass
        if ref_aln and seq in centroid_sequence_ids:
            continue
        member = int(seq.split('_')[0])
        ORF_ID = int(seq.split('_')[-1])
        isolate_name = isolate_list[member].replace(";",
                                                    "") + ";" + seq
        # generate DNA sequence
        ORFNodeVector = high_scoring_ORFs[member][ORF_ID]
        CDS = shd_arr[0].generate_sequence(ORFNodeVector[0], ORFNodeVector[1], overlap)

        output_sequences.append(
            SeqRecord(Seq(CDS), id=isolate_name, description=""))
        seq_no += 1
    # Put gene of interest sequences in a generator, with corrected isolate names
    output_sequences = (x for x in output_sequences)
    # set filename to gene name, if more than one sequence to be aligned
    if (not ref_aln and seq_no > 1) or (ref_aln and seq_no > 0):
        outname = temp_directory + node["name"] + ".fasta"
    else:
        # if number sequences is 0, do not output
        if seq_no > 0:
            # If only one sequence, output it to aligned directory and break
            outname = outdir + "/aligned_gene_sequences/" + node["name"] + ".aln.fas"
            if len(outname) >= 248:
                outname = outname[:248] + ".fasta"
            SeqIO.write(output_sequences, outname, 'fasta')
        return None, ref_outname

    # Write them to disk
    if len(outname) >= 248:
        outname = outname[:248] + ".fasta"
    SeqIO.write(output_sequences, outname, 'fasta')
    return outname, ref_outname


def generate_roary_gene_presence_absence(G, mems_to_isolates, orig_ids,
                                         ids_len_stop, output_dir):
    # arange isolates
    isolates = []
    mems_to_index = {}
    for i, mem in enumerate(mems_to_isolates):
        isolates.append(mems_to_isolates[mem])
        mems_to_index[str(mem)] = i

    # generate file
    with open(output_dir + "gene_presence_absence_roary.csv", 'w') as roary_csv_outfile, \
            open(output_dir + "gene_presence_absence.csv", 'w') as csv_outfile, \
            open(output_dir + "gene_presence_absence.Rtab", 'w') as Rtab_outfile:
        header = [
                     "Gene", "Non-unique Gene name", "Annotation", "No. isolates",
                     "No. sequences", "Avg sequences per isolate", "Genome Fragment",
                     "Order within Fragment", "Accessory Fragment",
                     "Accessory Order with Fragment", "QC", "Min group size nuc",
                     "Max group size nuc", "Avg group size nuc"
                 ] + isolates
        roary_csv_outfile.write(",".join(header) + "\n")
        csv_outfile.write(",".join(header[:3] + isolates) + "\n")
        Rtab_outfile.write("\t".join((["Gene"] + isolates)) + "\n")

        # Iterate through coponents writing out to file
        used_gene_names = set([""])
        unique_id_count = 0
        frag = 0
        entry_list = []
        entry_ext_list = []
        pres_abs_list = []
        entry_sizes = []
        entry_count = 0
        for component in nx.connected_components(G):
            frag += 1
            count = 0
            for node in component:
                count += 1
                len_mode = max(G.nodes[node]['lengths'],
                               key=G.nodes[node]['lengths'].count)
                name = '~~~'.join([
                    gn for gn in G.nodes[node]['annotation'].strip().strip(
                        ';').split(';') if gn != ''
                ])
                name = ''.join(e for e in name
                               if e.isalnum() or e in ["_", "~", "|"])
                name = name.replace("|", "-")
                if name not in used_gene_names:
                    entry = [name]
                    used_gene_names.add(name)
                    G.nodes[node]['name'] = name
                else:
                    G.nodes[node]['name'] = "group_" + str(unique_id_count)
                    entry = [G.nodes[node]['name']]
                    unique_id_count += 1
                entry.append(G.nodes[node]['annotation'])
                entry.append(G.nodes[node]['description'])
                entry.append(G.nodes[node]['size'])
                entry.append(len(G.nodes[node]['seqIDs']))
                entry.append((1.0 * len(G.nodes[node]['seqIDs'])) /
                             G.nodes[node]['size'])
                entry.append(frag)
                entry.append(count)
                entry += ["", "", ""]
                entry.append(np.min(G.nodes[node]['lengths']))
                entry.append(np.max(G.nodes[node]['lengths']))
                entry.append(np.mean(G.nodes[node]['lengths']))
                pres_abs = [""] * len(isolates)
                pres_abs_ext = [""] * len(isolates)
                entry_size = 0
                for seq in G.nodes[node]['seqIDs']:
                    sample_id = mems_to_index["_".join(seq.split("_")[:-2])]
                    if pres_abs[
                        sample_id] == "":  # ensures we only take the first one
                        if seq in orig_ids:
                            pres_abs[sample_id] = orig_ids[seq]
                            pres_abs_ext[sample_id] = orig_ids[seq]
                        else:
                            pres_abs[sample_id] = seq
                            pres_abs_ext[sample_id] = seq
                        entry_size += 1
                    else:
                        # this is similar to PIRATE output
                        if seq in orig_ids:
                            pres_abs[sample_id] += ";" + orig_ids[seq]
                            pres_abs_ext[sample_id] += ";" + orig_ids[seq]
                        else:
                            pres_abs[sample_id] += ";" + seq
                            pres_abs_ext[sample_id] += ";" + seq
                    if (abs(ids_len_stop[seq][0] - len_mode) /
                        len_mode) > (0.05 * len_mode):
                        pres_abs_ext[sample_id] += "_len"
                    if ids_len_stop[seq][1]:
                        pres_abs_ext[sample_id] += "_stop"

                entry += pres_abs
                entry_list.append(entry)
                entry_ext_list.append(entry[:3] + pres_abs_ext)
                pres_abs_list.append(pres_abs)
                entry_sizes.append((entry_size, entry_count))
                entry_count += 1

        # sort so that the most common genes are first (as in roary)
        entry_sizes = sorted(entry_sizes, reverse=True)
        for s, i in entry_sizes:
            roary_csv_outfile.write(",".join([str(e)
                                              for e in entry_list[i]]) + "\n")
            csv_outfile.write(",".join([str(e)
                                        for e in entry_ext_list[i]]) + "\n")
            Rtab_outfile.write(entry_list[i][0] + "\t")
            Rtab_outfile.write("\t".join(
                (["0" if e == "" else "1" for e in pres_abs_list[i]])) + "\n")

    return G


def generate_pan_genome_reference(G, output_dir, split_paralogs=False):
    # need to treat paralogs differently?
    centroids = set()
    records = []

    for node in G.nodes():
        if not split_paralogs and G.nodes[node]['centroid'][0] in centroids:
            continue
        records.append(
            SeqRecord(Seq(max(G.nodes[node]['dna'], key=lambda x: len(x))),
                      id=G.nodes[node]['name'],
                      description=""))
        for centroid in G.nodes[node]['centroid']:
            centroids.add(centroid)

    with open(output_dir + "pan_genome_reference.fa", 'w') as outfile:
        SeqIO.write(records, outfile, "fasta")

    return


def generate_common_struct_presence_absence(G,
                                            output_dir,
                                            mems_to_isolates,
                                            min_variant_support=2):
    # arange isolates
    isolates = []
    members = []
    for mem in mems_to_isolates:
        isolates.append(mems_to_isolates[mem])
        members.append(mem)

    struct_variants = {}
    for node in G.nodes():
        if G.degree[node] < 3: continue  # skip as linear
        for path in iter.combinations(G.edges(node), 2):
            in_both = (G[path[0][0]][path[0][1]]['members']
                       & G[path[1][0]][path[1][1]]['members'])
            if len(in_both) >= min_variant_support:
                struct_variants[(path[0][0], path[0][1], path[1][1])] = in_both

    header = []
    for variant in struct_variants:
        header.append("-".join([
            G.nodes[variant[1]]['name'], G.nodes[variant[0]]['name'],
            G.nodes[variant[2]]['name']
        ]))

    with open(output_dir + "struct_presence_absence.Rtab",
              'w') as Rtab_outfile:
        Rtab_outfile.write("\t".join((["Gene"] + isolates)) + "\n")
        for h, variant in zip(header, struct_variants):
            variant_calls = [h]
            for member in members:
                if member in struct_variants[variant]:
                    variant_calls.append("1")
                else:
                    variant_calls.append("0")
            Rtab_outfile.write("\t".join(variant_calls) + "\n")

    return


def generate_pan_genome_alignment(G, temp_dir, output_dir, threads, aligner,
                                  isolates, shd_arr_tup, high_scoring_ORFs, overlap, pool, ref_aln, verbose):
    unaligned_sequence_files = []
    unaligned_reference_files = []
    # Make a folder for the output alignments
    try:
        os.mkdir(output_dir + "aligned_gene_sequences")
    except FileExistsError:
        None

    # Multithread writing gene sequences to disk (temp directory) so aligners can find them
    for outname, ref_outname in pool.map(partial(output_dna_sequence, isolate_list=isolates,
                                                 temp_directory=temp_dir, outdir=output_dir, shd_arr_tup=shd_arr_tup,
                                                 high_scoring_ORFs=high_scoring_ORFs, overlap=overlap, ref_aln=ref_aln),
                                         G.nodes(data=True)):
        unaligned_sequence_files.append(outname)
        unaligned_reference_files.append(ref_outname)
    if ref_aln:
        # centroid files with paired sequence files
        ref_seq_pairs = [
            (temp_dir + unaligned_reference_files[i].split("/")[-1].split("_ref.")[0] + "_ref.aln.fas",
             unaligned_sequence_files[i])
            for i in range(0, len(unaligned_reference_files)) if unaligned_sequence_files[i] is not None
                                                                 and unaligned_reference_files[i] is not None]
        # reference files with no paired sequence files
        ref_seq_singles = [
            temp_dir + unaligned_reference_files[i].split("/")[-1].split("_ref.")[0] + "_ref.aln.fas"
            for i in range(0, len(unaligned_reference_files))
            if unaligned_sequence_files[i] is None and unaligned_reference_files[i] is not None]

    # remove single sequence files
    unaligned_sequence_files = filter(None, unaligned_sequence_files)
    unaligned_reference_files = filter(None, unaligned_reference_files)
    # Get Biopython command calls for each output gene sequences
    # check if ref_alignment being done
    if ref_aln:
        # conduct MSA on reference files, first using standard MSA
        commands = [
            get_alignment_commands(fastafile, None, aligner[:-4], threads)
            for fastafile in unaligned_reference_files if "_ref.aln.fas" not in fastafile
        ]
        if verbose: print("Aligning centroids...")
        multi_align_sequences(commands, temp_dir, threads, aligner[:-4])

        # move any centroid alignments which do not have associated sequence files
        for file in ref_seq_singles:
            os.rename(file, output_dir + "aligned_gene_sequences/" + file.split("/")[-1].split("_ref.")[0] + '.aln.fas')

        # repeat with reference-guided alignment
        commands = [
            get_alignment_commands(fastapair, output_dir, aligner, threads)
            for fastapair in ref_seq_pairs
        ]
        if verbose: print("Aligning remaining sequences...")
        multi_align_sequences(commands, output_dir + "aligned_gene_sequences/",
                              threads, aligner)
    else:
        commands = [
            get_alignment_commands(fastafile, output_dir, aligner, threads)
            for fastafile in unaligned_sequence_files
        ]
        # Run these commands in a multi-threaded way
        multi_align_sequences(commands, output_dir + "aligned_gene_sequences/",
                              threads, aligner)
    return


def get_unannotated_nodes(G):
    # Get the unannotated nodes based on bitscore
    unannotated_nodes = []
    for node in G.nodes(data=True):
        if float(G.nodes[node[0]]["bitscore"]) == 0:
            unannotated_nodes.append(node)
    return unannotated_nodes


def get_core_gene_nodes(G, threshold, num_isolates):
    # Get the core genes based on percent threshold
    core_nodes = []
    for node in G.nodes(data=True):
        if float(G.nodes[node[0]]["size"]) / float(num_isolates) > threshold:
            core_nodes.append(node)
    return core_nodes


def concatenate_core_genome_alignments(core_names, output_dir):
    alignments_dir = output_dir + "/aligned_gene_sequences/"
    # Open up each alignment that is assosciated with a core node
    alignment_filenames = os.listdir(alignments_dir)
    core_filenames = [
        x for x in alignment_filenames if x.split('.')[0] in core_names
    ]
    # Read in all these alignments
    gene_alignments = []
    isolates = set()
    for filename in core_filenames:
        gene_name = os.path.splitext(os.path.basename(filename))[0]
        alignment = AlignIO.read(alignments_dir + filename, 'fasta')
        gene_dict = {}
        for record in alignment:
            genome_id = record.id.split(";")[0]
            if genome_id in gene_dict:
                if str(record.seq).count("-") < str(
                        gene_dict[genome_id][1]).count("-"):
                    gene_dict[genome_id] = (record.id, record.seq)
            else:
                gene_dict[genome_id] = (record.id, record.seq)
            gene_length = len(record.seq)
            isolates.add(genome_id)
        gene_alignments.append((gene_name, gene_dict, gene_length))
    # Combine them
    isolate_aln = []
    for iso in isolates:
        seq = ""
        for gene in gene_alignments:
            if iso in gene[1]:
                seq += gene[1][iso][1]
            else:
                seq += "-" * gene[2]
        isolate_aln.append(SeqRecord(seq, id=iso, description=""))

    # Write out the two output files
    SeqIO.write(isolate_aln, output_dir + 'core_gene_alignment.aln', 'fasta')

    write_alignment_header(gene_alignments, output_dir)
    return core_filenames


def generate_core_genome_alignment(G, temp_dir, output_dir, threads, aligner,
                                   isolates, threshold, num_isolates, shd_arr_tup, high_scoring_ORFs,
                                   overlap, pool, ref_aln, verbose):
    unaligned_sequence_files = []
    unaligned_reference_files = []
    # Make a folder for the output alignments
    try:
        os.mkdir(output_dir + "aligned_gene_sequences")
    except FileExistsError:
        None

    # Get core nodes
    core_genes = get_core_gene_nodes(G, threshold, num_isolates)
    core_gene_names = [G.nodes[x[0]]["name"] for x in core_genes]

    # Multithread writing gene sequences to disk (temp directory) so aligners can find them
    for outname, ref_outname in pool.map(partial(output_dna_sequence, isolate_list=isolates,
                                                 temp_directory=temp_dir, outdir=output_dir, shd_arr_tup=shd_arr_tup,
                                                 high_scoring_ORFs=high_scoring_ORFs, overlap=overlap, ref_aln=ref_aln),
                                         core_genes):
        unaligned_sequence_files.append(outname)
        unaligned_reference_files.append(ref_outname)
    if ref_aln:
        # centroid files with paired sequence files
        ref_seq_pairs = [
            (temp_dir + unaligned_reference_files[i].split("/")[-1].split("_ref.")[0] + "_ref.aln.fas",
             unaligned_sequence_files[i])
            for i in range(0, len(unaligned_reference_files)) if unaligned_sequence_files[i] is not None
                                                                 and unaligned_reference_files[i] is not None]
        # reference files with no paired sequence files
        ref_seq_singles = [
            temp_dir + unaligned_reference_files[i].split("/")[-1].split("_ref.")[0] + "_ref.aln.fas"
            for i in range(0, len(unaligned_reference_files))
            if unaligned_sequence_files[i] is None and unaligned_reference_files[i] is not None]

    # remove single sequence files
    unaligned_sequence_files = filter(None, unaligned_sequence_files)
    unaligned_reference_files = filter(None, unaligned_reference_files)
    # Get Biopython command calls for each output gene sequences
    # check if ref_alignment being done
    if ref_aln:
        # conduct MSA on reference files, first using standard MSA
        commands = [
            get_alignment_commands(fastafile, None, aligner[:-4], threads)
            for fastafile in unaligned_reference_files if "_ref.aln.fas" not in fastafile
        ]
        if verbose: print("Aligning centroids...")
        multi_align_sequences(commands, temp_dir, threads, aligner[:-4])

        # move any centroid alignments which do not have associated sequence files
        for file in ref_seq_singles:
            os.rename(file, output_dir + "aligned_gene_sequences/" + file.split("/")[-1].split("_ref.")[0] + '.aln.fas')

        # repeat with reference-guided alignment
        commands = [
            get_alignment_commands(fastapair, output_dir, aligner, threads)
            for fastapair in ref_seq_pairs
        ]
        if verbose: print("Aligning remaining sequences...")
        multi_align_sequences(commands, output_dir + "aligned_gene_sequences/",
                              threads, aligner)
    else:
        # Get alignment commands
        commands = [
            get_alignment_commands(fastafile, output_dir, aligner, threads)
            for fastafile in unaligned_sequence_files
        ]
        # Run alignment commands
        multi_align_sequences(commands, output_dir + "aligned_gene_sequences/",
                              threads, aligner)
    # Concatenate them together to produce the two output files
    concatenate_core_genome_alignments(core_gene_names, output_dir)
    return

def generate_summary_stats(output_dir):
    with open(output_dir + "gene_presence_absence_roary.csv", 'r') as inhandle:
        gene_presence_absence = inhandle.read().splitlines()[1:]
    noSamples = len(gene_presence_absence[0].split(',')) - 14
    # Layout categories
    noCore = 0
    noSoftCore = 0
    noShell = 0
    noCloud = 0
    total_genes = 0
    # Iterate through GPA and summarise
    for gene in gene_presence_absence:
        proportion_present = float(gene.split(',')[3]) / noSamples * 100.0
        if proportion_present >= 99:
            noCore += 1
        elif proportion_present >= 95:
            noSoftCore += 1
        elif proportion_present >= 15:
            noShell += 1
        else:
            noCloud += 1
        total_genes += 1

    # write output
    with open(output_dir + "summary_statistics.txt", 'w') as outfile:
        output = ("Core genes\t(99% <= strains <= 100%)\t" + str(noCore) +
                  "\n" + "Soft core genes\t(95% <= strains < 99%)\t" +
                  str(noSoftCore) + "\n" +
                  "Shell genes\t(15% <= strains < 95%)\t" + str(noShell) +
                  "\n" + "Cloud genes\t(0% <= strains < 15%)\t" +
                  str(noCloud) + "\n" +
                  "Total genes\t(0% <= strains <= 100%)\t" + str(total_genes))
        outfile.write(output)

    return True
