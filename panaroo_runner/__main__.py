import os, sys
import tempfile
from Bio import SeqIO
import shutil
import networkx as nx
import argparse
import textwrap
import ast

# panaroo scripts
from panaroo.isvalid import *
from .generate_network import *
from panaroo.cdhit import check_cdhit_version
from panaroo.cdhit import run_cdhit
from .find_missing import *
# from panaroo.generate_alignments import check_aligner_install
from intbitset import intbitset
from ggCaller.shared_memory import *

# custom panaroo scripts
from .clean_network import *
from .generate_output import *
from .annotate import *


# debugging scripts
from .cdhit_align import *


class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            lines = []
            for l in text[2:].splitlines():
                if l == "":
                    lines += [""]
                else:
                    lines += textwrap.wrap(l, width=55)
            return lines
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)


def run_panaroo(pool, shd_arr_tup, high_scoring_ORFs, high_scoring_ORF_edges, cluster_id_list, cluster_dict, overlap,
                input_colours, output_dir, temp_dir, verbose, n_cpu, length_outlier_support_proportion, identity_cutoff,
                family_threshold, min_trailing_support, trailing_recursive, clean_edges, edge_support_threshold,
                merge_para, aln, alr, core, min_edge_support_sv, all_seq_in_graph, is_ref, write_idx, kmer, repeat,
                remove_by_consensus, search_radius, refind_prop_match, annotate, evalue, annotation_db, hmm_db):

    # load shared memory items
    existing_shm = shared_memory.SharedMemory(name=shd_arr_tup.name)
    shd_arr = np.ndarray(shd_arr_tup.shape, dtype=shd_arr_tup.dtype, buffer=existing_shm.buf)

    # check if reference-guided alignment specified
    ref_aln = False
    if "ref" in alr:
        ref_aln = True

    # Check cd-hit is installed
    check_cdhit_version()
    # Make sure aligner is installed if alignment requested
    if aln != None:
        if ref_aln:
            check_aligner_install(alr[:-4])
        else:
            check_aligner_install(alr)

    if verbose:
        print("Generating initial network...")

    # generate network from clusters and adjacency information
    G, centroid_contexts, seqid_to_centroid = generate_network(shd_arr[0], high_scoring_ORFs, high_scoring_ORF_edges,
                                                               cluster_id_list, cluster_dict, overlap, all_seq_in_graph)


    # merge paralogs
    if verbose:
        print("Processing paralogs...")
    G = collapse_paralogs(G, centroid_contexts, seqid_to_centroid, quiet=(not verbose))

    # write out pre-filter graph in GML format
    for node in G.nodes():
        G.nodes[node]['size'] = len(G.nodes[node]['members'])
        G.nodes[node]['genomeIDs'] = ";".join(
            [str(m) for m in G.nodes[node]['members']])
        G.nodes[node]['geneIDs'] = ";".join([m for m in G.nodes[node]['seqIDs']])
        G.nodes[node]['degrees'] = G.degree[node]
    for edge in G.edges():
        G.edges[edge[0], edge[1]]['genomeIDs'] = ";".join(
            [str(m) for m in G.edges[edge[0], edge[1]]['members']])
    nx.write_gml(G,
                 output_dir + "pre_filt_graph.gml",
                 stringizer=custom_stringizer)

    if verbose:
        print("collapse mistranslations...")

    # clean up translation errors
    G = collapse_families(G,
                          seqid_to_centroid=seqid_to_centroid,
                          outdir=temp_dir,
                          dna_error_threshold=0.98,
                          correct_mistranslations=True,
                          length_outlier_support_proportion=length_outlier_support_proportion,
                          n_cpu=n_cpu,
                          quiet=(not verbose))[0]

    if annotate:
        if verbose:
            print("annotating gene families...")

        # create directory for annotation
        annotation_temp_dir = os.path.join(temp_dir, "annotation")
        if not os.path.exists(annotation_temp_dir):
            os.mkdir(annotation_temp_dir)
        # make sure trailing forward slash is present
        annotation_temp_dir = os.path.join(annotation_temp_dir, "")

        # generate annotations
        G, annotation_list = iterative_annotation_search(G, annotation_temp_dir, annotation_db, hmm_db, evalue,
                                                         len(input_colours), n_cpu, pool)

    if verbose:
        print("collapse gene families...")

    # collapse gene families
    G, distances_bwtn_centroids, centroid_to_index = collapse_families(
        G,
        seqid_to_centroid=seqid_to_centroid,
        outdir=temp_dir,
        family_threshold=family_threshold,
        correct_mistranslations=False,
        length_outlier_support_proportion=length_outlier_support_proportion,
        n_cpu=n_cpu,
        quiet=(not verbose))

    if verbose:
        print("trimming contig ends...")

    # re-trim low support trailing ends
    G = trim_low_support_trailing_ends(G,
                                       min_support=min_trailing_support,
                                       max_recursive=trailing_recursive)

    if verbose:
        print("refinding genes...")

    # find genes that Prokka has missed
    G, high_scoring_ORFs = find_missing(G,
                                        shd_arr_tup,
                                        high_scoring_ORFs,
                                        is_ref=is_ref,
                                        write_idx=write_idx,
                                        kmer=kmer,
                                        repeat=repeat,
                                        isolate_names=input_colours,
                                        remove_by_consensus=remove_by_consensus,
                                        search_radius=search_radius,
                                        prop_match=refind_prop_match,
                                        pairwise_id_thresh=identity_cutoff,
                                        merge_id_thresh=max(0.8, family_threshold),
                                        pool=pool,
                                        n_cpu=n_cpu,
                                        verbose=verbose)

    # remove edges that are likely due to misassemblies (by consensus)

    # merge again in case refinding has resolved issues
    if verbose:
        print("collapse gene families with refound genes...")
    G = collapse_families(G,
                          seqid_to_centroid=seqid_to_centroid,
                          outdir=temp_dir,
                          family_threshold=family_threshold,
                          correct_mistranslations=False,
                          length_outlier_support_proportion=length_outlier_support_proportion,
                          n_cpu=n_cpu,
                          quiet=(not verbose),
                          distances_bwtn_centroids=distances_bwtn_centroids,
                          centroid_to_index=centroid_to_index)[0]

    if clean_edges:
        G = clean_misassembly_edges(
            G, edge_support_threshold=edge_support_threshold)

    # if requested merge paralogs
    if merge_para:
        G = merge_paralogs(G, seqid_to_centroid)

    # generate list of input isolate names
    isolate_names = [
        os.path.splitext(os.path.basename(x))[0] for x in input_colours
    ]
    G.graph['isolateNames'] = isolate_names
    mems_to_isolates = {}
    for i, iso in enumerate(isolate_names):
        mems_to_isolates[i] = iso

    if verbose:
        print("writing output...")

    # write out roary like gene_presence_absence.csv
    # get original annotation IDs, lengths and whether or
    # not an internal stop codon is present
    orig_ids = {}
    ids_len_stop = {}
    for node in G.nodes():
        for sid in G.nodes[node]['seqIDs']:
            orig_ids[sid] = sid
            mem = int(sid.split("_")[0])
            ORF_ID = int(sid.split("_")[-1])
            ORFNodeVector = high_scoring_ORFs[mem][ORF_ID]
            # determine if gene is refound. If it is, then determine if premature stop codon present
            if (ORF_ID < 0):
                ids_len_stop[sid] = (ORFNodeVector[2] / 3, ORFNodeVector[-2])
            else:
                ids_len_stop[sid] = (ORFNodeVector[2] / 3, False)
            if annotate:
                # here, add each sequence to its respective contig for each gff file.
                if sid in annotation_list[mem]:
                    annotation_list[mem][sid].append(ORFNodeVector[-1])
                else:
                    annotation_list[mem][sid] = [None, "hypothetical protein", 0, ORFNodeVector[-1]]

    # write roary output and summary stats file
    G = generate_roary_gene_presence_absence(G,
                                             mems_to_isolates=mems_to_isolates,
                                             orig_ids=orig_ids,
                                             ids_len_stop=ids_len_stop,
                                             output_dir=output_dir)
    # # Write out presence_absence summary
    # generate_summary_stats(output_dir=output_dir)

    # write pan genome reference fasta file
    generate_pan_genome_reference(G,
                                  output_dir=output_dir,
                                  split_paralogs=False)

    # write out common structural differences in a matrix format
    generate_common_struct_presence_absence(
        G,
        output_dir=output_dir,
        mems_to_isolates=mems_to_isolates,
        min_variant_support=min_edge_support_sv)

    # Write out core/pan-genome alignments
    # determine if reference-guided alignment being done
    if aln == "pan":
        if verbose: print("generating pan genome MSAs...")
        generate_pan_genome_alignment(G, temp_dir, output_dir, n_cpu,
                                      alr, isolate_names, shd_arr_tup,
                                      high_scoring_ORFs, overlap, pool, ref_aln, verbose)
        core_nodes = get_core_gene_nodes(G, core, len(input_colours))
        concatenate_core_genome_alignments(core_nodes, output_dir)
    elif aln == "core":
        if verbose: print("generating core genome MSAs...")
        generate_core_genome_alignment(G, temp_dir, output_dir,
                                       n_cpu, alr, isolate_names,
                                       core, len(input_colours), shd_arr_tup,
                                       high_scoring_ORFs, overlap, pool, ref_aln, verbose)

    # add helpful attributes and write out graph in GML format
    for node in G.nodes():
        G.nodes[node]['size'] = len(G.nodes[node]['members'])
        G.nodes[node]['centroid'] = ";".join(G.nodes[node]['centroid'])
        G.nodes[node]['dna'] = ";".join(conv_list(G.nodes[node]['dna']))
        G.nodes[node]['protein'] = ";".join(conv_list(
            G.nodes[node]['protein']))
        G.nodes[node]['genomeIDs'] = ";".join(
            [str(m) for m in G.nodes[node]['members']])
        G.nodes[node]['geneIDs'] = ";".join(G.nodes[node]['seqIDs'])
        G.nodes[node]['degrees'] = G.degree[node]
        G.nodes[node]['members'] = list(G.nodes[node]['members'])
        G.nodes[node]['seqIDs'] = list(G.nodes[node]['seqIDs'])

    for edge in G.edges():
        G.edges[edge[0], edge[1]]['genomeIDs'] = ";".join(
            [str(m) for m in G.edges[edge[0], edge[1]]['members']])
        G.edges[edge[0],
                edge[1]]['members'] = list(G.edges[edge[0],
                                                   edge[1]]['members'])

    nx.write_gml(G, output_dir + "final_graph.gml")

    # remove temporary directory
    shutil.rmtree(temp_dir)

    return

#
# if __name__ == '__main__':
#     main()
