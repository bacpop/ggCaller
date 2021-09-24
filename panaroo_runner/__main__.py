import os, sys
import tempfile
from Bio import SeqIO
import shutil
import networkx as nx
import argparse
import textwrap
import ast

from panaroo.isvalid import *
from panaroo.set_default_args import set_default_args
from panaroo.prokka import process_prokka_input
from .generate_network import generate_network
from panaroo.generate_output import *
from panaroo.clean_network import *
from panaroo.find_missing import find_missing
from panaroo.generate_alignments import check_aligner_install
from intbitset import intbitset

from panaroo.__init__ import __version__


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


def run_panaroo(DBG, high_scoring_ORFs, high_scoring_ORF_edges, cluster_id_list, cluster_dict, overlap, output_dir,
                verbose, all_seq_in_graph=True):
    # args = get_options(sys.argv[1:])
    # Check cd-hit is installed
    # check_cdhit_version()
    # Make sure aligner is installed if alignment requested
    # if args.aln != None:
    #     check_aligner_install(args.alr)

    # create directory if it isn't present already
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    # make sure trailing forward slash is present
    output_dir = os.path.join(output_dir, "")
    # Create temporary directory
    temp_dir = os.path.join(tempfile.mkdtemp(dir=output_dir), "")

    # check if input is a file containing filenames
    # if len(args.input_files) == 1:
    #     files = []
    #     with open(args.input_files[0], 'r') as infile:
    #         for line in infile:
    #             files.append(line.strip())
    #     args.input_files = files

    # if verbose:
    #     print("pre-processing gff3 files...")

    # convert input GFF3 files into summary files
    # process_prokka_input(args.input_files, output_dir,
    #                      args.filter_invalid, (not verbose),
    #                      args.n_cpu, args.table)

    # Cluster protein sequences using cdhit
    # cd_hit_out = output_dir + "combined_protein_cdhit_out.txt"
    # run_cdhit(input_file=output_dir + "combined_protein_CDS.fasta",
    #           output_file=cd_hit_out,
    #           id=args.id,
    #           s=args.len_dif_percent,
    #           quiet=(not verbose),
    #           n_cpu=args.n_cpu)

    # if verbose:
    print("Generating initial network...")

    # generate network from clusters and adjacency information
    G, centroid_contexts, seqid_to_centroid = generate_network(DBG, high_scoring_ORFs, high_scoring_ORF_edges,
                                                               cluster_id_list, cluster_dict, overlap, all_dna=False)

    test = 1

    ### commented from here
    # merge paralogs
    if verbose:
        print("Processing paralogs...")
    G = collapse_paralogs(G, centroid_contexts, quiet=(not verbose))

    # write out pre-filter graph in GML format
    for node in G.nodes():
        G.nodes[node]['size'] = len(G.nodes[node]['members'])
        G.nodes[node]['genomeIDs'] = ";".join(
            [str(m) for m in G.nodes[node]['members']])
        G.nodes[node]['geneIDs'] = ";".join([str(m) for m in G.nodes[node]['seqIDs']])
        G.nodes[node]['degrees'] = G.degree[node]
    for edge in G.edges():
        G.edges[edge[0], edge[1]]['genomeIDs'] = ";".join(
            [str(m) for m in G.edges[edge[0], edge[1]]['members']])
    nx.write_gml(G,
                 output_dir + "pre_filt_graph.gml",
                 stringizer=custom_stringizer)

    # if verbose:
    #     print("collapse mistranslations...")
    #
    # # clean up translation errors
    # G = collapse_families(G,
    #                       seqid_to_centroid=seqid_to_centroid,
    #                       outdir=temp_dir,
    #                       dna_error_threshold=0.98,
    #                       correct_mistranslations=True,
    #                       length_outlier_support_proportion=args.
    #                       length_outlier_support_proportion,
    #                       n_cpu=args.n_cpu,
    #                       quiet=(not verbose))[0]
    #
    # if verbose:
    #     print("collapse gene families...")
    #
    # # collapse gene families
    # G, distances_bwtn_centroids, centroid_to_index = collapse_families(
    #     G,
    #     seqid_to_centroid=seqid_to_centroid,
    #     outdir=temp_dir,
    #     family_threshold=args.family_threshold,
    #     correct_mistranslations=False,
    #     length_outlier_support_proportion=args.
    #         length_outlier_support_proportion,
    #     n_cpu=args.n_cpu,
    #     quiet=(not verbose))
    #
    # if verbose:
    #     print("trimming contig ends...")
    #
    # # re-trim low support trailing ends
    # G = trim_low_support_trailing_ends(G,
    #                                    min_support=args.min_trailing_support,
    #                                    max_recursive=args.trailing_recursive)
    #
    # if verbose:
    #     print("refinding genes...")
    #
    # # find genes that Prokka has missed
    # G = find_missing(G,
    #                  args.input_files,
    #                  dna_seq_file=output_dir + "combined_DNA_CDS.fasta",
    #                  prot_seq_file=output_dir +
    #                                "combined_protein_CDS.fasta",
    #                  gene_data_file=output_dir + "gene_data.csv",
    #                  remove_by_consensus=args.remove_by_consensus,
    #                  search_radius=args.search_radius,
    #                  prop_match=args.refind_prop_match,
    #                  pairwise_id_thresh=args.id,
    #                  merge_id_thresh=max(0.8, args.family_threshold),
    #                  n_cpu=args.n_cpu,
    #                  verbose=verbose)
    #
    # # remove edges that are likely due to misassemblies (by consensus)
    #
    # # merge again in case refinding has resolved issues
    # if verbose:
    #     print("collapse gene families with refound genes...")
    # G = collapse_families(G,
    #                       seqid_to_centroid=seqid_to_centroid,
    #                       outdir=temp_dir,
    #                       family_threshold=args.family_threshold,
    #                       correct_mistranslations=False,
    #                       length_outlier_support_proportion=args.
    #                       length_outlier_support_proportion,
    #                       n_cpu=args.n_cpu,
    #                       quiet=(not verbose),
    #                       distances_bwtn_centroids=distances_bwtn_centroids,
    #                       centroid_to_index=centroid_to_index)[0]
    #
    # if args.clean_edges:
    #     G = clean_misassembly_edges(
    #         G, edge_support_threshold=args.edge_support_threshold)
    #
    # # if requested merge paralogs
    # if args.merge_paralogs:
    #     G = merge_paralogs(G)
    #
    # isolate_names = [
    #     os.path.splitext(os.path.basename(x))[0] for x in args.input_files
    # ]
    # G.graph['isolateNames'] = isolate_names
    # mems_to_isolates = {}
    # for i, iso in enumerate(isolate_names):
    #     mems_to_isolates[i] = iso
    #
    # if verbose:
    #     print("writing output...")
    #
    # # write out roary like gene_presence_absence.csv
    # # get original annotaiton IDs, lengts and whether or
    # # not an internal stop codon is present
    # orig_ids = {}
    # ids_len_stop = {}
    # with open(output_dir + "gene_data.csv", 'r') as infile:
    #     next(infile)
    #     for line in infile:
    #         line = line.split(",")
    #         orig_ids[line[2]] = line[3]
    #         ids_len_stop[line[2]] = (len(line[4]), "*" in line[4][1:-3])
    #
    # G = generate_roary_gene_presence_absence(G,
    #                                          mems_to_isolates=mems_to_isolates,
    #                                          orig_ids=orig_ids,
    #                                          ids_len_stop=ids_len_stop,
    #                                          output_dir=output_dir)
    # # Write out presence_absence summary
    # generate_summary_stats(output_dir=output_dir)
    #
    # # write pan genome reference fasta file
    # generate_pan_genome_reference(G,
    #                               output_dir=output_dir,
    #                               split_paralogs=False)
    #
    # # write out common structural differences in a matrix format
    # generate_common_struct_presence_absence(
    #     G,
    #     output_dir=output_dir,
    #     mems_to_isolates=mems_to_isolates,
    #     min_variant_support=args.min_edge_support_sv)
    #
    # # add helpful attributes and write out graph in GML format
    # for node in G.nodes():
    #     G.nodes[node]['size'] = len(G.nodes[node]['members'])
    #     G.nodes[node]['centroid'] = ";".join(G.nodes[node]['centroid'])
    #     G.nodes[node]['dna'] = ";".join(conv_list(G.nodes[node]['dna']))
    #     G.nodes[node]['protein'] = ";".join(conv_list(
    #         G.nodes[node]['protein']))
    #     G.nodes[node]['genomeIDs'] = ";".join(
    #         [str(m) for m in G.nodes[node]['members']])
    #     G.nodes[node]['geneIDs'] = ";".join(G.nodes[node]['seqIDs'])
    #     G.nodes[node]['degrees'] = G.degree[node]
    #     G.nodes[node]['members'] = list(G.nodes[node]['members'])
    #     G.nodes[node]['seqIDs'] = list(G.nodes[node]['seqIDs'])
    #
    # for edge in G.edges():
    #     G.edges[edge[0], edge[1]]['genomeIDs'] = ";".join(
    #         [str(m) for m in G.edges[edge[0], edge[1]]['members']])
    #     G.edges[edge[0],
    #             edge[1]]['members'] = list(G.edges[edge[0],
    #                                                edge[1]]['members'])
    #
    # nx.write_gml(G, output_dir + "final_graph.gml")
    #
    # # Write out core/pan-genome alignments
    # if args.aln == "pan":
    #     if verbose: print("generating pan genome MSAs...")
    #     generate_pan_genome_alignment(G, temp_dir, output_dir, args.n_cpu,
    #                                   args.alr, isolate_names)
    #     core_nodes = get_core_gene_nodes(G, args.core, len(args.input_files))
    #     concatenate_core_genome_alignments(core_nodes, output_dir)
    # elif args.aln == "core":
    #     if verbose: print("generating core genome MSAs...")
    #     generate_core_genome_alignment(G, temp_dir, output_dir,
    #                                    args.n_cpu, args.alr, isolate_names,
    #                                    args.core, len(args.input_files))
    #
    # # remove temporary directory
    # shutil.rmtree(temp_dir)

    return


if __name__ == '__main__':
    main()
