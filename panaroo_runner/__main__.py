import os, sys
import tempfile
from Bio import SeqIO
import networkx as nx
import argparse
import textwrap
import ast
from collections import defaultdict
import _pickle as cPickle

# panaroo scripts
from panaroo.isvalid import *
from .generate_network import *
from panaroo.cdhit import check_cdhit_version, run_cdhit
from .find_missing import *

# custom panaroo scripts
from .clean_network import *
from .generate_output import *
from .annotate import *
from .generate_alignments import check_aligner_install

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
                merge_para, aln, alr, core, min_edge_support_sv, all_seq_in_graph, ref_list, write_idx, kmer, repeat,
                remove_by_consensus, search_radius, refind_prop_match, annotate, evalue, annotation_db, hmm_db,
                call_variants, ignore_pseduogenes, truncation_threshold, save_objects, refind):

    # load shared memory items
    existing_shm = shared_memory.SharedMemory(name=shd_arr_tup.name)
    shd_arr = np.ndarray(shd_arr_tup.shape, dtype=shd_arr_tup.dtype, buffer=existing_shm.buf)

    # check if reference-guided alignment specified
    ref_aln = False
    if alr == "ref":
        ref_aln = True

    # Check cd-hit is installed
    check_cdhit_version()
    # Make sure aligner and rapidnj is installed if alignment requested
    if aln != None:
        check_aligner_install()
        check_rapidnj_install()

    # check snp-sites installed
    if call_variants:
        check_snpsites_install()

    if verbose:
        print("Generating initial network...")

    # generate network from clusters and adjacency information
    G, centroid_contexts, seqid_to_centroid = generate_network(shd_arr[0], high_scoring_ORFs, high_scoring_ORF_edges,
                                                               cluster_id_list, cluster_dict, overlap, all_seq_in_graph)


    # check if G is empty before proceeding
    if G.number_of_nodes() == 0:
        print("No genes detected in graph...")
        return

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

    if annotate is not None:
        if verbose:
            print("annotating gene families...")

        # create directory for annotation
        annotation_temp_dir = os.path.join(temp_dir, "annotation")
        if not os.path.exists(annotation_temp_dir):
            os.mkdir(annotation_temp_dir)
        # make sure trailing forward slash is present
        annotation_temp_dir = os.path.join(annotation_temp_dir, "")

        # generate annotations
        G, high_scoring_ORFs = iterative_annotation_search(G, high_scoring_ORFs, annotation_temp_dir, annotation_db,
                                                           hmm_db, evalue, annotate, n_cpu, pool)

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

    if refind:
        if verbose:
            print("refinding genes...")

        # find genes that Prokka has missed
        G, high_scoring_ORFs = find_missing(G,
                                            shd_arr_tup,
                                            high_scoring_ORFs,
                                            kmer=kmer,
                                            repeat=repeat,
                                            isolate_names=input_colours,
                                            remove_by_consensus=remove_by_consensus,
                                            search_radius=search_radius,
                                            prop_match=refind_prop_match,
                                            pairwise_id_thresh=identity_cutoff,
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
        print("writing Roary output...")

    # write out roary like gene_presence_absence.csv
    # get original annotation IDs, lengths and whether or
    # not an internal stop codon is present
    orig_ids = {}
    ids_len_stop = {}
    contig_annotation = defaultdict(list)
    for node in G.nodes():
        length_centroid = G.nodes[node]['longCentroidID'][0]
        for sid in G.nodes[node]['seqIDs']:
            orig_ids[sid] = sid
            mem = int(sid.split("_")[0])
            ORF_ID = int(sid.split("_")[-1])
            ORF_info = high_scoring_ORFs[mem][ORF_ID]
            ORF_len = ORF_info[2]
            # determine if gene is refound. If it is, then determine if premature stop codon present
            if (ORF_ID < 0):
                ids_len_stop[sid] = (ORF_len / 3, ORF_info[3])
            else:
                ids_len_stop[sid] = (ORF_len / 3, False)
            if annotate is not None and ref_list[mem]:
                # annotated genes
                if len(ORF_info) == 7 or ORF_ID < 0:
                    # add each sequence to its respective contig for each gff file.
                    annotation = ORF_info[-1]
                else:
                    annotation = ("prediction", "hypothetical protein", 0, "hypothetical protein")

                # annotate potential pseudogene if fits criteria of length, premature stop codon or length not multiple of 3
                if ORF_len < (length_centroid * truncation_threshold) or (ORF_ID < 0 and (ORF_info[3] is True
                                                                                          or ORF_len % 3 != 0)):
                    description = annotation[-1] + ", potential psuedogene"
                    annotation = annotation[0:3] + (description,)
                contig_annotation[mem].append((ORF_ID, annotation))

    # write output GFF
    if annotate is not None and any(ref_list):
        if verbose:
            print("writing GFF files...")
        generate_GFF(shd_arr[0], high_scoring_ORFs, input_colours, isolate_names, contig_annotation, output_dir,
                     overlap, write_idx, ref_list, n_cpu)

    # write roary output and summary stats file
    G = generate_roary_gene_presence_absence(G,
                                             mems_to_isolates=mems_to_isolates,
                                             orig_ids=orig_ids,
                                             ids_len_stop=ids_len_stop,
                                             output_dir=output_dir,
                                             threads=n_cpu)

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

    if verbose:
        print("writing gene fasta...")
    print_ORF_calls(high_scoring_ORFs, os.path.join(output_dir, "gene_calls.fasta"),
                    input_colours, overlap, shd_arr[0], truncation_threshold, G)

    # Write out core/pan-genome alignments
    # determine if reference-guided alignment being done
    if aln == "pan":
        if verbose: print("generating pan genome MSAs...")
        generate_pan_genome_alignment(G, temp_dir, output_dir, n_cpu, isolate_names, shd_arr_tup,
                                      high_scoring_ORFs, overlap, ref_aln, call_variants, verbose,
                                      ignore_pseduogenes, truncation_threshold, pool)
        core_nodes = get_core_gene_nodes(G, core, len(input_colours))
        core_gene_names = [G.nodes[x[0]]["name"] for x in core_nodes]
        concatenate_core_genome_alignments(core_gene_names, output_dir, isolate_names, n_cpu)
    elif aln == "core":
        if verbose: print("generating core genome MSAs...")
        generate_core_genome_alignment(G, temp_dir, output_dir,
                                       n_cpu, isolate_names, core, len(input_colours), shd_arr_tup,
                                       high_scoring_ORFs, overlap, ref_aln, call_variants, verbose,
                                       ignore_pseduogenes, truncation_threshold, pool)


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

    if verbose:
        print("writing graph file...")
    nx.write_gml(G, output_dir + "final_graph.gml")

    if save_objects:
        # create directory if it isn't present already
        objects_dir = output_dir + "ggc_data"
        if not os.path.exists(objects_dir):
            os.mkdir(objects_dir)

        # make sure trailing forward slash is present
        objects_dir = os.path.join(objects_dir, "")

        # serialise graph object and high scoring ORFs to future reading
        shd_arr[0].data_out(objects_dir + "ggc_graph.dat")
        with open(objects_dir + "high_scoring_orfs.dat", "wb") as o:
            cPickle.dump(high_scoring_ORFs, o)

        # create index of all high_scoring_ORFs node_IDs
        node_index = defaultdict(list)
        for colour, gene_dict in high_scoring_ORFs.items():
            for ORF_ID, ORF_info in gene_dict.items():
                entry_ID = str(colour) + "_" + str(ORF_ID)
                for node in ORF_info[0]:
                    node_index[abs(node)].append(entry_ID)

        with open(objects_dir + "node_index.dat", "wb") as o:
            cPickle.dump(node_index, o)


    return

#
# if __name__ == '__main__':
#     main()
