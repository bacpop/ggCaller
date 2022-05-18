# imports
import argparse
from ggCaller.graph_traversal import *
import ggCaller_cpp
import shutil
from models.__main__ import *
from ggCaller.shared_memory import *
from Bio import Seq
from panaroo_runner.set_default_args import *
from panaroo_runner.__main__ import run_panaroo
from panaroo_runner.generate_output import print_ORF_calls
from panaroo_runner.annotate import check_diamond_install, check_HMMER_install, generate_HMMER_index, \
    generate_diamond_index
from collections import defaultdict
import ast
import tempfile
import json

def get_options():
    description = 'Generates ORFs from a Bifrost graph.'
    parser = argparse.ArgumentParser(description=description,
                                     prog='ggcaller')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--graph',
                    default=None,
                    help='Bifrost GFA file generated by Bifrost build. ')
    IO.add_argument('--colours',
                    default=None,
                    help='Bifrost colours file generated by Bifrost build.  ')
    IO.add_argument('--not-ref',
                    action="store_false",
                    default=True,
                    help='If using existing graph, was not graph built exclusively with assembled genomes.  '
                         '[Default = False] ')
    IO.add_argument('--refs',
                    default=None,
                    help='List of reference genomes (one file path per line). ')
    IO.add_argument('--reads',
                    default=None,
                    help='List of read files (one file path per line). ')
    IO.add_argument('--query',
                    default=None,
                    help='List of unitig sequences to query (either FASTA or one sequence per line) ')
    IO.add_argument('--codons',
                    default=None,
                    help='JSON file containing start and stop codon sequences. ')
    IO.add_argument('--kmer',
                    type=int,
                    default=31,
                    help='K-mer size used in Bifrost build (bp). '
                         '[Default = 31] ')
    IO.add_argument('--save',
                    action="store_true",
                    default=False,
                    help='Save graph objects for sequence querying. '
                         '[Default = False] ')
    IO.add_argument('--data',
                    default=None,
                    help='Directory containing data from previous ggCaller run generated via "--save" ')
    IO.add_argument(
        "--all-seq-in-graph",
        dest="all_seq_in_graph",
        help=("Retains all DNA sequence for each gene cluster in the Panaroo graph " +
              "output. Off by default as it uses a large amount of space."),
        action='store_true',
        default=False)
    IO.add_argument('--out',
                    default='ggCaller_output',
                    help='Output directory ')
    Settings = parser.add_argument_group('ggCaller traversal and gene-calling cut-off settings')
    Settings.add_argument('--max-path-length',
                          type=int,
                          default=20000,
                          help='Maximum path length during ORF finding (bp). '
                               '[Default = 20000] ')
    Settings.add_argument('--min-orf-length',
                          type=int,
                          default=90,
                          help='Minimum ORF length to return (bp). '
                               '[Default = 90] ')
    Settings.add_argument('--max-ORF-overlap',
                          type=int,
                          default=60,
                          help='Maximum overlap allowed between overlapping ORFs. '
                               '[Default = 60] ')
    Settings.add_argument('--min-path-score',
                          type=int,
                          default=100,
                          help='Minimum total Balrog score for a path of ORFs to be returned. '
                               '[Default = 100] ')
    Settings.add_argument('--min-orf-score',
                          type=int,
                          default=100,
                          help='Minimum individual Balrog score for an ORF to be returned. '
                               '[Default = 100] ')
    Settings.add_argument('--max-orf-orf-distance',
                          type=int,
                          default=10000,
                          help='Maximum distance for graph traversal during ORF connection (bp). '
                               '[Default = 10000] ')
    Settings.add_argument('--query-id',
                          type=float,
                          default=0.8,
                          help='Ratio of query-kmers to required to match in graph. '
                               '[Default = 0.8] ')
    Algorithm = parser.add_argument_group('Settings to avoid/include algorithms')
    Algorithm.add_argument('--no-filter',
                           action="store_true",
                           default=False,
                           help='Do not filter ORF calls using Balrog. Will return all ORF calls. '
                                '[Default = False] ')
    Algorithm.add_argument('--no-write-idx',
                           action="store_false",
                           default=True,
                           help='Do not write FMIndexes to file. '
                                '[Default = False] ')
    Algorithm.add_argument('--no-write-graph',
                           action="store_false",
                           default=True,
                           help='Do not write Bifrost GFA and colours to file. '
                                '[Default = False] ')
    Algorithm.add_argument('--repeat',
                           action="store_true",
                           default=False,
                           help='Enable traversal of nodes multiple times. '
                                '[Default = False] ')
    Algorithm.add_argument('--no-clustering',
                           action="store_true",
                           default=False,
                           help='Do not cluster ORFs. '
                                '[Default = False] ')
    Algorithm.add_argument('--no-refind',
                           dest="refind",
                           action="store_false",
                           default=True,
                           help='Do not refind uncalled genes '
                                '[Default = False] ')
    Clustering = parser.add_argument_group('Gene clustering options.')
    Clustering.add_argument('--identity-cutoff',
                            type=float,
                            default=0.98,
                            help='Minimum identity at amino acid level between two ORFs for clustering. '
                                 '[Default = 0.98] ')
    Clustering.add_argument('--len-diff-cutoff',
                            type=float,
                            default=0.98,
                            help='Minimum ratio of length between two ORFs for clustering.  '
                                 '[Default = 0.98] ')
    Clustering.add_argument(
        "--family-threshold",
        dest="family_threshold",
        help="protein family sequence identity threshold (default=0.7)",
        type=float)
    Clustering.add_argument("--merge-paralogs",
                            dest="merge_paralogs",
                            help="don't split paralogs",
                            action='store_true',
                            default=False)
    Panaroo_mode_opts = parser.add_argument_group('Panaroo run mode options')
    Panaroo_mode_opts.add_argument(
        "--clean-mode",
        dest="mode",
        help=
        ('''R|The stringency mode at which to run panaroo. Must be one of 'strict',\
    'moderate' or 'sensitive'. Each of these modes can be fine tuned using the\
     additional parameters in the 'Graph correction' section.

    strict: 
    Requires fairly strong evidence (present in  at least 5%% of genomes)\
     to keep likely contaminant genes. Will remove genes that are refound more often than\
     they were called originally.

    moderate: 
    Requires moderate evidence (present in  at least 1%% of genomes)\
     to keep likely contaminant genes. Keeps genes that are refound more often than\
     they were called originally.

    sensitive: 
    Does not delete any genes and only performes merge and refinding\
     operations. Useful if rare plasmids are of interest as these are often hard to\
     disguish from contamination. Results will likely include  higher number of\
     spurious annotations.'''),
        choices=['strict', 'moderate', 'sensitive'],
        required=False)

    Panaroo_annotation = parser.add_argument_group('Panaroo gene cluster annotation options')
    Panaroo_annotation.add_argument("--annotation",
                                    dest="annotate",
                                    help="Annotate genes using diamond (fast) or diamond and HMMscan (sensitive)."
                                         "If not specified, no annotation done",
                                    choices=["fast", "sensitive"],
                                    default=None)
    Panaroo_annotation.add_argument("--diamonddb",
                                    dest="annotation_db",
                                    help="Diamond database. Defaults are 'Bacteria' or 'Viruses'. Can also "
                                         "specify path to fasta file for custom database generation",
                                    default="Bacteria")
    Panaroo_annotation.add_argument("--hmmdb",
                                    dest="hmm_db",
                                    help="HMMER hmm profile file. Default is Uniprot HAMAP. Can also"
                                         "specify path to pre-built hmm profile file generated using hmmbuild",
                                    type=str,
                                    default="default")
    Panaroo_annotation.add_argument("--evalue",
                                    dest="evalue",
                                    help="Maximum e-value to return for DIAMOND and HMMER searches during annotation",
                                    default=0.001,
                                    type=float)
    Panaroo_annotation.add_argument("--truncation-threshold",
                                    dest="truncation_threshold",
                                    help="Sequences in a gene family cluster below this proportion of the length of the"
                                         "centroid will be annotated as 'potential pseudogene'",
                                    default=0.8,
                                    type=float)

    Panaroo_refind = parser.add_argument_group('Panaroo gene-refinding options')
    Panaroo_refind.add_argument(
        "--search-radius",
        dest="search_radius",
        help=("the distance in nucleotides surronding the " +
              "neighbour of an accessory gene in which to search for it"),
        default=5000,
        type=int)
    Panaroo_refind.add_argument(
        "--refind-prop-match",
        dest="refind_prop_match",
        help=("the proportion of an accessory gene that must " +
              "be found in order to consider it a match"),
        default=0.2,
        type=float)
    Panaroo_graph = parser.add_argument_group('Panaroo graph correction stringency options')
    Panaroo_graph.add_argument(
        "--remove-invalid-genes",
        dest="filter_invalid",
        action='store_true',
        default=False,
        help=("removes annotations that do not conform to the expected Prokka" +
              " format such as those including premature stop codons."))
    Panaroo_graph.add_argument(
        "--min-trailing-support",
        dest="min_trailing_support",
        help=("minimum cluster size to keep a gene called at the " +
              "end of a contig"),
        type=int)
    Panaroo_graph.add_argument(
        "--trailing-recursive",
        dest="trailing_recursive",
        help=("number of times to perform recursive trimming of low support " +
              "nodes near the end of contigs"),
        type=int)
    Panaroo_graph.add_argument(
        "--edge-support-threshold",
        dest="edge_support_threshold",
        help=("minimum support required to keep an edge that has been flagged" +
              " as a possible mis-assembly"),
        type=float)
    Panaroo_graph.add_argument(
        "--length-outlier-support-proportion",
        dest="length_outlier_support_proportion",
        help=
        ("proportion of genomes supporting a gene with a length more " +
         "than 1.5x outside the interquatile range for genes in the same cluster"
         +
         " (default=0.01). Genes failing this test will be re-annotated at the "
         + "shorter length"),
        type=float,
        default=0.1)
    Panaroo_graph.add_argument(
        "--remove-by-consensus",
        dest="remove_by_consensus",
        type=ast.literal_eval,
        choices=[True, False],
        help=
        ("if a gene is called in the same region with similar sequence a minority "
         + "of the time, remove it. One of 'True' or 'False'"),
        default=None)
    Panaroo_graph.add_argument(
        "--high-var-flag",
        dest="cycle_threshold_min",
        help=(
                "minimum number of nested cycles to call a highly variable gene " +
                "region (default = 5)."),
        type=int,
        default=5)
    Panaroo_graph.add_argument(
        "--min-edge-support-sv",
        dest="min_edge_support_sv",
        help=("minimum edge support required to call structural variants" +
              " in the presence/absence sv file"),
        type=int)
    Panaroo_graph.add_argument(
        "--no-clean-edges",
        dest="clean_edges",
        help=("Turn off edge filtering in the final output graph."),
        action='store_false',
        default=True)

    Panaroo_aln = parser.add_argument_group('Gene alignment options')
    Panaroo_aln.add_argument(
        "--alignment",
        dest="aln",
        help=("Output alignments of core genes or all genes. Options are" +
              " 'core' and 'pan'. Default: 'None'"),
        type=str,
        choices=['core', 'pan'],
        default=None)
    Panaroo_aln.add_argument(
        "--aligner",
        dest="alr",
        help=
        "Specify an aligner. Options:'ref' for reference-guided MSA and 'def' for default standard MSA",
        type=str,
        choices=['def', 'ref'],
        default="def")
    Panaroo_aln.add_argument("--core-threshold",
                             dest="core",
                             help="Core-genome sample threshold (default=0.95)",
                             type=float,
                             default=0.95)
    Panaroo_aln.add_argument("--no-variants",
                             dest="call_variants",
                             help="Do not call variants using SNP-sites after alignment.",
                             action='store_false',
                             default=True)
    Panaroo_aln.add_argument("--ignore-pseduogenes",
                             dest="ignore_pseduogenes",
                             help="Ignore ORFs annotated as 'potential pseudogenes' in alignment",
                             action='store_true',
                             default=False)

    # Other options
    Misc = parser.add_argument_group('Misc. options')
    Misc.add_argument("--quiet",
                      dest="verbose",
                      help="suppress additional output",
                      action='store_false',
                      default=True)
    Misc.add_argument('--threads',
                      type=int,
                      default=1,
                      help='Number of threads to use. '
                           '[Default = 1] ')

    return parser.parse_args()

def main():
    # parse command line arguments for ggCaller
    options = get_options()

    # determine if references/assemblies present
    ref_set = set()
    if options.refs is not None:
        with open(options.refs, "r") as f:
            for line in f.readlines():
                line = line.strip("\n")
                ref_set.add(line)

    if (options.refs is not None and options.reads is None) or (
            options.graph is not None and options.colours is not None
            and options.not_ref):
        is_ref = True
    else:
        is_ref = False

    # define start/stop codons
    if options.codons is not None:
        with open(options.codons, "r") as json_file:
            try:
                data = json.load(json_file)
                start_codons = data["codons"]["start"]
                stop_codons_for = data["codons"]["stop"]
                stop_codons_rev = [str((Seq(i)).reverse_complement()) for i in stop_codons_for]
            except:
                print("Please specify codons in the format shown in codons.json.")
                sys.exit(1)
    else:
        start_codons = ["ATG", "GTG", "TTG"]
        stop_codons_for = ["TAA", "TGA", "TAG"]
        stop_codons_rev = ["TTA", "TCA", "CTA"]

    # initialise graph
    graph = ggCaller_cpp.create_graph()

    # create directory if it isn't present already
    if not os.path.exists(options.out):
        os.mkdir(options.out)

    # make sure trailing forward slash is present
    output_dir = os.path.join(options.out, "")

    # if build graph specified, build graph and then call ORFs
    if (options.graph is not None) and (options.colours is not None) and (options.query is None):
        graph_tuple = graph.read(options.graph, options.colours, stop_codons_for, stop_codons_rev,
                                 options.threads, is_ref, ref_set)
    # query unitigs in previous saved ggc graph
    elif (options.graph is not None) and (options.colours is not None) and (options.refs is None) and \
            (options.query is not None):
        if options.data is None:
            print("Please specify a ggc_data directory from a previous ggCaller run.")
            sys.exit(1)
        search_graph(graph, options.graph, options.colours, options.query, options.data, output_dir, options.query_id,
                     options.threads)
        print("Finished.")
        sys.exit(0)
    # if refs file specified for building
    elif (options.graph is None) and (options.colours is None) and (options.refs is not None) and (
            options.reads is None) and (
            options.query is None):
        graph_tuple = graph.build(options.refs, options.kmer, stop_codons_for, stop_codons_rev,
                                  options.threads, True, options.no_write_graph, "NA", ref_set)
    # if reads file specified for building
    elif (options.graph is None) and (options.colours is None) and (options.refs is None) and (
            options.reads is not None) and (options.query is None):
        graph_tuple = graph.build(options.reads, options.kmer, stop_codons_for, stop_codons_rev,
                                  options.threads, False, options.no_write_graph, "NA", ref_set)
    # if both reads and refs file specified for building
    elif (options.graph is None) and (options.colours is None) and (options.refs is not None) and (
            options.reads is not None) and (options.query is None):
        graph_tuple = graph.build(options.refs, options.kmer, stop_codons_for, stop_codons_rev,
                                  options.threads, False, options.no_write_graph, options.reads, ref_set)
    else:
        print("Error: incorrect number of input files specified. Please only specify the below combinations:\n"
              "- Bifrost GFA and Bifrost colours file (with/without list of reference files)\n"
              "- Bifrost GFA, Bifrost colours file and list of query sequences\n"
              "- List of reference files\n"
              "- List of read files\n"
              "- A list of reference files and a list of read files.")
        sys.exit(1)

    # unpack ORF pair into overlap dictionary and list for gene scoring
    input_colours, nb_colours, overlap, ref_list = graph_tuple

    # download balrog and annotation files
    db_dir = download_db()

    # set rest of panaroo arguments
    options = set_default_args(options, nb_colours)
    annotation_db = options.annotation_db
    hmm_db = options.hmm_db

    if options.annotate is not None:
        # check diamond and HMMER are installed correctly
        check_diamond_install()
        check_HMMER_install()

        # unpack annotation database
        if annotation_db == "Bacteria" or annotation_db == "Viruses":
            db_id = annotation_db
            diamond_dir = os.path.join(db_dir, "diamond")
            annotation_db = os.path.join(diamond_dir, annotation_db)

            if not os.path.exists(annotation_db):
                print("Unzipping protein annotation file...")
                tar = tarfile.open(annotation_db + ".tar.gz", mode="r:gz")
                tar.extractall(diamond_dir)
                tar.close()

            annotation_db = os.path.join(annotation_db, db_id + ".dmnd")

        # if custom annotation database specified, then create diamond db if not present already
        else:
            annotation_db = os.path.abspath(annotation_db)
            if ".dmnd" not in annotation_db:
                print("Generating diamond index...")
                annotation_db = generate_diamond_index(annotation_db)

        # set-up hmm_db
        if hmm_db == "default":
            hmm_dir = os.path.join(db_dir, "hmm")
            hmm_db = os.path.join(hmm_dir, "HAMAP.hmm")
        else:
            hmm_db = os.path.abspath(hmm_db)

        if not os.path.exists(hmm_db + ".h3f") and options.annotate == "sensitive":
            print("Generating HMMER index...")
            generate_HMMER_index(hmm_db)

    # Create temporary directory
    temp_dir = os.path.join(tempfile.mkdtemp(dir=output_dir), "")

    # Create temp_file for cluster_map
    cluster_file = os.path.join(temp_dir, "cluster_map.dat")

    # load models models if required
    if not options.no_filter:
        print("Loading gene models...")
        ORF_model_file, TIS_model_file = load_balrog_models()

    else:
        ORF_model_file, TIS_model_file = "NA", "NA"

    gene_tuple = graph.findGenes(options.repeat, overlap, options.max_path_length,
                                 options.no_filter, stop_codons_for, start_codons, options.min_orf_length,
                                 options.max_ORF_overlap, input_colours, ORF_model_file,
                                 TIS_model_file, options.min_orf_score, options.min_path_score,
                                 options.max_orf_orf_distance, not options.no_clustering,
                                 options.identity_cutoff, options.len_diff_cutoff, options.threads, cluster_file)

    high_scoring_ORFs, high_scoring_ORF_edges = gene_tuple

    # generate ORF clusters
    if not options.no_clustering:
        with SharedMemoryManager() as smm:
            # generate shared numpy arrays
            total_arr = np.array([graph])
            array_shd, array_shd_tup = generate_shared_mem_array(total_arr, smm)
            with Pool(processes=options.threads) as pool:
                run_panaroo(pool, array_shd_tup, high_scoring_ORFs, high_scoring_ORF_edges,
                            cluster_file, overlap, input_colours, output_dir, temp_dir, options.verbose,
                            options.threads, options.length_outlier_support_proportion, options.identity_cutoff,
                            options.family_threshold, options.min_trailing_support, options.trailing_recursive,
                            options.clean_edges, options.edge_support_threshold, options.merge_paralogs, options.aln,
                            options.alr, options.core, options.min_edge_support_sv, options.all_seq_in_graph, ref_list,
                            options.no_write_idx, overlap + 1, options.repeat, options.remove_by_consensus,
                            options.search_radius, options.refind_prop_match, options.annotate, options.evalue,
                            annotation_db, hmm_db, options.call_variants, options.ignore_pseduogenes,
                            options.truncation_threshold, options.save, options.refind)

    else:
        print_ORF_calls(high_scoring_ORFs, os.path.join(output_dir, "gene_calls.fasta"),
                        input_colours, overlap, graph)

        if options.save:
            # create directory if it isn't present already
            objects_dir = output_dir + "ggc_data"
            if not os.path.exists(objects_dir):
                os.mkdir(objects_dir)

            # make sure trailing forward slash is present
            objects_dir = os.path.join(objects_dir, "")

            # serialise graph object and high scoring ORFs to future reading
            graph[0].data_out(objects_dir + "ggc_graph.dat")
            with open(objects_dir + "high_scoring_orfs.dat", "wb") as o:
                cPickle.dump(high_scoring_ORFs, o)

            # create index of all high_scoring_ORFs node_IDs
            node_index = defaultdict(list)
            for colour, gene_dict in high_scoring_ORFs.items():
                for ORF_ID, ORF_info in gene_dict.items():
                    entry_ID = str(colour) + "_" + str(ORF_ID)
                    for node in ORF_info[0]:
                        node_index[node].append(entry_ID)

            with open(objects_dir + "node_index.dat", "wb") as o:
                cPickle.dump(node_index, o)

    # remove temporary directory
    shutil.rmtree(temp_dir)

    print("Finished.")

    sys.exit(0)


if __name__ == '__main__':
    main()
