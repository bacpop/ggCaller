import os, sys
import _pickle as cPickle
from Bio import SeqIO
import networkx as nx
from collections import defaultdict
import ggCaller_cpp

# TODO update this to work with ggCaller ORF map outputs and panaroo graph which contains annotation
def search_graph(graph, graphfile, coloursfile, queryfile, objects_dir, output_dir, query_id, num_threads):
    # check if objects_dir present, if not exit
    ggc_data_dir = os.path.join(objects_dir, "ggc_data")
    objects_dir = os.path.join(objects_dir, "")

    if not os.path.exists(ggc_data_dir):
        print("Please specify a directory from a ggCaller run with '--save'")
        sys.exit(1)

    # read in and parse sequences in queryfile
    id_vec = []
    query_vec = []
    fasta_file = False

    with open(queryfile, "r") as f:
        first_line = f.readline()
        if ">" in first_line:
            fasta_file = True

    if fasta_file:
        fasta_sequences = SeqIO.parse(open(queryfile), 'fasta')
        for entry in fasta_sequences:
            query_vec.append(str(entry.seq))
            id_vec.append(str(entry.id))
    else:
        with open(queryfile, "r") as f:
            for line in f:
                if line[0] != ">":
                    query_vec.append(line.strip())

    # read in graph object and high_scoring ORFs and node_index
    graph.data_in(os.path.join(objects_dir, "ORF_dir", "kmer_array.dat"), graphfile, coloursfile, num_threads)

    # query the sequences in the graph
    print("Querying unitigs in graph...")
    input_colours, kmer, query_nodes = graph.search_graph(query_vec, query_id, num_threads)

    # parse isolate names
    isolate_names = [
        os.path.splitext(os.path.basename(x))[0] for x in input_colours
    ]

    # load in gene to DBG node mappings
    with open(os.path.join(ggc_data_dir, "node_index.dat"), "rb") as input_file:
        node_index = cPickle.load(input_file)

    # load in paths to ORF data
    with open(os.path.join(objects_dir, "ORF_dir", "ORF_file_paths.dat"), "rb") as input_file:
        ORF_file_paths = cPickle.load(input_file)

    # try to load panaroo graph, not present if panaroo not run
    G = None
    ids_len_stop = None
    try:
        G = nx.read_gml('t.gml')
        # load in gene to panaroo node mappings
        with open(ggc_data_dir + "ORF_to_node_map.dat", "rb") as input_file:
            ids_len_stop = cPickle.load(input_file)
    except:
        pass

    outfile = output_dir + "matched_queries.fasta"
    print("Matching overlapping ORFs...")


    ORF_dict = defaultdict(lambda: defaultdict(list))
    for i in range(len(query_nodes)):
        for node in query_nodes[i]:
            for ORF in node_index[node]:
                split_ID = ORF.split("_")
                colour = int(split_ID[0])
                ORF_ID = int(split_ID[-1])
                ORF_dict[colour][ORF_ID].append(i)

    with open(outfile, "w") as f:
        for colour, ORF_element in ORF_dict.items():
            # load ORF info
            ORF_map = ggCaller_cpp.read_ORF_file(ORF_file_paths[colour])
            fasta_ID = isolate_names[colour] + "_" + str(ORF_ID)
            
            # iterate over each ORF that has match
            for ORF_ID, aligned in ORF_element.items():
                ORF_info = ORF_map[ORF_ID]
                seq = graph.generate_sequence(ORF_info[0], ORF_info[1], kmer - 1)
                if fasta_file:
                    id_set = set([id_vec[i] for i in aligned])
                    queries = ";".join([i for i in id_set])
                else:
                    queries = ";".join([query_vec[i] for i in aligned])
                
                # get annotation for node
                annotation = "NA"
                if G is not None:
                    delim = "_0_" if ORF_ID >= 0 else "_refound_"
                    pan_ORF_id = str(colour) + delim + str(ORF_ID)
                    # if pan_ORF_id not in ORF_map means has been removed so remove it from ORF_map
                    node = ids_len_stop[pan_ORF_id][2]
                    annotation = G.nodes[node]["description"]

                # add annotation
                f.write(
                        ">" + fasta_ID + " ggcID=" + ORF + " QUERY=" + queries + " annotation=" + annotation + "\n" + seq + "\n")

    return
