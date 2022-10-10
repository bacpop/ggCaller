import os, sys
import _pickle as cPickle

def search_graph(graph, graphfile, coloursfile, queryfile, objects_dir, output_dir, query_id, num_threads):
    # check if objects_dir present, if not exit
    if not os.path.exists(objects_dir):
        print("Please specify a ggc_data directory")
        sys.exit(1)

    objects_dir = os.path.join(objects_dir, "")

    # read in and parse sequences in queryfile
    query_vec = []
    with open(queryfile, "r") as f:
        for line in f:
            if line[0] != ">":
                query_vec.append(line.strip())

    # read in graph object and high_scoring ORFs and node_index
    graph.data_in(objects_dir + "ggc_graph.dat", graphfile, coloursfile, num_threads)

    # query the sequences in the graph
    print("Querying unitigs in graph...")
    input_colours, kmer, query_nodes = graph.search_graph(query_vec, query_id, num_threads)

    # parse isolate names
    isolate_names = [
        os.path.splitext(os.path.basename(x))[0] for x in input_colours
    ]

    with open(objects_dir + "high_scoring_orfs.dat", "rb") as input_file:
        high_scoring_ORFs = cPickle.load(input_file)

    with open(objects_dir + "node_index.dat", "rb") as input_file:
        node_index = cPickle.load(input_file)

    outfile = output_dir + "matched_queries.fasta"
    print("Matching overlapping ORFs...")

    ORF_dict = {}
    for i in range(len(query_nodes)):
        for node in query_nodes[i]:
            for ORF_ID in node_index[node]:
                if ORF_ID not in ORF_dict:
                    ORF_dict[ORF_ID] = []
                ORF_dict[ORF_ID].append(i)

    with open(outfile, "w") as f:
        for ORF, aligned in ORF_dict.items():
            split_ID = ORF.split("_")
            colour = int(split_ID[0])
            ORF_ID = int(split_ID[-1])
            fasta_ID = isolate_names[colour] + "_" + str(ORF_ID)
            ORF_info = high_scoring_ORFs[colour][ORF_ID]
            seq = graph.generate_sequence(ORF_info[0], ORF_info[1], kmer - 1)
            queries = ";".join([query_vec[i] for i in aligned])

            # add annotation
            f.write(
                    ">" + fasta_ID + " ggcID=" + ORF + " QUERY=" + queries + " annotation=" + ORF_info[-1] + "\n" + seq + "\n")

    return
