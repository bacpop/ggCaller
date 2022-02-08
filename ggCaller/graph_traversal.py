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
    input_colours, kmer, query_nodes = graph.search_graph(graphfile, coloursfile, query_vec, query_id, num_threads)

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
    with open(outfile, "w") as f:
        for i in range(len(query_nodes)):
            query_set = set()
            for node in query_nodes[i]:
                query_set.update(node_index[node])
                for ORF in query_set:
                    split_ID = ORF.split("_")
                    colour = int(split_ID[0])
                    ORF_ID = int(split_ID[1])
                    fasta_ID = isolate_names[colour] + "_" + str(ORF_ID).zfill(5)
                    ORF_info = high_scoring_ORFs[colour][ORF_ID]
                    seq = graph.generate_sequence(ORF_info[0], ORF_info[1], kmer - 1)
                    # add annotation if available
                    if len(ORF_info) == 8 or ORF_ID < 0:
                        ORF_annotation = ORF_info[-1]
                        f.write(
                            ">" + fasta_ID + " " + ORF_annotation[-1] + " QUERY=" + query_vec[i] + "\n" + seq + "\n")
                    else:
                        f.write(">" + fasta_ID + " QUERY=" + query_vec[i] + "\n" + seq + "\n")

    return
