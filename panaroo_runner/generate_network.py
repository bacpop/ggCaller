from collections import defaultdict
import networkx as nx
from intbitset import intbitset
import ggCaller_cpp


def generate_network(DBG, overlap, ORF_file_paths, Edge_file_paths, cluster_file):
    # read in cluster_dict
    cluster_dict = ggCaller_cpp.read_cluster_file(cluster_file)

    # associate sequences with their clusters
    seq_to_cluster = {}
    seqid_to_centroid = {}
    cluster_centroids = {}
    cluster_members = defaultdict(list)
    cluster_centroid_data = {}
    cluster_id = 0

    # iterate over all ORFs in ORF_file_paths to determine centroids and add clusters
    for colour_ID, file_path in ORF_file_paths.items():
        ORF_map = ggCaller_cpp.read_ORF_file(file_path)

        for ORF_ID, ORFNodeVector in ORF_map.items():
            # identify centroids
            ORF_ID_str = str(colour_ID) + "_" + str(ORF_ID)

            if ORF_ID_str in cluster_dict:
                # access ORF information for centroid from high_scoring_ORFs
                pan_centroid_ID = str(colour_ID) + "_0_" + str(ORF_ID)

                seq = DBG.generate_sequence(ORFNodeVector[0], ORFNodeVector[1], overlap)
                current_hash = DBG.rb_hash(seq)

                # add information to cluster_centroid_data
                cluster_centroid_data[cluster_id] = {
                    'ORF_info': ORFNodeVector,
                    'hash': current_hash,
                }

                # append centroid to cluster
                cluster_centroids[cluster_id] = pan_centroid_ID

                # for each ORF, add ORF_id and cluster_id to respective dictionaries
                for ORF_ID_pair in cluster_dict[ORF_ID_str]:
                    # generate a panaroo sequence ID for current ORF
                    genome_id = ORF_ID_pair[0]
                    local_id = ORF_ID_pair[1]

                    pan_ORF_id = str(genome_id) + "_0_" + str(local_id)

                    # index sequences to clusters and the number of edges they have
                    seq_to_cluster[pan_ORF_id] = [cluster_id, 0]
                    cluster_members[cluster_id].append(pan_ORF_id)

                cluster_id += 1

    # clear cluster_dict
    cluster_dict.clear()

    # determine paralogs if required
    paralogs = set()
    for clust in cluster_members:
        # determine if paralogs present by checking colours for each entry in cluster_id_list
        genomes = [int(s.split("_")[0]) for s in cluster_members[clust]]
        if len(genomes) != len(set(genomes)):
            paralogs.add(clust)

    # build graph using adjacency information and optionally split paralogs
    G = nx.Graph()
    centroid_context = defaultdict(list)
    n_nodes = len(cluster_members)

    # hold which nodes hold genes from different colours
    colour_to_nodes = {}

    # iterating over each cluster, adding edges between clusters that each have a connection
    for current_cluster, ORF_members in cluster_members.items():
        # check if cluster is currently in graph
        add_cluster = True if not G.has_node(current_cluster) else False

        # check if current cluster contains paralogs
        has_paralogs = True if current_cluster in paralogs else False

        ORFNodeVector = cluster_centroid_data[current_cluster]['ORF_info']

        for ORF_id in ORF_members:
            # parse genome id and local ORF id
            parsed_id = ORF_id.split("_")
            genome_id = int(parsed_id[0])
            local_id = int(parsed_id[-1])

            # initialise cluster to add
            cluster_to_add = current_cluster

            if add_cluster:
                if has_paralogs:
                    seq = DBG.generate_sequence(ORFNodeVector[0], ORFNodeVector[1], overlap)
                    seq_hash = DBG.rb_hash(seq)
                    # create a new paralog
                    n_nodes += 1
                    cluster_to_add = n_nodes
                    centroid_context[
                        cluster_centroids[current_cluster]].append(
                        (cluster_to_add, genome_id, seq_hash))
                    # overwrite seq_to_cluster for current ORF
                    seq_to_cluster[ORF_id][0] = cluster_to_add
                G.add_node(
                    cluster_to_add,
                    CID=cluster_to_add,
                    size=1,
                    centroid=[cluster_centroids[current_cluster]],
                    maxLenId=0,
                    members=intbitset([genome_id]),
                    seqIDs=set([ORF_id]),
                    hasEnd=False,
                    ORF_info=[(ORFNodeVector[0],
                               ORFNodeVector[1])],
                    hash=[cluster_centroid_data[current_cluster]['hash']],
                    annotation='',
                    bitscore=0,
                    description='',
                    lengths=[
                        ORFNodeVector[2]
                    ],
                    paralog=has_paralogs,
                    mergedDNA=False)
                # if has_paralogs == true, then need to add new node for each sequence
                # otherwise can stop adding clusters
                if not has_paralogs:
                    add_cluster = False

            else:
                # check if ORF_id already added to the cluster
                if ORF_id not in G.nodes[current_cluster]['seqIDs']:
                    G.nodes[current_cluster]['size'] += 1
                    G.nodes[current_cluster]['members'].add(genome_id)
                    G.nodes[current_cluster]['seqIDs'].add(ORF_id)
            
            # add to colour_to_nodes
            if genome_id not in colour_to_nodes:
                colour_to_nodes[genome_id] = set()
            colour_to_nodes[genome_id].add(current_cluster)

    # iterative over each colour, adding edges between nodes with genes of that colour
    for genome_id, node_set in colour_to_nodes.items():
        ORF_edges = ggCaller_cpp.read_edge_file(Edge_file_paths[genome_id])
        for node in node_set:                        
            for index, ORF_id in enumerate(G.nodes[node]['seqIDs']):
                parsed_id = ORF_id.split("_")
                genome_id = int(parsed_id[0])
                local_id = int(parsed_id[-1])

                # if edge present between ORFs add, otherwise pass
                if local_id in ORF_edges:
                    edge_set = ORF_edges[local_id]
                else:
                    continue

                for neighbour in edge_set:
                    pan_neigbour_id = str(genome_id) + "_0_" + str(neighbour)

                    # ensure neighbour is high scoring ORF, otherwise ignore
                    if pan_neigbour_id in seq_to_cluster:
                        neighbour_cluster = seq_to_cluster[pan_neigbour_id][0]
                        # update edge count for ORF
                        seq_to_cluster[pan_neigbour_id][1] += 1
                        seq_to_cluster[ORF_id][1] += 1
                    else:
                        continue

                    # add edge between current ORF and neighbour
                    if G.has_edge(node, neighbour_cluster):
                        G[node][neighbour_cluster]['size'] += 1
                        G[node][neighbour_cluster]['members'].add(genome_id)
                    else:
                        G.add_edge(node,
                                neighbour_cluster,
                                size=1,
                                members=intbitset([genome_id]))

    # iterate over seq_to_cluster, identifying nodes that have ends
    for ORF_ID, cluster_list in seq_to_cluster.items():
            if cluster_list[1] < 2:
                G.nodes[cluster_list[0]]['hasEnd'] = True

    return G, centroid_context
