from collections import defaultdict
import networkx as nx
from intbitset import intbitset
import ggCaller_cpp


def generate_network(DBG, overlap, ORF_file_paths, Edge_file_paths, cluster_file):
    # read in cluster_dict
    # TODO save pair here that holds ORFs removed for low scores after centroid scored
    cluster_dict, ORFs_to_remove = ggCaller_cpp.read_cluster_file(cluster_file)

    # associate sequences with their clusters
    seq_to_cluster = {}
    ORF_length_map = {}
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

            # only hold lengths of genes that are not in a cluster
            if ORF_ID_str not in seq_to_cluster:
                seq = DBG.generate_sequence(ORFNodeVector[0], ORFNodeVector[1], overlap)
                current_hash = DBG.rb_hash(seq)
                ORF_length_map[ORF_ID_str] = (ORFNodeVector[2], current_hash)

            if ORF_ID_str in cluster_dict:
                # access ORF information for centroid from high_scoring_ORFs
                pan_centroid_ID = str(colour_ID) + "_0_" + str(ORF_ID)

                # make sure ORF wasn't removed after centroid scored
                if pan_centroid_ID in ORFs_to_remove:
                    continue

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

                    # make sure ORF wasn't removed after centroid scored
                    if pan_ORF_id in ORFs_to_remove:
                        continue

                    # only hold lengths of genes that are not in a cluster
                    if ORF_ID_str in ORF_length_map:
                        del ORF_length_map[ORF_ID_str]

                    # index sequences to clusters and the number of edges they have
                    seq_to_cluster[pan_ORF_id] = [cluster_id, 0]
                    cluster_members[cluster_id].append(pan_ORF_id)

                cluster_id += 1
                del cluster_dict[ORF_ID_str]

    # determine if any clusters have not been added as centroid was removed, add a new centroid
    for centroid, ORF_IDs in cluster_dict.items():
        current_centroid = ""
        current_length = 0
        current_hash = 0
        cluster_list = []
        
        # iterate through, if in ORF_length_map then is real ORF
        for ORF_ID_pair in ORF_IDs:
            genome_id = ORF_ID_pair[0]
            local_id = ORF_ID_pair[1]

            pan_ORF_id = str(genome_id) + "_0_" + str(local_id)

            if pan_ORF_id in ORF_length_map:
                new_centroid = False
                length, hash = ORF_length_map[pan_ORF_id]
                
                # add to seq_to_cluster for new ORFs
                seq_to_cluster[pan_ORF_id] = [cluster_id, 0]

                # assign centroid first on length, then hash, then genome index
                if length > current_length:
                    new_centroid = True
                elif length == current_length:
                    if hash < current_hash:
                        new_centroid = True
                    elif hash == current_hash:
                        centroid_genome_ID = int(current_centroid.split("_")[0])
                        if genome_id < centroid_genome_ID:
                            new_centroid = True
                
                # add to end if not centroid
                if new_centroid == True:
                    if current_centroid != "":
                        cluster_list.append(current_centroid)
                    current_centroid = pan_ORF_id
                    current_length = length
                    current_hash = hash
                else:
                    cluster_list.append(pan_ORF_id)

        # add only if other genes found in cluster
        if current_centroid != "":
            # ensure centroid at start of list
            cluster_list = [current_centroid] + cluster_list

            # add cluster member
            cluster_members[cluster_id] = cluster_list
            cluster_centroids[cluster_id] = current_centroid
            
            # index sequences to clusters and the number of edges they have
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
                
                # add to colour_to_nodes
                if genome_id not in colour_to_nodes:
                    colour_to_nodes[genome_id] = set()
                colour_to_nodes[genome_id].add(cluster_to_add)
            else:
                # check if ORF_id already added to the cluster
                if ORF_id not in G.nodes[current_cluster]['seqIDs']:
                    G.nodes[current_cluster]['size'] += 1
                    G.nodes[current_cluster]['members'].add(genome_id)
                    G.nodes[current_cluster]['seqIDs'].add(ORF_id)
            

    # iterative over each colour, adding edges between nodes with genes of that colour
    for colour_id, node_set in colour_to_nodes.items():
        ORF_edges = ggCaller_cpp.read_edge_file(Edge_file_paths[colour_id])
        for node in node_set:                        
            for index, ORF_id in enumerate(G.nodes[node]['seqIDs']):
                parsed_id = ORF_id.split("_")

                # ensure colours between ORF and file match
                genome_id = int(parsed_id[0])
                if genome_id != colour_id:
                    continue

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
