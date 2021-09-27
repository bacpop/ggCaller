from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter, defaultdict
import networkx as nx
from panaroo.clean_network import collapse_paralogs
import numpy as np
from scipy.sparse import csc_matrix, lil_matrix
from intbitset import intbitset


def generate_network(DBG, high_scoring_ORFs, high_scoring_ORF_edges,
                     cluster_id_list, cluster_dict, overlap, all_dna=False):
    # associate sequences with their clusters
    seq_to_cluster = {}
    seqid_to_centroid = {}
    cluster_centroids = {}
    cluster_members = defaultdict(list)
    # iterate over cluster_dict, parsing all sequences to clusters
    cluster_id = 0
    for cluster, ORF_list in cluster_dict.items():
        # append centroid to cluster
        cluster_centroids[cluster_id] = cluster
        # for each ORF, add ORF_id and cluster_id to respective dictionaries
        for ORF_id in ORF_list:
            seq_to_cluster[ORF_id] = cluster_id
            cluster_members[cluster_id].append(ORF_id)

            # also add population ID to high_scoring_ORFs for later reference
            genome_id = cluster_id_list[ORF_id][0]
            local_id = cluster_id_list[ORF_id][1]
            ORF_tuple = high_scoring_ORFs[genome_id][local_id]
            ORF_tuple = list(ORF_tuple)
            ORF_tuple[6] = ORF_id
            high_scoring_ORFs[genome_id][local_id] = ORF_tuple
        cluster_id += 1

    # determine paralogs if required
    paralogs = set()
    for clust in cluster_members:
        # determine if paralogs present by checking colours for each entry in cluster_id_list
        genomes = [cluster_id_list[s][0] for s in cluster_members[clust]]
        if len(genomes) != len(set(genomes)):
            paralogs.add(clust)

    # Load meta data such as sequence and annotation
    cluster_centroid_data = {}
    for cluster_id, centroid_id in cluster_centroids.items():
        # parse genome id and local ORF id
        genome_id = cluster_id_list[centroid_id][0]
        local_id = cluster_id_list[centroid_id][1]

        # access ORF information for centroid from high_scoring_ORFs
        ORFNodeVector = high_scoring_ORFs[genome_id][local_id]
        ORF_nodelist = ORFNodeVector[0]
        ORF_node_coords = ORFNodeVector[1]

        # generate dna sequence and protein sequence
        dna_sequence = DBG.generate_sequence(ORF_nodelist, ORF_node_coords, overlap)
        prot_sequence = str(Seq(dna_sequence).translate())

        # add information to cluster_centroid_data
        cluster_centroid_data[cluster_id] = {
            'prot_sequence': prot_sequence,
            'dna_sequence': dna_sequence,
            'annotation': '',
            'description': '',
        }

    # build graph using adjacency information and optionally split paralogs
    G = nx.Graph()
    centroid_context = defaultdict(list)
    n_nodes = len(cluster_members)

    # iterating over each cluster, adding edges between clusters that each have a connection
    for current_cluster, ORF_members in cluster_members.items():
        # check if cluster is currently in graph
        add_cluster = True if not G.has_node(current_cluster) else False

        # check if current cluster contains paralogs
        has_paralogs = True if current_cluster in paralogs else False

        for ORF_id in ORF_members:
            # parse genome id and local ORF id
            genome_id = cluster_id_list[ORF_id][0]
            local_id = cluster_id_list[ORF_id][1]

            # map ORF ID to centroid ID
            seqid_to_centroid[str(ORF_id)] = cluster_centroids[current_cluster]

            # parse neighbour information for current ORF
            edge_set = high_scoring_ORF_edges[genome_id][local_id]

            # check if current ORF is end of contig
            has_end = True if len(edge_set) == 0 else False

            # initialise cluster to add
            cluster_to_add = current_cluster

            # # initialise paralog node
            # paralog_node = None

            if add_cluster:
                if has_paralogs:
                    # create a new paralog
                    n_nodes += 1
                    cluster_to_add = n_nodes
                    centroid_context[
                        cluster_centroids[current_cluster]].append(
                        [cluster_to_add, genome_id])
                G.add_node(
                    cluster_to_add,
                    size=1,
                    centroid=[cluster_centroids[current_cluster]],
                    maxLenId=0,
                    members=intbitset([genome_id]),
                    seqIDs=set([str(ORF_id)]),
                    hasEnd=has_end,
                    protein=[
                        cluster_centroid_data[current_cluster]
                        ['prot_sequence']
                    ],
                    dna=[
                        cluster_centroid_data[current_cluster]
                        ['dna_sequence']
                    ],
                    annotation=cluster_centroid_data[current_cluster]
                    ['annotation'],
                    description=cluster_centroid_data[current_cluster]
                    ['description'],
                    lengths=[
                        len(cluster_centroid_data[current_cluster]
                            ['dna_sequence'])
                    ],
                    longCentroidID=(len(
                        cluster_centroid_data[current_cluster]
                        ['dna_sequence']),
                                    cluster_centroids[current_cluster]),
                    paralog=has_paralogs,
                    mergedDNA=False)
                add_cluster = False

            else:
                # check if ORF_id already added to the cluster
                if str(ORF_id) not in G.nodes[current_cluster]['seqIDs']:
                    G.nodes[current_cluster]['size'] += 1
                    G.nodes[current_cluster]['members'].add(genome_id)
                    G.nodes[current_cluster]['seqIDs'].add(str(ORF_id))
                    if G.nodes[current_cluster]['hasEnd'] == False: G.nodes[current_cluster]['hasEnd'] = has_end
                    G.nodes[current_cluster]['lengths'].append(
                        len(cluster_centroid_data[current_cluster]
                            ['dna_sequence']))
                    if all_dna:
                        G.nodes[current_cluster]['dna'] += [
                            cluster_centroid_data[current_cluster]
                            ['dna_sequence']
                        ]

            # # check if paralog exists, if so add paralog node
            # if has_paralogs:
            #     # create a new paralog
            #     n_nodes += 1
            #     paralog_node = n_nodes
            #     centroid_context[
            #         cluster_centroids[current_cluster]].append(
            #         [paralog_node, genome_id])
            #     G.add_node(
            #         paralog_node,
            #         size=1,
            #         centroid=[cluster_centroids[current_cluster]],
            #         maxLenId=0,
            #         members=intbitset([genome_id]),
            #         seqIDs=set([ORF_id]),
            #         hasEnd=has_end,
            #         protein=[
            #             cluster_centroid_data[current_cluster]['prot_sequence']
            #         ],
            #         dna=[
            #             cluster_centroid_data[current_cluster]['dna_sequence']
            #         ],
            #         annotation=cluster_centroid_data[current_cluster]
            #         ['annotation'],
            #         description=cluster_centroid_data[current_cluster]
            #         ['description'],
            #         lengths=[
            #             len(cluster_centroid_data[current_cluster]
            #                 ['dna_sequence'])
            #         ],
            #         longCentroidID=(len(cluster_centroid_data[current_cluster]
            #                             ['dna_sequence']),
            #                         cluster_centroids[current_cluster]),
            #         paralog=has_paralogs,
            #         mergedDNA=False)

            # iterate through edge_set, adding nodes and then adding edges if required
            for neighbour in edge_set:
                # parse neighbour information from high_scoring_ORFs
                neighbour_id = high_scoring_ORFs[genome_id][neighbour][6]
                # neighbour_genome_id = cluster_id_list[neighbour][0]
                # neighbour_local_id = cluster_id_list[neighbour][1]
                neighbour_cluster = seq_to_cluster[neighbour_id]

                # check that neighbour exists, if not add to graph
                add_neighbour = True if not G.has_node(neighbour_cluster) else False

                # check if current cluster contains paralogs
                neighbour_has_paralogs = True if neighbour_cluster in paralogs else False

                # initialise cluster to add
                neighbour_cluster_to_add = neighbour_cluster

                # # initialise neighbour paralog node
                # neighbour_paralog_node = None

                # # parse neighbour information for current ORF
                neighbour_edge_set = high_scoring_ORF_edges[genome_id][neighbour]

                # check if current ORF is end of contig
                neighbour_has_end = True if len(neighbour_edge_set) == 0 else False

                # add neighbour if not present in graph already
                if add_neighbour:
                    if neighbour_has_paralogs:
                        # create a new paralog
                        n_nodes += 1
                        neighbour_cluster_to_add = n_nodes
                        centroid_context[
                            cluster_centroids[neighbour_cluster]].append(
                            [neighbour_cluster_to_add, genome_id])
                    G.add_node(
                        neighbour_cluster_to_add,
                        size=1,
                        centroid=[cluster_centroids[neighbour_cluster]],
                        maxLenId=0,
                        members=intbitset([genome_id]),
                        seqIDs=set([str(neighbour_id)]),
                        hasEnd=neighbour_has_end,
                        protein=[
                            cluster_centroid_data[neighbour_cluster]
                            ['prot_sequence']
                        ],
                        dna=[
                            cluster_centroid_data[neighbour_cluster]
                            ['dna_sequence']
                        ],
                        annotation=cluster_centroid_data[neighbour_cluster]
                        ['annotation'],
                        description=cluster_centroid_data[neighbour_cluster]
                        ['description'],
                        lengths=[
                            len(cluster_centroid_data[neighbour_cluster]
                                ['dna_sequence'])
                        ],
                        longCentroidID=(len(
                            cluster_centroid_data[neighbour_cluster]
                            ['dna_sequence']),
                                        cluster_centroids[neighbour_cluster]),
                        paralog=neighbour_has_paralogs,
                        mergedDNA=False)
                else:
                    # check if neighbour_id already added to the cluster
                    if str(neighbour_id) not in G.nodes[neighbour_cluster]['seqIDs']:
                        G.nodes[neighbour_cluster]['size'] += 1
                        G.nodes[neighbour_cluster]['members'].add(genome_id)
                        G.nodes[neighbour_cluster]['seqIDs'].add(str(neighbour_id))
                        if G.nodes[neighbour_cluster]['hasEnd'] == False:
                            G.nodes[neighbour_cluster]['hasEnd'] = neighbour_has_end
                        G.nodes[neighbour_cluster]['lengths'].append(
                            len(cluster_centroid_data[neighbour_cluster]
                                ['dna_sequence']))
                        if all_dna:
                            G.nodes[neighbour_cluster]['dna'] += [
                                cluster_centroid_data[neighbour_cluster]
                                ['dna_sequence']
                            ]

                # # check if paralog exists, if so add paralog node
                # if neighbour_has_paralogs:
                #     # create a new paralog
                #     n_nodes += 1
                #     neighbour_paralog_node = n_nodes
                #     centroid_context[
                #         cluster_centroids[neighbour_cluster]].append(
                #         [paralog_node, neighbour_genome_id])
                #     G.add_node(
                #         neighbour_paralog_node,
                #         size=1,
                #         centroid=[cluster_centroids[neighbour_cluster]],
                #         maxLenId=0,
                #         members=intbitset([neighbour_genome_id]),
                #         seqIDs=set([neighbour]),
                #         hasEnd=neighbour_has_end,
                #         protein=[
                #             cluster_centroid_data[neighbour_cluster]['prot_sequence']
                #         ],
                #         dna=[
                #             cluster_centroid_data[neighbour_cluster]['dna_sequence']
                #         ],
                #         annotation=cluster_centroid_data[neighbour_cluster]
                #         ['annotation'],
                #         description=cluster_centroid_data[neighbour_cluster]
                #         ['description'],
                #         lengths=[
                #             len(cluster_centroid_data[neighbour_cluster]
                #                 ['dna_sequence'])
                #         ],
                #         longCentroidID=(len(cluster_centroid_data[neighbour_cluster]
                #                             ['dna_sequence']),
                #                         cluster_centroids[neighbour_cluster]),
                #         paralog=neighbour_has_paralogs,
                #         mergedDNA=False)

                # add edge between current ORF and neighbour
                if G.has_edge(cluster_to_add, neighbour_cluster_to_add):
                    G[cluster_to_add][neighbour_cluster_to_add]['size'] += 1
                    G[cluster_to_add][neighbour_cluster_to_add]['members'].add(genome_id)
                else:
                    G.add_edge(cluster_to_add,
                               neighbour_cluster_to_add,
                               size=1,
                               members=intbitset([genome_id]))

                # # add edge between current ORF paralog and neighbour
                # if paralog_node is not None:
                #     if G.has_edge(paralog_node, neighbour):
                #         G[paralog_node][neighbour]['size'] += 1
                #         G[paralog_node][neighbour]['members'].add(genome_id)
                #     else:
                #         G.add_edge(paralog_node,
                #                    neighbour,
                #                    size=1,
                #                    members=intbitset([genome_id]))
                #
                # # add edge between current ORF and neighbour paralog
                # if neighbour_paralog_node is not None:
                #     if G.has_edge(current_cluster, neighbour):
                #         G[current_cluster][neighbour_paralog_node]['size'] += 1
                #         G[current_cluster][neighbour_paralog_node]['members'].add(genome_id)
                #     else:
                #         G.add_edge(current_cluster,
                #                    neighbour_paralog_node,
                #                    size=1,
                #                    members=intbitset([genome_id]))

                # currently don't add edge between current cluster paralog and neighbour paralog.

    return G, centroid_context, seqid_to_centroid
