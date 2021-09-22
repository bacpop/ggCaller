from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter, defaultdict
import networkx as nx
from panaroo.clean_network import collapse_paralogs
import numpy as np
from scipy.sparse import csc_matrix, lil_matrix
from intbitset import intbitset


def add_neighbour(G, edge_set, cluster_id_list, seq_to_cluster, high_scoring_ORF_edges, paralogs, cluster_centroids,
                  cluster_centroid_data, n_nodes, centroid_context):
    for neighbour in edge_set:
        # parse neighbour information
        neighbour_genome_id = cluster_id_list[neighbour][0]
        neighbour_local_id = cluster_id_list[neighbour][1]
        neighbour_cluster = seq_to_cluster[neighbour]

        # check that neighbour exists, if not add to graph
        add_neighbour = True if not G.has_node(neighbour_cluster) else False

        # check if current cluster contains paralogs
        neighbour_has_paralogs = True if neighbour_cluster in paralogs else False

        # parse neighbour information for current ORF
        neighbour_edge_set = high_scoring_ORF_edges[neighbour_genome_id][neighbour_local_id]

        # check if current ORF is end of contig
        neighbour_has_end = True if len(neighbour_edge_set) < 2 else False

        if add_neighbour:
            G.add_node(
                neighbour_cluster,
                size=1,
                centroid=[cluster_centroids[neighbour_cluster]],
                maxLenId=0,
                members=intbitset([neighbour_genome_id]),
                seqIDs=set([neighbour]),
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

            # check if paralog exists, if so add paralog node
            if neighbour_has_paralogs:
                # create a new paralog
                n_nodes += 1
                paralog_node = n_nodes
                centroid_context[
                    cluster_centroids[neighbour_cluster]].append(
                    [paralog_node, neighbour_genome_id])
                G.add_node(
                    paralog_node,
                    size=1,
                    centroid=[cluster_centroids[neighbour_cluster]],
                    maxLenId=0,
                    members=intbitset([neighbour_genome_id]),
                    seqIDs=set([neighbour]),
                    hasEnd=neighbour_has_end,
                    protein=[
                        cluster_centroid_data[neighbour_cluster]['prot_sequence']
                    ],
                    dna=[
                        cluster_centroid_data[neighbour_cluster]['dna_sequence']
                    ],
                    annotation=cluster_centroid_data[neighbour_cluster]
                    ['annotation'],
                    description=cluster_centroid_data[neighbour_cluster]
                    ['description'],
                    lengths=[
                        len(cluster_centroid_data[neighbour_cluster]
                            ['dna_sequence'])
                    ],
                    longCentroidID=(len(cluster_centroid_data[neighbour_cluster]
                                        ['dna_sequence']),
                                    cluster_centroids[neighbour_cluster]),
                    paralog=neighbour_has_paralogs,
                    mergedDNA=False)


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

    # create a dictionary of indexes for paralog context
    centroid_index = {}
    for i, centroid in enumerate(cluster_centroids.values()):
        centroid_index[centroid] = i
    n_centroids = len(centroid_index.keys())

    # build graph using adjacency information and optionally split paralogs
    G = nx.Graph()
    centroid_context = defaultdict(list)
    n_nodes = len(cluster_members)
    # temp_nodes = []
    # prev = None

    # iterating over each cluster, adding edges between clusters that each have a connection
    for current_cluster, ORF_members in cluster_members.items():
        # check if cluster is currently in graph
        add_node = True if not G.has_node(current_cluster) else False

        # check if current cluster contains paralogs
        has_paralogs = True if current_cluster in paralogs else False

        for ORF_id in ORF_members:
            # parse genome id and local ORF id
            genome_id = cluster_id_list[ORF_id][0]
            local_id = cluster_id_list[ORF_id][1]

            # map ORF ID to centroid ID
            seqid_to_centroid[ORF_id] = cluster_centroids[current_cluster]

            # parse neighbour information for current ORF
            edge_set = high_scoring_ORF_edges[genome_id][local_id]

            # check if current ORF is end of contig
            has_end = True if len(edge_set) < 2 else False

            # initialise paralog node
            paralog_node = None

            if add_node:
                G.add_node(
                    current_cluster,
                    size=1,
                    centroid=[cluster_centroids[current_cluster]],
                    maxLenId=0,
                    members=intbitset([genome_id]),
                    seqIDs=set([ORF_id]),
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
                add_node = False

                # check if paralog exists, if so add paralog node
                if has_paralogs:
                    # create a new paralog
                    n_nodes += 1
                    paralog_node = n_nodes
                    centroid_context[
                        cluster_centroids[current_cluster]].append(
                        [paralog_node, genome_id])
                    G.add_node(
                        paralog_node,
                        size=1,
                        centroid=[cluster_centroids[current_cluster]],
                        maxLenId=0,
                        members=intbitset([genome_id]),
                        seqIDs=set([ORF_id]),
                        hasEnd=has_end,
                        protein=[
                            cluster_centroid_data[current_cluster]['prot_sequence']
                        ],
                        dna=[
                            cluster_centroid_data[current_cluster]['dna_sequence']
                        ],
                        annotation=cluster_centroid_data[current_cluster]
                        ['annotation'],
                        description=cluster_centroid_data[current_cluster]
                        ['description'],
                        lengths=[
                            len(cluster_centroid_data[current_cluster]
                                ['dna_sequence'])
                        ],
                        longCentroidID=(len(cluster_centroid_data[current_cluster]
                                            ['dna_sequence']),
                                        cluster_centroids[current_cluster]),
                        paralog=has_paralogs,
                        mergedDNA=False)

            # add edges between node and all neighbours in neighbour set
            for neighbour in edge_set:
                # iterate through edge_set, adding nodes and then adding edges if required
                for neighbour in edge_set:
                    # parse neighbour information
                    neighbour_genome_id = cluster_id_list[neighbour][0]
                    neighbour_local_id = cluster_id_list[neighbour][1]
                    neighbour_cluster = seq_to_cluster[neighbour]

                    # check that neighbour exists, if not add to graph
                    add_neighbour = True if not G.has_node(neighbour_cluster) else False

                    # check if current cluster contains paralogs
                    neighbour_has_paralogs = True if neighbour_cluster in paralogs else False

                    # parse neighbour information for current ORF
                    neighbour_edge_set = high_scoring_ORF_edges[neighbour_genome_id][neighbour_local_id]

                    # check if current ORF is end of contig
                    neighbour_has_end = True if len(neighbour_edge_set) < 2 else False

                    if add_neighbour:
                        G.add_node(
                            neighbour_cluster,
                            size=1,
                            centroid=[cluster_centroids[neighbour_cluster]],
                            maxLenId=0,
                            members=intbitset([neighbour_genome_id]),
                            seqIDs=set([neighbour]),
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

                        # check if paralog exists, if so add paralog node
                        if neighbour_has_paralogs:
                            # create a new paralog
                            n_nodes += 1
                            paralog_node = n_nodes
                            centroid_context[
                                cluster_centroids[neighbour_cluster]].append(
                                [paralog_node, neighbour_genome_id])
                            G.add_node(
                                paralog_node,
                                size=1,
                                centroid=[cluster_centroids[neighbour_cluster]],
                                maxLenId=0,
                                members=intbitset([neighbour_genome_id]),
                                seqIDs=set([neighbour]),
                                hasEnd=neighbour_has_end,
                                protein=[
                                    cluster_centroid_data[neighbour_cluster]['prot_sequence']
                                ],
                                dna=[
                                    cluster_centroid_data[neighbour_cluster]['dna_sequence']
                                ],
                                annotation=cluster_centroid_data[neighbour_cluster]
                                ['annotation'],
                                description=cluster_centroid_data[neighbour_cluster]
                                ['description'],
                                lengths=[
                                    len(cluster_centroid_data[neighbour_cluster]
                                        ['dna_sequence'])
                                ],
                                longCentroidID=(len(cluster_centroid_data[neighbour_cluster]
                                                    ['dna_sequence']),
                                                cluster_centroids[neighbour_cluster]),
                                paralog=neighbour_has_paralogs,
                                mergedDNA=False)

            # need to think about connecting each current node and paralog with each neighbour node and paralog.
            # also need to add a paralog for every genome a paralog is found, not just one per cluster!

            # check if ORF has more than 2 neighbours. If not, at end of contig
            if len(neighbour_set) < 2:
                # we're at the start of a contig
                if G.has_node(current_cluster) and (current_cluster not in paralogs):
                    G.nodes[current_cluster]['size'] += 1
                    G.nodes[current_cluster]['members'].add(genome_id)
                    G.nodes[current_cluster]['seqIDs'].add(ORF_id)
                    G.nodes[current_cluster]['hasEnd'] = True
                    G.nodes[current_cluster]['lengths'].append(
                        len(cluster_centroid_data[current_cluster]['dna_sequence']))
                    if all_dna:
                        G.nodes[current_cluster]['dna'] += [
                            cluster_centroid_data[current_cluster]['dna_sequence']
                        ]
                # need to iterate over neightbour_set, adding all connections that have not already been added to graph
                else:
                    if current_cluster in paralogs:
                        # create a new paralog
                        n_nodes += 1
                        paralog_node = n_nodes
                        temp_nodes.append(paralog_node)
                        centroid_context[
                            cluster_centroids[current_cluster]].append(
                            [paralog_node, genome_id])
                    # add non paralog node
                    G.add_node(
                        paralog_node,
                        size=1,
                        centroid=[cluster_centroids[current_cluster]],
                        maxLenId=0,
                        members=intbitset([genome_id]),
                        seqIDs=set([ORF_id]),
                        hasEnd=True,
                        protein=[
                            cluster_centroid_data[current_cluster]['prot_sequence']
                        ],
                        dna=[
                            cluster_centroid_data[current_cluster]['dna_sequence']
                        ],
                        annotation=cluster_centroid_data[current_cluster]
                        ['annotation'],
                        description=cluster_centroid_data[current_cluster]
                        ['description'],
                        lengths=[
                            len(cluster_centroid_data[current_cluster]
                                ['dna_sequence'])
                        ],
                        longCentroidID=(len(cluster_centroid_data[current_cluster]
                                            ['dna_sequence']),
                                        cluster_centroids[current_cluster]),
                        paralog=(current_cluster in paralogs),
                        mergedDNA=False)
            else:
                is_paralog = current_cluster in paralogs
                if is_paralog:
                    # create a new paralog
                    n_nodes += 1
                    paralog_node = n_nodes
                    temp_nodes.append(paralog_node)
                    centroid_context[cluster_centroids[current_cluster]].append(
                        [paralog_node, genome_id])
                    G.add_node(
                        paralog_node,
                        size=1,
                        centroid=[cluster_centroids[current_cluster]],
                        maxLenId=0,
                        members=intbitset([genome_id]),
                        seqIDs=set([ORF_id]),
                        hasEnd=False,
                        protein=[
                            cluster_centroid_data[current_cluster]['prot_sequence']
                        ],
                        dna=[
                            cluster_centroid_data[current_cluster]['dna_sequence']
                        ],
                        annotation=cluster_centroid_data[current_cluster]
                        ['annotation'],
                        description=cluster_centroid_data[current_cluster]
                        ['description'],
                        lengths=[
                            len(cluster_centroid_data[current_cluster]
                                ['dna_sequence'])
                        ],
                        longCentroidID=(len(cluster_centroid_data[current_cluster]
                                            ['dna_sequence']),
                                        cluster_centroids[current_cluster]),
                        paralog=True,
                        mergedDNA=False)
                    # add edge between nodes

                    G.add_edge(prev,
                               neighbour,
                               size=1,
                               members=intbitset([genome_id]))
                else:
                    if not G.has_node(current_cluster):
                        # we need to add the gene in
                        G.add_node(
                            current_cluster,
                            size=1,
                            centroid=[cluster_centroids[current_cluster]],
                            maxLenId=0,
                            members=intbitset([genome_id]),
                            seqIDs=set([ORF_id]),
                            hasEnd=False,
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
                            paralog=is_paralog,
                            mergedDNA=False)
                        # add edge between nodes
                        G.add_edge(prev,
                                   current_cluster,
                                   size=1,
                                   members=intbitset([genome_id]))
                    else:
                        G.nodes[current_cluster]['size'] += 1
                        G.nodes[current_cluster]['members'].add(genome_id)
                        G.nodes[current_cluster]['seqIDs'].add(id)
                        G.nodes[current_cluster]['lengths'].append(
                            len(cluster_centroid_data[current_cluster]
                                ['dna_sequence']))
                        if all_dna:
                            G.nodes[current_cluster]['dna'] += [
                                cluster_centroid_data[current_cluster]
                                ['dna_sequence']
                            ]
                        if G.has_edge(prev, current_cluster):
                            G[prev][current_cluster]['size'] += 1
                            G[prev][current_cluster]['members'].add(genome_id)
                        else:
                            G.add_edge(prev,
                                       current_cluster,
                                       size=1,
                                       members=intbitset([genome_id]))
                    prev = current_cluster

    return G, centroid_context, seqid_to_centroid
