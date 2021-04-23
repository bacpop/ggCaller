import sys
import graph_tool.all as gt
from balrog.__main__ import *
import ggCaller_cpp
import collections
import numpy as np

try:
    from multiprocessing import shared_memory
    from multiprocessing.managers import SharedMemoryManager

    NumpyShared = collections.namedtuple('NumpyShared', ('name', 'shape', 'dtype'))
except ImportError as e:
    sys.stderr.write("This version of ggCaller requires python v3.8 or higher\n")
    sys.exit(1)


# @profile
def traverse_components(component, tc, component_list, edge_weights, minimum_path_score):
    # initilise high scoring ORF set to return
    high_scoring_ORFs = set()

    # generate subgraph view
    u = gt.GraphView(tc, vfilt=component_list == component)

    # initialise list of high scoring ORFs and their associated score for single component
    high_scoring_ORFs_temp = []

    high_score_temp = 0

    # iterate over edges, determine which are source and sink nodes for connected components
    vertices = u.get_vertices()
    in_degs = u.get_in_degrees(u.get_vertices())
    out_degs = u.get_out_degrees((u.get_vertices()))

    # calculate start nodes and add to list, if in-degree is 0
    start_vertices = [vertices[i] for i in range(len(in_degs)) if in_degs[i] == 0]

    # calculate end nodes and add to list, if in-degree is 0
    end_vertices = [vertices[i] for i in range(len(out_degs)) if out_degs[i] == 0]

    # iterate over start and stop vertices
    for start in start_vertices:
        # add the score of the first node
        start_score = u.vertex_properties["score"][start]

        # check if start vertex is lone node
        if start in end_vertices:
            vertex_list = [start]
            # replace high_score_ORFs with new high scoring ORF path
            if start_score > high_score_temp:
                high_score_temp = start_score
                high_scoring_ORFs_temp = vertex_list
        else:
            for end in end_vertices:
                score = start_score
                vertex_list, edge_list = gt.shortest_path(u, start, end, weights=edge_weights,
                                                          negative_weights=True, dag=True)
                # ensure a path has been found. If not, pass.
                if edge_list:
                    for e in edge_list:
                        # add the inverse of the score to get the true score
                        score -= edge_weights[e]

                    # replace high_score_ORFs with new high scoring ORF path
                    if score > high_score_temp:
                        high_score_temp = score
                        high_scoring_ORFs_temp = vertex_list

    # for highest scoring path, see if greater than cut-off, if so, add high scoring ORFs to set
    if high_score_temp >= minimum_path_score:
        for node in high_scoring_ORFs_temp:
            high_scoring_ORFs.add(u.vertex_properties["ID"][node])

    return high_scoring_ORFs


#@profile
def call_true_genes(ORF_score_dict, ORF_overlap_dict, minimum_path_score):
    # initilise high scoring ORF set to return
    high_scoring_ORFs_all = set()

    # create empty, directed graph
    g = gt.Graph()

    # create a dictionaries assign each ORF an index in the graph to vertices
    ORF_index = {}

    # add vertexes to graph, store ORF information in ORF_index
    for ORF in ORF_score_dict.keys():
        v = g.add_vertex()
        ORF_index[ORF] = g.vertex_index[v]

    # add edges and edge weights between connected ORFs using ORF_overlap_dict. ORF1 is sink, ORF2 is source
    for ORF1, overlap_dict in ORF_overlap_dict.items():
        # check that ORF1 nodes exist in graph, may have been removed due to low score
        if ORF1 in ORF_index:
            for ORF2 in overlap_dict.keys():
                # check that ORF2 nodes exist in graph as before
                if ORF2 in ORF_index:
                    # add new edge between the two ORFs, where ORF2 is the source and ORF1 is the sink
                    e = g.add_edge(g.vertex(ORF_index[ORF2]), g.vertex(ORF_index[ORF1]))

    # determine if cycles present. If so, break them by removing edge before repeated node and re-test
    cycle = True
    try:
        circuit = next(gt.all_circuits(g))
    except StopIteration:
        cycle = False

    while cycle:
        end_cycle = circuit[-1]
        start_cycle = circuit[0]
        e = g.edge(end_cycle, start_cycle)
        g.remove_edge(e)
        try:
            circuit = next(gt.all_circuits(g))
        except StopIteration:
            break

    # generate a transative closure of the graph to add all directed edges and add vertex properties
    tc = gt.transitive_closure(g)

    # clear original graph
    g.clear()

    # get label components of tc
    components = gt.label_components(tc, directed=False)[0].a

    # create vertex property map to store node IDs and scores
    vertex_ID = tc.new_vertex_property("int")
    vertex_score = tc.new_vertex_property("double")
    for ORF_ID, index in ORF_index.items():
        vertex_ID[tc.vertex(index)] = ORF_ID
        vertex_score[tc.vertex(index)] = ORF_score_dict[ORF_ID]

    # add vertex IDs and scores as internal property of graph
    tc.vertex_properties["ID"] = vertex_ID
    tc.vertex_properties["score"] = vertex_score

    # create edge weights property
    edge_weights = tc.new_edge_property("double")

    # create incompatible list for edge iteration
    incomp_edges = []

    # iterate over edges and assign weights, and identify invalid overlaps to remove
    for e in tc.edges():
        # parse sink and source nodes, order same as before
        ORF1 = tc.vertex_properties["ID"][e.target()]
        ORF2 = tc.vertex_properties["ID"][e.source()]

        # parse sink ORF score calculated by Balrog
        ORF_score = ORF_score_dict[ORF1]

        # check if edges are present by overlap detection. If not, set edge weight as ORF1 score
        if ORF1 not in ORF_overlap_dict:
            edge_weights[e] = -(ORF_score)
            continue
        elif ORF2 not in ORF_overlap_dict[ORF1]:
            edge_weights[e] = -(ORF_score)
            continue
        # parse overlap info, calculate edge weight
        overlap_type, abs_overlap = ORF_overlap_dict[ORF1][ORF2]

        # set penalty to 0 for "n" or "w" or "i" overlaps
        penalty = 0
        if overlap_type == "i":
            # add incompatible ORFs as tuple, added as source (1) and sink(2)
            incomp_edges.append(e)
        elif overlap_type == "u":
            penalty = unidirectional_penalty_per_base * abs_overlap
        elif overlap_type == "c":
            penalty = convergent_penalty_per_base * abs_overlap
        elif overlap_type == "d":
            penalty = divergent_penalty_per_base * abs_overlap
        # set ORF score to negative so that shortest algorithm can be applied
        edge_weights[e] = -(ORF_score - penalty)

    # add edge_weights as internal property
    tc.edge_properties["weights"] = edge_weights

    # iterate over incompatible edges to remove them
    for e in incomp_edges:
        tc.remove_edge(e)

    # iterate over components, find highest scoring path within component with multiprocessing to determine geniest path through components
    for component in set(components):
        high_scoring_ORFs = traverse_components(component, tc, components, edge_weights, minimum_path_score)
        high_scoring_ORFs_all.update(high_scoring_ORFs)

    return high_scoring_ORFs_all


def run_calculate_ORFs(node_set_tuple, graph_vector, repeat, overlap, max_path_length, is_ref, no_filter,
                       stop_codons_for, start_codons, min_ORF_length, max_ORF_overlap, minimum_ORF_score,
                       minimum_path_score, write_idx, input_colours, aa_kmer_set, model, model_tis):
    # unpack tuple
    colour_ID, node_set = node_set_tuple

    # print("Started analysing: " + str(colour_ID))

    # load shared memory items
    graph_vector_shm = shared_memory.SharedMemory(name=graph_vector.name)
    graph_vector = np.ndarray(graph_vector.shape, dtype=graph_vector.dtype, buffer=graph_vector_shm.buf)

    # aa_kmer_set_shm = shared_memory.SharedMemory(name=aa_kmer_set.name)
    # aa_kmer_set = np.ndarray(aa_kmer_set.shape, dtype=aa_kmer_set.dtype, buffer=aa_kmer_set_shm.buf)

    model_shm = shared_memory.SharedMemory(name=model.name)
    model = np.ndarray(model.shape, dtype=model.dtype, buffer=model_shm.buf)

    model_tis_shm = shared_memory.SharedMemory(name=model_tis.name)
    model_tis = np.ndarray(model_tis.shape, dtype=model_tis.dtype, buffer=model_tis_shm.buf)

    # parse data from np_arrays
    model_obj = model[0]
    model_tis_obj = model_tis[0]
    graph_vector_list = graph_vector.tolist()

    # determine all ORFs in Bifrost graph
    ORF_overlap_dict, ORF_vector = ggCaller_cpp.calculate_ORFs(graph_vector_list, colour_ID, node_set, repeat,
                                                               overlap, max_path_length, is_ref, no_filter,
                                                               stop_codons_for, start_codons, min_ORF_length,
                                                               max_ORF_overlap, write_idx, input_colours[colour_ID])

    # if no filter specified, just copy ORF_vector to true_genes
    if no_filter:
        true_genes = ORF_vector
    else:
        # calculate scores for genes
        ORF_score_dict = score_genes(ORF_vector, graph_vector, minimum_ORF_score, overlap, model_obj, model_tis_obj,
                                     aa_kmer_set)

        # determine highest scoring genes
        high_scoring_ORFs = call_true_genes(ORF_score_dict, ORF_overlap_dict, minimum_path_score)

        # initiate true genes list
        true_genes = [None] * len(high_scoring_ORFs)

        for index, ORF_id in enumerate(high_scoring_ORFs):
            # add only high scoring ORFs to true_genes
            true_genes[index] = ORF_vector[ORF_id]

    #print("Finished analysing: " + str(colour_ID))

    return colour_ID, true_genes


def generate_shared_mem_array(in_array, smm):
    """Generates a shared memory representation of a numpy array"""
    array_raw = smm.SharedMemory(size=in_array.nbytes)
    array_shared = np.ndarray(in_array.shape, dtype=in_array.dtype, buffer=array_raw.buf)
    array_shared[:] = in_array[:]
    array_shared = NumpyShared(name=array_raw.name, shape=in_array.shape, dtype=in_array.dtype)
    return (array_shared)
