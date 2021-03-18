import graph_tool.all as gt
from balrog.__main__ import *

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


def call_true_genes(colour_ORF_tuple, minimum_path_score, ORF_score_dict, ORF_overlap_dict):
    colour, ORF_ID_list = colour_ORF_tuple

    # initilise high scoring ORF set to return
    high_scoring_ORFs_all = set()

    # create empty, directed graph
    g = gt.Graph()

    # create a dictionaries assign each ORF an index in the graph to vertices
    ORF_index = {}

    # create new graph item with node labels
    # vertex_ID = g.new_vertex_property("int")

    # add vertexes to graph, store ORF information in ORF_index
    for ORF in ORF_ID_list:
        # check if ORF is present in ORF_score_dict, if not, score too low so do not add
        if ORF in ORF_score_dict:
            v = g.add_vertex()
            # vertex_ID[v] = ORF
            ORF_index[ORF] = g.vertex_index[v]

    # add vertex sequences to graph
    # g.vertex_properties["ID"] = vertex_ID

    # add edges and edge weights between connected ORFs using ORF_overlap_dict. ORF1 is sink, ORF2 is source
    if ORF_overlap_dict[colour]:
        for ORF1, overlap_dict in ORF_overlap_dict[colour].items():
            # check that ORF1 nodes exist in graph, may have been removed due to low score
            if ORF1 in ORF_index:
                for ORF2 in overlap_dict.keys():
                    # check that ORF2 nodes exist in graph as before
                    if (ORF2 in ORF_index):
                        # add new edge between the two ORFs, where ORF2 is the source and ORF1 is the sink
                        e = g.add_edge(g.vertex(ORF_index[ORF2]), g.vertex(ORF_index[ORF1]))

    # generate a transative closure of the graph to add all directed edges and add vertex properties
    tc = gt.transitive_closure(g)

    # clear original graph
    g.clear()

    # get label components of tc
    components = gt.label_components(tc, directed=False)[0].a

    # component_ids = set(component_assignments.a)

    # tc.vertex_properties["component"] = component_assignments

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
        if ORF1 not in ORF_overlap_dict[colour]:
            edge_weights[e] = -(ORF_score)
            continue
        elif ORF2 not in ORF_overlap_dict[colour][ORF1]:
            edge_weights[e] = -(ORF_score)
            continue
        # parse overlap info, calculate edge weight
        overlap_type, abs_overlap = ORF_overlap_dict[colour][ORF1][ORF2]

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

    # create graph view
    #tc.save(colour + "_1threads_tc.graphml")
    #g.save(colour + "_1threads_original.graphml")

    # for debugging if negative cycles encountered.
    # for c in gt.all_circuits(g):
    #    print(c)

    # iterate over components, find highest scoring path within component with multiprocessing to determine geniest path through components
    for component in set(components):
        high_scoring_ORFs = traverse_components(component, tc, components, edge_weights, minimum_path_score)
        high_scoring_ORFs_all.update(high_scoring_ORFs)

    return colour, high_scoring_ORFs_all


def update_colour(colour1, colour2):
    updated_colour = ""
    for i in range(0, len(colour1)):
        if colour1[i] == "1" or colour2[i] == "1":
            updated_colour += "1"
        else:
            updated_colour += "0"
    return updated_colour