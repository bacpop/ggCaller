"""
Python scripts that, if they are compatible, eliminate and compact nodes and edges
"""

import logging

GRAPH_LOGGER = logging.getLogger(__name__)

def tuple_to_string(node):
    """Given a tuple of a node id and orientation
    :return tuple element concatenete by |:
    """
    node_id, node_orn = node
    return node_id+"|"+node_orn

def reverse_and_complement(string):
    """Given a genomic string
    :return the reverese&complement of the string: If is specify
    :return the same string: If is not specify (*)
    """
    reverse_dict = dict([('A', 'T'), ('T', 'A'), ('C', 'G'), ('G', 'C'), ('*', '*')])
    complement_string = ''.join([reverse_dict[c] for c in string])
    return complement_string[::-1]

def reverse_strand(strand):
    """Given a strand
    :return the opposite strand: if is specify
    :return None: if strand is not defined
    """
    dict_inverted = dict([('+', '-'), ('-', '+'), (None, None)])
    return dict_inverted[strand]

def update_graph(gfa_, keep_id, remove_id, new_seq, overlap, orn):
    """Given the id of the nodes in the edge to compact, the new sequence,
    orientation and overlap
    calcolate the new sequence lenght, the edges IN or OUT from the node
    to remove are updated with the node to keep and the node with
    remove_id is removed.
    """
    gfa_.node()[keep_id]['sequence'] = new_seq
    if not new_seq == '*':
        gfa_.node()[keep_id]['slen'] = len(gfa_.node(keep_id)['sequence'])
    else:
        if gfa_.node()[keep_id]['slen'] and gfa_.node()[remove_id]['slen']:
            gfa_.node()[keep_id]['slen'] += gfa_.node()[remove_id]['slen'] - overlap
        else:
            gfa_.node()[keep_id]['slen'] = None

        if 'fu' in gfa_.node()[remove_id]:
            remove_fu = gfa_.node()[remove_id].get('fu').lstrip('Z:')
        else:
            remove_fu = remove_id

        if 'fu' in gfa_.node()[keep_id]:
            gfa_.node()[keep_id]['fu'] += '_'+remove_fu
        else:
            gfa_.node()[keep_id]['fu'] = 'Z:'+keep_id+'_'+remove_fu

    #fix gfa_.node()[keep_id]['option']

    remove_edge_dict = dict()
    data_update_edges = gfa_.edge()[remove_id]
    for node in data_update_edges:
        for edge_id in data_update_edges[node]:
            if data_update_edges[node][edge_id]['from_node'] == remove_id:
                parse_from_node = keep_id
                parse_from_orn = data_update_edges[node][edge_id]['from_orn']
                if not parse_from_orn:
                    parse_to_node = data_update_edges[node][edge_id]['to_node']
                    parse_to_orn = data_update_edges[node][edge_id]['to_orn']
                    parse_from_positions = data_update_edges[node][edge_id]['from_positions']
                    parse_to_positions = data_update_edges[node][edge_id]['to_positions']
                    keep_slen = int(gfa_.node()[keep_id]['slen'])
                    remove_slen = int(gfa_.node()[remove_id]['slen'])
                    if orn == '-+':
                        parse_to_orn = reverse_strand(parse_to_orn)
                        parse_from_positions = (str(remove_slen-int(parse_from_positions[1])), \
                            str(remove_slen-int(parse_from_positions[0])))
                    elif orn == '+-':
                        parse_to_orn = reverse_strand(parse_to_orn)
                        parse_from_positions = (str(keep_slen-int(parse_from_positions[1])), \
                            str(keep_slen-int(parse_from_positions[0])))
                    elif orn == '++':
                        parse_from_positions = (str(keep_slen-remove_slen+ \
                        	int(parse_from_positions[0])), \
                            str(keep_slen-remove_slen+int(parse_from_positions[1])))
                else:
                    if orn == '-+' or orn == '+-':
                        parse_from_orn = reverse_strand(parse_from_orn)
                    parse_to_node = data_update_edges[node][edge_id]['to_node']
                    parse_to_orn = data_update_edges[node][edge_id]['to_orn']
            else:
                parse_to_node = keep_id
                parse_to_orn = data_update_edges[node][edge_id]['to_orn']
                if orn == '-+' or orn == '+-':
                    parse_to_orn = reverse_strand(parse_to_orn)
                parse_from_node = data_update_edges[node][edge_id]['from_node']
                parse_from_orn = data_update_edges[node][edge_id]['from_orn']
            parse_overlap = data_update_edges[node][edge_id]['alignment']

            if parse_from_orn:
                parse_new_edge = 'L\t'+parse_from_node+'\t'+parse_from_orn+'\t'+ \
                    parse_to_node+'\t'+parse_to_orn+'\t'+parse_overlap
            else:
                parse_new_edge = 'F\t'+parse_from_node+'\t'+parse_to_node+parse_to_orn+'\t'+ \
                    parse_from_positions[0]+'\t'+parse_from_positions[1]+'\t'+ \
                    parse_to_positions[0]+'\t'+parse_to_positions[1]+'\t'+parse_overlap[0]

            gfa_.add_edge(parse_new_edge)
            remove_edge_dict[edge_id] = True

    for remove_edge_id in remove_edge_dict:
        gfa_.remove_edge(remove_edge_id)
    gfa_.remove_node(remove_id)

def compact_sequence(gfa_, from_node, to_node):
    """Given the id and orientation of the nodes in the edge to compact
    get the sequence and CIGAR information and
    make the new sequence string.
    :return the new string sequence, the aligment/overlap, and orientation.
    """
    from_id, from_orn = from_node
    to_id, to_orn = to_node

    edges = gfa_._search_edge_by_nodes((from_id, to_id))
    edge = None
    for edge in edges:
        overlap = int(edges[edge]['alignment'].rstrip('M'))
    gfa_.remove_edge(edge)

    from_seq = gfa_.node(from_id)['sequence']
    to_seq = gfa_.node(to_id)['sequence']

    if from_orn == '-' and to_orn == '-':
        if from_seq == '*' or to_seq == '*':
            new_seq = '*'
        else:
            new_seq = to_seq+from_seq[overlap:]
        return new_seq, overlap, '--'
    elif from_orn == '+' and to_orn == '+':
        if from_seq == '*' or to_seq == '*':
            new_seq = '*'
        else:
            new_seq = from_seq+to_seq[overlap:]
        return new_seq, overlap, '++'
    elif from_orn == '-':
        if from_seq == '*' or to_seq == '*':
            new_seq = '*'
        else:
            new_seq = reverse_and_complement(to_seq)[:-overlap]+from_seq
        return new_seq, overlap, '-+'
    else:
        if from_seq == '*' or to_seq == '*':
            new_seq = '*'
        else:
            new_seq = from_seq+reverse_and_complement(to_seq)[overlap:]
        return new_seq, overlap, '+-'

def update_dictionaries_by_nodes(nodes, count_dictionaries, orn):
    """Given the nodes to compact and the dictionaries
    the edges and nodes are update or deleted in the dictionaries.
    """
    keep_node, remove_node = nodes
    keep_id, keep_orn = keep_node
    remove_id, remove_orn = remove_node
    from_dict, to_dict = count_dictionaries

    from_dict[keep_id][keep_orn] = []
    to_dict[remove_id][remove_orn] = []

    if from_dict.get(remove_id):
        for orientation in from_dict.get(remove_id):
            if from_dict.get(remove_id).get(orientation):
                links = from_dict.get(remove_id).get(orientation)
                from_dict[remove_id][orientation] = []
                if orn == '+-' or orn == '-+':
                    orientation = reverse_strand(orientation)
                if not from_dict.get(keep_id):
                    from_dict[keep_id] = {}
                if from_dict.get(keep_id).get(orientation):
                    for link in links:
                        from_dict[keep_id][orientation].append(link)
                else:
                    from_dict[keep_id][orientation] = links

    if to_dict.get(remove_id):
        for orientation in to_dict.get(remove_id):
            if to_dict.get(remove_id).get(orientation):
                links = to_dict.get(remove_id).get(orientation)
                to_dict[remove_id][orientation] = []
                if orn == '+-' or orn == '-+':
                    orientation = reverse_strand(orientation)
                if not to_dict.get(keep_id):
                    to_dict[keep_id] = {}
                if to_dict.get(keep_id).get(orientation):
                    for link in links:
                        to_dict[keep_id][orientation].append(link)
                else:
                    to_dict[keep_id][orientation] = links

def compression_graph_by_nodes(gfa_):
    """Create dictionaries with edges's id and with nodes involved.
    By the dictionaries determine what edge can be removed. The ends
    of segment(node) with only 1 edge IN or OUT are compatible with compression.
    Passing data to the other functions update the graph and the dictionaries.
    """
    eid_dict = dict()
    from_dict = dict()
    to_dict = dict()

    data_edges = gfa_.edge()
    for node1 in data_edges:
        for node2 in data_edges[node1]:
            for eid in data_edges[node1][node2]:
                if not eid_dict.get(eid):
                    eid_dict[eid] = True

                    from_id = data_edges[node1][node2][eid]['from_node']
                    from_orn = data_edges[node1][node2][eid]['from_orn']
                    to_id = data_edges[node1][node2][eid]['to_node']
                    to_orn = data_edges[node1][node2][eid]['to_orn']

                    if from_dict.get(from_id):
                        from_dict[from_id][from_orn].append(tuple_to_string((to_id, to_orn)))
                    else:
                        from_dict[from_id] = {}
                        from_dict[from_id][from_orn] = [tuple_to_string((to_id, to_orn))]
                        from_dict[from_id][reverse_strand(from_orn)] = []

                    if to_dict.get(to_id):
                        to_dict[to_id][to_orn].append(tuple_to_string((from_id, from_orn)))
                    else:
                        to_dict[to_id] = {}
                        to_dict[to_id][to_orn] = [tuple_to_string((from_id, from_orn))]
                        to_dict[to_id][reverse_strand(to_orn)] = []

    edge_no_text = 0
    count_edge_compacted = 0
    for from_id, orientation in from_dict.items():
        for from_orn in orientation:
            if len(orientation.get(from_orn)) == 1:
                to_node = orientation.get(from_orn)
                to_id, to_orn = to_node[0].split('|')
                if gfa_.node(from_id) and gfa_.node(to_id):
                    if len(to_dict[to_id][to_orn]) == 1:
                        if to_dict.get(from_id):
                            if to_dict.get(from_id).get(reverse_strand(from_orn)):
                                break
                        if from_dict.get(to_id):
                            if from_dict.get(to_id).get(reverse_strand(to_orn)):
                                break
                        from_node = (from_id, from_orn)
                        to_node = (to_id, to_orn)
                        new_seq, overlap, orn = compact_sequence(gfa_, from_node, to_node)
                        update_graph(gfa_, from_node[0], to_node[0], new_seq, overlap, orn)
                        update_dictionaries_by_nodes((from_node, to_node), \
                            (from_dict, to_dict), orn)
                        count_edge_compacted += 1
                else:
                    edge_no_text += 1

    GRAPH_LOGGER.debug('%s edges has been compacted of total amount of %s', \
    	count_edge_compacted, len(eid_dict)-edge_no_text)

    return count_edge_compacted

def update_dictionaries_by_edges(eid_dict, count_dictionaries, eid, orn):
    """Given the list with the index of edge to compact and orientation
    the edges and nodes are update or deleted in the list.
    """
    keep_node = (eid_dict[eid]['from_id'], eid_dict[eid]['from_orn'])
    keep_id, _ = keep_node
    remove_node = (eid_dict[eid]['to_id'], eid_dict[eid]['to_orn'])
    remove_id, _ = remove_node
    from_dict, to_dict = count_dictionaries

    eid_dict[eid] = False
    from_dict[tuple_to_string(keep_node)] -= 1
    to_dict[tuple_to_string(remove_node)] -= 1

    for edge in eid_dict:
        if eid_dict[edge]:
            from_node = (eid_dict[edge]['from_id'], eid_dict[edge]['from_orn'])
            to_node = (eid_dict[edge]['to_id'], eid_dict[edge]['to_orn'])

            if eid_dict[edge]['from_id'] == remove_id:
                from_dict[tuple_to_string(from_node)] -= 1
                eid_dict[edge]['from_id'] = keep_id
                if orn == '-+' or orn == '+-':
                    eid_dict[edge]['from_orn'] = reverse_strand(eid_dict[edge]['from_orn'])
                    if not eid_dict[edge]['from_orn']:
                        to_dict[tuple_to_string(to_node)] -= 1
                        eid_dict[edge]['to_orn'] = reverse_strand(eid_dict[edge]['to_orn'])
                        to_node = (eid_dict[edge]['to_id'], eid_dict[edge]['to_orn'])
                        if to_dict.get(tuple_to_string(to_node)):
                            to_dict[tuple_to_string(to_node)] += 1
                        else:
                            to_dict[tuple_to_string(to_node)] = 1
                from_node = (eid_dict[edge]['from_id'], eid_dict[edge]['from_orn'])
                if from_dict.get(tuple_to_string(from_node)):
                    from_dict[tuple_to_string(from_node)] += 1
                else:
                    from_dict[tuple_to_string(from_node)] = 1

            if eid_dict[edge]['to_id'] == remove_id:
                to_dict[tuple_to_string(to_node)] -= 1
                eid_dict[edge]['to_id'] = keep_id
                if orn == '-+' or orn == '+-':
                    eid_dict[edge]['to_orn'] = reverse_strand(eid_dict[edge]['to_orn'])
                to_node = (eid_dict[edge]['to_id'], eid_dict[edge]['to_orn'])
                if to_dict.get(tuple_to_string(to_node)):
                    to_dict[tuple_to_string(to_node)] += 1
                else:
                    to_dict[tuple_to_string(to_node)] = 1

def compression_graph_by_edges(gfa_):
    """Create dictionaries with edges's id and with nodes involved.
    By the dictionaries determine what edge can be removed. The ends
    of segment(node) with only 1 edge IN or OUT are compatible with compression.
    Passing data to the other functions update the graph and the dictionaries.
    """
    eid_dict = dict()
    from_dict = dict()
    to_dict = dict()

    data_edges = gfa_.edge()
    for node1 in data_edges:
        for node2 in data_edges[node1]:
            for eid in data_edges[node1][node2]:
                if not eid_dict.get(eid):
                    eid_dict[eid] = {}

                    from_node = (data_edges[node1][node2][eid]['from_node'],\
                        data_edges[node1][node2][eid]['from_orn'])
                    if from_dict.get(tuple_to_string(from_node)):
                        from_dict[tuple_to_string(from_node)] += 1
                    else:
                        from_dict[tuple_to_string(from_node)] = 1

                    to_node = (data_edges[node1][node2][eid]['to_node'],\
                        data_edges[node1][node2][eid]['to_orn'])
                    if to_dict.get(tuple_to_string(to_node)):
                        to_dict[tuple_to_string(to_node)] += 1
                    else:
                        to_dict[tuple_to_string(to_node)] = 1

                    eid_dict[eid]['from_id'] = from_node[0]
                    eid_dict[eid]['from_orn'] = from_node[1]
                    eid_dict[eid]['to_id'] = to_node[0]
                    eid_dict[eid]['to_orn'] = to_node[1]

    count_edge_compacted = 0
    edge_no_text = 0

    for eid in eid_dict:
        from_node = (eid_dict[eid]['from_id'], eid_dict[eid]['from_orn'])
        to_node = (eid_dict[eid]['to_id'], eid_dict[eid]['to_orn'])
        if gfa_.node(from_node[0]) and gfa_.node(to_node[0]):
            if from_node[1] and from_dict.get(tuple_to_string(from_node)) == 1 and \
                to_dict.get(tuple_to_string(to_node)) == 1:
                inverted_from = (from_node[0], reverse_strand(from_node[1]))
                inverted_to = (to_node[0], reverse_strand(to_node[1]))
                if (not to_dict.get(tuple_to_string(inverted_from))) and \
                    (not from_dict.get(tuple_to_string(inverted_to))):
                    new_seq, overlap, orn = compact_sequence(gfa_, from_node, to_node)
                    update_graph(gfa_, from_node[0], to_node[0], new_seq, overlap, orn)
                    update_dictionaries_by_edges(eid_dict, (from_dict, to_dict), eid, orn)
                    count_edge_compacted += 1
        else:
            edge_no_text += 1

    GRAPH_LOGGER.debug('%s edges has been compacted of total amount of %s', \
    	count_edge_compacted, len(eid_dict)-edge_no_text)
