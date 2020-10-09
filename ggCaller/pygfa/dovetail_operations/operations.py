from ggCaller.pygfa.dovetail_operations.components.connected import dovetails_nodes_connected_components
from ggCaller.pygfa.dovetail_operations.components.connected import dovetails_nodes_connected_component
from ggCaller.pygfa.dovetail_operations.components.connected import dovetails_connected_components_subgraphs

from ggCaller.pygfa.dovetail_operations.components.biconnected import dovetails_articulation_points
from ggCaller.pygfa.dovetail_operations.linear_paths import *
from ggCaller.pygfa.dovetail_operations.simple_paths import *

from ..operations import nodes_connected_component

def dovetails_remove_small_components(gfa_, min_length):
    """Remove all the connected components where
    the sequences length is less than min_length.

    Find all the connected components nodes,
    for each component obtain the sum of the
    sequences length.
    If length is less than the given length remove the connected
    component nodes.

    :param min_length: An integer describing the required length
        to keep a connected component.

    :note:
       When connected components are computed only dovetail overlaps
      edges are considered.
    """
    if min_length < 0:
        raise ValueError("min_length must be >= 0")
    for conn_comp in dovetails_nodes_connected_components(gfa_):
        length = 0
        for nid in conn_comp:
            node_ = gfa_.nodes(identifier = nid)
            try:
                length += node_['slen']
            except (TypeError, KeyError):
                pass
        if length < min_length:
            for nid in conn_comp:
                gfa_.remove_node(nid)

def dovetails_remove_dead_ends(\
                      gfa_, \
                      min_length, \
                      safe_remove=False):
    """Remove all the nodes where its right
    degree and its left degree are the following (0,0), (1,0), (1,0)
    and the length of the sequence is less than the given length.
    The node to remove mustn't split its connected component
    in two.

    :param min_length:
    :param consider_sequence: If set try to get the sequence length
        where length field is not defined.
    :param safe_remove: If set the operation doesn't remove nodes
        where is not possible to obtain the length value.

    :note:
        Using the right and left degree, only dovetails overlaps
        are considered.
    """
    if min_length < 0:
        raise ValueError("min_length must be >= 0")

    art_points = set(dovetails_articulation_points(gfa_))
    to_remove = set()
    for nid, node_ in gfa_.nodes_iter(data=True):
        left_deg = gfa_.left_degree(nid)
        right_deg = gfa_.right_degree(nid)
        if (left_deg, right_deg) in [(0,0), (0,1), (1,0)] \
          and nid not in art_points:
            try:
                length = node_['slen']
                if length is None:
                   length = 0
                if length < min_length:
                    to_remove.add(nid)
            except KeyError:
                if not safe_remove:
                    to_remove.add(nid)

    for nid in to_remove:
        gfa_.remove_node(nid)
