from networkx.algorithms.components.connected import node_connected_component as nx_node_connected_component
from networkx.algorithms.components.connected import connected_components as nx_connected_components

import pygfa.gfa # required for GFAError (gives error otherwise)

def nodes_connected_component(gfa_, nid):
    """Return the connected component
    belonging to the given node.

    :param nid: The id of the node to find the reachable nodes.
    """
    if nid not in gfa_:
        raise pygfa.gfa.GFAError("The source node is not in the graph.")
    return nx_node_connected_component(\
                        gfa_._graph, nid)

def nodes_connected_components(gfa_):
    """Return a generator of sets with nodes of each weakly
    connected component in the graph.
    """
    return nx_connected_components(gfa_._graph)
