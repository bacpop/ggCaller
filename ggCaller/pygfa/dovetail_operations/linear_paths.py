"""
Module that contain operation to find linearh paths
in a GFA graph.
"""
def dovetails_linear_path(gfa_, node_, keys=False):
    """Return the oriented edges involved in a linear path
    where that contain the given node.
    """
    return gfa_.dovetails_linear_path_iter(node_, keys=keys)

def dovetails_linear_paths(gfa_, components=False, keys=False):
    seen = set(gfa_.nodes())
    while len(seen):
        source = seen.pop() # source is removed
        nodes = set(gfa_.dovetails_linear_path_traverse_nodes_iter(source))
        lin_path = list(dovetails_linear_path(gfa_, source, keys=keys))
        if len(lin_path) > 0:
            yield lin_path
        seen -= nodes # remove nodes in linear path
