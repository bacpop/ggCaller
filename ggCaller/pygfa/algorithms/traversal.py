def dfs_edges(gfa_, selector, source=None, keys=False, **args):
    """Custom dfs_edges to select custom edges
    while traversing.

    :param keys: If set return the keys of the edges of the dfs tree."""
    if source is None:
        # produce edges for all components
        nodes = gfa_.nodes()
    else:
        # produce edges for components with source
        nodes = [source]
    visited = set()
    for start in nodes:
        if start in visited:
            continue
        visited.add(start)
        stack = [(start, iter(selector(start, keys=keys, **args)))]
        while stack:
            parent, children = stack[-1]
            try:
                child = next(children)
                node_to = child[1]
                if node_to not in visited:
                    yield child
                    visited.add(node_to) # add the to node to the visited set
                    stack.append((node_to, iter(selector(node_to, keys=keys, **args))))
            except StopIteration:
                stack.pop()
