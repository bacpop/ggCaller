"""
Algorithms to find biconnected components in the graph considering
only dovetails overlaps.

Adapted using the networkx biconnected module:
networkx/networkx/algorithms/components/biconnected.py
"""

#    Copyright (C) 2011-2013 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.

def _dovetails_biconnected_dfs(gfa_, components=True):
    visited = set()
    for start in gfa_.nodes():
        if start in visited:
            continue
        discovery = {start:0} # "time" of first discovery of node during search
        low = {start:0}
        root_children = 0
        visited.add(start)
        edge_stack = []
        # stack = [(start, start, iter(G[start]))] # networkx (gfa_ should be G)
        stack = [(start, start, iter(gfa_.dovetails_neighbors(start)))] # PyGFA
        while stack:
            grandparent, parent, children = stack[-1]
            try:
                child = next(children)
                if grandparent == child:
                    continue
                if child in visited:
                    if discovery[child] <= discovery[parent]: # back edge
                        low[parent] = min(low[parent],discovery[child])
                        if components:
                            edge_stack.append((parent,child))
                else:
                    low[child] = discovery[child] = len(discovery)
                    visited.add(child)
                    # stack.append((parent, child, iter(G[child]))) # networkx
                    stack.append((parent, child, iter(gfa_.dovetails_neighbors(child)))) # PyGFA
                    if components:
                        edge_stack.append((parent,child))
            except StopIteration:
                stack.pop()
                if len(stack) > 1:
                    if low[parent] >= discovery[grandparent]:
                        if components:
                            ind = edge_stack.index((grandparent,parent))
                            yield edge_stack[ind:]
                            edge_stack=edge_stack[:ind]
                        else:
                            yield grandparent
                    low[grandparent] = min(low[parent], low[grandparent])
                elif stack: # length 1 so grandparent is root
                    root_children += 1
                    if components:
                        ind = edge_stack.index((grandparent,parent))
                        yield edge_stack[ind:]
        if not components:
            # root node is articulation point if it has more than 1 child
            if root_children > 1:
                yield start

def dovetails_articulation_points(gfa_):
    """Redefinition of articulation point
    for dovetails connected components.

    An articulation point or cut vertex is any node whose removal
    (along with all its incident edges) increases the number ofconnected
    components of a graph. An undirected connected graph without
    articulation points is biconnected. Articulation points belong to
    more than one biconnected component of a graph.
    """
    return _dovetails_biconnected_dfs(gfa_, components=False)
