"""
Iterators used by the GFA graph.
This iterators work considering only edges
representing dovetails overlaps.
"""

# Note that "self" will be referred to the GFA class, keep this in mind
# while reading this class.

from networkx.exception import NetworkXError
from pygfa.algorithms.traversal import dfs_edges

class DovetailIterator:
    def dovetails_iter(self, nbunch=None, keys=False, data=False):
        """Return an iterator on edges that describe
        dovetail overlaps with the given node.

        :notes:
            It seems that networkx edges_iter keeps track of
            edges already seen, so the edge (u,v) is in the
            results but edge (v,u) is not.
        """
        for from_node, to_node, key, edge_ in self.edges_iter(nbunch, keys=True, data=True):
            if 'is_dovetail' in edge_ and edge_['is_dovetail']:
                if data is True:
                    yield (from_node, to_node, key, edge_) if keys \
                      else (from_node, to_node, edge_)
                else:
                    yield (from_node, to_node, key) if keys \
                      else (from_node, to_node)

    def dovetails_nbunch_iter(self, nbunch=None):
        """Return an iterator checking that the given nbunch nodes
        are in the graphs.
        Consider only nodes involved into a dovetail overlap.
        """
        dovetails_nodes = set()
        for from_, to_ in self.dovetails_iter():
            dovetails_nodes.add(from_)
            dovetails_nodes.add(to_)

        if nbunch is None: # include all nodes
            bunch = iter(dovetails_nodes)
        elif nbunch in dovetails_nodes: #if nbunch is a single node
            bunch = iter([nbunch])
        else:
            def bunch_iter(nlist, adj):
                try:
                    for n in nlist:
                        if n in adj:
                            yield n
                # it seems impossible to raise an exception
                # here, since dovetails_iter takes care of
                # nodes checking...
                # It's a possible dead code?
                except TypeError as e:
                    message = e.args[0]
                    # capture error for non-sequence/iterator nbunch.
                    if 'iter' in message:
                        raise NetworkXError( \
                            "nbunch is not a node or a sequence of nodes.")
                    # capture error for unhashable node.
                    elif 'hashable' in message:
                        raise NetworkXError( \
                            "Node %s in the sequence nbunch is not a valid node."%n)
                    else:
                        raise
            bunch = bunch_iter(nbunch, dovetails_nodes)
        return bunch

    def dovetails_neighbors_iter(self, nbunch=None, keys=False, data=False):
        """Return an iterator over neighbors nodes considering
        all nodes in nbunch as source node.

        :notes:
            This method is used to check right and left links among
            sequences, so from_node is needed.
            If only to_node in neighborhood are need, consider using
            `dovetails_neighbors'.
        """
        for from_node, to_node, key, edge_ in self.dovetails_iter(nbunch, keys=True, data=True):
            if data is True:
                yield (from_node, to_node, key, edge_) if keys \
                  else (from_node, to_node, edge_)
            else:
                yield (from_node, to_node, key) if keys \
                  else (from_node, to_node)

    def dovetails_neighbors(self, nbunch=None):
        """Return a list of all the right and left segments
        of the given nodes.
        """
        return list(to_ for from_, to_ in self.dovetails_neighbors_iter(nbunch))

    def right_end_iter(self, nbunch, keys=False, data=False):
        """Return an iterator over dovetail edges where
        nodes id  right-segment end is taken into account
        in the overlap
        """
        try:
            if nbunch is None:
                nids = set(self.nodes())
            elif isinstance(nbunch, str):
                raise TypeError
            else:
                nids = set(nbunch)
        except TypeError:
            nids = set()
            nids.add(nbunch)
        for nid in nids:
            for from_node, to_node, key, edge_ in self.dovetails_neighbors_iter(nid, keys=True, data=True):
                if nid == edge_["from_node"] \
                  and edge_["from_segment_end"] == "R":
                    if data is True:
                        yield (from_node, to_node, key, edge_) if keys \
                          else (from_node, to_node, edge_)
                    else:
                        yield (from_node, to_node, key) if keys \
                          else (from_node, to_node)
                if nid == edge_["to_node"] \
                  and edge_["to_segment_end"] == "R":
                    if data is True:
                        yield (from_node, to_node, key, edge_) if keys \
                          else (from_node, to_node, edge_)
                    else:
                        yield (from_node, to_node, key) if keys \
                          else (from_node, to_node)

    def right(self, nbunch=None):
        """Return all the nodes connected to the right
        end of the given node sequence.
        """
        return list(to_ for from_, to_ in self.right_end_iter(nbunch))

    def right_degree_iter(self, nbunch=None):
        return ((x, len(self.right(x))) for x in self._graph.nbunch_iter(nbunch))

    def right_degree(self, nbunch=None):
        if nbunch in self:      # return a single node
            return next(self.right_degree_iter(nbunch))[1]
        else:           # return a dict
            return dict(self.right_degree_iter(nbunch))

    def left_end_iter(self, nbunch=None, keys=False, data=False):
        """Return an iterator over dovetail edges where
        left segment-end of the nodes ids given are taken into account
        in the overlap
        """
        try:
            if nbunch is None:
                nids = set(self.nodes())
            elif isinstance(nbunch, str):
                raise TypeError
            else:
                nids = set(nbunch)
        except TypeError:
            nids = set()
            nids.add(nbunch)

        for nid in nids:
            for from_node, to_node, key, edge_ in self.dovetails_neighbors_iter(nid, keys=True, data=True):

                if nid == edge_["from_node"] \
                  and edge_["from_segment_end"] == "L":
                    if data is True:
                        yield (from_node, to_node, key, edge_) if keys \
                          else (from_node, to_node, edge_)
                    else:
                        yield (from_node, to_node, key) if keys \
                          else (from_node, to_node)

                if nid == edge_["to_node"] \
                  and edge_["to_segment_end"] == "L":
                    if data is True:
                        yield (from_node, to_node, key, edge_) if keys \
                          else (from_node, to_node, edge_)
                    else:
                        yield (from_node, to_node, key) if keys \
                          else (from_node, to_node)

    def left(self, nbunch=None):
        """Return all the nodes connected to the left
        end of the given node sequence.
        """
        return list(to_ for from_, to_ in self.left_end_iter(nbunch))

    def left_degree_iter(self, nbunch=None):
        return ((x, len(self.left(x))) for x in self._graph.nbunch_iter(nbunch))

    def left_degree(self, nbunch=None):
        if nbunch in self:      # return a single node
            return next(self.left_degree_iter(nbunch))[1]
        else:           # return a dict
            return dict(self.left_degree_iter(nbunch))

    def dovetails_linear_path_traverse_nodes_iter(self, source):
        """Traverse all nodes adjacent to source node where
        the right degree and left degree of each node
        is 1.

        :param source: One of the node in the linear path.
            It doesn't matter if it's one of the end of
            the linear path.

        :notes:
            If the source node it's not one of the end node of the path,
            the result is not an iterator over the ordered node
            of the linear path, but an iterator where nodes are returned
            by their distance from the source node.

            The code is the same as networkx _plain_bfs.
        """
        if source not in self:
            return ()
        seen = set()
        nextlevel = {source}
        while nextlevel:
            thislevel = nextlevel
            nextlevel = set()
            for v in thislevel:
                if v not in seen \
                  and self.right_degree(v) <= 1 \
                  and self.left_degree(v) <= 1:
                    yield v
                    seen.add(v)
                    nextlevel.update(self.right(v))
                    nextlevel.update(self.left(v))

    def dovetails_linear_path_traverse_edges_iter(self, source, keys=False):
        """Traverse all nodes adjacent to source node where
        the right degree and left degree of each node
        is 1.

        :param source: One of the node in the linear path.
            It doesn't matter if it's one of the end of
            the linear path.

        :notes:
            If the source node it's not one of the end node of the path,
            the result is not an iterator over the ordered node
            of the linear path, but an iterator where nodes are returned
            by their distance from the source node.

            If the source node is an isolated node, then this method
            returns an empty list (no edge is found). Use
            dovetails_linear_path_traverse_nodes_iter instead.

            Same code as _plain_bfs_dovetails_with_edges function.
        """
        is_linear = lambda gfa_, node_: gfa_.right_degree(node_) <= 1 \
          and gfa_.left_degree(node_) <= 1

        if not is_linear(self, source):
            return iter([])

        seen = set()
        seen.add(source)
        queue = list()
        queue.append(source)
        while len(queue):
            from_node = queue.pop()
            for from_, to_, key in self.dovetails_neighbors_iter(from_node, \
                                                            keys=True):
                if to_ not in seen:
                  if is_linear(self, to_):
                    yield (from_, to_, key) if keys \
                      else (from_, to_)
                    seen.add(to_)
                    queue.append(to_)

    def dovetails_linear_path_iter(self, source, keys=False):
        """Return an iterator over the linear path
        whose source node belongs to, starting from one end
        of the path to another.

        :param source: One of the node in the linear path.
        """
        path_nodes = list(self.dovetails_linear_path_traverse_nodes_iter(source))
        if path_nodes == []:
            return iter([])
        # find one of the end node of the path
        index = 0
        while index < len(path_nodes):
            node_ = path_nodes[index]
            for adj in self.dovetails_neighbors(node_):
                if adj not in path_nodes: # found one of the path end node
                    return self.dovetails_linear_path_traverse_edges_iter(node_, keys=keys)
            index += 1
        # here we are in the situation where a path has been found
        # (we checked for path_nodes == [] before), but
        # we didn't find a path extreme... so
        # we have a circular path!
        # Just take the source as starting point.
        return dfs_edges(self, self.dovetails_iter, source, keys=keys)
        
        
