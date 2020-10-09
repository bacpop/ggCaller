"""
GFA representation through a networkx MulitGraph.

The dovetail operations are available thanks to
the dovetail_operation.Iterator class, that considers only
dovetail overlaps edges.

:TODO:
    * Rewrite pprint method.
"""
import logging
import copy
import re
import os
import warnings

import networkx as nx
from networkx.classes.function import all_neighbors as nx_all_neighbors

from ggCaller.pygfa.graph_element.parser import header, segment, link, containment, path
from ggCaller.pygfa.graph_element.parser import edge, gap, fragment, group, line
from ggCaller.pygfa.graph_element import node, edge as ge, subgraph as sg
from ggCaller.pygfa.serializer import gfa1_serializer as gs1, gfa2_serializer as gs2

from ggCaller.pygfa.dovetail_operations.iterator import DovetailIterator

from ggCaller.pygfa.graph_operations.compression import compression_graph_by_nodes, compression_graph_by_edges
from ggCaller.pygfa.graph_operations.overlap_consistency import check_overlap
#from benchmark.extract_subgraph import extract_subgraph

GRAPH_LOGGER = logging.getLogger(__name__)

class InvalidSearchParameters(Exception):
    pass
class InvalidElementError(Exception):
    pass
class GFAError(Exception):
    pass


class Element:
    """Represent the types of graph a GFA graph object can have.
    """
    NODE = 0
    EDGE = 1
    SUBGRAPH = 2

def _index(obj, other):
    """Given an object O and a list
    of objects L check that exist an object O'
    in the list such that O == O'.
    
    :return True: If O' exists.
    :return: The position of O' in the list.
    """
    found = False
    index = 0
    max_len = len(other)
    while not found and index < max_len:
        if obj == other[index]:
            found = True
        else:
            index += 1
    return found, index

class GFA(DovetailIterator):
    """GFA will use a networkx MultiGraph as structure to contain
    the elements of the specification.
    GFA graphs directly accept only instances coming from the
    graph_elements package, but can contains any kind of data
    undirectly by accessing the `_graph` attribute.
    """

    def __init__(self, base_graph=None, is_rGFA = None):
        """Creates a GFA graph.

        If :param base_graph: is not `None` use the graph provided
        as base for the graph.

        A `virtual id` is assigned to edges(graph edges) that don't
        have an id.
        Their id will be `virtual_#` where `#` will be given
        by `next_virtual_id`.

        :param base graph: An instance of a networkx.MultiGraph.
        """
        if base_graph != None and not isinstance(base_graph, nx.MultiGraph):
            raise GFAError("{0} ".format(type(base_graph)) \
                            + "cannot be used as base " \
                            + "graph, "\
                            + "use networkx.MultiGraph instead.")
        self._graph = nx.MultiGraph(base_graph)
        self._subgraphs = {}
        self._next_virtual_id = 0 if base_graph is None else \
                                self._find_max_virtual_id()
        self._is_rGFA = is_rGFA

    def __contains__(self, id_):
        try:
            if self._graph.has_node(id_) :
                return True
            edge_keys = (key for from_node in self._graph.adj \
                            for to_node in self._graph.adj[from_node] \
                                for key in self._graph.adj[from_node][to_node])
            if id_ in edge_keys:
                return True
            if id_ in self._subgraphs:
                return True
            return False
        except TypeError:
            return False


    def clear(self):
        """Clear all GFA object elements.

        Call networkx `clear` method, reset the virtual id counter and
        delete all the subgraphs.
        """
        self._graph.clear()
        self._next_virtual_id = 0
        self._subgraphs = {}


    def _get_virtual_id(self, increment=True):
        """Return the next virtual id value available.

        :param increment: If set to False, the virtual id is not
            incremented. Useful mainly in interactive mode.
        """
        key = self._next_virtual_id
        if increment:
            self._next_virtual_id += 1
        return key


    def _find_max_virtual_id(self):
        """Traverse the graph to find the greatest virtual id value.
        """
        # nodes cannot have a virtual_id, so don't search inside them
        virtual_rxp = "^virtual_([\d]+)$"
        regexp = re.compile(virtual_rxp)
        virtual_keys = [0]

        for from_, to_, key in self.edges_iter(keys=True):
            match = regexp.fullmatch(key)
            if match:
                virtual_keys.append(int(match.group(1)))

        for key, data in self._subgraphs.items():
            match = regexp.fullmatch(key)
            if match:
                virtual_keys.append(int(match.group(1)))

        return max(virtual_keys)


    def nodes(self, data=False, with_sequence=False, identifier=None):
        """Return a list of the nodes in the graph.

        :param with_sequence: If set return only nodes with
            a `sequence` property.
        """
        #return list(self.nodes_iter(data=data, with_sequence=with_sequence))
        if not identifier is None :
            if self._graph.has_node(identifier):
                return self._graph.nodes(data=data)[identifier]
            else:
                return;

        if with_sequence is True:
            return list(self.nodes_iter(data=data, with_sequence=with_sequence))

        return self._graph.nodes(data=data)

    def edges(self,identifier = None, adj_dict = False, **kwargs):
        """Return all the edges in the graph."""
        #return list(self._graph.edges(**kwargs))
        
        if not identifier is None:
            if isinstance(identifier, tuple):
                return self._search_edge_by_nodes(identifier)
            else:
                return self._search_edge_by_key(identifier)

        if adj_dict is True:
            return self._graph.adj

        return self._graph.edges(**kwargs)

    def subgraphs(self, identifier=None):
        """An interface to access to the subgraphs inside
        the GFA object.

        If `identifier` is `None` all the graph Subgraph objects are
        returned.
        """
        if identifier is None:
            return self._subgraphs
        else:
            if identifier in self._subgraphs:
                return self._subgraphs[identifier]

    def _search_edge_by_key(self, edge_key):
        from_node, to_node = self._get_edge_end_nodes(edge_key)
        if (from_node, to_node) != (None, None):
            return self._graph.get_edge_data(from_node, to_node, edge_key)
        return None


    def _search_edge_by_nodes(self, nodes):
        """Search for edge and edges providing end nodes.

        If given a tuple with from_node and to_node return all the edges
        between the two nodes.

        If a third element is present in the tuple return the exact edge
        between the two nodes with the key specified by the third element.
        If no match is found return `None`.

        :returns list of dictionary: If `nodes` is a two element tuple.
        :returns dictionary: Otherwise.
        """
        if len(nodes) < 2:
            raise InvalidSearchParameters("At least two values are required.")
        from_node = nodes[0]
        to_node = nodes[1]
        try:
            if len(nodes) > 2:
                key = nodes[2]
                return self._graph.get_edge_data(from_node, to_node, key)
            return self._graph.get_edge_data(from_node, to_node)
        except:
            return None


    def _get_edge_end_nodes(self, edge_key):
        """Given an edge key return a tuple that contains
        the end nodes for that edge.
        """
        for from_node, to_node, key in self.edges_iter(keys=True):
            if key == edge_key:
                return from_node, to_node
        return None, None


    def get(self, key):
        """Return the element pointed by the specified key."""
        if self._graph.has_node(key):
            return self.nodes(data = True, identifier = key)
        if key in self._subgraphs:
            return self._subgraphs[key]
        edge_ = self._search_edge_by_key(key)
        if not edge_ is None:
            return edge_


    def as_graph_element(self, key):
        """Given a key of an existing node, edge or subgraph, return
        its equivalent graph element object.
        """
        element = self.get(key)
        if element is None:
            raise InvalidElementError(\
                    "No graph element has the given key: {0}".format(key))

        # Subgraph objects don't need to be converted
        if sg.is_subgraph(element):
            return copy.deepcopy(element)

        tmp_list = copy.deepcopy(element)
        try:
            if 'nid' in element:
                tmp_list.pop('nid')
                tmp_list.pop('sequence')
                tmp_list.pop('slen')
                return node.Node(\
                                element['nid'], \
                                element['sequence'], \
                                element['slen'], \
                                opt_fields=tmp_list)
            if 'eid' in element:
                tmp_list.pop('eid')
                tmp_list.pop('from_node')
                tmp_list.pop('from_orn')
                tmp_list.pop('to_node')
                tmp_list.pop('to_orn')
                tmp_list.pop('from_positions')
                tmp_list.pop('to_positions')
                tmp_list.pop('alignment')
                tmp_list.pop('variance')
                tmp_list.pop('distance')
                edge_ = ge.Edge(\
                                element['eid'], \
                                element['from_node'], element['from_orn'], \
                                element['to_node'], element['to_orn'], \
                                element['from_positions'], \
                                element['to_positions'], \
                                element['alignment'], element['distance'], \
                                element['variance'], \
                                opt_fields=tmp_list, \
                                is_dovetail=element['is_dovetail'] \
                                )
                return edge_
        except KeyError:
            return None


    def add_graph_element(self, element):
        """Add a graph element -Node, Edge or Subgraph- object to
        the graph."""
        if isinstance(element, node.Node):
            self.add_node(element)
        elif isinstance(element, ge.Edge):
            self.add_edge(element)
        elif isinstance(element, sg.Subgraph):
            self.add_subgraph(element)


    def add_node(self, new_node, safe=False):
        """Add a graph_element Node to the GFA graph
        using the node id as key.

        Its sequence and sequence length will be individual attributes
        on the graph and all the remainders optional field will be stored
        individually as node data.

        :param new_node: A graph_element.Node object or a string
            that can represent a node (such as the Segment line).
        :param safe: If set check if the given identifier has already
            been added to the graph, and in that case raise
            an exception
        """
        if isinstance(new_node, str) and new_node[0] == "S":
            if segment.is_segmentv1(new_node):
                new_node = node.Node.from_line(\
                    segment.SegmentV1.from_string(new_node.strip()))
            else:
                new_node = node.Node.from_line(\
                    segment.SegmentV2.from_string(new_node.strip()))

        if not node.is_node(new_node):
            raise node.InvalidNodeError("The object given is not a node.")

        if safe and new_node.nid in self:
            raise GFAError("An element with the same id already exists.")

        if self._is_rGFA == True:
            if not self.check_rGFA_node(new_node):
                raise node.InvalidNodeError("{0}".format(new_node.nid)+\
                                                " cannot be a rGFA node")

        self._graph.add_node(\
                              new_node.nid, \
                              nid=new_node.nid, \
                              sequence=new_node.sequence, \
                              slen=new_node.slen, \
                              **new_node.opt_fields)
        return True


    def remove_node(self, nid):
        """Remove a node with nid as its node id.

        Edges containing nid as end node will be automatically
        deleted.

        :param nid: The id belonging to the node to delete.
        :raise InvalidNodeError: If `nid` doesn't point to any node.
        """
        try:
            self._graph.remove_node(nid)
        except:
            raise node.InvalidNodeError("{0} doesn't point".format(nid) \
                                        + " to any node in the graph.")


    def nodes_iter(self, data=False, with_sequence=False):
        """Return an iterator over nodes in the graph.

        :para with_sequence: If set return only nodes with
            a sequence property.
        """
        if with_sequence is True:
            if data is True:
                return iter((nid, data_) \
                         for nid, data_ in self._graph.nodes(data = True) \
                            if 'sequence' in data_)
            else:
                return iter(nid \
                         for nid, data_ in self._graph.nodes(data = True) \
                            if 'sequence' in data_)
        else:
            return iter(list(self._graph.nodes(data = data)))


    def nbunch_iter(self, nbunch=None):
        """Return an iterator of nodes contained in nbunch that are
        also in the graph.

        Interface to the networkx method.
        """
        return self._graph.nbunch_iter(nbunch=nbunch)


    def add_edge(self, new_edge, safe=False):
        """Add a graph_element Edge or a networkx edge to the GFA
        graph using  the edge id as key.

        If its id is `*` or `None` the edge will be given a
        **virtual_id**, in either case the original edge id will
        be preserved as edge attribute.

        All edge attributes will be stored as netwrorkx edge
        attributes and  all the remainders optional field will be stored
        individually as edge data.
        """
        if isinstance(new_edge, str):
            if new_edge[0] == 'L':
                new_edge = ge.Edge.from_line(\
                    link.Link.from_string(new_edge.strip()))
            elif new_edge[0] == 'C':
                new_edge = ge.Edge.from_line(\
                    containment.Containment.from_string(new_edge.strip()))
            elif new_edge[0] == 'E':
                new_edge = ge.Edge.from_line(\
                    edge.Edge.from_string(new_edge.strip()))
            elif new_edge[0] == 'G':
                new_edge = ge.Edge.from_line(\
                    gap.Gap.from_string(new_edge.strip()))
            elif new_edge[0] == 'F':
                new_edge = ge.Edge.from_line(\
                    fragment.Fragment.from_string(new_edge.strip()))
            else:
                raise ge.InvalidEdgeError(\
                    "The string given doesn't represent a GFA line that could" \
                    + " be represented as an edge,\n" \
                    + "given: {0}".format(new_edge))

        if not ge.is_edge(new_edge):
            raise ge.InvalidEdgeError("The object is not a valid edge.")

        key = new_edge.eid
        if new_edge.eid is None or new_edge.eid == '*':
            key = "virtual_{0}".format(self._get_virtual_id())

        if safe:
            edge_exists = key in self
            node1_exists = new_edge.from_node in self
            node2_exists = new_edge.to_node in self
            if edge_exists:
                raise GFAError("An element with the same id already exists.")
            if not (node1_exists and \
                       node2_exists):
                raise GFAError("From/To node are not already in the graph.")

        if self._is_rGFA:

            if not self.check_rGFA_edge(new_edge):
                #tmp_field_type,tmp_SR_max = self.max_SR(new_edge.from_node,\
                #                                            new_edge.to_node)
                #tmp_opt_field = line.OptField("SR", tmp_SR_max, tmp_field_type)
                tmp_fields = re.split(":", self.max_SR(new_edge.from_node,new_edge.to_node))[1:]
                tmp_opt_field = line.OptField("SR",tmp_fields[1],tmp_fields[0])
                new_edge.opt_fields.update({"SR" : tmp_opt_field})
                warnings.warn("{0} + {1}".format(new_edge.from_node,new_edge.to_node) +\
                    ": SR has been set to the maximum value between the SR of the adjacent nodes")

        self._graph.add_edge( \
                               new_edge.from_node, new_edge.to_node, key=key, \
                               eid=new_edge.eid, \
                               from_node=new_edge.from_node, \
                               from_orn=new_edge.from_orn, \
                               to_node=new_edge.to_node, \
                               to_orn=new_edge.to_orn, \
                               from_positions=new_edge.from_positions, \
                               to_positions=new_edge.to_positions, \
                               alignment=new_edge.alignment, \
                               distance=new_edge.distance, \
                               variance=new_edge.variance, \
                               is_dovetail=new_edge.is_dovetail, \
                               from_segment_end=new_edge.from_segment_end, \
                               to_segment_end=new_edge.to_segment_end, \
                               **new_edge.opt_fields \
                               )

    def max_SR(self, from_node, to_node):
        #return SR is the greater between the input nodes
        # SR:field_type:value
        sr_1 = re.split(":", \
            str(self.as_graph_element(from_node).opt_fields['SR']))[2]
        sr_2 = re.split(":", \
            str(self.as_graph_element(to_node).opt_fields['SR']))[2]
        tmp_max = max(sr_1,sr_2)
        tmp_field_type = re.split(":", str(self.as_graph_element(from_node).opt_fields['SR']))[1]

        #return tmp_field_type, str(tmp_max)
   
        return str.join(":", ("SR", tmp_field_type, str(tmp_max)))
    
    def check_rGFA(self, force = False):
        #returns true if the graph is rGFA, otherwise False
        #if self._is_rGFA == None or force == True:
        if force == True:
            if self.check_rGFA_nodes(self.nodes())\
                 and self.check_rGFA_edges(self.edges()):

                return True

            return False

        return self._is_rGFA


    #def check_rGFA(self, force = False):

        #if self._is_rGFA == None or force == True:
    #    if force == True:
    #        if self.check_rGFA_nodes(self.node())\
    #             and self.check_rGFA_edges(self.edges()):
                
    #            self._is_rGFA = True

    #        else:
    #            self._is_rGFA = False
        
    #    return self._is_rGFA


    def check_rGFA_nodes(self, nodes):
        #check that the nodes are suitable for an rGFA graph
        #if isinstance(nodes,list)
        #else isinstance(nodes, dict or view)
        if isinstance(nodes,list):
            for node in nodes:
                if  self.check_rGFA_node(self.nodes(identifier = node)):
                    continue
                else:

                    return False
        else:
            for node in nodes:
                if self.check_rGFA_node(nodes[node]):
                    continue
                else:

                    return False

        return True

    def check_rGFA_node(self,node):
        if isinstance(node, dict):
            if "SN" in node and \
                "SO" in node and \
                "SR" in node:

                return True
        else:
            if "SN" in node.opt_fields and \
                "SO" in node.opt_fields and \
                "SR" in node.opt_fields:

                return True
        
        return False

    def check_rGFA_edges(self, edges):
        #check that the edges are suitable for an rGFA graph
        for edge in edges:
            if not self.check_rGFA_edge(edge):
          
                return False
        
        return True

    def check_rGFA_edge(self, edge):

        if isinstance(edge, tuple):
            tmp_libr = self._search_edge_by_nodes(edge)
            if "SR" in tmp_libr[next(iter(tmp_libr))]:
                
                return True
        
        else:
            if "SR" in edge.opt_fields:

                return True

        return False 



    def remove_edge(self, identifier):
        """Remove an edge or all edges identified by an id
        or by a tuple with end node, respectively.

        * If `identifier` is a two elements tuple remove all the
            all the edges between the two nodes.

        * If `identifier` is a three elements tuple remove the edge
            specified by the third element of the tuple with end nodes
            given by the first two elements of the tuple itself.

        * If `identifier` is not a tuple, treat it as it should be
            an edge id.

        :raise InvalidEdgeError: If `identifier` is not in the cases
            described above.
        """
        try:
            if isinstance(identifier, tuple):
                if len(identifier) == 2:
                    self.remove_edges(identifier[0], identifier[1])
                else:
                    self._graph.remove_edge(identifier[0], \
                                            identifier[1], \
                                            identifier[2])
            else:
                from_node, to_node = self._get_edge_end_nodes(identifier)
                self._graph.remove_edge(from_node, \
                                        to_node, \
                                        identifier)
        except nx.NetworkXError as nxe:
            raise ge.InvalidEdgeError(nxe)


    def remove_edges(self, from_node, to_node):
        """Remove all the direct edges between the two nodes given.

        Call iteratively remove_edge (remove a not specified edge
        from `from_node` and `to_node`) for n-times where n is
        the number of edges between the given nodes,
        removing all the edges indeed.
        """
        num_edges = len(self.edges(identifier = (from_node, to_node)))
        for edge_ in range(0, num_edges):
            self._graph.remove_edge(from_node, to_node)


    def edges_iter(self, nbunch=None, data=False, keys=False, default=None):
        """Interface to networx edges iterator."""
        return iter(self._graph.edges(nbunch=nbunch,\
            data=data, \
            keys=keys, \
            default=default))


    def add_subgraph(self, subgraph, safe=False):
        """Add a Subgraph object to the graph.

        The object is not altered in any way.
        A deepcopy of the object given is attached to the graph.
        """
        if isinstance(subgraph, str):
            if subgraph[0] == "P":
                subgraph = sg.Subgraph.from_line(\
                    path.Path.from_string(subgraph))
            elif subgraph[0] == "O":
                subgraph = sg.Subgraph.from_line(\
                    group.OGroup.from_string(subgraph))
            elif subgraph[0] == "U":
                subgraph = sg.Subgraph.from_line(\
                    group.UGroup.from_string(subgraph))
            else:
                raise sg.InvalidSubgraphError(\
                    "The string given cannot be represented as a subgraph,\n" \
                    + "given: {0}".format(subgraph))
        if not sg.is_subgraph(subgraph):
            raise sg.InvalidSubgraphError("The object given is not a subgraph.")

        key = subgraph.sub_id
        if key == '*':
            key = "virtual_{0}".format(self._get_virtual_id())
        if safe and key in self:
            raise GFAError("An element with the same id already exists.")
        self._subgraphs[key] = copy.deepcopy(subgraph)


    def remove_subgraph(self, subgraph_id):
        """Remove the Subgraph object identified by the given id.
        """
        try:
            del(self._subgraphs[subgraph_id])
        except:
            raise sg.InvalidSubgraphError("The given id doesn't " \
                                         + " identify any subgraph.")


    def subgraphs_iter(self, data=False):
        """Return an iterator over subgraphs elements
        in the GFA graph.
        """
        if data is True:
            return iter(self._subgraphs.items())
        else:
            return iter(self._subgraphs)

    def get_subgraph(self, sub_key):
        """Return a GFA subgraph from the parent graph.

        Return a new GFA graph structure with the nodes,
        edges and subgraphs specified in the elements attributes
        of the subgraph object pointed by the id.

        The returned GFA is *independent* from the original object.

        :param sub_key: The id of a subgraph present in the GFA graph.
        :returns None: if the subgraph id doesn't exist.
        """
        if not sub_key in self._subgraphs:
            raise sg.InvalidSubgraphError(\
                "There is no subgraph pointed by this key.")
        subgraph = self._subgraphs[sub_key]
        sub_gfa = GFA()
        for id_, orn in subgraph.elements.items():
            # creating a new GFA graph and the add method,
            # the virtual id are recomputed
            sub_gfa.add_graph_element(self.as_graph_element(id_))
        return sub_gfa


    def subgraph(self, nbunch, copy=True):
        """Given a bunch of nodes return a graph with
        all the given nodes and the edges between them.

        The returne object is not a GFA Graph, but a
        MultiGraph. To create a new GFA graph, just
        use the GFA initializer an give the subgraph to it.

        Interface to the networkx subgraph method.
        Given a collection of nodes return a subgraph with the nodes
        given and all the edges between each pair of nodes.

        :param nbunch: The nodes.
        :param copy: If set to True return a copy of the subgraph.
        """
        subgraph_ = self._graph.subgraph(nbunch)
        if copy:
            return subgraph_.copy()
        return subgraph_


    def dovetails_subgraph(self, nbunch=None, copy=True):
        """Given a collection of nodes return a subgraph with the nodes
        given and all the edges between each pair of nodes.
        Only dovetails overlaps are considered.
        """
        bunch = list(self.nbunch_iter(nbunch))
        # create new graph and copy subgraph into it
        H = self._graph.__class__()
        # add node and attribute dictionaries
        H.add_nodes_from((n, self._graph.nodes[n]) for n in bunch)

        # add edge and attribute dictionaries

        """for n, nbrs in self._graph.adj.items():
            if n in bunch:
                for nbr, keydict in nbrs.items():
                    if nbr in bunch:
                        for key, d in keydict.items():
                            if d['is_dovetail'] is True:
                                H.add_edges_from([(n, nbr, key, d)])
        """ 
        H.add_edges_from((n, nbr, key, d) 
            for n, nbrs in self._graph.adj.items() if n in bunch 
            for nbr, keydict in nbrs.items() if nbr in bunch 
            for key, d in keydict.items() if d['is_dovetail'] is True)

        H.graph = self._graph.graph
        if copy is True:
            return H.copy()
        return H



    """def dovetails_subgraph(self, nbunch=None, copy=True):
        Given a collection of nodes return a subgraph with the nodes
        given and all the edges between each pair of nodes.
        Only dovetails overlaps are considered.

        bunch = list(self.nbunch_iter(nbunch))
        # create new graph and copy subgraph into it
        H = self._graph.__class__()
        # add node and attribute dictionaries
        H.add_nodes_from((n, self._graph.nodes[n]) for n in bunch)

        #BRUTTO
        H_adj = dict(H.adj)
        for tmp in H_adj:
            H_adj[tmp] = dict(H_adj[tmp])

        # filter edges based on is_dovetail property
        self_adj = self._graph.adjlist_inner_dict_factory()
        for from_node in self._graph.adj:
            self_adj[from_node] = self._graph.adjlist_inner_dict_factory()
            for to_node in self._graph.adj[from_node]:
                self_adj[from_node][to_node] = self._graph.adjlist_inner_dict_factory()
                for edge_ in self._graph.adj[from_node][to_node]:
                    if self._graph.adj[from_node][to_node][edge_]['is_dovetail'] is True:
                        self_adj[from_node][to_node][edge_] = \
                          self._graph.adj[from_node][to_node][edge_]

        # add nodes and edges (undirected method)
        for n in H:
            Hnbrs = H.adjlist_inner_dict_factory()
            H_adj[n] = Hnbrs
            for nbr, edgedict in self_adj[n].items():
                if nbr in H_adj:
                    # add both representations of edge: n-nbr and nbr-n
                    # they share the same edgedict
                    ed = edgedict.copy()
                    Hnbrs[nbr] = ed
                    H_adj[nbr][n] = ed

        H.graph = self._graph.graph
        if copy is True:
            return H.copy()
        return H
    """

    def neighbors(self, nid):
        """Return all the nodes id of the nodes connected to
        the given node.

        Return all the predecessors and successors of the
        given source node.

        :params nid: The id of the selected node
        """
        if self.nodes(identifier = nid) is None:
            raise GFAError("The source node is not in the graph.")
        return list(nx_all_neighbors(self._graph, nid))

    def search(self, \
               comparator, \
               limit_type=None):
        """Perform a query applying the comparator on each graph element.
        """
        if limit_type == Element.NODE:
            return self.search_on_nodes(comparator)

        elif limit_type == Element.EDGE:
            return self.search_on_edges(comparator)

        elif limit_type == Element.SUBGRAPH:
            return self.search_on_subgraph(comparator)

        retval = []
        retval.extend(self.search_on_nodes(comparator))
        retval.extend(self.search_on_edges(comparator))
        retval.extend(self.search_on_subgraph(comparator))
        return retval


    def search_on_nodes(self, comparator):
        retval = []
        for key, data in self.nodes_iter(data=True):
            try:
                if comparator(data):
                    retval.append(key)
            except KeyError: pass
        return retval


    def search_on_edges(self, comparator):
        retval = []
        for u, v, key, data in self.edges_iter(data=True, keys=True):
            try:
                if comparator(data):
                    retval.append(key)
            except KeyError: pass
        return retval


    def search_on_subgraph(self, comparator):
        retval = []
        for key, data in self._subgraphs.items():
            data = data.as_dict()
            try:
                if comparator(data):
                    retval.append(key)
            except KeyError: pass
        return retval


    def from_string(self, string):
        """Add a GFA string to the graph once it has been
        converted.

        :TODO:
            Maybe this could be used instead of checking for line type
            in the add_xxx methods...
        """
        lines = re.split("\n", string)
        for line_ in lines:
            line_ = line_.strip()
            if len(line_) < 1:
                continue
            if line_[0] == 'S':
                if segment.is_segmentv1(line_):
                    self.add_graph_element(\
                        node.Node.from_line(\
                            segment.SegmentV1.from_string(line_)))
                else:
                    self.add_graph_element(\
                        node.Node.from_line(\
                            segment.SegmentV2.from_string(line_)))
            elif line_[0] == 'L':
                self.add_graph_element(\
                    ge.Edge.from_line(\
                        link.Link.from_string(line_)))
            elif line_[0] == 'C':
                self.add_graph_element(\
                    ge.Edge.from_line(\
                        containment.Containment.from_string(line_)))
            elif line_[0] == 'E':
                self.add_graph_element(\
                    ge.Edge.from_line(\
                        edge.Edge.from_string(line_)))
            elif line_[0] == 'G':
                self.add_graph_element(\
                    ge.Edge.from_line(\
                        gap.Gap.from_string(line_)))
            elif line_[0] == 'F':
                self.add_graph_element(\
                    ge.Edge.from_line(\
                        fragment.Fragment.from_string(line_)))
            elif line_[0] == 'P':
                self.add_graph_element(\
                    sg.Subgraph.from_line(\
                        path.Path.from_string(line_)))
            elif line_[0] == 'O':
                self.add_graph_element(\
                    sg.Subgraph.from_line(\
                        group.OGroup.from_string(line_)))
            elif line_[0] == 'U':
                self.add_graph_element(\
                    sg.Subgraph.from_line(\
                        group.UGroup.from_string(line_)))


    # This method has been checked manually
    @classmethod
    def from_file(cls, filepath, is_rGFA = None): # pragma: no cover
        """Parse the given file and return a GFA object.
        """
        pygfa_ = GFA(is_rGFA=is_rGFA)
        file_handler = open(filepath)
        file_content = file_handler.read()
        file_handler.close()
        pygfa_.from_string(file_content)
        return pygfa_


    def pprint(self): # pragma: no cover
        """A basic pretty print function for nodes and edges.
        """
        string = "\nGRAPH:\nNodes: [\n"
        for node, datas in self.nodes_iter(data=True):
            string += str(node) + "\t: {"
            for name, data in datas.items():
                string += str(name) + ": " + str(data) + "\t"
            string += "}\n"
        string += "]\n"

        string += "\nEdges: [\n"
        for from_node, to_node, key, datas in self.edges_iter( \
                                                            keys=True,\
                                                            data=True):
            string += str(key) + "\t: {"
            for name, data in datas.items():
                string += str(name) + ": " + str(data) + "\t"
            string += "}\n"
        string += "]\n"

        string += "\nSubgraphs: [\n"
        for key, data in self._subgraphs.items():
            string += str(key) + "\t: {" + str(data) +  "}\n"

        string += "]\n"
        return string


    def dump(self, gfa_version=1, out=None):
        try:
            dump_ = ""
            if gfa_version == 1:
                dump_ = gs1.serialize_gfa(self)
            elif gfa_version == 2:
                dump_ = gs2.serialize_gfa(self)
            else:
                raise ValueError("Invalid GFA output version.")
            if out is None:
                return dump_

            with open(out, 'w') as out_file:
                out_file.write(dump_)
        except EnvironmentError as env_error:
            GRAPH_LOGGER.error(repr(env_error))


    def _make_edge_lut(self):
        """Return a lookup table that associate each edge id with a best
        match unique id dependent on edge information (from_node, to_node).

        All the edges between a pair of nodes will have the same alias,
        so a single alias will collect a list of real id.
        False positive will also be place inside this lists.

        If both from_node and to_node are not specified the id
        is placed into a set for further processing.
        """
        virtual_rxp = "^virtual_([\d]+)$"
        regexp = re.compile(virtual_rxp)
        edge_lut = {}
        pure_virtuals = []
        for from_node, to_node, edge_, data_ in self.edges_iter(keys=True, data=True):
            from_data = self.nodes(identifier = from_node)
            to_data = self.nodes(identifier = to_node)
            match = regexp.fullmatch(edge_)
            if match is not None:
                from_sequence = ""
                to_sequence = ""

                if 'sequence' in from_data \
                  and from_data['sequence'] not in (None, '*'):
                    from_sequence = from_data['sequence']
                if 'sequence' in to_data and \
                  to_data['sequence'] not in (None, '*'):
                    to_sequence = to_data['sequence']

                if not from_sequence and not to_sequence:
                    pure_virtuals.append(edge_)
                else:
                    alias = from_sequence + to_sequence
                    if alias not in edge_lut:
                        edge_lut[alias] = []
                    edge_lut[alias].append(edge_)
            else:
                edge_lut[edge_] = [edge_]
        return edge_lut, pure_virtuals


    def _make_edge_table(self):
        """Create a table to list each edge id to
        its end nodes.

        A networkx Multigraph edge is identified by its end nodes and by
        a key/id. To access a specific edge one should do:
        >>> graph.edge[from_node][to_node][edge_key]

        This way, to identify an edge by only its key, a search
        operation that could take :math:`O(n)` is required.

        This function instead write end nodes and key in the opposite
        order, so that a search by key operation could be done
        in constant time, but requires :math:`O(n)` space.
        """
        edges = {}
        for from_node, to_node, key in self.edges_iter(keys=True):
            edges[key] = (from_node, to_node)
        return edges

    def _look_for_edge(self, key, edge_table):
        from_node, to_node = edge_table[key]
        return self._graph.get_edge_data(from_node, to_node, key)

    def __eq__(self, other):
        """
        :TODO:
            * make a lut for subgraphs (try to think for a way to write
              _make_edge_lut in a reusable way...
        """
        try:
            # Nodes must be defined, so there is no reason to
            # create a LUT
            for nid, node_ in self.nodes_iter(data=True):
                if node_ != other.nodes(identifier = nid):
                    return False

            self_edge_table = self._make_edge_table()
            other_edge_table = other._make_edge_table()
            self_lut, self_edge_virtuals = self._make_edge_lut()
            other_lut, other_edge_virtuals = other._make_edge_lut()
            for alias, list_ids in self_lut.items():
                while len(list_ids):
                    id_ = list_ids.pop()
                    found = False
                    index = 0
                    edge_ = self._look_for_edge(id_, self_edge_table)
                    while not found and index < len(other_lut[alias]):
                        other_id = other_lut[alias][index]
                        if edge_ == other._look_for_edge(\
                                            other_id, \
                                            other_edge_table):
                            found = True
                        else:
                            index += 1
                    if not found:
                        return False
                    # if is found remove it from list
                    # to speed up next searches.
                    other_lut[alias].pop(index)
                # if other_lut has other ids attached to that alias, then
                # graphs are not equals
                #if not len(other_lut[alias]):
                #    return False

            for edge_ in self_edge_virtuals:
                found, index = _index(edge_, other_edge_virtuals)
                if not found:
                    return False
                other_edge_virtuals.pop(index)

            # I think it's difficult to have lots of subgraphs
            # If I am wrong a subgraphs lut will be made and the comparison
            # should be nearly linear in time
            self_subgraphs = [sub.as_dict() for sub in self.subgraphs().values()]
            other_subgraphs = [sub.as_dict() for sub in other.subgraphs().values()]
            for sub_ in self_subgraphs:
                found, index = _index(sub_, other_subgraphs)
                if not found:
                    return False
                other_subgraphs.pop(index)

        except (AttributeError, KeyError) as e:
            return False
        return True

    def __neq__(self, other):
        return not self == other

    def compression(self, type_compression='by_nodes'):
        if type_compression == 'by_edges':
            compression_graph_by_edges(self)
        else:
            count_edge_compacted = compression_graph_by_nodes(self)
            while not count_edge_compacted == 0:
                count_edge_compacted = compression_graph_by_nodes(self)

    def overlap_consistency(self, external_file=None):
        FOLDER, _ = os.path.split(__file__)
        return check_overlap(self, FOLDER.rstrip('pygfa'), external_file)

#    def subgraphs_extractor(self, n_source, distance):
#        extract_subgraph(self, n_source, distance)

if __name__ == '__main__': #pragma: no cover
    pass
