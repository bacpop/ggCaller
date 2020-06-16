"""
GFA2 Serializer for nodes, edge, Subgraphs and networkx graphs.

Can serialize either one of the object from the group mentioned
before or from a dictionary with equivalent key.
"""

import copy
import logging

import networkx as nx

from pygfa.graph_element.parser import field_validator as fv
from pygfa.serializer import utils

class GFA2SerializationError(Exception): pass

serializer_logger = logging.getLogger(__name__)

DEFAULT_IDENTIFIER = "no identifier given."

SEGMENT_FIELDS = [\
    fv.GFA2_ID, \
    fv.GFA2_INT, \
    fv.GFA2_SEQUENCE]

EDGE_FIELDS = [\
    fv.GFA2_OPTIONAL_ID, \
    fv.GFA2_REFERENCE, \
    fv.GFA2_REFERENCE, \
    fv.GFA2_POSITION, \
    fv.GFA2_POSITION, \
    fv.GFA2_POSITION, \
    fv.GFA2_POSITION, \
    fv.GFA2_ALIGNMENT]

FRAGMENT_FIELDS = [\
    fv.GFA2_ID, \
    fv.GFA2_REFERENCE, \
    fv.GFA2_POSITION, \
    fv.GFA2_POSITION, \
    fv.GFA2_POSITION, \
    fv.GFA2_POSITION, \
    fv.GFA2_ALIGNMENT]

GAP_FIELDS = [\
    fv.GFA2_OPTIONAL_ID, \
    fv.GFA2_REFERENCE, \
    fv.GFA2_REFERENCE, \
    fv.GFA2_INT, \
    fv.GFA2_OPTIONAL_INT]

UGROUP_FIELDS = [fv.GFA2_OPTIONAL_ID, fv.GFA2_IDS]
OGROUP_FIELDS = [fv.GFA2_OPTIONAL_ID, fv.GFA2_REFERENCES]
################################################################################
# NODE SERIALIZER
################################################################################
def serialize_node(node_, identifier=DEFAULT_IDENTIFIER):
    """Serialize to the GFA2 specification a graph_element Node or a
    dictionary that has the same informations.

    If sequence length is undefined (for example, after parsing
    a GFA1 Sequence line) a sequence length of 0 is automatically
    added in the serialization process.

    :param node: A Graph Element Node or a dictionary
    :identifier: If set help gaining useful debug information.
    :returns "": If the object cannot be serialized to GFA.
    """
    identifier = utils._check_identifier(identifier)
    try:
        if isinstance(node_, dict):
            node_dict = copy.deepcopy(node_)
            # do not modify node_dict since it's not a copy
            node_length = node_['slen']
            if node_length is None:
                node_length = 0

            # 'slen' has been seitched to node_length, but
            # now 'slen' must be removed
            node_dict.pop('slen')
            defined_fields = [ \
                                node_dict.pop('nid'), \
                                node_length, \
                                node_dict.pop('sequence') \
                             ]
            fields = ["S"]
            fields.append(str(node_['nid']))
            fields.append(str(node_length))
            fields.append(str(node_['sequence']))
            fields.extend(utils._serialize_opt_fields(node_dict))
        else:
            # do not modify node_ since it's not a copy
            node_length = node_.slen
            if node_length is None:
                node_length = 0
            defined_fields = [ \
                               node_.nid, \
                               node_.sequence, \
                               node_length \
                             ]
            fields = ["S"]
            fields.append(str(node_.nid))
            fields.append(str(node_length))
            fields.append(str(node_.sequence))
            fields.extend(utils._serialize_opt_fields(node_.opt_fields))

        if not utils. _are_fields_defined(defined_fields) or \
           not utils._check_fields(fields[1:], SEGMENT_FIELDS):
            raise GFA2SerializationError("Required node elements " \
                                        + "missing or invalid.")

        return str.join("\t", fields)
    except(AttributeError, KeyError, GFA2SerializationError) as e:
        serializer_logger.debug(utils._format_exception(identifier, e))
        return ""

################################################################################
# EDGE SERIALIZER
################################################################################
def serialize_edge(edge_, identifier=DEFAULT_IDENTIFIER):
    """Converts to a GFA2 line the given edge.
    """
    identifier = utils._check_identifier(identifier)
    try:
        if isinstance(edge_, dict):
            if edge_['eid'] is None: # edge_ is a fragment
                return _serialize_to_fragment(edge_, identifier)
            if edge_['distance'] != None or \
              edge_['variance'] != None: # edge_ is a gap
                return _serialize_to_gap(edge_, identifier)
            return _serialize_to_edge(edge_, identifier)
        else:
            if edge_.eid is None: # edge_ is a fragment
                return _serialize_to_fragment(edge_, identifier)
            if edge_.distance != None or \
              edge_.variance != None: # edge_ is a gap
                return _serialize_to_gap(edge_, identifier)
            return _serialize_to_edge(edge_)

    except (KeyError, AttributeError) as e:
        serializer_logger.debug(utils._format_exception(identifier, e))
        return ""


def _serialize_to_fragment(fragment_, identifier=DEFAULT_IDENTIFIER):
    identifier = utils._check_identifier(identifier)
    try:
        if isinstance(fragment_, dict):

            fragment_dict = copy.deepcopy(fragment_)
            utils._remove_common_edge_fields(fragment_dict)
            defined_fields = [\
                                fragment_['from_node'], \
                                fragment_['to_node'], \
                                fragment_['to_orn'], \
                                fragment_['from_positions'][0], \
                                fragment_['from_positions'][1], \
                                fragment_['to_positions'][0], \
                                fragment_['to_positions'][1], \
                                fragment_['alignment'] \
                            ]
            fields = ["F"]
            fields.append(str(fragment_['from_node']))
            fields.append(str(fragment_['to_node']) + str(fragment_['to_orn']))
            fields.append(str(fragment_['from_positions'][0]))
            fields.append(str(fragment_['from_positions'][1]))
            fields.append(str(fragment_['to_positions'][0]))
            fields.append(str(fragment_['to_positions'][1]))
            fields.append(str(fragment_['alignment']))
            fields.extend(utils._serialize_opt_fields(fragment_dict))
        else:
            defined_fields = [\
                                fragment_.from_node, \
                                fragment_.to_node, \
                                fragment_.to_orn, \
                                fragment_.from_positions[0], \
                                fragment_.from_positions[1], \
                                fragment_.to_positions[0], \
                                fragment_.to_positions[1], \
                                fragment_.alignment \
                             ]
            fields = ["F"]
            fields.append(str(fragment_.from_node))
            fields.append(str(fragment_.to_node) + str(fragment_.to_orn))
            fields.append(str(fragment_.from_positions[0]))
            fields.append(str(fragment_.from_positions[1]))
            fields.append(str(fragment_.to_positions[0]))
            fields.append(str(fragment_.to_positions[1]))
            fields.append(str(fragment_.alignment))
            fields.extend(utils._serialize_opt_fields(fragment_.opt_fields))

        if not utils. _are_fields_defined(defined_fields) or \
           not utils._check_fields(fields[1:], FRAGMENT_FIELDS):
            raise GFA2SerializationError("Required Fragment elements " \
                                        + "missing or invalid.")

        return str.join("\t", fields)

    except(KeyError, AttributeError, GFA2SerializationError) as e:
        serializer_logger.debug(utils._format_exception(identifier, e))
        return ""


def _serialize_to_gap(gap_, identifier=DEFAULT_IDENTIFIER):
    identifier = utils._check_identifier(identifier)
    try:
        if isinstance(gap_, dict):
            gap_dict = copy.deepcopy(gap_)
            utils._remove_common_edge_fields(gap_dict)
            defined_fields = [\
                                gap_['eid'], \
                                gap_['from_node'], \
                                gap_['from_orn'], \
                                gap_['to_node'], \
                                gap_['to_orn'], \
                                gap_['distance'], \
                                gap_['variance'] \
                            ]
            fields = ["G"]
            fields.append(str(gap_['eid']))
            fields.append(str(gap_['from_node']) + str(gap_['from_orn']))
            fields.append(str(gap_['to_node']) + str(gap_['to_orn']))
            fields.append(str(gap_['distance']))
            fields.append(str(gap_['variance']))

            fields.extend(utils._serialize_opt_fields(gap_dict))
            return str.join("\t", fields)
        else:
            defined_fields = [\
                                gap_.eid, \
                                gap_.from_node, \
                                gap_.from_orn, \
                                gap_.to_node, \
                                gap_.to_orn, \
                                gap_.distance, \
                                gap_.variance \
                            ]
            fields = ["G"]
            fields.append(str(gap_.eid))
            fields.append(str(gap_.from_node) + str(gap_.from_orn))
            fields.append(str(gap_.to_node) + str(gap_.to_orn))
            fields.append(str(gap_.distance))
            fields.append(str(gap_.variance))
            fields.extend(utils._serialize_opt_fields(gap_.opt_fields))

        if not utils. _are_fields_defined(defined_fields) or \
           not utils._check_fields(fields[1:], GAP_FIELDS):
            raise GFA2SerializationError("Required Gap elements " \
                                        + "missing or invalid.")
        return str.join("\t", fields)
    except(AttributeError, KeyError, GFA2SerializationError) as e:
        serializer_logger.debug(utils._format_exception(identifier, e))
        return ""


def _serialize_to_edge(edge_, identifier=DEFAULT_IDENTIFIER):
    identifier = utils._check_identifier(identifier)
    try:
        if isinstance(edge_, dict):

            edge_dict = copy.deepcopy(edge_)
            utils._remove_common_edge_fields(edge_dict)
            defined_fields = [ \
                                edge_['eid'], \
                                edge_['from_node'], \
                                edge_['from_orn'], \
                                edge_['to_node'], \
                                edge_['to_orn'], \
                                edge_['from_positions'][0], \
                                edge_['from_positions'][1], \
                                edge_['to_positions'][0], \
                                edge_['to_positions'][1], \
                                edge_['alignment'] \
                             ]
            fields = ["E"]
            fields.append(str(edge_['eid']))
            fields.append(str(edge_['from_node']) + str(edge_['from_orn']))
            fields.append(str(edge_['to_node']) + str(edge_['to_orn']))
            fields.append(str(edge_['from_positions'][0]))
            fields.append(str(edge_['from_positions'][1]))
            fields.append(str(edge_['to_positions'][0]))
            fields.append(str(edge_['to_positions'][1]))
            fields.append(str(edge_['alignment']))
            fields.extend(utils._serialize_opt_fields(edge_dict))
        else:
            defined_fields = [ \
                                edge_.eid, \
                                edge_.from_node, \
                                edge_.from_orn, \
                                edge_.to_node, \
                                edge_.to_orn, \
                                edge_.from_positions[0], \
                                edge_.from_positions[1], \
                                edge_.to_positions[0], \
                                edge_.to_positions[1], \
                                edge_.alignment \
                             ]
            fields = ["E"]
            fields.append(str(edge_.eid))
            fields.append(str(edge_.from_node) + str(edge_.from_orn))
            fields.append(str(edge_.to_node) + str(edge_.to_orn))
            fields.append(str(edge_.from_positions[0]))
            fields.append(str(edge_.from_positions[1]))
            fields.append(str(edge_.to_positions[0]))
            fields.append(str(edge_.to_positions[1]))
            fields.append(str(edge_.alignment))
            fields.extend(utils._serialize_opt_fields(edge_.opt_fields))

        if not utils. _are_fields_defined(defined_fields) or \
           not utils._check_fields(fields[1:], EDGE_FIELDS):
            raise GFA2SerializationError("Required Edge elements " \
                                        + "missing or invalid.")

        return str.join("\t", fields)
    except(KeyError, AttributeError, GFA2SerializationError) as e:
        serializer_logger.debug(utils._format_exception(identifier, e))
        return ""

################################################################################
# SUBGRAPH SERIALIZER
################################################################################
def are_elements_oriented(subgraph_elements):
    """Check wheter all the elements of a subgraph have
    an orientation value `[+/-]`.
    """
    for id_, orientation in subgraph_elements.items():
        if orientation is None:
            return False
    return True

def _serialize_subgraph_elements(subgraph_elements, gfa_=None):
    """Serialize the elements belonging to a subgraph.

    Check if the orientation is provided for each element of the
    subgraph.

    :param subgraph_elements: The elements of a Subgraph.

    TODO
    ----
      Refactor list comprehension to function or cycle.

    """
    return str.join(" ", \
                    [str(id) \
                         + ((str(orientation)) if orientation != None \
                        else "") \
                            for id, orientation in subgraph_elements.items()])


def serialize_subgraph(subgraph_, identifier=DEFAULT_IDENTIFIER, gfa_=None):
    """Serialize a Subgraph object or an equivalent dictionary.

    :returns "": If subgraph cannot be serialized.

    :TODO:
        Check with `gfa` for OGroup in UGroup.
        See GFA2 spec.
    """
    identifier = utils._check_identifier(identifier)
    try:
        if isinstance(subgraph_, dict):
            subgraph_dict = copy.deepcopy(subgraph_)
            defined_fields = [\
                                subgraph_dict.pop('sub_id'), \
                                subgraph_dict.pop('elements') \
                             ]
            fields = ["O"] if are_elements_oriented(\
                                        subgraph_['elements']) else \
                     ["U"]
            fields.append(str(subgraph_['sub_id']))
            fields.append(_serialize_subgraph_elements(\
                                        subgraph_['elements'], gfa_))
            if 'overlaps' in subgraph_:
                subgraph_dict.pop('overlaps')
            fields.extend(utils._serialize_opt_fields(subgraph_dict))
        else:
            opt_fields = copy.deepcopy(subgraph_.opt_fields)
            defined_fields = [\
                                subgraph_.sub_id, \
                                subgraph_.elements \
                             ]
            fields = ["O"] if are_elements_oriented(subgraph_.elements) else \
                     ["U"]
            fields.append(str(subgraph_.sub_id))
            fields.append(_serialize_subgraph_elements(subgraph_.elements, gfa_))
            if 'overlaps' in subgraph_.opt_fields:
                opt_fields.pop('overlaps')
            fields.extend(utils._serialize_opt_fields(subgraph_.opt_fields))

        group_fields = OGROUP_FIELDS if fields[0] == "O" else \
                       UGROUP_FIELDS
        if not utils. _are_fields_defined(defined_fields) or \
           not utils._check_fields(fields[1:], group_fields):
            raise GFA2SerializationError("Required Subgraph elements " \
                                        + "missing or invalid.")

        return str.join("\t", fields)
    except(KeyError, ValueError, AttributeError, GFA2SerializationError) as e:
        serializer_logger.debug(utils._format_exception(identifier, e))
        return ""

################################################################################
# SERIALIZE GRAPH
################################################################################
def serialize_graph(graph, write_header=True):
    """Serialize a networkx.MultiGraph or a derivative object.

    :param graph: A networkx.MultiGraph instance.
    :param write_header: If set to True put a GFA2 header as first
        line.
    """
    if not isinstance(graph, nx.MultiGraph):
        raise ValueError("The object to serialize must be an instance" \
                        +" of a networkx.MultiGraph.")

    if write_header:
        string_serialize = "H\tVN:Z:2.0\n"

    for node_id, node_ in graph.nodes_iter(data=True):
        node_serialize = serialize_node(node_, node_id)
        if len(node_serialize) > 0:
            string_serialize += node_serialize + "\n"

    for from_node, to_node, key in graph.edges_iter(keys=True):
        edge_serialize = serialize_edge(graph.edge[from_node][to_node][key], key)
        if len(edge_serialize) > 0:
            string_serialize += edge_serialize + "\n"

    return string_serialize


def serialize_gfa(gfa_):
    """Serialize a GFA object into a GFA2 file.

    TODO:
        maybe process  the header fields here
    """
    gfa_serialize = serialize_graph(gfa_._graph, write_header=True)

    for sub_id, subgraph_ in gfa_.subgraphs().items():
        subgraph_serialize = serialize_subgraph(subgraph_, sub_id)
        if len(subgraph_serialize) > 0:
            gfa_serialize += subgraph_serialize + "\n"
    return gfa_serialize


if __name__ == '__main__': # pragma: no cover
    pass
