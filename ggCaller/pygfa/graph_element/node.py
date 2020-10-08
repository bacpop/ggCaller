import copy

from pygfa.graph_element.parser import segment
from pygfa.graph_element.parser import line, field_validator as fv

class InvalidNodeError(Exception):
    pass

def is_node(obj):
    """Check wheter the given object is a Node object.

    :param obj: Any Python object.
    :returns True: If obj can be treated as a Node object.
    """
    try:
        return obj.nid != None and \
          obj.sequence != None and \
          hasattr(obj, 'slen') and \
          hasattr(obj, 'opt_fields')
    except: return False


class Node:
    """A Node object that abstract the GFA1 and GFA2 Sequence concepts.

    GFA graphs will operate on Nodes, by adding them directly to their
    structures.

    Node accepts elements  (ids, sequences, lengths and so on) from
    the more tolerant of the two specification.
    So, a sequence will be accepted if and only if is a valid GFA2
    sequence, since GFA2 sequence is more tolerant than GFA1 sequence.
    """

    def __init__(self, node_id, sequence, length, opt_fields={}):
        """Construct a Node given an id, a sequence and a length.

        :param node_id: A node id given in the form of a string.
        :param sequence: A GFA1 or a GFA2 sequence.
        :param length: The length of the sequence. Can be `None`.
        :param opt_fields: A dictionary of Field or OptField objects.

        :raises InvalidNodeError: If node_id or sequence are undefined
            or not valid as the specification states.
        """
        if not isinstance(node_id, str) or node_id == '*':
            raise InvalidNodeError(\
                "A Node has always a defined id of type string, " \
                + "given {0} of type {1}".format(node_id, \
                                                type(node_id)))

        # checks sequence validation against GFA2 sequence specification
        # that is more tolerant than the GFA1 one
        if not(isinstance(sequence, str) and \
          fv.is_valid(sequence, fv.GFA2_SEQUENCE)):
            raise InvalidNodeError(\
                "A sequence must be of type string and must be a " \
                 + "valid GFA2 sequence," \
                 + "given '{0}' of type {1}".format(sequence, type(sequence)))

        if not( \
                (isinstance(length, int) and int(length) >= 0) or \
                length is None \
               ):
            raise InvalidNodeError(\
                "Sequence length must be a number >= 0, " + \
                "given {0} of type {1}".format(length, type(length)))
        self._nid = node_id
        self._sequence = sequence
        self._slen = length
        self._opt_fields = {}
        for key, field in opt_fields.items():
            if line.is_field(field):
                self._opt_fields[key] = copy.deepcopy(field)


    @property
    def nid(self):
        return self._nid

    @property
    def sequence(self):
        return self._sequence

    @property
    def slen(self):
        return self._slen

    @property
    def opt_fields(self):
        return self._opt_fields

    @classmethod
    def from_line(cls, segment_line):
        """Given a Segment Line construct a Node from it.

        If segment_line is a GFA1 Segment segment_line then the sequence length
        taken into account will be the value of the optional
        field `LN` if specified in the line fields.

        :param segment_line: A valid Segment Segment_line.
        :raises InvalidSegment_lineError: If the given segment_line
             is not valid.
        """
        try:
            fields = copy.deepcopy(segment_line.fields)
            if segment.is_segmentv1(segment_line):
                fields.pop('name')
                fields.pop('sequence')


                length = None
                if segment_line.fields['sequence'].value != "*":
                    length = len(segment_line.fields['sequence'].value)
                if 'LN' in segment_line.fields:
                    length = segment_line.fields['LN'].value
                    fields.pop('LN')

                return Node( \
                            segment_line.fields['name'].value, \
                            segment_line.fields['sequence'].value, \
                            length, \
                            fields)
            else:
                fields.pop('sid')
                fields.pop('sequence')
                fields.pop('slen')
                return Node( \
                            segment_line.fields['sid'].value, \
                            segment_line.fields['sequence'].value, \
                            segment_line.fields['slen'].value, \
                            fields)
        except(KeyError, AttributeError):
            raise line.InvalidLineError("The given line cannot be "\
                                        + "a Node.")


    def __eq__(self, other):
        try:
            if self.nid != other.nid or \
              self.sequence != other.sequence or \
              self.slen != other.slen:
                return False

            for key, item in self.opt_fields.items():
                if not key in other.opt_fields or \
                  not self.opt_fields[key] == other.opt_fields[key]:
                    return False

        except: return False
        return True


    def __neq__(self, other):
        return not self == other


    def __str__(self): # pragma: no cover
        fields = ("nid", "sequence", "slen", "opt_fields")
        opt_fields = []
        if len(self.opt_fields) > 0:
            opt_fields = str.join(\
                                   ",\t", \
                                   [str(field) for key, field in self.opt_fields.items()])
        values = (str(self.nid), \
                  str(self.sequence), \
                  str(self.slen), \
                   "{" + str(opt_fields) + "}" \
                  )
        assoc = [str.join(" : ", pair) for pair in zip(fields, values)]
        return str.join(",\t", assoc)


if __name__ == '__main__': # pragma: no cover
    pass
