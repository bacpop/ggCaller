import re

from pygfa.graph_element.parser import line
from pygfa.graph_element.parser import field_validator as fv

class Edge(line.Line):

    def __init__(self):
        super().__init__('E')

    REQUIRED_FIELDS = { \
    'eid' : fv.GFA2_OPTIONAL_ID, \
    'sid1' : fv.GFA2_REFERENCE, \
    'sid2' : fv.GFA2_REFERENCE, \
    'beg1' : fv.GFA2_POSITION, \
    'end1' : fv.GFA2_POSITION, \
    'beg2' : fv.GFA2_POSITION, \
    'end2' : fv.GFA2_POSITION, \
    'alignment' : fv.GFA2_ALIGNMENT \
    }

    @classmethod
    def from_string(cls, string):
        """Extract the Edge fields from the string.

        The string can contains the E character at the begin or can
        only contains the fields of the Edge directly.
        """
        if len(string.split()) == 0:
            raise line.InvalidLineError("Cannot parse the empty string.")
        fields = re.split('\t', string)
        efields = []
        if fields[0] == 'E':
            fields = fields[1:]

        if len(fields) < len(cls.REQUIRED_FIELDS):
            raise line.InvalidLineError("The minimum number of field for "
                                        + "Edge line is not reached.")


        edge = Edge()

        eid_f = fv.validate(fields[0], cls.REQUIRED_FIELDS['eid'])
        efields.append(line.Field('eid', eid_f))

        sid1_f = fv.validate(fields[1], cls.REQUIRED_FIELDS['sid1'])
        efields.append(line.Field('sid1', sid1_f))

        sid2_f = fv.validate(fields[2], cls.REQUIRED_FIELDS['sid2'])
        efields.append(line.Field('sid2', sid2_f))

        beg1_f = fv.validate(fields[3], cls.REQUIRED_FIELDS['beg1'])
        efields.append(line.Field('beg1', beg1_f))

        end1_f = fv.validate(fields[4], cls.REQUIRED_FIELDS['end1'])
        efields.append(line.Field('end1', end1_f))

        beg2_f = fv.validate(fields[5], cls.REQUIRED_FIELDS['beg2'])
        efields.append(line.Field('beg2', beg2_f))

        end2_f = fv.validate(fields[6], cls.REQUIRED_FIELDS['end2'])
        efields.append(line.Field('end2', end2_f))

        alignment_f = fv.validate(fields[7], cls.REQUIRED_FIELDS['alignment'])
        efields.append(line.Field('alignment', alignment_f))

        for field in fields[8:]:
            efields.append(line.OptField.from_string(field))

        for field in efields:
            edge.add_field(field)

        return edge

if __name__ == '__main__': # pragma: no cover
    pass
