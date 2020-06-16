import re

from pygfa.graph_element.parser import line, field_validator as fv

class Containment(line.Line):

    def __init__(self):
        super().__init__('C')

    REQUIRED_FIELDS = { \
    'from' : fv.GFA1_NAME, \
    'from_orn' : fv.GFA1_ORIENTATION, \
    'to': fv.GFA1_NAME, \
    'to_orn' : fv.GFA1_ORIENTATION, \
    'pos' : fv.GFA1_INT, \
    'overlap' : fv.GFA1_CIGAR \
    }

    PREDEFINED_OPTFIELDS = { \
    'NM' : fv.TYPE_i, \
    'RC' : fv.TYPE_i, \
    'ID' : fv.TYPE_Z \
    }

    @classmethod
    def from_string(cls, string):
        """Extract the containment fields from the string.

        The string can contains the C character at the begin or can
        only contains the fields of the containment directly.
        """
        if len(string.split()) == 0:
            raise line.InvalidLineError("Cannot parse the empty string.")
        fields = re.split('\t', string)
        cfields = []
        if fields[0] == 'C':
            fields = fields[1:] #skip the first field(the C)

        if len(fields) < len(cls.REQUIRED_FIELDS):
            raise line.InvalidLineError("The minimum number of field for "
                                        + "Containment line is not reached.")

        containment = Containment()

        from_name = fv.validate(fields[0], cls.REQUIRED_FIELDS['from'])
        from_orn = fv.validate(fields[1], cls.REQUIRED_FIELDS['from_orn'])
        to_name = fv.validate(fields[2], cls.REQUIRED_FIELDS['to'])
        to_orn = fv.validate(fields[3], cls.REQUIRED_FIELDS['to_orn'])
        pos = fv.validate(fields[4], cls.REQUIRED_FIELDS['pos'])
        overlap = fv.validate(fields[5], cls.REQUIRED_FIELDS['overlap'])

        cfields.append(line.Field('from', from_name))
        cfields.append(line.Field('from_orn', from_orn))
        cfields.append(line.Field('to', to_name))
        cfields.append(line.Field('to_orn', to_orn))
        cfields.append(line.Field('pos', pos))
        cfields.append(line.Field('overlap', overlap))

        for field in fields[6:]:
            cfields.append(line.OptField.from_string(field))

        for field in cfields:
            containment.add_field(field)

        return containment

if __name__ == '__main__': # pragma: no cover
    pass
