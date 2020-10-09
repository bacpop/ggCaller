import re

from ggCaller.pygfa.graph_element.parser import line, field_validator as fv

class Path(line.Line):

    def __init__(self):
        super().__init__('P')

    REQUIRED_FIELDS = { \
    'path_name' : fv.GFA1_NAME, \
    'seqs_names' : fv.GFA1_NAMES, \
    'overlaps': fv.GFA1_CIGARS \
    }

    PREDEFINED_OPTFIELDS = {}

    @classmethod
    def from_string(cls, string):
        """Extract the path fields from the string.

        The string can contains the P character at the begin or can
        just contains the fields of the path directly.
        """
        if len(string.split()) == 0:
            raise line.InvalidLineError("Cannot parse the empty string.")
        fields = re.split('\t', string)
        pfields = []
        if fields[0] == 'P':
            fields = fields[1:]

        if len(fields) < len(cls.REQUIRED_FIELDS):
            raise line.InvalidLineError("The minimum number of field for "
                                        + "Path line is not reached.")
        path = Path()
        path_name = fv.validate(fields[0], cls.REQUIRED_FIELDS['path_name'])
        sequences_names = [fv.validate(label, \
                            cls.REQUIRED_FIELDS['seqs_names']) \
                            for label in fields[1].split(",")
                          ]

        overlaps = fv.validate(fields[2], cls.REQUIRED_FIELDS['overlaps'])

        pfields.append(line.Field('path_name', path_name))
        pfields.append(line.Field('seqs_names', sequences_names))
        pfields.append(line.Field('overlaps', overlaps))

        for field in fields[3:]:
            pfields.append(line.OptField.from_string(field))

        for field in pfields:
            path.add_field(field)

        return path

if __name__ == '__main__': # pragma: no cover
    pass
