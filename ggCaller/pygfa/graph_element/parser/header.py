import re

from ggCaller.pygfa.graph_element.parser import line, field_validator as fv

class Header(line.Line):

    PREDEFINED_OPTFIELDS = \
      { \
        'VN' : fv.TYPE_Z, \
        'TS' : fv.TYPE_i \
      }

    def __init__(self):
        super().__init__('H')

    @classmethod
    def from_string(cls, string):
        """Extract the header fields from the string.

        The string can contains the H character at the begin or can
        only contains the fields of the header directly.
        """
        if len(string.split()) == 0:
            raise line.InvalidLineError("Cannot parse the empty string.")
        fields = re.split('\t', string)
        hfields = []
        if fields[0] == 'H':
            fields = fields[1:]

        # header lines have no required fields.
        header = Header()
        for field in fields:
            hfields.append(line.OptField.from_string(field))

        for field in hfields:
            header.add_field(field)

        return header

if __name__ == '__main__': # pragma: no cover
    pass
