import re

from pygfa.graph_element.parser import line, field_validator as fv

class Fragment(line.Line):

    def __init__(self):
        super().__init__('F')

    REQUIRED_FIELDS = { \
    'sid' : fv.GFA2_ID, \
    'external' : fv.GFA2_REFERENCE, \
    'sbeg' : fv.GFA2_POSITION, \
    'send' : fv.GFA2_POSITION, \
    'fbeg' : fv.GFA2_POSITION, \
    'fend' : fv.GFA2_POSITION, \
    'alignment' : fv.GFA2_ALIGNMENT \
    }

    @classmethod
    def from_string(cls, string):
        """Extract the fragment fields from the string.

        The string can contains the F character at the begin or can
        only contains the fields of the fragment directly.
        """
        if len(string.split()) == 0:
            raise line.InvalidLineError("Cannot parse the empty string.")
        fields = re.split('\t', string)
        ffields = []
        if fields[0] == 'F':
            fields = fields[1:]

        if len(fields) < len(cls.REQUIRED_FIELDS):
            raise line.InvalidLineError("The minimum number of field for "
                                        + "Fragment line is not reached.")

        fragment = Fragment()
        sid_f = fv.validate(fields[0], cls.REQUIRED_FIELDS['sid'])
        ffields.append(line.Field('sid', sid_f))

        external_f = fv.validate(fields[1], cls.REQUIRED_FIELDS['external'])
        ffields.append(line.Field('external', external_f))

        sbeg_f = fv.validate(fields[2], cls.REQUIRED_FIELDS['sbeg'])
        ffields.append(line.Field('sbeg', sbeg_f))

        send_f = fv.validate(fields[3], cls.REQUIRED_FIELDS['send'])
        ffields.append(line.Field('send', send_f))

        fbeg_f = fv.validate(fields[4], cls.REQUIRED_FIELDS['fbeg'])
        ffields.append(line.Field('fbeg', fbeg_f))

        fend_f = fv.validate(fields[5], cls.REQUIRED_FIELDS['fend'])
        ffields.append(line.Field('fend', fend_f))

        alignment_f = fv.validate(fields[6], cls.REQUIRED_FIELDS['alignment'])
        ffields.append(line.Field('alignment', alignment_f))

        for field in fields[7:]:
            ffields.append(line.OptField.from_string(field))

        for field in ffields:
            fragment.add_field(field)

        return fragment

if __name__ == '__main__': # pragma: no cover
    pass
