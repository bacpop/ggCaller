"""
Field validation module to check each field string against GFA1
and GFA2 specification.
"""
import re

class InvalidFieldError(Exception):
    """Exception raised when an invalid field is provided."""

class UnknownDataTypeError(Exception):
    """Exception raised when the datatype provided is not in
    the `DATASTRING_VALIDATION_REGEXP` dictionary.
    """

class FormatError(Exception):
    """Exception raised when a wrong type of object is given
    to the validator.
    """

TYPE_A = 'A'
TYPE_i = 'i'
TYPE_f = 'f'
TYPE_Z = 'Z'
JSON = 'J'
HEX_BYTE_ARRAY = 'H'
DEC_ARRAY = 'B'

GFA1_NAME = 'lbl'
GFA1_NAMES = 'lbs'
GFA1_ORIENTATION = 'orn'
GFA1_SEQUENCE = 'seq'
GFA1_CIGAR = 'cig'
GFA1_CIGARS = 'cgs'
GFA1_INT = 'pos'

GFA2_ID = 'id'
GFA2_IDS = 'ids'
GFA2_REFERENCE = 'ref'
GFA2_REFERENCES = 'rfs'
GFA2_INT = 'int'
GFA2_TRACE = 'trc'
GFA2_ALIGNMENT = 'aln'
GFA2_POSITION = 'pos2'
GFA2_CIGAR = 'cig2'
GFA2_SEQUENCE = 'seq2'
GFA2_OPTIONAL_INT = 'oint'
GFA2_OPTIONAL_ID = 'oid'


# These are the types of value a field can assume.
# These are the same as the ones in rgfa, I've extended
# the list to support GFA2.
#
# GFA2: 'id', 'ids', 'ref', 'rfs', 'cig2', 'oid'(opt_id),
# 'trc', 'aln', 'pos2', 'seq2', 'int', 'oint'

DATASTRING_VALIDATION_REGEXP = \
  {\
  TYPE_A : "^[!-~]", \
  # any printable character
  #
  TYPE_i : "^[-+]?[0-9]+$", \
  # Signed integer
  #
  TYPE_f : "^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$", \
  # Single-precision floating number
  #
  TYPE_Z : "^[ !-~]+$", \
  # Printable string, including space
  #
  JSON : "^[ !-~]+$",  \
  # JSON, excluding new-line and tab characters
  #
  HEX_BYTE_ARRAY : "^[0-9A-F]+$", \
  # Byte array in the Hex format
  #
  DEC_ARRAY : "^[cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+$", \
  # Integer or numeric array
  #
  GFA1_NAME : "^[!-)+-<>-~][!-~]*$", \
  # segment/path label(segment name)
  #
  GFA1_ORIENTATION : "^\+|-$", \
  #segment orientation
  #
  ###
  #'lbs' : "^[!-)+-<>-~][!-~]*[+-](,[!-)+-<>-~][!-~]*[+-])+$",
  # multiple labels with orientations, comma-sep
  #
  # Changed according to issue 59, since the comma is accepted by [!-~],
  # it's not possible to make a clear regexp for an array of labels,
  # so the implementation has been modified to reflect
  # this behaviour, splitting the labels and checking them one by one
  # with the new lbs regexp beyond.
  #
  GFA1_NAMES : "^[!-)+-<>-~][!-~]*[+-]$", \
  GFA1_SEQUENCE : "^\*$|^[A-Za-z=.]+$", \
  # nucleotide sequence(segment sequence)
  #
  GFA1_INT : "^[0-9]*$", \
  # positive integer(CLAIM ISSUE HERE, MOVE TO -> int)
  #
  GFA1_CIGAR : "^(\*|(([0-9]+[MIDNSHPX=])+))$", # CIGAR string \
  GFA1_CIGARS : "^(\*|(([0-9]+[MIDNSHPX=])+))(,(\*|(([0-9]+[MIDNSHPX=])+)))*$", \
  # multiple CIGARs, comma-sep \
  #
  'cmt' : ".*", \
  # content of comment line, everything is allowed
  #
  GFA2_ID : "^[!-~]+$", \
  # it's the lbl for GFA2
  #
  GFA2_IDS : "^[!-~]+([ ][!-~]+)*$", \
  # it's the lbs for GFA2
  #
  GFA2_REFERENCE : "^[!-~]+[+-]$", \
  GFA2_REFERENCES : "^[!-~]+[+-]([ ][!-~]+[+-])*$", \
  # array of references
  #
  GFA2_INT : "^[0-9]+$", \
  # GFA1 has pos to describe any positive integer,
  # but pos accept the empty string, while 'int' doesn't
  #
  GFA2_TRACE : "^[0-9]+(,[0-9]+)*$", \
  GFA2_ALIGNMENT : "^\*$|^[0-9]+(,[0-9]+)*$|^([0-9]+[MDIP])+$", \
  GFA2_POSITION : "^[0-9]+\$?$", \
  # pos2 represent a position in GFA2, it's similar in NO WAY
  # to pos which represent a positive integer in GFA1
  #
  GFA2_CIGAR : "^([0-9]+[MDIP])+$", \
  # CIGAR string for GFA2
  #
  GFA2_SEQUENCE : "^\*$|^[!-~]+$", # seq2 is a GFA2 sequence,
  # it's more flexible than GFA1 seq
  #
  GFA2_OPTIONAL_ID : "^\*$|^[!-~]+$", \
  # optional id  for GFA2
  #
  GFA2_OPTIONAL_INT : "^\*$|^[0-9]+$" \
  # optional int
  #
  }


def is_valid(string, datatype):
    """Check if the string respects the datatype.

    :param datatype: The type of data corresponding to the string.
    :returns: True if the string respect the type defined by the datatype.
    :raises UnknownDataTypeError: If the datatype is not presents in
        `DATASTRING_VALIDATION_REGEXP`.
    :raises UnknownFormatError: If string is not python string.

    :TODO:
        Fix exception reference in the documentation.
    """
    if not isinstance(string, str):
        raise FormatError("A string must be given to validate it, " \
                         + "given:{0}".format(string))
    if not datatype in DATASTRING_VALIDATION_REGEXP:
        raise UnknownDataTypeError(\
                                    "Invalid field datatype," + \
                                    "given: {0}".format(datatype) \
                                  )
    regexp = DATASTRING_VALIDATION_REGEXP[datatype]
    if not re.fullmatch(regexp, string):
        return False
    return True


def is_dazzler_trace(string):
    return is_valid(string, GFA2_TRACE)

def is_gfa1_cigar(string):
    """Check if the given string is a valid CIGAR string
    as defined in the GFA1 specification.
    """
    return string != "*" and is_valid(string, GFA1_CIGAR)


def is_gfa2_cigar(string):
    """Check if the given string is a valid CIGAR string
    as defined in the GFA2 specification.
    """
    return string != "*" and is_valid(string, GFA2_CIGAR)


def validate(string, datatype):
    """Return a value from the given string with the type closer to the
    one it's represented.
    """
    if not is_valid(string, datatype):
        raise InvalidFieldError("The string cannot be validated within " \
                               + "its datatype,\n" \
                               + "given string : " \
                               + "{0}\ndatatype: {1}.".format(string, \
                                                              datatype))
    if datatype in (TYPE_i,):
        return int(string)
    elif datatype in(GFA1_INT, GFA2_INT):
        # fullmatch grants that we have a string whose int value is >= 0

        # position = int(string)
        # if position < 0:
        #     raise Exception("Position must be >= 0.")
        return int(string)

    elif datatype in (GFA2_OPTIONAL_INT, ):
        if string == "*":
            return string
        return int(string)

    elif datatype in (TYPE_f, ):
        return float(string)

    elif datatype in (GFA1_CIGARS, ):
        return string.split(",")
    elif datatype in (GFA2_ALIGNMENT, ):
        # string is either * or a trace or a cigar
        if string == "*":
            return string
        elif is_valid(string, GFA2_CIGAR):
            return validate(string, GFA2_CIGAR)
        return validate(string, GFA2_TRACE)

    elif datatype in (JSON, ):
        return string # TODO: ask if the json must be manipulated
    elif datatype in(GFA2_IDS, GFA2_REFERENCES):
        return string.split()
    else:
        # 'orn', 'A', 'Z', 'seq', 'lbl', 'cig', 'cig2', 'H',
        # 'B', 'trc', 'id', 'ref', pos2', 'seq2', 'oid', 'lbs'
        return string


if __name__ == '__main__': # pragma: no cover
    pass
