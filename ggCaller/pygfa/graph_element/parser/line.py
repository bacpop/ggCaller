import re

from pygfa.graph_element.parser import field_validator as fv

class InvalidLineError(Exception):
    """Exception raised when making a Line object from a string.
    The number of fields gained by splittin the string
    must be equal to or great than the number of required field
    ecluding the optional first field indicating the type of the line.
    """

# support for duck typing
def is_field(field):
    """Check if the given object is a valid field

    A field is valid if it has at least a name and a value
    attribute/property.
    """
    for attr in('name', 'value'):
        if not hasattr(field, attr):
            return False
        if field.name is None or field.value is None:
            return False
    if not isinstance(field.name, str):
        return False
    return True


def is_optfield(field):
    """Check if the given object is an optfield

    A field is an optfield if it's a field with name that
    match a given expression and its type is defined.
    """
    return is_field(field) and \
      re.fullmatch('[A-Za-z0-9]' * 2, field.name) and \
      hasattr(field, '_type') and \
      field.type != None




class Line:
    """
    A generic Line, it's unlikely that it will be directly instantiated
    (but could be done so).
    Its subclasses should be used instead.

    It's possible to instatiate a Line to save a custom line in a gfa file.
    """

    REQUIRED_FIELDS = {}
    PREDEFINED_OPTFIELDS = {}
    # this will contain tha name of the required optfield for
    # each kind of line and the ralative type of value the value of
    # the field must contains


    def __init__(self, line_type=None):
        self._fields = {}
        self._type = line_type

    @classmethod
    def is_valid(cls, line_):
        """Check if the line is valid.

        Defining the method here allows to have automatically validated
        all the line of the specifications.
        """
        # use polymorphism to get the type and the required fields of a
        # specific kind of line.
        instance = cls()
        try:
            if line_.type != cls().type:
                return False
            for required_field in instance.REQUIRED_FIELDS:
                if not required_field in line_.fields:
                    return False
            return True
        except (AttributeError, KeyError):
            return False


    @classmethod
    def get_static_fields(cls):
        keys = []
        values = []
        for key, value in cls.REQUIRED_FIELDS.items():
            keys.append(key)
            values.append(value)

        for key, value in cls.PREDEFINED_OPTFIELDS.items():
            keys.append(key)
            values.append(value)
        return dict(zip(keys, values))


    @property
    def type(self):
        return self._type


    @property
    def fields(self):
        return self._fields


    def add_field(self, field):
        """Add a field to the line.

        It's possible to add a Field if an only if its name
        is in the `REQUIRED_FIELDS` dictionary. Otherwise
        the field will be considered as an optional field and
        an InvalidFieldError will be raised.

        :param field: The field to add to the line
        :raises InvalidFieldError: If a 'name' and a 'value' attributes are
            not found or the field has already been added.

        :note:
            If you want to add a Field for a custom Line object be
            sure to add its name to the REQUIRED_FIELDS dictionary
            for that particular Line subclass.
        """
        if not(is_field(field) or is_optfield(field)):
            raise  fv.InvalidFieldError("A valid field must be attached")

        if field.name in self.fields:
            raise ValueError(\
                    "This field is already been added, field name: '{0}'.".format(field.name))

        if field.name in self.REQUIRED_FIELDS:
            self._fields[field.name] = field
        else: # here we are appending an optfield
            if not is_optfield(field):
                raise fv.InvalidFieldError(\
                    "Cannot add an invalid OptField.")
            self._fields[field.name] = field
        return True


    def remove_field(self, field):
        """
        If the field is contained in the line it gets removed.
        Otherwise it does nothing, without raising any exception.
        """
        field_name = field
        if is_field(field):
            field_name = field.name

        if field_name in self.fields:
            self.fields.pop(field_name)


    @classmethod
    def from_string(cls, string): # pragma: no cover
        raise NotImplementedError

    def __eq__(self, other):
        try:
            return self.type == other.type \
              and self.fields == other.fields
        except:
            return False


    def __neq__(self, other):
        return not self == other


    def __str__(self): # pragma: no cover
        tmp_str = "line_type: {0}, fields: [".format(str(self.type))
        field_strings = []

        for field in self.fields:
            field_strings.append(str(field))

        tmp_str += str.join(", ", field_strings) + "]"
        return tmp_str




class Field:
    """This class represent any required field.

    The type of field is bound to the field name.
    """
    def __init__(self, name, value):
        self._name = name
        self._value = value

    @property
    def name(self):
        return self._name

    @property
    def value(self):
        return self._value


    def __eq__(self, other):
        try:
            return self.name == other.name and \
              self.value == other.value
        except:
            return False

    def __neq__(self, other):
        return not self == other

    def __str__(self): # pragma: no cover
        return str.join(":", (self.name, str(self.value)))




class OptField(Field):
    """An Optional field of the form `TAG:TYPE:VALUE`, where:
    TAG match [A-Za-z0-9][A-Za-z0-9]
    TYPE match [AiZfJHB]
    """
    def __init__(self, name, value, field_type):
        if not re.fullmatch('[A-Za-z0-9]' * 2, name):
            raise ValueError("Invalid optfield name, given '{0}'".format(name))

        if not re.fullmatch("^[ABHJZif]$", field_type):
            raise ValueError("Invalid type for an optional field.")

        self._name = name
        self._type = field_type
        self._value = fv.validate(value, field_type)


    @property
    def type(self):
        return self._type


    @classmethod
    def from_string(cls, string):
        """Create an OptField with a given string.
        """
        groups = re.split(":", string.strip())
        if len(groups) != 3:
            raise ValueError(\
                    "OptField must have a name, a type and a value," \
                    + " given{0}".format(string))

        optfield = OptField(groups[0], groups[2], groups[1])
        return optfield


    def __eq__(self, other):
        try:
            return self.name == other.name and \
              self.value == other.value and \
              self.type == other.type
        except:
            return False


    def __neq__(self, other):
        return not self == other


    def __str__(self): # pragma: no cover
        return str.join(":", (self.name, self.type, str(self.value)))

if __name__ == '__main__': # pragma: no cover
    pass
