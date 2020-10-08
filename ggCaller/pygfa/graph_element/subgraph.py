import copy
import collections

from pygfa.graph_element.parser import line


class InvalidSubgraphError(Exception): pass


def is_subgraph(obj):
    try:
        return obj.sub_id != None and \
          obj.elements != None and \
          hasattr(obj, 'opt_fields')
    except: return False

class Subgraph:

    def __init__(self, graph_id, elements, opt_fields={}):
        """Create a Subgraph object.

        Subgraphs can be created from the following lines:
        * Path
        * OGroup
        * UGroup

        :param graph_id: The id of the subgraph. It can be '*', but
            cannot be None.
        :param elements: An ordered dictionary of node id and orientations.
        :param opt_fields: A dictionary of Fields or OptFields.

        :note:
            In case of Path line, the `overlaps` won't be represented as
            an explicit Subgraph attribute, but it will be added
            into the opt_fields dictionary with the name `overlaps`.

            This choice has been done since the overlaps can be computed
            with the alignment of each edge in the Subgraph.
        """
        if not isinstance(graph_id, str):
            raise InvalidSubgraphError("A subgraph has always an id " \
                                       + "of type string, " \
                                       +"given {0} of type {1}".format( \
                                                                graph_id, \
                                                                type(graph_id)))
        if not isinstance(elements, dict):
            raise InvalidSubgraphError(\
                    "A dictionary of elements id:orientation is required.")
        self._sub_id = graph_id
        self._elements = copy.deepcopy(elements)
        self._opt_fields = {}
        for key, field in opt_fields.items():
            if line.is_field(field):
                self._opt_fields[key] = copy.deepcopy(field)

    def is_path(self):
        for element, orn in self.elements.items():
            if orn is None:
                return False
        return True

    @property
    def sub_id(self):
        return self._sub_id

    @property
    def elements(self):
        return self._elements

    @property
    def opt_fields(self):
        return self._opt_fields

    @classmethod
    def from_line(cls, line_):
        try:
            fields = copy.deepcopy(line_.fields)
            if line_.type == 'P':
                fields.pop('path_name')
                fields.pop('seqs_names')
                names = collections.OrderedDict(\
                            (ref[0:-1], ref[-1:]) \
                                for ref in line_.fields['seqs_names'].value)
                return Subgraph( \
                                line_.fields['path_name'].value, \
                                names, \
                                fields)
            if line_.type == 'O':
                fields.pop('oid')
                fields.pop('references')
                refs = collections.OrderedDict(\
                            (ref[0:-1], ref[-1:]) \
                                for ref in line_.fields['references'].value)
                return Subgraph( \
                                line_.fields['oid'].value, \
                                refs, \
                                fields)
            if line_.type == 'U':
                fields.pop('uid')
                fields.pop('ids')
                ids = collections.OrderedDict(\
                            (id, None) \
                                for id in line_.fields['ids'].value)
                return Subgraph( \
                                line_.fields['uid'].value, \
                                ids, \
                                fields)
        except(KeyError, AttributeError):
            raise line.InvalidLineError("The given line cannot be "\
                                        + "a Subgraph.")


    def as_dict(self):
        """Turn the Subgraph into a dictionary,

        Put all fields and the optional fields into a dictionary.
        """
        retval = {}
        retval['sub_id'] = self.sub_id
        retval['elements'] = self.elements
        for key, value in self.opt_fields.items():
            retval[key] = value
        return retval


    def __eq__(self, other):
        try:
            if self.sub_id != other.sub_id \
              or self.elements != other.elements:
                return False
            for key, field in self.opt_fields.items():
                if not key in other.opt_fields \
                  or field != other.opt_fields[key]:
                    return False
            return True
        except: return False


    def __neq__(self, other):
        return not self == other


    def __str__(self): # pragma: no cover
        fields = ("sub_id", "elements", "opt_fields")
        opt_fields_ = []
        if len(self.opt_fields) > 0:
            opt_fields_ = str.join(\
                                   ",\t", \
                                   [str(field) for key, field in self.opt_fields.items()])
        elements_ = str.join(\
                        "\t", \
                        [id+(orn if orn != None else "") \
                             for id, orn in self.elements.items()])

        values = [self.sub_id, \
                  elements_, \
                  "{" + str(opt_fields_) + "}" \
                 ]
        assoc = [str.join(" : ", pair) for pair in zip(fields, values)]
        return str.join(",\t", assoc)


if __name__ == '__main__': # pragma: no cover
    pass
