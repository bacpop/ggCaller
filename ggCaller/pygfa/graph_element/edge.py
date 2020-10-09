import copy

from ggCaller.pygfa.graph_element.parser import line

class InvalidEdgeError(Exception):
    pass

def is_edge(obj):
    try:
        return \
          hasattr(obj, 'eid') and \
          obj.from_node != None and \
          obj.to_node != None and \
          hasattr(obj, 'from_positions') and \
          hasattr(obj, 'from_orn') and \
          hasattr(obj, 'to_positions') and \
          hasattr(obj, 'to_orn') and \
          hasattr(obj, 'alignment')
    except: return False


class Edge:

    def __init__(self, \
                  edge_id, \
                  from_node, from_orn, \
                  to_node, to_orn, \
                  from_positions, \
                  to_positions, \
                  alignment, \
                  distance=None, \
                  variance=None, \
                  opt_fields={}, \
                  is_dovetail=False):
        """Construct an Edge from a Line subclass that represent can be
        represented as an Edge.

        Line subclasses that can be valid Edge objects are:
        * Edge
        * Fragment
        * Gap
        * Link
        * Containment

        The Edge class try to represent the majority of fields
        this lines have in common.

        Fields that aren't fundamental (such as the 'pos' field in the
        Containment line) will be put in the opt_fields
        dictionary.

        :param eid: The id of the edge, can be '*' or None. If None
            that edge is considered to be a Fragment line in future
            operations.
        :param from_node: The node where the edge comes from.
        :param from_orn: The orientation of the source node sequence.
        :param to_node: The node where the edge ends to.
        :param to_orn: The orientation of the end node sequence.
        :param from_positions: A tuple composed by 2 string values
            that represent the positions where the alignment occurs
            in the source node.
        :param from_positions: A tuple composed by 2 string values
            that represent the positions where the alignment occurs
            in the source end node.
        :param alignment: A string that represents the alignment between
            the nodes specified.
        :param distance: The distance between two nodes indicated in a
            GFA2 Gap line.
        :param variance: The variance values indicated in a Gap line.
        :param opt_fields: A dictionary of OptFields.

        :raises InvalidEdgeError: If the line cannot be represented as
            Edge.

        :note:
            Eventhough the field is called opt_fields, actually it can
            contains also Field object like the 'pos' field
            of the Containment line. This will happen for line fields
            that haven't been represented as an Edge atrtibute.
        """
        if not(isinstance(from_positions, tuple) and \
            len(from_positions) == 2):
            raise InvalidEdgeError("Invalid from_node tuple: given: {0}".format(\
                str(from_positions)))
        if not(isinstance(to_positions, tuple) and \
            len(to_positions) == 2):
            raise InvalidEdgeError("Invalid to_position tuple: given: {0}".format(\
                str(to_positions)))
        self._eid = edge_id
        self._from_node = from_node
        self._from_orn = from_orn
        self._to_node = to_node
        self._to_orn = to_orn
        self._from_positions = copy.deepcopy(from_positions)
        self._to_positions = copy.deepcopy(to_positions)
        self._alignment = alignment

        self._distance = distance
        self._variance = variance

        self._opt_fields = {}
        for key, field in opt_fields.items():
            if line.is_field(field):
                self._opt_fields[key] = copy.deepcopy(field)

        self._is_dovetail = is_dovetail
        self._from_segment_end = None
        self._to_segment_end = None
        if self._is_dovetail:
            self._set_segments_end()


    @property
    def eid(self):
        return self._eid

    @property
    def from_node(self):
        return self._from_node

    @property
    def from_orn(self):
        return self._from_orn

    @property
    def to_orn(self):
        return self._to_orn

    @property
    def to_node(self):
        return self._to_node

    @property
    def from_positions(self):
        return self._from_positions

    @property
    def to_positions(self):
        return self._to_positions

    @property
    def alignment(self):
        return self._alignment

    @property
    def distance(self):
        return self._distance

    @property
    def variance(self):
        return self._variance

    @property
    def opt_fields(self):
        return self._opt_fields

    @property
    def is_dovetail(self):
        return self._is_dovetail

    @property
    def from_segment_end(self):
        return self._from_segment_end

    @property
    def to_segment_end(self):
        return self._to_segment_end

    def _set_segments_end(self):
        """Set the segments ends considerend by the
        edge, only for dovetail overlaps.
        Do nothing otherwise.
        """
        if not self.is_dovetail:
            return

        # check if it's a GFA1 Link 
        if self.from_positions == (None, None) \
          or self.to_positions == (None, None):
            self._set_segments_end_link()
        else:
            self._set_segments_end_edge()


    def _set_segments_end_link(self):
        if self.from_orn == "+":
            self._from_segment_end = "R"
        else:
            self._from_segment_end = "L"
        if self.to_orn == "+":
            self._to_segment_end = "L"
        else:
            self._to_segment_end = "R"


    def _set_segments_end_edge(self):
        beg1, end1 = self.from_positions
        beg2, end2 = self.to_positions

        # the dovetail is from the end of from_node to
        # the beginning of to_node, just like a GFA1 Link
        if beg2 == "0":
            self._set_segments_end_link()
        else: # dovetail between end of to_node and begin of from_node
            if self.from_orn == "+":
                self._from_segment_end = "L"
            else:
                self._from_segment_end = "R"
            if self.to_orn == "+":
                self._to_segment_end = "R"
            else:
                self._to_segment_end = "L"


    @classmethod
    def from_line(cls, line_):
        try:
            fields = copy.deepcopy(line_.fields)
            if line_.type == 'L':
                if 'ID' in line_.fields:
                    fields.pop('ID')
                fields.pop('from')
                fields.pop('from_orn')
                fields.pop('to')
                fields.pop('to_orn')
                fields.pop('overlap')
                return Edge( \
                    '*' if 'ID' not in line_.fields else \
                        line_.fields['ID'].value, \
                    line_.fields['from'].value, \
                    line_.fields['from_orn'].value, \
                    line_.fields['to'].value, \
                    line_.fields['to_orn'].value, \
                    (None, None), \
                    (None, None), \
                    line_.fields['overlap'].value, \
                    opt_fields=fields, \
                    is_dovetail=True)

            if line_.type == 'C':
                if 'ID' in line_.fields:
                    fields.pop('ID')
                fields.pop('from')
                fields.pop('from_orn')
                fields.pop('to')
                fields.pop('to_orn')
                fields.pop('overlap')
                return Edge( \
                    '*' if 'ID' not in line_.fields else \
                    line_.fields['ID'].value, \
                    line_.fields['from'].value, \
                    line_.fields['from_orn'].value, \
                    line_.fields['to'].value, \
                    line_.fields['to_orn'].value, \
                   (None, None), \
                   (None, None), \
                    line_.fields['overlap'].value,\
                    opt_fields=fields)

            if line_.type == 'F':
                fields.pop('sid')
                fields.pop('external')
                fields.pop('sbeg')
                fields.pop('send')
                fields.pop('fbeg')
                fields.pop('fend')
                fields.pop('alignment')
                return Edge( \
                    None, \
                    line_.fields['sid'].value, None, \
                    line_.fields['external'].value[0:-1], \
                    line_.fields['external'].value[-1:], \
                    (line_.fields['sbeg'].value, \
                        line_.fields['send'].value), \
                    (line_.fields['fbeg'].value, line_.fields['fend'].value), \
                    line_.fields['alignment'].value, \
                    opt_fields=fields)

            if line_.type == 'E':
                fields.pop('eid')
                fields.pop('sid1')
                fields.pop('sid2')
                fields.pop('beg1')
                fields.pop('end1')
                fields.pop('beg2')
                fields.pop('end2')
                fields.pop('alignment')

                beg1 = line_.fields['beg1'].value
                end1 = line_.fields['end1'].value
                beg2 = line_.fields['beg2'].value
                end2 = line_.fields['end2'].value
                is_dovetail_ = False
                #if (beg1 == "0" or end1[-1:] == "$") \
                #  and \
                #   (beg2 == "0" or end2[-1:] == "$"):
                #    is_dovetail_ = True

                if (beg1 == "0" and end2[-1:] == "$") \
                  or \
                   (beg2 == "0" and end1[-1:] == "$"):
                    is_dovetail_ = True

                    
                return Edge( \
                    line_.fields['eid'].value, \
                    line_.fields['sid1'].value[0:-1], \
                    line_.fields['sid1'].value[-1:], \
                    line_.fields['sid2'].value[0:-1], \
                    line_.fields['sid2'].value[-1:], \
                    (line_.fields['beg1'].value, line_.fields['end1'].value), \
                    (line_.fields['beg2'].value, line_.fields['end2'].value), \
                    line_.fields['alignment'].value, \
                    opt_fields=fields, \
                    is_dovetail=is_dovetail_)

            if line_.type == 'G':
                fields.pop('gid')
                fields.pop('sid1')
                fields.pop('sid2')
                fields.pop('distance')
                fields.pop('variance')
                return Edge( \
                    line_.fields['gid'].value, \
                    line_.fields['sid1'].value[0:-1], \
                    line_.fields['sid1'].value[-1:], \
                    line_.fields['sid2'].value[0:-1], \
                    line_.fields['sid2'].value[-1:], \
                    (None, None), \
                    (None, None), \
                    None, \
                    line_.fields['distance'].value, \
                    line_.fields['variance'].value, \
                    opt_fields=fields)

        except(KeyError, AttributeError):
            raise line.InvalidLineError("The given line cannot be "\
                                        + "an Edge.")


    def __eq__(self, other):
        try:
            if self.eid != other.eid or \
              self.from_node != other.from_node or \
              self.from_orn != other.from_orn or \
              self.to_node != other.to_node or \
              self.to_orn != other.to_orn or \
              self.from_positions != other.from_positions or \
              self.to_positions != other.to_positions or \
              self.alignment != other.alignment or \
              self.distance != other.distance or \
              self.variance != other.variance:
                return False
            for key, field in self.opt_fields.items():
                if not key in other.opt_fields or \
                  self.opt_fields[key] != other.opt_fields[key]:
                    return False
        except: return False
        return True


    def __neq__(self, other):
        return not self == other


    def __str__(self): # pragma: no cover
        fields = ("eid", "from_node", "from_orn", "to_node", "to_orn", \
                  "from_positions", "to_positions", "alignment", \
                  "distance", "variance", "opt_fields")

        opt_fields = []
        if len(self.opt_fields) > 0:
            opt_fields = str.join( \
                                   ",\t", \
                                   [str(field) for key, field in self.opt_fields.items()])
        values = (str(self.eid), \
                  str(self.from_node), str(self.from_orn), \
                  str(self.to_node), str(self.to_orn), \
                  str(self.from_positions), str(self.to_positions), \
                  str(self.alignment), str(self.distance), \
                  str(self.variance), "{" + str(opt_fields) + "}")
        assoc = [str.join(" : ", pair) for pair in zip(fields, values)]
        return str.join(",\t", assoc)


if __name__ == '__main__': # pragma: no cover
    pass
