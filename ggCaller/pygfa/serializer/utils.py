from ggCaller.pygfa.graph_element.parser import line, field_validator as fv

SERIALIZATION_ERROR_MESSAGGE = "Couldn't serialize object identified by: "

def _format_exception(identifier, exception):
    return SERIALIZATION_ERROR_MESSAGGE + identifier \
      + "\n\t" + repr(exception)

def _remove_common_edge_fields(edge_dict):
    edge_dict.pop('eid')
    edge_dict.pop('from_node')
    edge_dict.pop('from_orn')
    edge_dict.pop('to_node')
    edge_dict.pop('to_orn')
    edge_dict.pop('from_positions')
    edge_dict.pop('to_positions')
    edge_dict.pop('alignment')
    edge_dict.pop('distance')
    edge_dict.pop('variance')

def _serialize_opt_fields(opt_fields):
    fields = []
    for key, opt_field in opt_fields.items():
        if line.is_optfield(opt_field):
            fields.append(str(opt_field))
    return fields

def _are_fields_defined(fields):
    try:
        for field in fields:
            if field is None:
                return False
    except:
        return False
    return True

def _check_fields(fields, required_fields):
    """Check if each field has the correct format as
    stated from the specification.
    """
    try:
        for field in range(0, len(required_fields)):
            if not fv.is_valid(fields[field], required_fields[field]):
                return False
        return True
    except:
        return False

def _check_identifier(identifier):
    if not isinstance(identifier, str):
        identifier = "'{0}' - id of type {1}.".format(\
                                            str(identifier), \
                                            type(identifier) \
                                            )
    return identifier
