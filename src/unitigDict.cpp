// ggCaller header
#include "unitigDict.h"
#include "indexing.h"

// add codon, copy semantics
void unitigDict::add_codon (const bool& full, const bool& forward, const int& frame, const uint8_t& array) {
    if (full)
    {
        _full_codon[forward][frame] = array;

        // update forward/reverse stop codon presence
        if (forward && !_forward_stop_defined)
        {
            if (_full_codon[forward][frame] != 0)
            {
                _forward_stop = true;
            }
            _forward_stop_defined = true;
        } else if (!forward && !_reverse_stop_defined)
        {
            if (_full_codon[forward][frame] != 0)
            {
                _reverse_stop = true;
            }
            _reverse_stop_defined = true;
        }
    } else {
        _part_codon[forward][frame] = array;
    }
}

// add size, copy semantics
void unitigDict::add_size(const size_t& full_len, const size_t& part_len) {
    std::pair pair(full_len, part_len);
    _unitig_size = pair;
}

// add size, move semantics
void unitigDict::add_size(size_t& full_len, size_t& part_len) {
    _unitig_size = make_pair(full_len, part_len);
}

// check if head and tail colours are the same
void unitigDict::_check_head_tail_equal() {
    if (_unitig_head_colour == _unitig_tail_colour)
    {
        _head_tail_colours_equal = true;
        _unitig_full_colour = _unitig_head_colour;
    }
    else
    {
        _head_tail_colours_equal = false;
        _end_contig = true;
        _unitig_full_colour = bool_or(_unitig_head_colour, _unitig_tail_colour);
    }
}

uint8_t unitigDict::get_codon_arr (bool full, bool forward, bool frame) const {
    if (full)
    {
        return _full_codon.at(forward).at(frame);
    } else
    {
        return _part_codon.at(forward).at(frame);
    }
}

std::vector<uint8_t> unitigDict::get_codon_dict (bool full, bool forward) const {
    if (full)
    {
        return _full_codon.at(forward);
    } else
    {
        return _part_codon.at(forward);
    }
}
