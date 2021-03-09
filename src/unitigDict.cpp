// ggCaller header
#include "ggCaller_classes.h"

// add codon, copy semantics
void unitigDict::add_codon (const bool& full, const bool& forward, const int& frame, const uint8_t& array) {
    if (full)
    {
        full_codon[forward][frame] = array;

        // update forward/reverse stop codon presence
        if (forward && !forward_stop_defined)
        {
            if (full_codon[forward][frame] != 0)
            {
                forward_stop = true;
            }
            forward_stop_defined = true;
        } else if (!forward && !reverse_stop_defined)
        {
            if (full_codon[forward][frame] != 0)
            {
                reverse_stop = true;
            }
            reverse_stop_defined = true;
        }
    } else {
        part_codon[forward][frame] = array;
    }
}

// add codon, move semantics
void unitigDict::add_codon (const bool& full, const bool& forward, const int& frame, uint8_t& array) {
    if (full)
    {
        full_codon[forward][frame] = std::move(array);

        // update forward/reverse stop codon presence
        if (forward && !forward_stop_defined)
        {
            if (full_codon[forward][frame] != 0)
            {
                forward_stop = true;
            }
            forward_stop_defined = true;
        } else if (!forward && !reverse_stop_defined)
        {
            if (full_codon[forward][frame] != 0)
            {
                reverse_stop = true;
            }
            reverse_stop_defined = true;
        }
    } else {
        part_codon[forward][frame] = std::move(array);
    }
}

// add size, copy semantics
void unitigDict::add_size(const size_t& full_len, const size_t& part_len) {
    std::pair pair(full_len, part_len);
    unitig_size = pair;
}

// add size, move semantics
void unitigDict::add_size(size_t& full_len, size_t& part_len) {
    unitig_size = make_pair(full_len, part_len);
}
