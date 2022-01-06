// ggCaller header
#include "unitigDict.h"
#include "indexing.h"

std::mutex mtx2;

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
        _unitig_full_colour = _unitig_head_colour;
        _unitig_full_colour |= _unitig_tail_colour;
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

void unitigDict::add_neighbour_colour (bool strand, int neighbour_ID, size_t colour_ID)
{
    // search for neighbour in _neighbours, add colour_ID to colour set
    auto it = find_if(_neighbours[strand].begin( ), _neighbours[strand].end( ), [=](auto item){return get< 0 >(item) == neighbour_ID;});
    if (it != _neighbours[strand].end())
    {
        mtx2.lock();
        get<2>(*it).insert(colour_ID);
        mtx2.unlock();
    }
}

// non-member function for generating unitig sequence
std::string unitig_seq(const int& node_id,
                       const GraphVector& graph_vector,
                       const ColoredCDBG<>& ccdbg)
{
    std::string head_kmer = graph_vector.at(abs(node_id) - 1).head_kmer();

    const Kmer km = Kmer(head_kmer.c_str());

    auto um = ccdbg.find(km, true);

    std::string seq;

    if (node_id < 0)
    {
        seq = reverse_complement(um.referenceUnitigToString());
    } else
    {
        seq = um.referenceUnitigToString();
    }

    return seq;
}

// non-member function for generating sequence from DBG
std::string generate_sequence_nm(const std::vector<int>& nodelist,
                                 const std::vector<indexPair>& node_coords,
                                 const size_t& overlap,
                                 const GraphVector& graph_vector,
                                 const ColoredCDBG<>& ccdbg)
{
    std::string sequence;
    for (size_t i = 0; i < nodelist.size(); i++)
    {
        // parse information
        const auto& id = nodelist[i];
        const auto& coords = node_coords[i];

        // initialise sequence items
        std::string substring;
        const std::string seq = unitig_seq(id, graph_vector, ccdbg);

        if (sequence.empty())
        {
            // get node_seq_len, add one as zero indexed
            int node_seq_len = (std::get<1>(coords) - std::get<0>(coords)) + 1;
            substring = seq.substr(std::get<0>(coords), node_seq_len);
        } else
        {
            // get node_seq_len, add one as zero indexed
            int node_seq_len = (std::get<1>(coords) - overlap) + 1;
            // need to account for overlap, if overlap is greater than the end of the node, sequence already accounted for
            if (node_seq_len > 0)
            {
                substring = seq.substr(overlap, node_seq_len);
            } else
            {
                break;
            }
        }
        sequence += substring;
    }
    return sequence;
}