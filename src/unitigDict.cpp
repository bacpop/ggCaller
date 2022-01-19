// ggCaller header
#include "unitigDict.h"
#include "indexing.h"

std::mutex mtx2;

// Clear method for CompactedDBG
void MyUnitigMap::clear(const UnitigMap<MyUnitigMap>& um_dest){

}

// Clear method for CompactedDBG
void MyUnitigMap::clear(const UnitigColorMap<MyUnitigMap>& um_dest){
    return;
}

// Concatenation method for ColoredCDBG
void MyUnitigMap::concat(const UnitigColorMap<MyUnitigMap>& um_dest, const UnitigColorMap<MyUnitigMap>& um_src){
    return;
}

// Concatenation method for ColoredCDBG
void MyUnitigMap::concat(const UnitigMap<MyUnitigMap>& um_dest, const UnitigMap<MyUnitigMap>& um_src){
    return;
}

// Extraction method for ColoredCDBG
void MyUnitigMap::extract(const UnitigColors* uc_dest, const UnitigColorMap<MyUnitigMap>& um_src, const bool last_extraction){
    return;
}

// Extraction method for ColoredCDBG
void MyUnitigMap::extract(const UnitigMap<MyUnitigMap>& um_src, bool last_extraction){
    return;
}

// add codon, copy semantics
void MyUnitigMap::add_codon (const bool& full, const bool& forward, const int& frame, const uint8_t& array) {
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

//// add size, copy semantics
//void MyUnitigMap::add_size(const size_t& full_len, const size_t& part_len) {
//    std::pair pair(full_len, part_len);
//    _unitig_size = pair;
//}
//
//// add size, move semantics
//void MyUnitigMap::add_size(size_t& full_len, size_t& part_len) {
//    _unitig_size = make_pair(full_len, part_len);
//}

// check if head and tail colours are the same
void MyUnitigMap::_check_head_tail_equal() {
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

uint8_t MyUnitigMap::get_codon_arr (bool full, bool forward, bool frame) const {
    if (full)
    {
        return _full_codon.at(forward).at(frame);
    } else
    {
        return _part_codon.at(forward).at(frame);
    }
}

std::vector<uint8_t> MyUnitigMap::get_codon_dict (bool full, bool forward) const {
    if (full)
    {
        return _full_codon.at(forward);
    } else
    {
        return _part_codon.at(forward);
    }
}

void MyUnitigMap::add_neighbour_colour (bool strand, int neighbour_ID, size_t colour_ID)
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

// function to get unitig data from an ID number, returning pointer to data
std::pair<UnitigColorMap<MyUnitigMap>, MyUnitigMap*> get_um_data (ColoredCDBG<MyUnitigMap>& ccdbg,
                                                                  const std::vector<Kmer>& head_kmer_arr,
                                                                  const int id)
{
    const Kmer head_kmer = head_kmer_arr.at(abs(id) - 1);
    auto um = ccdbg.find(head_kmer, true);

    DataAccessor<MyUnitigMap>* da = um.getData();
    MyUnitigMap* um_data = da->getData(um);

    return {um, um_data};
}


// const qualified get_um_data, returning copy of data
std::pair<const_UnitigMap<DataAccessor<MyUnitigMap>, DataStorage<MyUnitigMap>>, MyUnitigMap*> get_um_data (const ColoredCDBG<MyUnitigMap>& ccdbg,
                                                                                                           const std::vector<Kmer>& head_kmer_arr,
                                                                                                           const int id)
{
    const Kmer head_kmer = head_kmer_arr.at(abs(id) - 1);
    auto um = ccdbg.find(head_kmer, true);

    const DataAccessor<MyUnitigMap>* da = um.getData();
    const MyUnitigMap* unitig_map = da->getData(um);

    return {um, const_cast<MyUnitigMap *>(unitig_map)};
}


// function to get unitig data from an head_kmer, returning pointer to data
std::pair<UnitigColorMap<MyUnitigMap>, MyUnitigMap*> get_um_data (ColoredCDBG<MyUnitigMap>& ccdbg,
                                                                  const Kmer& head_kmer)
{
    auto um = ccdbg.find(head_kmer, true);

    DataAccessor<MyUnitigMap>* da = um.getData();
    MyUnitigMap* um_data = da->getData(um);

    return {um, um_data};
}

// const qualified get_um_data, returning copy of data
std::pair<const_UnitigMap<DataAccessor<MyUnitigMap>, DataStorage<MyUnitigMap>>, MyUnitigMap*> get_um_data (const ColoredCDBG<MyUnitigMap>& ccdbg,
                                                                                                           const Kmer& head_kmer)
{
    auto um = ccdbg.find(head_kmer, true);

    const DataAccessor<MyUnitigMap>* da = um.getData();
    const MyUnitigMap* unitig_map = da->getData(um);

    return {um, const_cast<MyUnitigMap *>(unitig_map)};
}

// retrieve id from k-mer mapping
size_t get_id (const ColoredCDBG<MyUnitigMap>& ccdbg,
               const Kmer& head_kmer)
{
    auto um = ccdbg.find(head_kmer, true);
    auto da = um.getData(); // Get DataAccessor from unitig
    const MyUnitigMap* data = da->getData(um); // Get boolean from DataAccessor

    return data->id;
}

// non-member function for generating sequence from DBG
std::string generate_sequence_nm(const std::vector<int>& nodelist,
                                 const std::vector<indexPair>& node_coords,
                                 const size_t& overlap,
                                 const ColoredCDBG<MyUnitigMap>& ccdbg,
                                 const std::vector<Kmer>& head_kmer_arr)
{
    std::string sequence;
    for (size_t i = 0; i < nodelist.size(); i++)
    {
        // parse information
        const auto& id = nodelist[i];
        const auto& coords = node_coords[i];

        // initialise sequence items
        std::string substring;

        // get a reference to the unitig map object
        auto um_pair = get_um_data(ccdbg, head_kmer_arr, id);
        auto& um = um_pair.first;

        // reverse sequence if strand is negative
        std::string seq;
        if (id >= 0)
        {
            seq = um.referenceUnitigToString();
        } else
        {
            seq = reverse_complement(um.referenceUnitigToString());
        }

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