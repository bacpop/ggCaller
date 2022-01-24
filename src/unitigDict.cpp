// ggCaller header
#include "unitigDict.h"
#include "indexing.h"

std::mutex mtx2;


void MyUnitigMap::clear(const UnitigColorMap<MyUnitigMap>& um_dest)
{
    if (!um_dest.isEmpty && (um_dest.getGraph() != nullptr)){

        DataStorage<MyUnitigMap>* ds = um_dest.getGraph()->getData();

        if (ds != nullptr){

            ds->remove(um_dest);

            da_id = 0;
        }
    }
}

// Concatenation method for ColoredCDBG
void MyUnitigMap::concat(const UnitigColorMap<MyUnitigMap>& um_dest, const UnitigColorMap<MyUnitigMap>& um_src){
    return;
}

// Extraction method for ColoredCDBG
void MyUnitigMap::extract(const UnitigColors* uc_dest, const UnitigColorMap<MyUnitigMap>& um_src, const bool last_extraction){
    return;
}

string MyUnitigMap::serialize(const const_UnitigColorMap<MyUnitigMap>& um_src) const {
    return "";
}

// add codon, copy semantics
void MyUnitigMap::add_codon (const bool forward, const std::bitset<3>& array) {
    // skip if no stop detected
    if (array.none())
    {
        return;
    }

    // set first position of bits
    size_t a = 0;

    // if reverse, then set to second half of the bitset (3 for full, 9 for part)
    if (!forward)
    {
        a = 3;
    }

    // add to a for the remaining bits
    size_t b = a + 1;
    size_t c = a + 2;

    _full_codon[a] = array[0];
    _full_codon[b] = array[1];
    _full_codon[c] = array[2];

    // update forward/reverse stop codon presence
    if (forward && !_forward_stop)
    {
        _forward_stop = true;
    } else if (!forward && !_reverse_stop)
    {
        _reverse_stop = true;
    }
}

void MyUnitigMap::add_codon (const bool forward, const std::bitset<9>& array) {
    // skip if no stop detected
    if (array.none())
    {
        return;
    }

    // set first position of bits
    size_t a = 0;

    // if reverse, then set to second half of the bitset (3 for full, 9 for part)
    if (!forward)
    {
        a = 9;
    }

    // add to a for the remaining bits
    size_t b = a + 1;
    size_t c = a + 2;
    size_t d = a + 3;
    size_t e = a + 4;
    size_t f = a + 5;
    size_t g = a + 6;
    size_t h = a + 7;
    size_t i = a + 8;

    _part_codon[a] = array[0];
    _part_codon[b] = array[1];
    _part_codon[c] = array[2];
    _part_codon[d] = array[3];
    _part_codon[e] = array[4];
    _part_codon[f] = array[5];
    _part_codon[g] = array[6];
    _part_codon[h] = array[7];
    _part_codon[i] = array[8];
}

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

    // clear the head and tail colours
    _unitig_head_colour.clear();
    _unitig_tail_colour.clear();
}

std::bitset<3> MyUnitigMap::get_codon_arr (const bool full, const bool forward, const int frame) const {
    // set first position of bits
    size_t a = 0;

    // if reverse, then set to second half of the bitset (3 for full, 9 for part)
    if (!forward)
    {
        if (full)
        {
            a = 3;
        } else
        {
            a = 9;
        }
    }

    // determine the frame of the bits
    a += 3 * frame;

    // add to a for the remaining bits
    size_t b = a + 1;
    size_t c = a + 2;

    std::bitset<3> array;

    if (full)
    {
        array[0] = _full_codon[a];
        array[1] = _full_codon[b];
        array[2] = _full_codon[c];

    } else
    {
        array[0] = _part_codon[a];
        array[1] = _part_codon[b];
        array[2] = _part_codon[c];
    }

    return array;
}

//std::vector<uint8_t> MyUnitigMap::get_codon_dict (bool full, bool forward) const {
//    if (full)
//    {
//        return _full_codon.at(forward);
//    } else
//    {
//        return _part_codon.at(forward);
//    }
//}

//void MyUnitigMap::add_neighbour_colour (bool strand, int neighbour_ID, size_t colour_ID)
//{
//    // search for neighbour in _neighbours, add colour_ID to colour set
//    auto it = find_if(_neighbours[strand].begin( ), _neighbours[strand].end( ), [=](auto item){return std::get< 0 >(item) == neighbour_ID;});
//    if (it != _neighbours[strand].end())
//    {
//        mtx2.lock();
//        std::get<2>(*it).insert(colour_ID);
//        mtx2.unlock();
//    }
//}

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

//// retrieve id from k-mer mapping
//size_t get_id (const ColoredCDBG<MyUnitigMap>& ccdbg,
//               const Kmer& head_kmer)
//{
//    auto um = ccdbg.find(head_kmer, true);
//    auto da = um.getData(); // Get DataAccessor from unitig
//    const MyUnitigMap* data = da->getData(um); // Get boolean from DataAccessor
//
//    return data->get_id();
//}

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

        const std::string kmer = head_kmer_arr[abs(id) - 1].toString();

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