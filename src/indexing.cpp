// ggCaller header
#include "ggCaller_classes.h"

// define mutex for safe addition to robinhood_maps
std::mutex mtx1;

ColoredCDBG<> buildGraph (const std::string& infile_1,
                          const std::string& infile_2,
                          const bool is_ref,
                          const int kmer,
                          const int threads,
                          const bool verb,
                          const bool write_graph,
                          const std::string& output_prefix)
{
    std::ifstream infile1(infile_1);
    std::ifstream infile2(infile_2);
    CCDBG_Build_opt opt;


    opt.k = kmer;
    opt.nb_threads = threads;
    opt.verbose = verb;
    opt.prefixFilenameOut = output_prefix;

    std::string filename;
    if (is_ref && (infile_2 == "NA")) {
        while (std::getline(infile1, filename))
        {
            opt.filename_ref_in.push_back(filename);
        }
    } else if (!is_ref && (infile_2 == "NA"))
    {
        while (std::getline(infile1, filename))
        {
            opt.filename_seq_in.push_back(filename);
        }
    } else {
        while (std::getline(infile1, filename))
        {
            opt.filename_ref_in.push_back(filename);
        }
        while (std::getline(infile2, filename))
        {
            opt.filename_seq_in.push_back(filename);
        }
    }

    ColoredCDBG<> ccdbg(opt.k);
    ccdbg.buildGraph(opt);
    ccdbg.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose);
    ccdbg.buildColors(opt);

    if (write_graph)
    {
        ccdbg.write(opt.prefixFilenameOut, opt.nb_threads, opt.verbose);
    }

    return ccdbg;
}

std::vector<std::size_t> findIndex(const std::string& seq,
                                   const std::string& subseq,
                                   const int start_index,
                                   const int overlap,
                                   const bool reverse)
{
    // search in forward direction from start index forwards
    if (!reverse)
    {
        std::vector<std::size_t> index_list;
        size_t pos = seq.find(subseq, start_index);
        while(pos != string::npos)
        {
            index_list.push_back(pos - overlap);
            pos = seq.find(subseq,pos+1);
        }
        return index_list;
    } else // search in reverse direction from start index backwards
    {
        std::vector<std::size_t> index_list;
        size_t pos = seq.rfind(subseq, start_index);
        while(pos != string::npos && pos != 0)
        {
            index_list.push_back(start_index - pos);
            pos = seq.rfind(subseq,pos - 1);
        }
        return index_list;
    }
}

// calculate stop codon frames for each reading frame
uint8_t calculateFrame_binary (const std::vector<std::size_t>& index_list)
{
    uint8_t codon_binary = 0;

    for (const auto& index : index_list)
    {
        uint8_t modulus = index % 3;
        codon_binary |= 1UL << modulus;
    }


    return codon_binary;
}

// switch codon reading frames depending on partial reading frame used
uint8_t switchFrame_binary (const uint8_t& binary_array, const int& frame)
{
    uint8_t codon_binary = binary_array;

    bool bit0 = codon_binary & (1 << 0);
    bool bit1 = codon_binary & (1 << 1);
    bool bit2 = codon_binary & (1 << 2);

    if (frame == 1)
    {
        codon_binary ^= (-bit2 ^ codon_binary) & (1UL << 0);
        codon_binary ^= (-bit0 ^ codon_binary) & (1UL << 1);
        codon_binary ^= (-bit1 ^ codon_binary) & (1UL << 2);

    } else if (frame == 2)
    {
        codon_binary ^= (-bit1 ^ codon_binary) & (1UL << 0);
        codon_binary ^= (-bit2 ^ codon_binary) & (1UL << 1);
        codon_binary ^= (-bit0 ^ codon_binary) & (1UL << 2);
    }
    bit0 = codon_binary & (1 << 0);
    bit1 = codon_binary & (1 << 1);
    bit2 = codon_binary & (1 << 2);

    return codon_binary;
}

template <class T, class U, bool is_const>
std::vector<bool> generate_colours(const UnitigMap<DataAccessor<T>, DataStorage<U>, is_const> unitig,
                                   const size_t& nb_colours,
                                   const size_t position)
{
    // get colours information for unitig
    const auto colourset = unitig.getData()->getUnitigColors(unitig);
    std::vector<bool> colours_arr(nb_colours, 0);

    // initialise a iterator, will only determine colours of single kmer, have to do for head and tail as may be different
    UnitigColors::const_iterator it_uc = colourset->begin(unitig);
    UnitigColors::const_iterator it_uc_end = colourset->end();

    // iterate over specified kmer colours, update colours array with presence of those colours
    for (it_uc; it_uc != it_uc_end; it_uc++)
    {
        if (it_uc.getKmerPosition() == position)
        {
            colours_arr[it_uc.getColorID()] = 1;
        } else if (it_uc.getKmerPosition() > position)
        {
           break;
        }
    }
    return colours_arr;
}

template<class T>
std::vector<std::pair<std::string, bool>> get_neighbours (const T& neighbour_iterator)
{
    std::vector<std::pair<std::string, bool>> neighbour_vector;
    for (const auto um : neighbour_iterator)
    {
        std::string head_kmer = um.getUnitigHead().toString();
        bool unitig_strand = um.strand;
        auto head_kmer_pair = std::make_pair(head_kmer, unitig_strand);
        neighbour_vector.push_back(std::move(head_kmer_pair));
    }
    return neighbour_vector;
}

template <class T, class U, bool is_const>
unitigDict analyse_unitigs_binary (const ColoredCDBG<>& ccdbg,
                                    UnitigMap<DataAccessor<T>, DataStorage<U>, is_const> um,
                                    const std::vector<std::string>& codon_for,
                                    const std::vector<std::string>& codon_rev,
                                    const int& kmer,
                                    const size_t& nb_colours)
{
    // initialise unitig_map
    unitigDict unitig_map;

    // generate string from unitig
    const std::string unitig = um.referenceUnitigToString();
    const size_t unitig_len = unitig.size();

    // add unitig sequence to unitig_map
    unitig_map.unitig_seq = unitig;

    // get head kmer for unitig, add to kmer dictionary
    const Kmer head_kmer_binary = um.getUnitigHead();
    const std::string head_kmer = head_kmer_binary.toString();
    unitig_map.add_head(head_kmer);

    // calculate unitig sizes
    unitig_map.add_size(unitig_len, unitig_len - (kmer-1));

    // calculate colours for unitig
    std::vector<bool> unitig_colours = generate_colours(um, nb_colours, 0);
    unitig_map.unitig_head_colour = std::move(unitig_colours);

    // generate colours for tail also
    const Kmer tail_kmer = um.getUnitigTail();
    const auto um_tail = ccdbg.find(tail_kmer, true);
    const size_t tail_pos = unitig_len - kmer;
    unitig_colours = generate_colours(um_tail, nb_colours, tail_pos);
    unitig_map.unitig_tail_colour = std::move(unitig_colours);

    // check if head and tail colours are equal
    if (unitig_map.unitig_head_colour != unitig_map.unitig_tail_colour)
    {
        unitig_map.head_tail_colours_equal = false;
        // calculate the full assortment of colours for a unitig
        unitig_map.unitig_full_colour = std::move(add_colours_array(unitig_map.unitig_head_colour, unitig_map.unitig_tail_colour));
    } else
    {
        // if colours match, just set to head colour
        unitig_map.unitig_full_colour = unitig_map.unitig_head_colour;
    }

    // generate successor head kmers (need to flip sign to get successors in reverse strand)
    unitig_map.succ_heads = std::move(get_neighbours(um.getSuccessors()));
    um.strand = !um.strand;
    unitig_map.pred_heads = std::move(get_neighbours(um.getSuccessors()));


    // analyse codon presence/absence
    //initialise codon search vectors
    std::vector<std::size_t> full_indices_pos;
    std::vector<std::size_t> part_indices_pos;
    std::vector<std::size_t> full_indices_neg;
    std::vector<std::size_t> part_indices_neg;

    std::vector<std::size_t> found_indices;

    // find the full and part indices of each of the codons in the unitig in positive strand
    for (const auto& codon : codon_for)
    {
        // full positive strand
        found_indices = findIndex(unitig, codon, 0, 0, false);
        full_indices_pos.insert(full_indices_pos.end(), make_move_iterator(found_indices.begin()), make_move_iterator(found_indices.end()));
        found_indices.clear();

        // part positive strand, ensure that unitig length is sufficient to have at least 1 codon length once overlap negated
        if (unitig_len > kmer + 2) {
            found_indices = findIndex(unitig, codon, kmer - 1, kmer - 1,false);
            part_indices_pos.insert(part_indices_pos.end(), make_move_iterator(found_indices.begin()), make_move_iterator(found_indices.end()));
            found_indices.clear();
        }
    }

    // find the full and part indices of each of the codons in the unitig in negative strand
    for (const auto& codon : codon_rev)
    {
        // full negative strand

        // issue here when actually find a match
        found_indices = findIndex(unitig, codon, unitig_len - 3, 0, true);
        full_indices_neg.insert(full_indices_neg.end(), make_move_iterator(found_indices.begin()), make_move_iterator(found_indices.end()));
        found_indices.clear();

        // part negative strand, ensure that unitig length is sufficient to have at least 1 codon length once overlap negated
        if (unitig_len > kmer + 2)
        {
            found_indices = findIndex(unitig, codon, (unitig_len - 3) - (kmer - 1), 0, true);
            part_indices_neg.insert(part_indices_neg.end(), make_move_iterator(found_indices.begin()), make_move_iterator(found_indices.end()));
            found_indices.clear();
        }
    }

    // insert frame mapping into unitig dictionaries, using head-kmer as the identifier
    // full unitig dictionary
    const uint8_t full_binary_pos = calculateFrame_binary(full_indices_pos);
    const uint8_t full_binary_neg = calculateFrame_binary(full_indices_neg);

    unitig_map.add_codon(true, true, 0, full_binary_pos);
    unitig_map.add_codon(true, false, 0, full_binary_neg);

    // part unitig dictionary
    const uint8_t part_binary_pos = calculateFrame_binary(part_indices_pos);
    const uint8_t part_binary_neg = calculateFrame_binary(part_indices_neg);

    unitig_map.add_codon(false, true, 0, part_binary_pos);
    unitig_map.add_codon(false, true, 1, switchFrame_binary(part_binary_pos, 1));
    unitig_map.add_codon(false, true, 2, switchFrame_binary(part_binary_pos, 2));

    unitig_map.add_codon(false, false, 0, part_binary_neg);
    unitig_map.add_codon(false, false, 1, switchFrame_binary(part_binary_neg, 1));
    unitig_map.add_codon(false, false, 2, switchFrame_binary(part_binary_neg, 2));

    return unitig_map;
}

GraphPair index_graph(const ColoredCDBG<>& ccdbg,
                       const std::vector<std::string>& stop_codons_for,
                       const std::vector<std::string>& stop_codons_rev,
                       const int kmer,
                       const size_t nb_colours)
{
    // get all head kmers for parrellelisation
    std::vector<Kmer> head_kmer_arr;
    for (const auto& um : ccdbg)
    {
        head_kmer_arr.push_back(um.getUnitigHead());
    }

    // structures for results
    UnitigVector graph_vector(head_kmer_arr.size());
    NodeColourVector node_colour_vector(nb_colours);
    robin_hood::unordered_map<std::string, size_t> head_kmer_map;

    // run unitig indexing in parallel
    size_t unitig_id = 1;
    #pragma omp parallel
    {
        UnitigVector graph_vector_private;
        NodeColourMap node_colour_vector_private(nb_colours);
        robin_hood::unordered_map<std::string, size_t> head_kmer_map_private;
        #pragma omp for nowait
        for (auto it = head_kmer_arr.begin(); it < head_kmer_arr.end(); it++)
        {
            // convert Kmer defined in *it to unitig
            auto unitig = ccdbg.find(*it, true);

            // generate results per unitig
            unitigDict unitig_map = std::move(analyse_unitigs_binary(ccdbg, unitig, stop_codons_for, stop_codons_rev, kmer, nb_colours));
            #pragma omp atomic capture
            unitig_map.unitig_id = unitig_id++;

            // add to node_colour_map_private
            for (size_t i = 0; i < unitig_map.unitig_full_colour.size(); i++)
            {
                if (unitig_map.unitig_full_colour.at(i))
                {
                    node_colour_vector_private[i].push_back(unitig_map.unitig_id)
                }
            }

            // add head_kmer and unitig id to map
            head_kmer_map_private[unitig_map.head_kmer] = unitig_map.unitig_id;

            // add unitig to graph_vector, minus 1 as zero based
            graph_vector[unitig_map.unitig_id - 1] = std::move(unitig_map));
        }
        #pragma omp critical
        {
            head_kmer_map.insert(head_kmer_map_private.begin(), head_kmer_map_private.end());

            // update node_colour_vector with calculated colours
            for (int i = 0; i < node_colour_vector_private.size(); i++)
            {
                node_colour_vector[i].insert(node_colour_vector[i].end(), std::make_move_iterator(node_colour_vector_private[i].begin()), std::make_move_iterator(node_colour_vector_private[i].end()));
            }
        }
    }
    // iterate over entries, determine correct successors/predecessors (i.e. have correct colours)
    size_t graph_map_size = graph_map.size() + 1;
    #pragma omp parallel
    {
        #pragma omp for nowait
        for (size_t i = 1; i < graph_map_size; i++)
        {
            // get a copy of the unitig_Dict object
            mtx1.lock();
            auto unitig_map = graph_map.at(i);
            mtx1.unlock();

            // initialise neighbours vectors
            unitig_map.neighbours.insert(std::pair<bool, NeighbourVector>(true, NeighbourVector()));
            unitig_map.neighbours.insert(std::pair<bool, NeighbourVector>(false, NeighbourVector()));

            // iterate over the connected nodes in successors
            for (const auto& succ : unitig_map.succ_heads)
            {
                // get numerical successor ID
                const size_t succ_id = head_kmer_map.at(succ.first);

                // get copy of successor unitig_ID
                mtx1.lock();
                auto adj_unitig_map = graph_map.at(succ_id);
                mtx1.unlock();

                // ensure colours are viable between the current and neigbouring unitig. Base this off head/tail colours depending on orientation
                // for current untiig, should negate tail, as this portion will overlap with next unitig
                std::vector<bool> colours;
                if (succ.second)
                {
                    colours = std::move(negate_colours_array(unitig_map.unitig_tail_colour, adj_unitig_map.unitig_head_colour));
                } else
                {
                    colours = std::move(negate_colours_array(unitig_map.unitig_tail_colour, adj_unitig_map.unitig_tail_colour));
                }

                // calculate sum_colours
                int sum_colours = accumulate(colours.begin(), colours.end(), 0);

                // if colours are viable, add successor information to current unitig
                if (sum_colours != 0)
                {
                    //generate integer value of successor ID, if negative strand ID will be negative etc.
                    int succ_id_int = (succ.second) ? succ_id : succ_id * -1;
                    std::pair<int, std::unordered_map<int, uint8_t>> neighbour (succ_id_int, adj_unitig_map.part_codon.at(succ.second));
                    unitig_map.neighbours[true].push_back(std::move(neighbour));
                }
            }

            // clear succ_heads
            unitig_map.succ_heads.clear();
            unitig_map.succ_heads.shrink_to_fit();

            // iterate over the connected nodes in predecessors
            for (const auto& pred : unitig_map.pred_heads)
            {
                // get numerical successor ID
                const size_t pred_id = head_kmer_map.at(pred.first);

                // get copy of successor unitig_ID
                mtx1.lock();
                auto adj_unitig_map = graph_map.at(pred_id);
                mtx1.unlock();

                // ensure colours are viable between the current and neigbouring unitig. Base this off head/tail colours depending on orientation
                // for current untiig, should negate head, as this portion will overlap with next unitig
                std::vector<bool> colours;
                if (pred.second)
                {
                    colours = std::move(negate_colours_array(unitig_map.unitig_head_colour, adj_unitig_map.unitig_head_colour));
                } else
                {
                    colours = std::move(negate_colours_array(unitig_map.unitig_head_colour, adj_unitig_map.unitig_tail_colour));
                }

                // calculate sum_colours
                int sum_colours = accumulate(colours.begin(), colours.end(), 0);

                // if colours are viable, add successor information to current unitig
                if (sum_colours != 0)
                {
                    //generate integer value of successor ID, if negative strand ID will be negative etc.
                    int pred_id_int = (pred.second) ? pred_id : pred_id * -1;
                    std::pair<int, std::unordered_map<int, uint8_t>> neighbour (pred_id_int, adj_unitig_map.part_codon.at(pred.second));
                    unitig_map.neighbours[false].push_back(std::move(neighbour));
                }
            }

            // clear pred_heads
            unitig_map.pred_heads.clear();
            unitig_map.pred_heads.shrink_to_fit();

            // if there are no successors/predecessors or untig colours change, means that unitig is first/last in sequence. Therefore,
            // update full codon array to be 3 frame stop index to enable complete traversal
            if (unitig_map.neighbours.at(true).empty() || unitig_map.neighbours.at(false).empty() || !unitig_map.head_tail_colours_equal)
            {
                // set first three bits to 1
                uint8_t full_binary = 7;


                // update full codon indexes
                unitig_map.add_codon(true, true, 0, full_binary);
                unitig_map.add_codon(true, false, 0, full_binary);

                unitig_map.forward_stop = true;
                unitig_map.reverse_stop = true;

                // set unitig to end of contig
                unitig_map.end_contig = true;
            }

            // update graph_map with new entry
            mtx1.lock();
            graph_map[i] = std::move(unitig_map);
            mtx1.unlock();
        }
    }

    const auto graph_pair = std::make_pair(graph_vector, node_colour_vector);
    return graph_pair;
}