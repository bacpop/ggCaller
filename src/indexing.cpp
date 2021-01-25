// ggCaller header
#include "ggCaller_classes.h"

ColoredCDBG<> buildGraph (const std::string& infile_1,
                          const std::string& infile_2,
                          const bool& is_ref,
                          const int kmer,
                          const int threads,
                          const bool verb,
                          const bool& write_graph,
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
                                   const size_t& nb_colours)
{
    // get colours information for unitig
    const auto colourset = unitig.getData()->getUnitigColors(unitig);
    std::vector<bool> colours_arr(nb_colours, 0);

    // initialise a iterator, will only iterate over colours of head kmer (this will be same for whole unitig)
    UnitigColors::const_iterator it_uc = colourset->begin(unitig);
    UnitigColors::const_iterator it_uc_end = colourset->end();

    // iterate over kmer head position only, update colours array with presence of those colours
    for (it_uc; it_uc != it_uc_end; it_uc++)
    {
        if (it_uc.getKmerPosition() == 0)
        {
            colours_arr[it_uc.getColorID()] = 1;
        }
    }
    return colours_arr;
}

template <class T, class U, bool is_const>
unitigDict analyse_unitigs_binary (const ColoredCDBG<>& ccdbg,
                            const UnitigMap<DataAccessor<T>, DataStorage<U>, is_const> um,
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

    // get head kmer for unitig, add to kmer dictionary
    const Kmer head_kmer_binary = um.getUnitigHead();
    const std::string head_kmer = head_kmer_binary.toString();
    unitig_map.add_head(head_kmer);

    //cout << "Calculating unitig lengths" << endl;

    // calculate unitig sizes
    unitig_map.add_size(unitig_len, unitig_len - (kmer-1));

    //cout << "Calculating colours" << endl;

    // calculate colours for unitig
    std::vector<bool> unitig_colours = generate_colours(um, nb_colours);
    unitig_map.add_colour(unitig_colours);

    // analyse codon presence/absence
    //initialise codon search vectors
    std::vector<std::size_t> full_indices_pos;
    std::vector<std::size_t> part_indices_pos;
    std::vector<std::size_t> full_indices_neg;
    std::vector<std::size_t> part_indices_neg;

    std::vector<std::size_t> found_indices;

    // find the full and part indices of each of the codons in the unitig in positive strand
    //cout << "Searching positive strands..." << endl;
    for (const auto& codon : codon_for)
    {
        //cout << "Identifying codon presence: " << codon << endl;
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
        //cout << "Finished." << endl;
    }

    //cout << "Searching negative strands..." << endl;
    // find the full and part indices of each of the codons in the unitig in negative strand
    for (const auto& codon : codon_rev)
    {
        //cout << "Identifying codon presence: " << codon << endl;
        // full negative strand
        //cout << "Full unitig analysis..." << endl;

        // issue here when actually find a match
        found_indices = findIndex(unitig, codon, unitig_len - 3, 0, true);
        full_indices_neg.insert(full_indices_neg.end(), make_move_iterator(found_indices.begin()), make_move_iterator(found_indices.end()));
        found_indices.clear();

        //cout << "Part unitig analysis..." << endl;
        // part negative strand, ensure that unitig length is sufficient to have at least 1 codon length once overlap negated
        if (unitig_len > kmer + 2)
        {
            found_indices = findIndex(unitig, codon, (unitig_len - 3) - (kmer - 1), 0, true);
            part_indices_neg.insert(part_indices_neg.end(), make_move_iterator(found_indices.begin()), make_move_iterator(found_indices.end()));
            found_indices.clear();
        }
        //cout << "Finished." << endl;
    }

    // insert frame mapping into unitig dictionaries, using head-kmer as the identifier
    // full unitig dictionary
    const uint8_t full_binary_pos = calculateFrame_binary(full_indices_pos);
    const uint8_t full_binary_neg = calculateFrame_binary(full_indices_neg);

    unitig_map.add_codon(true, true, 0, full_binary_pos);
    unitig_map.add_codon(true, true, 1, switchFrame_binary(full_binary_pos, 1));
    unitig_map.add_codon(true, true, 2, switchFrame_binary(full_binary_pos, 2));

    unitig_map.add_codon(true, false, 0, full_binary_neg);
    unitig_map.add_codon(true, false, 1, switchFrame_binary(full_binary_neg, 1));
    unitig_map.add_codon(true, false, 2, switchFrame_binary(full_binary_neg, 2));

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

GraphTuple index_graph(const ColoredCDBG<>& ccdbg,
                       const std::vector<std::string>& stop_codons_for,
                       const std::vector<std::string>& stop_codons_rev,
                       const int& kmer,
                       const size_t& nb_colours)
{
    // get all head kmers for parrellelisation
    std::vector<Kmer> head_kmer_arr;
    for (const auto um : ccdbg)
    {
        head_kmer_arr.push_back(um.getUnitigHead());
    }

    // structures for results
    robin_hood::unordered_map<std::string, unitigDict> graph_map;
    std::vector<std::string> stop_list_for;
    std::vector<std::string> stop_list_rev;

    // run unitig indexing in parallel
    size_t unitig_id = 0;
    #pragma omp parallel
    {
        robin_hood::unordered_map<std::string, unitigDict> graph_map_private;
        std::vector<std::string> stop_list_for_private;
        std::vector<std::string> stop_list_rev_private;
        #pragma omp for nowait
        for (auto it = head_kmer_arr.begin(); it < head_kmer_arr.end(); it++)
        {
            // convert Kmer defined in *it to unitig
            auto unitig = ccdbg.find(*it, true);
            auto unitig_str = unitig.referenceUnitigToString();

            // generate results per unitig
            unitigDict unitig_map = std::move(analyse_unitigs_binary(ccdbg, unitig, stop_codons_for, stop_codons_rev, kmer, nb_colours));
            #pragma omp atomic capture
            unitig_map.unitig_id = unitig_id++;

            // add results to private maps and vectors
            if (unitig_map.forward_stop)
            {
                stop_list_for_private.push_back(unitig_map.head_kmer);
            }
            if (unitig_map.reverse_stop)
            {
                stop_list_rev_private.push_back(unitig_map.head_kmer);
            }
            graph_map_private[unitig_map.head_kmer] = std::move(unitig_map);
        }
        #pragma omp critical
        {
            graph_map.insert(graph_map_private.begin(), graph_map_private.end());
            stop_list_for.insert(stop_list_for.end(), std::make_move_iterator(stop_list_for_private.begin()), std::make_move_iterator(stop_list_for_private.end()));
            stop_list_rev.insert(stop_list_rev.end(), std::make_move_iterator(stop_list_rev_private.begin()), std::make_move_iterator(stop_list_rev_private.end()));
        }
    }

    const auto graph_tuple = std::make_tuple(graph_map, stop_list_for, stop_list_rev);
    return graph_tuple;
}