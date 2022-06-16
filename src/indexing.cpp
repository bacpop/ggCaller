// ggCaller header
#include "indexing.h"

ColoredCDBG<MyUnitigMap> buildGraph (const std::string& infile_1,
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

    ColoredCDBG<MyUnitigMap> ccdbg(opt.k);
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
std::bitset<3> calculateFrame_binary_full (const std::vector<std::size_t>& index_list)
{
    std::bitset<3> binary_array;

    for (const auto& index : index_list)
    {
        size_t modulus = index % 3;
        binary_array[modulus] = 1;
    }

    return binary_array;
}

std::bitset<9> calculateFrame_binary_part (const std::vector<std::size_t>& index_list)
{
    std::bitset<3> binary_array = calculateFrame_binary_full(index_list);

    std::bitset<9> full_binary;

    bool bit0 = binary_array[0];
    bool bit1 = binary_array[1];
    bool bit2 = binary_array[2];

    // go through binary_array, setting for each frame
    full_binary[0] = bit0;
    full_binary[1] = bit1;
    full_binary[2] = bit2;
    full_binary[3] = bit2;
    full_binary[4] = bit0;
    full_binary[5] = bit1;
    full_binary[6] = bit1;
    full_binary[7] = bit2;
    full_binary[8] = bit0;

    return full_binary;
}

template <class T, class U, bool is_const>
boost::dynamic_bitset<> generate_colours(const UnitigMap<DataAccessor<T>, DataStorage<U>, is_const> unitig,
                                           const size_t nb_colours,
                                           const size_t position)
{
    // get colours information for unitig
    const auto colourset = unitig.getData()->getUnitigColors(unitig);
    boost::dynamic_bitset<> colours_arr(nb_colours);

    // initialise a iterator, will only determine colours of single kmer, have to do for head and tail as may be different
    UnitigColors::const_iterator it_uc = colourset->begin(unitig);
    UnitigColors::const_iterator it_uc_end = colourset->end();

    // iterate over specified kmer colours, update colours array with presence of those colours
    for (it_uc; it_uc != it_uc_end; it_uc++)
    {
        if (it_uc.getKmerPosition() == position)
        {
            // check here, seems to be issue where colours are greater than nb_colours for given kmer
            const auto col = it_uc.getColorID();
            if (col < nb_colours)
            {
                colours_arr[col] = 1;
            }
        } else if (it_uc.getKmerPosition() > position)
        {
           break;
        }
    }

    return colours_arr;
}

template<class T>
std::vector<std::pair<Kmer, bool>> get_neighbours (const T& neighbour_iterator)
{
    std::vector<std::pair<Kmer, bool>> neighbour_vector;
    for (const auto um : neighbour_iterator)
    {
        Kmer head_kmer = um.getUnitigHead();
        bool unitig_strand = um.strand;
        auto head_kmer_pair = std::make_pair(head_kmer, unitig_strand);
        neighbour_vector.push_back(std::move(head_kmer_pair));
    }
    return neighbour_vector;
}

template <class T, class U, bool is_const>
void analyse_unitigs_binary (ColoredCDBG<MyUnitigMap>& ccdbg,
                            UnitigMap<DataAccessor<T>, DataStorage<U>, is_const> um,
                            const std::vector<std::string>& codon_for,
                            const std::vector<std::string>& codon_rev,
                            const int& kmer,
                            const size_t& nb_colours)
{
    // initialise unitig_dict
    DataAccessor<MyUnitigMap>* da = um.getData();
    MyUnitigMap* um_data = da->getData(um);

    // generate string from unitig
    const std::string unitig = um.referenceUnitigToString();
    const size_t unitig_len = unitig.size();

    // get head kmer for unitig, add to kmer dictionary
    const Kmer head_kmer_binary = um.getUnitigHead();
    const std::string head_kmer = head_kmer_binary.toString();


    // calculate colours for unitig
    boost::dynamic_bitset<> full_unitig_colour = generate_colours(um, nb_colours, 0);

    // generate colours for tail also
    const Kmer tail_kmer = um.getUnitigTail();
    const auto um_tail = ccdbg.find(tail_kmer, true);
    const size_t tail_pos = unitig_len - kmer;
    full_unitig_colour |= generate_colours(um_tail, nb_colours, tail_pos);
    um_data->add_full_colour(full_unitig_colour);

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
    const auto full_binary_pos = calculateFrame_binary_full(full_indices_pos);
    const auto full_binary_neg = calculateFrame_binary_full(full_indices_neg);

    um_data->add_codon(true, full_binary_pos);
    um_data->add_codon(false, full_binary_neg);

    // part unitig dictionary
    const auto part_binary_pos = calculateFrame_binary_part(part_indices_pos);
    const auto part_binary_neg = calculateFrame_binary_part(part_indices_neg);

    um_data->add_codon(true, part_binary_pos);
    um_data->add_codon(false, part_binary_neg);
}

NodeColourVector index_graph(std::vector<Kmer>& head_kmer_arr,
                             ColoredCDBG<MyUnitigMap>& ccdbg,
                             const std::vector<std::string>& stop_codons_for,
                             const std::vector<std::string>& stop_codons_rev,
                             const int kmer,
                             const size_t nb_colours,
                             const std::vector<std::string>& input_colours,
                             const boost::dynamic_bitset<>& ref_set)
{
    // get all head kmers and add IDs
    size_t unitig_id = 1;
    for (auto& um : ccdbg)
    {
        head_kmer_arr.push_back(um.getUnitigHead());
        DataAccessor<MyUnitigMap>* da = um.getData();
        MyUnitigMap* um_data = da->getData(um);
        um_data->set_id(unitig_id++);
    }

    NodeColourVector node_colour_vector(nb_colours);

    // run unitig indexing in parallel
    #pragma omp parallel
    {
        NodeColourVector node_colour_vector_private(nb_colours);
        std::unordered_map<std::string, size_t> head_kmer_map_private;
        #pragma omp for
        for (auto it = head_kmer_arr.begin(); it < head_kmer_arr.end(); it++)
        {
            // convert Kmer defined in *it to unitig
            UnitigColorMap<MyUnitigMap> um = ccdbg.find(*it, true);

            // generate results per unitig
            analyse_unitigs_binary(ccdbg, um, stop_codons_for, stop_codons_rev, kmer, nb_colours);
            DataAccessor<MyUnitigMap>* da = um.getData();
            MyUnitigMap* um_data = da->getData(um);

            // add to node_colour_map_private
            for (size_t i = 0; i < um_data->full_colour().size(); i++)
            {
                if (um_data->full_colour()[i])
                {
                    node_colour_vector_private[i].push_back(um_data->get_id());
                }
            }
        }
        #pragma omp critical
        {
            // update node_colour_vector with calculated colours
            for (int i = 0; i < node_colour_vector_private.size(); i++)
            {
                node_colour_vector[i].insert(node_colour_vector[i].end(), make_move_iterator(node_colour_vector_private[i].begin()), make_move_iterator(node_colour_vector_private[i].end()));
            }
        }
    }

    // return node_colour vector
    return node_colour_vector;
}