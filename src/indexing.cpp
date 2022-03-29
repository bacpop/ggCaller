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
    boost::dynamic_bitset<> unitig_colours = generate_colours(um, nb_colours, 0);
    um_data->add_head_colour(unitig_colours);

    // generate colours for tail also
    const Kmer tail_kmer = um.getUnitigTail();
    const auto um_tail = ccdbg.find(tail_kmer, true);
    const size_t tail_pos = unitig_len - kmer;
    unitig_colours = generate_colours(um_tail, nb_colours, tail_pos);
    um_data->add_tail_colour(unitig_colours);
    
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

void update_neighbour_index(ColoredCDBG<MyUnitigMap>& ccdbg,
                            const std::vector<Kmer>& head_kmer_arr)
{
    // iterate over entries, determine correct successors/predecessors (i.e. have correct colours)
    size_t head_kmer_arr_size = head_kmer_arr.size();
    #pragma omp parallel
    {
        #pragma omp for nowait
        for (int i = 1; i < head_kmer_arr_size + 1; i++)
        {
            // get a reference to the unitig map object
            auto um_pair = get_um_data(ccdbg, head_kmer_arr, i);
            auto& um = um_pair.first;
            auto& um_data = um_pair.second;

            // keep track of the number of neighbours
            size_t num_succ = 0;
            size_t num_pred = 0;

            // get copy of unitig_dict.unitig_full_colour to determine whether unitig gains a new colour not present in
            // gains a new colour not present in predecessors/successors
            auto full_colours = um_data->full_colour();

            // iterate over the connected nodes in successors
            for (const auto& succ : get_neighbours(um.getSuccessors()))
            {
                // get copy of successor unitig
                const auto adj_um_pair = get_um_data(ccdbg, succ.first);
                const auto& adj_um_data = adj_um_pair.second;

                // ensure colours are viable between the current and neigbouring unitig. Base this off head/tail colours depending on orientation
                // for current untiig, should negate tail, as this portion will overlap with next unitig
                boost::dynamic_bitset<> colours;
                if (succ.second)
                {
                    colours = um_data->tail_colour();
                    colours &= adj_um_data->head_colour();
                } else
                {
                    colours = um_data->tail_colour();
                    colours &= adj_um_data->tail_colour();
                }

                // negate adjacent full colours from full colours
                auto adj_unitig_dict_full = adj_um_data->full_colour();
                adj_unitig_dict_full.flip();
                full_colours &= adj_unitig_dict_full;

                // if colours are viable, add successor information to current unitig
                if (!colours.none())
                {
                    num_succ++;
                }
            }

            // determine if unitig has new colour not found in successors, if so set end_contig as true
            if (!full_colours.none())
            {
                um_data->set_end_contig(true);
            }

            // reset full_colours to repeat with predecessors
            full_colours = um_data->full_colour();

            // iterate over the connected nodes in predecessors, flipping sign to get successors in reverse strand
            um.strand = !um.strand;
            for (const auto& pred : get_neighbours(um.getSuccessors()))
            {
                // get copy of predecessor unitig_ID
                const auto adj_um_pair = get_um_data(ccdbg, pred.first);
                const auto& adj_um_data = adj_um_pair.second;

                // ensure colours are viable between the current and neigbouring unitig. Base this off head/tail colours depending on orientation
                // for current untiig, should negate head, as this portion will overlap with next unitig
                boost::dynamic_bitset<> colours;
                if (pred.second)
                {
                    colours = um_data->head_colour();
                    colours &= adj_um_data->head_colour();
                } else
                {
                    colours = um_data->head_colour();
                    colours &= adj_um_data->tail_colour();
                }

                // negate adjacent full colours from full colours
                auto adj_unitig_dict_full = adj_um_data->full_colour();
                adj_unitig_dict_full.flip();
                full_colours &= adj_unitig_dict_full;

                // if colours are viable, add successor information to current unitig
                if (!colours.none())
                {
                    num_pred++;
                }
            }

            // determine if unitig has new colour not found in predecessors, if so set end_contig as true
            if (!full_colours.none())
            {
                um_data->set_end_contig(true);
            }

            // if there are no successors/predecessors or unitig colours change, means that unitig is first/last in sequence. Therefore,
            // update full codon array to be 3 frame stop index to enable complete traversal
            if (num_succ == 0 || num_pred == 0 || !um_data->head_tail_colours_equal() || um_data->end_contig())
            {
                // set all bits to 1
                std::bitset<3> full_binary;
                full_binary.set();

                // update full codon indexes
                um_data->add_codon(true, full_binary);
                um_data->add_codon(false, full_binary);

                um_data->set_forward_stop(true);
                um_data->set_reverse_stop(true);

                // set unitig to end of contig
                um_data->set_end_contig(true);
            }
        }
    }
}

void calculate_genome_paths(const std::vector<Kmer>& head_kmer_arr,
                            ColoredCDBG<MyUnitigMap>& ccdbg,
                            const std::string& fasta_file,
                            const int kmer,
                            const int colour_ID)
{
    // generate the index
    fm_index_coll ref_index;

    // create fm index file name
    std::string idx_file_name = fasta_file + ".fmp";

    // initialise string of nodes for FM-index generation
    std::string genome_path;

    // open the file handler
    gzFile fp = gzopen(fasta_file.c_str(), "r");

    if (fp == 0) {
        perror("fopen");
        exit(1);
    }
    // initialize seq
    kseq_t *seq = kseq_init(fp);

    // for setting end contig
    std::bitset<3> full_binary;
    full_binary.set();

    // read sequence
    size_t contig_ID = 1;
    int l;
    while ((l = kseq_read(seq)) >= 0) {
        std::string entry = seq->seq.s;

        const size_t num_kmers = entry.length() - kmer + 1;

        // roll through the sequence, generating k-mers and querying them in graph
        if (num_kmers > 0) {
            // initialise variables for contig, add initial delimeter
            std::string contig_path = ",";

            const char *query_str = entry.c_str();

            // create int to identify if head and tail kmers in unitig have been traversed
            int prev_head = 0;

            for (KmerIterator it_km(query_str), it_km_end; it_km != it_km_end; ++it_km)
            {
                auto um = ccdbg.find(it_km->first);

                // if found, add to FM-index string
                if (!um.isEmpty) {
                    int strand = um.strand ? 1 : -1;

                    DataAccessor<MyUnitigMap>* da = um.getData();
                    MyUnitigMap* um_data = da->getData(um);

                    // look for the head in the graph, determine node ID and add to genome_path
                    int node_ID = um_data->get_id() * strand;

                    if (prev_head != node_ID) {
                        // if at start of contig i.e. prev_head==0, set as end_contig
                        if (!prev_head)
                        {
                            if (!um_data->end_contig())
                            {
                                // update full codon indexes
                                um_data->add_codon(true, full_binary);
                                um_data->add_codon(false, full_binary);

                                um_data->set_forward_stop(true);
                                um_data->set_reverse_stop(true);

                                // set unitig to end of contig
                                um_data->set_end_contig(true);
                            }
                        }

                        // set prev_head to new node
                        prev_head = node_ID;

                        // clear unitig colours
                        um_data->clear_colours();

                        // add new node
                        const std::string node_entry = std::to_string(node_ID);
                        contig_path += node_entry + ",";
                    }
                }
            }

            // map to last entry and assign end-contig
            auto um_pair = get_um_data(ccdbg, head_kmer_arr, prev_head);
            auto& um_data = um_pair.second;

            if (!um_data->end_contig()) {
                // update full codon indexes
                um_data->add_codon(true, full_binary);
                um_data->add_codon(false, full_binary);

                um_data->set_forward_stop(true);
                um_data->set_reverse_stop(true);

                // set unitig to end of contig
                um_data->set_end_contig(true);
            }

            // add delimiter between contigs
            genome_path += contig_path;
            genome_path += ";";
            contig_ID++;
        }
    }

    // destroy seq and fp objects
    kseq_destroy(seq);
    gzclose(fp);

    std::string outfile_name = fasta_file + "_node_ids.txt";
    ofstream outfile;
    outfile.open(outfile_name);

    outfile << genome_path;
    outfile.close();


    sdsl::construct_im(ref_index, genome_path, 1); // generate index
    store_to_file(ref_index, idx_file_name); // save it
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
        #pragma omp for nowait
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

    // determine if end contigs present if any of the input colours are not references
    if (ref_set.count() != ref_set.size())
    {
        update_neighbour_index(ccdbg, head_kmer_arr);
    }

    // generate index of references within ref_set
    std::vector<size_t> ref_index;
    for (int i = 0; i < input_colours.size(); i++)
    {
        if ((bool)ref_set[i])
        {
            ref_index.push_back(i);
        }
    }

    // generate FM-indexes of all reference fasta in node-space
    if (ref_set.any())
    {
        cout << "Mapping contigs to graph..." << endl;
        #pragma omp parallel
        {
            #pragma omp for nowait
            for (int i = 0; i < ref_index.size(); i++)
            {
                calculate_genome_paths(head_kmer_arr, ccdbg, input_colours[ref_index[i]], kmer, ref_index[i]);
            }
        }
    }

    // return node_colour vector
    return node_colour_vector;
}