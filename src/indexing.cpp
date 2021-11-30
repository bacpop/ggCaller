// ggCaller header
#include "indexing.h"

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
boost::dynamic_bitset<> generate_colours(const UnitigMap<DataAccessor<T>, DataStorage<U>, is_const> unitig,
                                           const size_t& nb_colours,
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
    // initialise unitig_dict
    unitigDict unitig_dict;

    // generate string from unitig
    const std::string unitig = um.referenceUnitigToString();
    const size_t unitig_len = unitig.size();

    // add unitig sequence to unitig_dict
    unitig_dict.add_seq(unitig);

    // get head kmer for unitig, add to kmer dictionary
    const Kmer head_kmer_binary = um.getUnitigHead();
    const std::string head_kmer = head_kmer_binary.toString();
    unitig_dict.add_head(head_kmer);

    // calculate unitig sizes
    unitig_dict.add_size(unitig_len, unitig_len - (kmer-1));

    // calculate colours for unitig
    boost::dynamic_bitset<> unitig_colours = generate_colours(um, nb_colours, 0);
    unitig_dict.add_head_colour(std::move(unitig_colours));

    // generate colours for tail also
    const Kmer tail_kmer = um.getUnitigTail();
    const auto um_tail = ccdbg.find(tail_kmer, true);
    const size_t tail_pos = unitig_len - kmer;
    unitig_colours = generate_colours(um_tail, nb_colours, tail_pos);
    unitig_dict.add_tail_colour(std::move(unitig_colours));
    
    // generate successor head kmers (need to flip sign to get successors in reverse strand)
    unitig_dict.set_succ(std::move(get_neighbours(um.getSuccessors())));
    um.strand = !um.strand;
    unitig_dict.set_pred(std::move(get_neighbours(um.getSuccessors())));


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

    unitig_dict.add_codon(true, true, 0, full_binary_pos);
    unitig_dict.add_codon(true, false, 0, full_binary_neg);

    // part unitig dictionary
    const uint8_t part_binary_pos = calculateFrame_binary(part_indices_pos);
    const uint8_t part_binary_neg = calculateFrame_binary(part_indices_neg);

    unitig_dict.add_codon(false, true, 0, part_binary_pos);
    unitig_dict.add_codon(false, true, 1, switchFrame_binary(part_binary_pos, 1));
    unitig_dict.add_codon(false, true, 2, switchFrame_binary(part_binary_pos, 2));

    unitig_dict.add_codon(false, false, 0, part_binary_neg);
    unitig_dict.add_codon(false, false, 1, switchFrame_binary(part_binary_neg, 1));
    unitig_dict.add_codon(false, false, 2, switchFrame_binary(part_binary_neg, 2));

    return unitig_dict;
}

void update_neighbour_index(GraphVector& graph_vector,
                            const robin_hood::unordered_map<std::string, size_t>& head_kmer_map)
{
    // iterate over entries, determine correct successors/predecessors (i.e. have correct colours)
    size_t graph_vector_size = graph_vector.size();
    #pragma omp parallel
    {
        #pragma omp for nowait
        for (size_t i = 0; i < graph_vector_size; i++)
        {
            // get a copy of the unitig_Dict object
            mtx1.lock();
            auto unitig_dict = graph_vector.at(i);
            mtx1.unlock();


            // get copy of unitig_dict.unitig_full_colour to determine whether unitig gains a new colour not present in
            // gains a new colour not present in predecessors/successors
            auto full_colours = unitig_dict.full_colour();

            // iterate over the connected nodes in successors
            for (const auto& succ : unitig_dict.get_succs())
            {
                // get numerical successor ID
                const size_t succ_id = head_kmer_map.at(succ.first);

                // get copy of successor unitig_ID
                mtx1.lock();
                auto adj_unitig_dict = graph_vector.at(succ_id - 1);
                mtx1.unlock();

                // ensure colours are viable between the current and neigbouring unitig. Base this off head/tail colours depending on orientation
                // for current untiig, should negate tail, as this portion will overlap with next unitig
                boost::dynamic_bitset<> colours;
                if (succ.second)
                {
                    colours = unitig_dict.tail_colour();
                    colours &= adj_unitig_dict.head_colour();
                } else
                {
                    colours = unitig_dict.tail_colour();
                    colours &= adj_unitig_dict.tail_colour();
                }

                // negate adjacent full colours from full colours
                auto adj_unitig_dict_full = adj_unitig_dict.full_colour();
                adj_unitig_dict_full.flip();
                full_colours &= adj_unitig_dict_full;

                // if colours are viable, add successor information to current unitig
                if (!colours.none())
                {
                    //generate integer value of successor ID, if negative strand ID will be negative etc.
                    int succ_id_int = (succ.second) ? succ_id : succ_id * -1;
                    std::pair<int, std::vector<uint8_t>> neighbour (succ_id_int, adj_unitig_dict.get_codon_dict(false, succ.second));
                    unitig_dict.add_neighbour(true, neighbour);
                }
            }

            // determine if unitig has new colour not found in successors, if so set end_contig as true
            if (!full_colours.none())
            {
                unitig_dict.set_end_contig(true);
            }

            // reset full_colours to repeat with predecessors
            full_colours = unitig_dict.full_colour();

            // clear succ_heads
            unitig_dict.clear_succ();

            // iterate over the connected nodes in predecessors
            for (const auto& pred : unitig_dict.get_preds())
            {
                // get numerical successor ID
                const size_t pred_id = head_kmer_map.at(pred.first);

                // get copy of successor unitig_ID
                mtx1.lock();
                auto adj_unitig_dict = graph_vector.at(pred_id - 1);
                mtx1.unlock();

                // ensure colours are viable between the current and neigbouring unitig. Base this off head/tail colours depending on orientation
                // for current untiig, should negate head, as this portion will overlap with next unitig
                boost::dynamic_bitset<> colours;
                if (pred.second)
                {
                    //colours = std::move(bool_and(unitig_dict.head_colour(), adj_unitig_dict.head_colour()));
                    colours = unitig_dict.head_colour();
                    colours &= adj_unitig_dict.head_colour();
                } else
                {
                    //colours = std::move(bool_and(unitig_dict.head_colour(), adj_unitig_dict.tail_colour()));
                    colours = unitig_dict.head_colour();
                    colours &= adj_unitig_dict.tail_colour();
                }

                auto adj_unitig_dict_full = adj_unitig_dict.full_colour();
                adj_unitig_dict_full.flip();
                full_colours &= adj_unitig_dict_full;

                // if colours are viable, add successor information to current unitig
                if (!colours.none())
                {
                    //generate integer value of successor ID, if negative strand ID will be negative etc.
                    int pred_id_int = (pred.second) ? pred_id : pred_id * -1;
                    std::pair<int, std::vector<uint8_t>> neighbour (pred_id_int, adj_unitig_dict.get_codon_dict(false, pred.second));
                    unitig_dict.add_neighbour(false, neighbour);
                }
            }

            // determine if unitig has new colour not found in predecessors, if so set end_contig as true
            if (!full_colours.none())
            {
                unitig_dict.set_end_contig(true);
            }

            // clear pred_heads
            unitig_dict.clear_pred();

            // if there are no successors/predecessors or unitig colours change, means that unitig is first/last in sequence. Therefore,
            // update full codon array to be 3 frame stop index to enable complete traversal
            if (unitig_dict.get_neighbours(true).empty() || unitig_dict.get_neighbours(false).empty() || !unitig_dict.head_tail_colours_equal() || unitig_dict.end_contig())
            {
                // set first three bits to 1
                uint8_t full_binary = 7;

                // update full codon indexes
                unitig_dict.add_codon(true, true, 0, full_binary);
                unitig_dict.add_codon(true, false, 0, full_binary);

                unitig_dict.set_forward_stop(true);
                unitig_dict.set_reverse_stop(true);

                // set unitig to end of contig
                unitig_dict.set_end_contig(true);
            }

            // update graph_vector with new entry
            mtx1.lock();
            graph_vector[i] = std::move(unitig_dict);
            mtx1.unlock();
        }
    }
}

NodeContigMapping calculate_genome_paths(const robin_hood::unordered_map<std::string, size_t>& head_kmer_map,
                                        const ColoredCDBG<>& ccdbg,
                                        const std::string& fasta_file,
                                        const int& kmer)
{
    // generate the index
    fm_index_coll ref_index;

    // create fm index file name
    std::string idx_file_name = fasta_file + ".fm";

    // initialise string of nodes for FM-index generation
    std::string genome_path;
    NodeContigMapping node_contig_mappings;

    // open the file handler
    gzFile fp = gzopen(fasta_file.c_str(), "r");

    if (fp == 0) {
        perror("fopen");
        exit(1);
    }
    // initialize seq
    kseq_t *seq = kseq_init(fp);

    // read sequence
    size_t contig_ID = 1;
    int l;
    while ((l = kseq_read(seq)) >= 0) {
        std::string entry = seq->seq.s;

        const size_t num_kmers = entry.length() - kmer + 1;

        // roll through the sequence, generating k-mers and querying them in graph
        if (num_kmers > 0) {
            // initialise variables for contig
            std::string contig_path;

            const char *query_str = entry.c_str();

            // create int to identify if head and tail kmers in unitig have been traversed
            int prev_head = 0;
            int kmer_counter = 0;

            size_t kmer_index = 1;
            for (KmerIterator it_km(query_str), it_km_end; it_km != it_km_end; ++it_km)
            {
                auto um = ccdbg.find(it_km->first);

                // if found, add to FM-index string
                if (!um.isEmpty) {
                    std::string head_kmer = um.getUnitigHead().toString();
                    int strand = um.strand ? 1 : -1;

                    // look for the head in the graph, determine node ID and add to genome_path
                    int node_ID = head_kmer_map.at(head_kmer) * strand;

                    if (prev_head != node_ID) {
                        if (prev_head != 0)
                        {
                            // if moving to new node, need to update previous node distance
                            std::get<3>(node_contig_mappings.back().second) = std::get<2>(node_contig_mappings.back().second) + kmer_counter;
                        }
                        prev_head = node_ID;

                        // add new node
                        const std::string node_entry = std::to_string(node_ID);
                        contig_path += node_entry + ",";
                        node_contig_mappings.push_back({abs(node_ID) - 1, {contig_ID, kmer_index, um.dist, 0, um.strand}});

                        // reset kmer counter
                        kmer_counter = 0;
                    }

                    // if at last kmer, need to determine position within unitig
                    if (kmer_index == num_kmers)
                    {
                        std::get<3>(node_contig_mappings.back().second) = std::get<2>(node_contig_mappings.back().second) + kmer_counter;
                    }
                }
                kmer_counter++;
                kmer_index++;
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

    sdsl::construct_im(ref_index, genome_path, 1); // generate index
    store_to_file(ref_index, idx_file_name); // save it

    return node_contig_mappings;
}

NodeColourVector index_graph(GraphVector& graph_vector,
                              robin_hood::unordered_map<std::string, size_t>& head_kmer_map,
                              const ColoredCDBG<>& ccdbg,
                              const std::vector<std::string>& stop_codons_for,
                              const std::vector<std::string>& stop_codons_rev,
                              const int kmer,
                              const size_t nb_colours,
                              const bool is_ref,
                              const std::vector<std::string>& input_colours)
{
    // get all head kmers
    std::vector<Kmer> head_kmer_arr;
    for (const auto& um : ccdbg)
    {
        head_kmer_arr.push_back(um.getUnitigHead());
    }

    // structures for results
    GraphVector graph_vector_private(head_kmer_arr.size());

    NodeColourVector node_colour_vector(nb_colours);

    // run unitig indexing in parallel
    size_t unitig_id = 1;
    #pragma omp parallel
    {
        NodeColourVector node_colour_vector_private(nb_colours);
        robin_hood::unordered_map<std::string, size_t> head_kmer_map_private;
        #pragma omp for nowait
        for (auto it = head_kmer_arr.begin(); it < head_kmer_arr.end(); it++)
        {
            // convert Kmer defined in *it to unitig
            auto unitig = ccdbg.find(*it, true);

            // generate results per unitig
            unitigDict unitig_dict = std::move(analyse_unitigs_binary(ccdbg, unitig, stop_codons_for, stop_codons_rev, kmer, nb_colours));
            #pragma omp atomic capture
            unitig_dict.id = unitig_id++;

            // add to node_colour_map_private
            for (size_t i = 0; i < unitig_dict.full_colour().size(); i++)
            {
                if (unitig_dict.full_colour()[i])
                {
                    node_colour_vector_private[i].push_back(unitig_dict.id);
                }
            }

            // add head_kmer and unitig id to map
            head_kmer_map_private[unitig_dict.head_kmer()] = unitig_dict.id;

            // add unitig to graph_vector, minus 1 as zero based
            graph_vector_private[unitig_dict.id - 1] = std::move(unitig_dict);
        }
        #pragma omp critical
        {
            head_kmer_map.insert(head_kmer_map_private.begin(), head_kmer_map_private.end());

            // update node_colour_vector with calculated colours
            for (int i = 0; i < node_colour_vector_private.size(); i++)
            {
                node_colour_vector[i].insert(node_colour_vector[i].end(), make_move_iterator(node_colour_vector_private[i].begin()), make_move_iterator(node_colour_vector_private[i].end()));
            }
        }
    }
    // update neighbour index in place within graph_vector
    update_neighbour_index(graph_vector_private, head_kmer_map);

    // assign the graph vector to the graph _GraphVector
    graph_vector = std::move(graph_vector_private);

    // generate FM-indexes of all fastas in node-space
    std::vector<NodeContigMapping> colour_contig_mappings(nb_colours);
    if (is_ref)
    {
        for (int i = 0; i < nb_colours; i++)
        {
            colour_contig_mappings[i] = std::move(calculate_genome_paths(head_kmer_map, ccdbg, input_colours[i], kmer));
        }
    }

    // add contig locations to graph_vector
    for (size_t i = 0; i < nb_colours; i++)
    {
        const auto& contig_mappings = colour_contig_mappings.at(i);
        for (const auto& node_entry : contig_mappings)
        {
            graph_vector[node_entry.first].add_contig_coords(i, node_entry.second);
        }
    }

    // return node_colour vector
    return node_colour_vector;
}