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
                                   const bool reverse)
{
    std::vector<std::size_t> index_list;
    size_t pos = seq.find(subseq, start_index);

    while(pos != string::npos)
    {
        index_list.push_back(pos);
        pos = seq.find(subseq,pos+1);
    }

    if (reverse)
    {
        for (auto& entry : index_list)
        {
            // reverse by adding 3 (two to start from start of stop codon, 1 to make zero indexed)
            entry = seq.size() - (entry + 3);
        }
    }
    return index_list;
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
                            size_t& num_stops,
                            size_t& num_codons,
                            const std::vector<std::string>& stop_codon_for,
                            const std::vector<std::string>& stop_codon_rev,
                            const std::vector<std::string>& start_codon_for,
                            const std::vector<std::string>& start_codon_rev,
                            const int& kmer,
                            const size_t& nb_colours,
                            tbb::concurrent_unordered_map<std::string, tbb::concurrent_unordered_set<int>>& start_freq_set,
                            const int& aa_kmer)
{
    // initialise unitig_dict
    DataAccessor<MyUnitigMap>* da = um.getData();
    MyUnitigMap* um_data = da->getData(um);

    // generate string from unitig
    const std::string unitig = um.referenceUnitigToString();
    const size_t unitig_len = unitig.size();

    // number of 3-mers is sequence length negate 2, multiplied by two for forward and reverse strand
    #pragma omp atomic
    num_codons += (unitig_len - 2) * 2;

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

    // use overlap - 2, taking into account unitigs where stop codon is split at end
    const int overlap = kmer - 3;

    // find the full and part indices of each of the codons in the unitig in positive strand
    for (const auto& codon : stop_codon_for)
    {
        // full positive strand
        std::vector<std::size_t> found_indices = findIndex(unitig, codon, 0, false);

        #pragma omp atomic
        num_stops += found_indices.size();

        full_indices_pos.insert(full_indices_pos.end(), found_indices.begin(), found_indices.end());

        // part positive strand, take off overlap - 2, as may be case that codon in upstream node is split at end of unitig
        for (const auto& entry : found_indices)
        {
            if (entry >= (overlap))
            {
                part_indices_pos.push_back(entry - (overlap));
            }
        }
    }

    // find the full and part indices of each of the codons in the unitig in negative strand

    // get this working, now working without rfind, need to determine correct coorindates for partial unitigs
    for (const auto& codon : stop_codon_rev)
    {
        // full negative strand
        std::vector<std::size_t> found_indices = findIndex(unitig, codon, 0, true);

        #pragma omp atomic
        num_stops += found_indices.size();

        full_indices_neg.insert(full_indices_neg.end(), found_indices.begin(), found_indices.end());

        // part negative strand, taking into account potential for stop codon to be cut off from previous unitig
        for (const auto& entry : found_indices)
        {
            if (entry >= (overlap))
            {
                part_indices_neg.push_back(entry - (overlap));
            }
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

    // look for start codons forward
    for (const auto& codon : start_codon_for)
    {
        // full positive strand
        std::vector<std::size_t> found_indices = findIndex(unitig, codon, 0, false);
        // if no start codons found, pass
        if (found_indices.size() == 0)
        {
            continue;
        }

        // pull out start codon positions
        for (const auto& pos : found_indices)
        {
            if (unitig.size() - pos >= kmer)
            {
                std::string start_site_DNA = unitig.substr(pos, kmer);
                std::string start_site_AA = (translate(start_site_DNA)).aa();

                if (start_site_AA.find('*') != std::string::npos)
                {
                    continue;
                }

                const int num_kmers = start_site_AA.size() - aa_kmer;

                std::vector<std::string> AA_kmers(num_kmers);

                for (int kmer_index = 0; kmer_index < num_kmers; ++kmer_index)
                {
                    AA_kmers[kmer_index] = get_kmers(start_site_AA, kmer_index, aa_kmer);
                }

                // add colours to start_freq_set
                for (int i = 0; i < nb_colours; i++)
                {
                    if ((bool)full_unitig_colour[i])
                    {
                        for (const auto& entry_aa : AA_kmers)
                        {
                            start_freq_set[entry_aa].insert(i);
                        }
                    }
                }
            }
        }
    }

    // reverse complement unitig
    const std::string rev_unitig = reverse_complement(unitig);

    // look for start codons forward
    for (const auto& codon : start_codon_rev)
    {
        // full positive strand
        std::vector<std::size_t> found_indices = findIndex(unitig, codon, 0, true);
        // if no start codons found, pass
        if (found_indices.size() == 0)
        {
            continue;
        }

        // pull out start codon positions
        for (const auto& pos : found_indices)
        {
            if (unitig.size() - pos >= kmer)
            {
                std::string start_site_DNA = rev_unitig.substr(pos, kmer);
                std::string start_site_AA = (translate(start_site_DNA)).aa();

                if (start_site_AA.find('*') != std::string::npos)
                {
                    continue;
                }

                const int num_kmers = start_site_AA.size() - aa_kmer;

                std::vector<std::string> AA_kmers(num_kmers);

                for (int kmer_index = 0; kmer_index < num_kmers; ++kmer_index)
                {
                    AA_kmers[kmer_index] = get_kmers(start_site_AA, kmer_index, aa_kmer);
                }

                // add colours to start_freq_set
                for (int i = 0; i < nb_colours; i++)
                {
                    if ((bool)full_unitig_colour[i])
                    {
                        for (const auto& entry_aa : AA_kmers)
                        {
                            start_freq_set[entry_aa].insert(i);
                        }
                    }
                }
            }
        }
    }
}

void calculate_genome_paths(const std::vector<Kmer>& head_kmer_arr,
                            ColoredCDBG<MyUnitigMap>& ccdbg,
                            const std::string& fasta_file,
                            const int kmer,
                            const int colour_ID,
                            const size_t nb_colours)
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

    // read sequence
    size_t contig_ID = 1;
    int l;
    while ((l = kseq_read(seq)) >= 0) {
        std::string entry = seq->seq.s;

        const int num_kmers = entry.length() - kmer + 1;

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
                            um_data->set_end_contig(colour_ID, nb_colours);
                        }

                        // set prev_head to new node
                        prev_head = node_ID;

                        // add new node
                        const std::string node_entry = std::to_string(node_ID);
                        contig_path += node_entry + ",";
                    }
                }
            }

            // map to last entry and assign end-contig
            auto um_pair = get_um_data(ccdbg, head_kmer_arr, prev_head);
            auto& um_data = um_pair.second;

            um_data->set_end_contig(colour_ID, nb_colours);

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
}

NodeColourVector index_graph(std::vector<Kmer>& head_kmer_arr,
                             ColoredCDBG<MyUnitigMap>& ccdbg,
                             float& stop_codon_freq,
                             const std::vector<std::string>& stop_codons_for,
                             const std::vector<std::string>& stop_codons_rev,
                             const std::vector<std::string>& start_codons_for,
                             const std::vector<std::string>& start_codons_rev,
                             const int kmer,
                             const size_t nb_colours,
                             const std::vector<std::string>& input_colours,
                             const boost::dynamic_bitset<>& ref_set,
                             robin_hood::unordered_map<std::string, size_t>& start_freq)
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

    // determine number of stops and codons
    size_t num_stops = 0;
    size_t num_codons = 0;

    // determine sharing of start positions
    tbb::concurrent_unordered_map<std::string, tbb::concurrent_unordered_set<int>> start_freq_set;

    // get amino acid k-mer length. Make it so that mutation in middle of string won't affect all k-mers.
    const int aa_kmer = std::round((float) kmer / (float) 6) - 1;

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
            analyse_unitigs_binary(ccdbg, um, num_stops, num_codons, stop_codons_for, stop_codons_rev, start_codons_for, start_codons_rev, kmer, nb_colours, start_freq_set, aa_kmer);
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

    // determine average stop codon frequency per codon
    stop_codon_freq = (float)num_stops / (float)num_codons;

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
            #pragma omp for
            for (int i = 0; i < ref_index.size(); i++)
            {
                calculate_genome_paths(head_kmer_arr, ccdbg, input_colours[ref_index[i]], kmer, ref_index[i], nb_colours);
            }
        }
    }

    // go through start_freq_set, determine coverage of start
    for (const auto& entry : start_freq_set)
    {
        start_freq[entry.first] = entry.second.size();
    }

    // return node_colour vector
    return node_colour_vector;
}