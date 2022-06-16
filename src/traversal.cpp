// ggCaller header
#include "traversal.h"
PathVector iter_nodes_binary (const ColoredCDBG<MyUnitigMap>& ccdbg,
                              const std::vector<Kmer>& head_kmer_arr,
                              const NodeTuple& head_node_tuple,
                              const size_t current_colour,
                              const size_t length_max,
                              const size_t overlap,
                              const bool repeat,
                              const bool is_ref,
                              const boost::dynamic_bitset<>& ref_set,
                              const fm_index_coll& fm_idx)
{
    // generate path list, vector for path and the stack
    PathVector path_list;
    std::vector<int> node_vector;
    NodeStack node_stack;

    // create first item in stack
    node_stack.push(head_node_tuple);

    while(!node_stack.empty())
    {
        // pop node in stack
        auto node_tuple = node_stack.top();
        node_stack.pop();

        // unpack tuple
        const size_t & pos_idx = std::get<0>(node_tuple);
        const int & node_id = std::get<1>(node_tuple);
        const auto & codon_arr = std::get<2>(node_tuple);
        const boost::dynamic_bitset<>& colour_arr = std::get<3>(node_tuple);
        const size_t & path_length = std::get<4>(node_tuple);

        // slice path, unless at first node
        if (pos_idx != 0)
        {
            node_vector = std::vector<int> (node_vector.begin(), node_vector.begin() + pos_idx);
        }

        // add node to path
        node_vector.push_back(node_id);

        // get length of vector for new pos_idx
        const size_t new_pos_idx = node_vector.size();

        // get unitig data
        auto um_pair = get_um_data(ccdbg, head_kmer_arr, node_id);
        auto& um = um_pair.first;
        auto& um_data = um_pair.second;

        // determine strand of unitig
        const bool strand = (node_id >= 0) ? true : false;

        // reverse the strand of the unitig
        if (!strand)
        {
            um.strand = !um.strand;
        }

        // iterate over neighbours, recurring through incomplete paths
        for (auto& neighbour_um : um.getSuccessors())
        {
            auto neighbour_da = neighbour_um.getData();
            auto neighbour_um_data = neighbour_da->getData(neighbour_um);
            const bool neighbour_strand = neighbour_um.strand;

            // calculate colours array
            auto updated_colours_arr = colour_arr;
            auto neighbour_colour = neighbour_um_data->full_colour();
            updated_colours_arr &= neighbour_colour;

            // determine if neighbour is in same colour as iteration, if not pass
            if (!(bool)updated_colours_arr[current_colour])
            {
                continue;
            }

            // parse neighbour information. Frame is next stop codon, with first dictating orientation and second the stop codon index
            const int neighbour_id = (neighbour_strand) ? neighbour_um_data->get_id() : neighbour_um_data->get_id() * -1;

            // check against fm-idx every other node, pass if not present
            if (is_ref && node_vector.size() % 2 == 1)
            {
                std::vector<int> check_vector = node_vector;
                check_vector.push_back(neighbour_id);
                std::pair<bool, bool> present = path_search(check_vector, fm_idx);
                if (!present.first)
                {
                    continue;
                }
            } else if (!is_ref && !repeat)
            {
                // if using reads, check if unitig has already been traversed, and pass if repeat not specified
                const bool is_in = std::find(node_vector.begin(), node_vector.end(), neighbour_id) != node_vector.end();
                if (is_in)
                {
                    continue;
                }
            }


            // if not is_ref, check that unitig is shared in at least one other colour
//            if (!is_ref)
//            {
//                if (neighbour_colour.count() < 2)
//                {
//                    continue;
//                }
//            }

            // check path length, if too long continue
            const size_t updated_path_length = path_length + (neighbour_um.size - overlap);

            if (updated_path_length > length_max)
            {
                // if path is only 2 nodes in length return as will contain very large sequences
                if (new_pos_idx == 1)
                {
                    // check if end node is reference sequence
                    if (!is_ref)
                    {
                        auto ref_check = neighbour_colour;
                        ref_check &= ref_set;
                        if (ref_check.none())
                        {
                            continue;
                        }
                    }

                    // create temporary path to account for reaching end of contig
                    std::vector<int> return_path = node_vector;
                    return_path.push_back(neighbour_id);

                    // update path_list to enable this to be returned
                    path_list.push_back(std::move(return_path));
                }
                continue;
            }

            // calculate modulus for path and updated cached path reading frame
            int modulus = path_length % 3;
            std::bitset<3> updated_codon_arr = (codon_arr & ~(neighbour_um_data->get_codon_arr(false, neighbour_strand, modulus)));

            // return path if end of contig found or if stop indexes paired
            if (neighbour_um_data->end_contig(current_colour) || updated_codon_arr.none())
            {
                // check if end node is reference sequence
                if (!is_ref)
                {
                    auto ref_check = neighbour_colour;
                    ref_check &= ref_set;
                    if (ref_check.none())
                    {
                        continue;
                    }
                }

                // create temporary path to account for reaching end of contig
                std::vector<int> return_path = node_vector;
                return_path.push_back(neighbour_id);

                // update path_list to enable this to be returned
                path_list.push_back(std::move(return_path));

                // if node is end_contig, continue to add to path
                if (neighbour_um_data->end_contig(current_colour) && !updated_codon_arr.none())
                {
                    // if no previous conditions are satisfied, prepare tuple for stack
                    NodeTuple new_node_tuple(new_pos_idx, neighbour_id, updated_codon_arr, updated_colours_arr, updated_path_length);

                    // add to stack
                    node_stack.push(new_node_tuple);
                }

                // move onto next neighbour
                continue;
            }

            // if no previous conditions are satisfied, prepare tuple for stack
            NodeTuple new_node_tuple(new_pos_idx, neighbour_id, updated_codon_arr, updated_colours_arr, updated_path_length);

            // add to stack
            node_stack.push(new_node_tuple);
        }
    }
    return path_list;
}

PathVector iter_nodes_path (const ColoredCDBG<MyUnitigMap>& ccdbg,
                            const std::vector<Kmer>& head_kmer_arr,
                            const NodeTuple& head_node_tuple,
                            const size_t current_colour,
                            const size_t length_max,
                            const size_t overlap,
                            const fm_index_coll& fm_idx,
                            const robin_hood::unordered_map<int, robin_hood::unordered_set<int>>& edge_map,
                            const bool ref_strand)
{
    // generate path list, vector for path and the stack
    PathVector path_list;
    std::vector<int> node_vector;
    NodeStack node_stack;

    // create first item in stack
    node_stack.push(head_node_tuple);

    while(!node_stack.empty())
    {
        // pop node in stack
        auto node_tuple = node_stack.top();
        node_stack.pop();

        // unpack tuple
        const size_t & pos_idx = std::get<0>(node_tuple);
        const int & node_id = std::get<1>(node_tuple);
        const auto & codon_arr = std::get<2>(node_tuple);
        const boost::dynamic_bitset<>& colour_arr = std::get<3>(node_tuple);
        const size_t & path_length = std::get<4>(node_tuple);

        // slice path, unless at first node
        if (pos_idx != 0)
        {
            node_vector = std::vector<int> (node_vector.begin(), node_vector.begin() + pos_idx);
        }

        // add node to path
        node_vector.push_back(node_id);

        // get length of vector for new pos_idx
        const size_t new_pos_idx = node_vector.size();

        // get unitig data
        auto um_pair = get_um_data(ccdbg, head_kmer_arr, node_id);
        auto& um = um_pair.first;
        auto& um_data = um_pair.second;

        // determine strand of unitig
        const bool strand = (node_id >= 0) ? true : false;

        // reverse the strand of the unitig
        if (!strand)
        {
            um.strand = !um.strand;
        }

        // iterate over neighbours, recurring through incomplete paths
        const auto& edge_it = edge_map.find(node_id);
        // if not found, treat as end contig
        if (edge_it == edge_map.end())
        {
            path_list.push_back(node_vector);
            continue;
        }
        for (auto& neighbour_id : edge_it->second)
        {
            auto neighbour_da_pair = get_um_data(ccdbg, head_kmer_arr, neighbour_id);
            auto& neighbour_um = um_pair.first;
            auto& neighbour_um_data = um_pair.second;
            const bool neighbour_strand = (neighbour_id >= 0) ? true : false;

            // calculate colours array
            auto updated_colours_arr = colour_arr;
            auto neighbour_colour = neighbour_um_data->full_colour();
            updated_colours_arr &= neighbour_colour;

            // determine if neighbour is in same colour as iteration, if not pass
            if (!(bool)updated_colours_arr[current_colour])
            {
                continue;
            }

            // check if real path
            {
                std::vector<int> check_vector = node_vector;
                check_vector.push_back(neighbour_id);
                bool present = path_search_strand(check_vector, fm_idx, ref_strand);
                if (!present)
                {
                    continue;
                }
            }

            // check path length, if too long continue
            const size_t updated_path_length = path_length + (neighbour_um.size - overlap);

            if (updated_path_length > length_max)
            {
                // if path is only 2 nodes in length return as will contain very large sequences
                if (new_pos_idx == 1)
                {
                    // create temporary path to account for reaching end of contig
                    std::vector<int> return_path = node_vector;
                    return_path.push_back(neighbour_id);

                    // update path_list to enable this to be returned
                    path_list.push_back(std::move(return_path));
                }
                continue;
            }

            // calculate modulus for path and updated cached path reading frame
            int modulus = path_length % 3;
            std::bitset<3> updated_codon_arr = (codon_arr & ~(neighbour_um_data->get_codon_arr(false, neighbour_strand, modulus)));

            // return path if end of contig found or if stop indexes paired
            if (neighbour_um_data->end_contig(current_colour) || updated_codon_arr.none())
            {
                // create temporary path to account for reaching end of contig
                std::vector<int> return_path = node_vector;
                return_path.push_back(neighbour_id);

                // update path_list to enable this to be returned
                path_list.push_back(std::move(return_path));

                // if node is end_contig, continue to add to path
                if (neighbour_um_data->end_contig(current_colour) && !updated_codon_arr.none())
                {
                    // if no previous conditions are satisfied, prepare tuple for stack
                    NodeTuple new_node_tuple(new_pos_idx, neighbour_id, updated_codon_arr, updated_colours_arr, updated_path_length);

                    // add to stack
                    node_stack.push(new_node_tuple);
                }

                // move onto next neighbour
                continue;
            }

            // if no previous conditions are satisfied, prepare tuple for stack
            NodeTuple new_node_tuple(new_pos_idx, neighbour_id, updated_codon_arr, updated_colours_arr, updated_path_length);

            // add to stack
            node_stack.push(new_node_tuple);
        }
    }
    return path_list;
}

ORFNodeRobMap calculate_genome_paths(const std::vector<Kmer>& head_kmer_arr,
                                    ColoredCDBG<MyUnitigMap>& ccdbg,
                                    const std::string& fasta_file,
                                    const int kmer,
                                    const int colour_ID,
                                    const size_t nb_colours,
                                    const size_t& max_path_length,
                                    const std::vector<std::string>& stop_codons_for,
                                    const std::vector<std::string>& start_codons_for,
                                    const size_t min_ORF_length,
                                    torch::jit::script::Module& TIS_model,
                                    const double& minimum_ORF_score,
                                    const bool no_filter,
                                    tbb::concurrent_unordered_map<size_t, float>& all_TIS_scores)
{
    // initialise structures to hold ORFs
    ORFNodeMap ORF_node_map;
    std::unordered_set<size_t> hashes_to_remove;

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
    while ((l = kseq_read(seq)) >= 0)
    {
        // generate map to hold true edges from contig
        robin_hood::unordered_map<int, robin_hood::unordered_set<int>> edge_map;

        // generate set to hold nodes to start traversal from
        robin_hood::unordered_set<int> forward_stop_nodes;
        robin_hood::unordered_set<int> reverse_stop_nodes;

        // generate string from contig
        std::string entry = seq->seq.s;

        const int num_kmers = entry.length() - kmer + 1;

        // roll through the sequence, generating k-mers and querying them in graph
        if (num_kmers > 0) {
            // initialise small local fm-index
            fm_index_coll local_index;

            // scope for contig path
            {
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

                        if (um_data->forward_stop())
                        {
                            forward_stop_nodes.insert(node_ID);
                        }

                        if (um_data->reverse_stop())
                        {
                            reverse_stop_nodes.insert(node_ID * -1);
                        }

                        if (prev_head != node_ID) {
                            // if at start of contig i.e. prev_head==0, set as end_contig
                            if (!prev_head)
                            {
                                um_data->set_end_contig(colour_ID, nb_colours);
                                forward_stop_nodes.insert(node_ID);
                            } else
                            {
                                edge_map[prev_head].insert(node_ID);
                                edge_map[node_ID * -1].insert(prev_head * -1);
                            }

                            // set prev_head to new node
                            prev_head = node_ID;

                            // add new node
                            const std::string node_entry = std::to_string(node_ID);
                            contig_path += node_entry + ",";
                        }
                    }
                }

                // add delimiter between contigs
                genome_path += contig_path;
                genome_path += ";";
                contig_ID++;

                // map to last entry and assign end-contig
                auto um_pair = get_um_data(ccdbg, head_kmer_arr, prev_head);
                auto& um_data = um_pair.second;

                um_data->set_end_contig(colour_ID, nb_colours);
                reverse_stop_nodes.insert(prev_head * -1);

                // construct small fm-index for ORF calling
                sdsl::construct_im(local_index, contig_path, 1);
            }

            // set for any end contigs
            std::bitset<3> full_binary;
            full_binary.set();

            // traverse graph in forward direction
            for (const auto& head_id : forward_stop_nodes)
            {
                // get unitig data
                const auto um_pair = get_um_data(ccdbg, head_kmer_arr, head_id);
                const auto& um = um_pair.first;
                const auto& um_data = um_pair.second;

                // gather unitig information from graph_vector
                std::bitset<3> codon_arr = um_data->get_codon_arr(true, true, 0);
                const size_t unitig_len = um.size;
                const boost::dynamic_bitset<> colour_arr = um_data->full_colour();

                // if end contig, set to full array
                if (um_data->end_contig(colour_ID))
                {
                    codon_arr = full_binary;
                }

                // generate node tuple for iteration
                NodeTuple head_node_tuple(0, head_id, codon_arr, colour_arr, unitig_len);

                // recur paths
                PathVector unitig_complete_paths = iter_nodes_path(ccdbg, head_kmer_arr, head_node_tuple, colour_ID, max_path_length, kmer -1, local_index, edge_map, true);

                // iterate over paths, calling ORFs
                if (!unitig_complete_paths.empty())
                {
                    // iterate over all_paths
                    for (int i = 0; i < unitig_complete_paths.size(); i++)
                    {
                        // generate all ORFs within the path for start and stop codon pairs
                        generate_ORFs(colour_ID, ORF_node_map, hashes_to_remove, ccdbg, head_kmer_arr, stop_codons_for, start_codons_for, unitig_complete_paths[i], kmer -1, min_ORF_length,
                                      true, local_index, TIS_model, minimum_ORF_score, no_filter, all_TIS_scores);
                    }
                }
            }

            // traverse graph in forward and reverse direction
            for (const auto& head_id : reverse_stop_nodes)
            {
                // get unitig data
                const auto um_pair = get_um_data(ccdbg, head_kmer_arr, head_id);
                const auto& um = um_pair.first;
                const auto& um_data = um_pair.second;

                // gather unitig information from graph_vector
                std::bitset<3> codon_arr = um_data->get_codon_arr(true, true, 0);
                const size_t unitig_len = um.size;
                const boost::dynamic_bitset<> colour_arr = um_data->full_colour();

                // if end contig, set to full array
                if (um_data->end_contig(colour_ID))
                {
                    codon_arr = full_binary;
                }

                // generate node tuple for iteration
                NodeTuple head_node_tuple(0, head_id, codon_arr, colour_arr, unitig_len);

                // recur paths
                PathVector unitig_complete_paths = iter_nodes_path(ccdbg, head_kmer_arr, head_node_tuple, colour_ID, max_path_length, kmer -1, local_index, edge_map, false);

                // iterate over paths, calling ORFs
                if (!unitig_complete_paths.empty())
                {
                    // iterate over all_paths
                    for (int i = 0; i < unitig_complete_paths.size(); i++)
                    {
                        // generate all ORFs within the path for start and stop codon pairs
                        generate_ORFs(colour_ID, ORF_node_map, hashes_to_remove, ccdbg, head_kmer_arr, stop_codons_for, start_codons_for, unitig_complete_paths[i], kmer -1, min_ORF_length,
                                      true, local_index, TIS_model, minimum_ORF_score, no_filter, all_TIS_scores);
                    }
                }
            }
        }
    }

    // destroy seq and fp objects
    kseq_destroy(seq);
    gzclose(fp);

    sdsl::construct_im(ref_index, genome_path, 1); // generate index
    store_to_file(ref_index, idx_file_name); // save it

    // remove hashes to remove from ORF_hash
    for (const auto& hash : hashes_to_remove)
    {
        auto it = ORF_node_map.find(hash);
        if (it != ORF_node_map.end())
        {
            ORF_node_map.erase(it);
        }
    }

    // generate empty pos_strand_map
    NodeStrandMap pos_strand_map;

    // group colours of ORFs together
    ORFNodeRobMap ORF_map = std::move(sort_ORF_indexes(ORF_node_map, pos_strand_map, ccdbg, head_kmer_arr, true));

    return ORF_map;
}

ORFNodeRobMap traverse_graph(const ColoredCDBG<MyUnitigMap>& ccdbg,
                         const std::vector<Kmer>& head_kmer_arr,
                         const size_t colour_ID,
                         const std::vector<size_t>& node_ids,
                         const bool repeat,
                         const size_t max_path_length,
                         const size_t overlap,
                         const bool is_ref,
                         const boost::dynamic_bitset<>& ref_set,
                         const fm_index_coll& fm_idx,
                         const std::vector<std::string>& stop_codons_for,
                         const std::vector<std::string>& start_codons_for,
                         const size_t min_ORF_length,
                         torch::jit::script::Module& TIS_model,
                         const double& minimum_ORF_score,
                         const bool no_filter,
                         tbb::concurrent_unordered_map<size_t, float>& all_TIS_scores)
{
    //initialise ORF_nodes_paths to add ORF sequences to
    ORFNodeMap ORF_node_map;
    std::unordered_set<size_t> hashes_to_remove;

    // set for any end contigs
    std::bitset<3> full_binary;
    full_binary.set();

    // traverse nodes in forward direction
    for (const auto& node_id : node_ids)
    {
        // get unitig data
        const auto um_pair = get_um_data(ccdbg, head_kmer_arr, node_id);
        const auto& um = um_pair.first;
        const auto& um_data = um_pair.second;

        // check if stop codons present. If not, pass
        if (!um_data->forward_stop() && !um_data->end_contig(colour_ID))
        {
            continue;
        }

        // generate integer version of unitig_id for recursion
        const int head_id = (int) node_id;

        // gather unitig information from graph_vector
        std::bitset<3> codon_arr = um_data->get_codon_arr(true, true, 0);
        const size_t unitig_len = um.size;
        const boost::dynamic_bitset<> colour_arr = um_data->full_colour();

        // if end contig, set to full array
        if (um_data->end_contig(colour_ID))
        {
            codon_arr = full_binary;
        }

        // check if end node is reference sequence
        if (!is_ref)
        {
            auto ref_check = colour_arr;
            ref_check &= ref_set;
            if (ref_check.none())
            {
                continue;
            }
        }

        // generate node tuple for iteration
        NodeTuple head_node_tuple(0, head_id, codon_arr, colour_arr, unitig_len);

        // recur paths
        PathVector unitig_complete_paths = iter_nodes_binary(ccdbg, head_kmer_arr, head_node_tuple, colour_ID, max_path_length, overlap, repeat, is_ref, ref_set, fm_idx);

        // iterate over paths, calling ORFs
        if (!unitig_complete_paths.empty())
        {
            // iterate over all_paths
            for (int i = 0; i < unitig_complete_paths.size(); i++)
            {
                // generate all ORFs within the path for start and stop codon pairs
                generate_ORFs(colour_ID, ORF_node_map, hashes_to_remove, ccdbg, head_kmer_arr, stop_codons_for, start_codons_for, unitig_complete_paths[i], overlap, min_ORF_length, is_ref, fm_idx, TIS_model, minimum_ORF_score, no_filter, all_TIS_scores);
            }
        }
    }

    // traverse nodes in reverse direction
    for (const auto& node_id : node_ids)
    {
        // get unitig data
        const auto um_pair = get_um_data(ccdbg, head_kmer_arr, node_id);
        const auto& um = um_pair.first;
        const auto& um_data = um_pair.second;

        // check if stop codons present. If not, pass
        if (!um_data->reverse_stop() && !um_data->end_contig(colour_ID))
        {
            continue;
        }

        // generate integer version of unitig_id for recursion, multiplied by -1 to indicate reversal
        const int head_id = ((int) node_id) * -1;

        // gather unitig information from graph_vector
        std::bitset<3> codon_arr = um_data->get_codon_arr(true, false, 0);
        const size_t unitig_len = um.size;
        const boost::dynamic_bitset<> colour_arr = um_data->full_colour();

        // if end contig, set to full array
        if (um_data->end_contig(colour_ID))
        {
            codon_arr = full_binary;
        }

        // check if end node is reference sequence
        if (!is_ref)
        {
            auto ref_check = colour_arr;
            ref_check &= ref_set;
            if (ref_check.none())
            {
                continue;
            }
        }

        // generate node tuple for iteration
        NodeTuple head_node_tuple(0, head_id, codon_arr, colour_arr, unitig_len);

        // recur paths
        PathVector unitig_complete_paths = iter_nodes_binary(ccdbg, head_kmer_arr, head_node_tuple, colour_ID, max_path_length, overlap, repeat, is_ref, ref_set, fm_idx);

        // iterate over paths, calling ORFs
        if (!unitig_complete_paths.empty())
        {
            // iterate over all_paths
            for (int i = 0; i < unitig_complete_paths.size(); i++)
            {
                // generate all ORFs within the path for start and stop codon pairs
                generate_ORFs(colour_ID, ORF_node_map, hashes_to_remove, ccdbg, head_kmer_arr, stop_codons_for, start_codons_for, unitig_complete_paths[i], overlap, min_ORF_length, is_ref, fm_idx, TIS_model, minimum_ORF_score, no_filter, all_TIS_scores);
            }
        }
    }

    // remove hashes to remove from ORF_hash
    for (const auto& hash : hashes_to_remove)
    {
        auto it = ORF_node_map.find(hash);
        if (it != ORF_node_map.end())
        {
            ORF_node_map.erase(it);
        }
    }

    // generate pos_strand_map to determine relative strands of each node for each colour if not given by alignment of contigs
    NodeStrandMap pos_strand_map;
    if (!is_ref)
    {
        pos_strand_map = std::move(calculate_pos_strand(ccdbg, head_kmer_arr, ORF_node_map));
    }

    // group colours of ORFs together
    ORFNodeRobMap ORF_map = std::move(sort_ORF_indexes(ORF_node_map, pos_strand_map, ccdbg, head_kmer_arr, is_ref));

    return ORF_map;
}