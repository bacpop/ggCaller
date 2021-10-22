#include "gene_refinding.h"

void assign_seq(const GraphVector& graph_vector,
                const PathVector& unitig_complete_paths,
                const int kmer,
                const bool is_ref,
                const fm_index_coll& fm_idx,
                std::string& stream_seq)
{
    // iterate over all the paths, determine which is the longest and real.
    // If multiple, choose from that with the lowest hash for first kmer
    for (const auto& unitig_path : unitig_complete_paths)
    {
        // initilise path_sequence
        std::string path_sequence;

        // generate the path sequence
        for (const auto &node : unitig_path)
        {
            std::string unitig_seq;
            std::string new_sequence = path_sequence;

            // parse out information from node integer value
            bool strand = (node >= 0) ? true : false;

            if (strand) {
                unitig_seq = graph_vector.at(abs(node) - 1).seq();
            } else {
                unitig_seq = reverse_complement(graph_vector.at(abs(node) - 1).seq());
            }

            if (new_sequence.empty())
            {
                new_sequence = unitig_seq;
            } else
            {
                new_sequence.append(unitig_seq.begin() + kmer - 1, unitig_seq.end());
            }

            // check against FMindex
            // check new_sequence is real if is_ref
            if (is_ref)
            {
                const bool present = check_colours(new_sequence, fm_idx);

                // check if real sequence, if not pass on the ORF, move to next highest
                if (!present)
                {
                    break;
                }
            }

            // assign new_sequence to path_sequence
            path_sequence = new_sequence;
        }

        // compare the path_sequence to current stream_seq
        if (path_sequence.size() > stream_seq.size())
        {
            stream_seq = path_sequence;
        } else if (path_sequence.size() == stream_seq.size() && path_sequence != stream_seq)
        {
            // if equal size, get the hash of the last kmer in each and assign to highest
            const size_t path_hash = hasher{}(path_sequence.substr(path_sequence.size() - kmer, kmer));
            const size_t stream_hash = hasher{}(stream_seq.substr(stream_seq.size() - kmer, kmer));

            if (path_hash > stream_hash)
            {
                stream_seq = path_sequence;
            }
        }
    }
}

PathVector iter_nodes_length (const GraphVector& graph_vector,
                              const NodeTuple& head_node_tuple,
                              const size_t& current_colour,
                              const size_t& radius,
                              const bool& repeat)
{
    // generate path list, vector for path and the stack
    PathVector path_list;
    std::vector<int> node_vector;
    NodeStack node_stack;

    // create node set for identification of repeats
    std::unordered_set<int> node_set;
    node_set.insert(std::get<1>(head_node_tuple));

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
        const sdsl::bit_vector & colour_arr = std::get<3>(node_tuple);
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

        // get unitig_dict entry in graph_vector
        const auto& node_dict = graph_vector.at(abs(node_id) - 1);

        // determine strand of unitig
        const bool strand = (node_id >= 0) ? true : false;

        // iterate over neighbours, recurring through incomplete paths
        for (const auto& neighbour : node_dict.get_neighbours(strand))
        {
            // parse neighbour information. Frame is next stop codon, with first dictating orientation and second the stop codon index
            const auto& neighbour_id = neighbour.first;
            const auto& frame = neighbour.second;

            // check if unitig has already been traversed, and pass if repeat not specified
            const bool is_in = node_set.find(neighbour_id) != node_set.end();
            if (!repeat && is_in)
            {
                continue;
            }

            // get reference to unitig_dict object for neighbour
            const auto& neighbour_dict = graph_vector.at(abs(neighbour_id) - 1);

            // calculate colours array
            auto updated_colours_arr = colour_arr;
            updated_colours_arr &= neighbour_dict.full_colour();

            // determine if neighbour is in same colour as iteration, if not return and pass
            if (!updated_colours_arr[current_colour])
            {
                continue;
            }

            // check path length, if too long continue
            const size_t updated_path_length = path_length + neighbour_dict.size().second;

            // if adding new node pushes path length over radius or neighbour is end contig, add neighbour and return
            if (updated_path_length >= radius || neighbour_dict.end_contig())
            {
                // create temporary path to account for reaching end of contig
                std::vector<int> return_path = node_vector;
                return_path.push_back(neighbour_id);

                // update path_list to enable this to be returned
                path_list.push_back(std::move(return_path));
                continue;
            }

            // if no previous conditions are satisfied, prepare tuple for stack
            NodeTuple new_node_tuple(new_pos_idx, neighbour_id, 0, updated_colours_arr, updated_path_length);

            // add to stack
            node_stack.push(new_node_tuple);

            // add node to node_set
            node_set.insert(neighbour_id);
        }
    }
    return path_list;
}

std::pair<std::string, std::string> traverse_outward(const GraphVector& graph_vector,
                                                     const size_t& colour_ID,
                                                     const ORFNodeVector& ORF_info,
                                                     const size_t& radius,
                                                     bool is_ref,
                                                     const bool write_idx,
                                                     const int kmer,
                                                     const std::string& FM_fasta_file,
                                                     const bool repeat)
{
    // initialise upstream and downstream strings
    std::string upstream_seq;
    std::string downstream_seq;

    // if no FM_fasta_file specified, cannot generate FM Index
    if (FM_fasta_file == "NA")
    {
        is_ref = false;
    }

    fm_index_coll fm_idx;
    if (is_ref)
    {
        fm_idx = index_fasta(FM_fasta_file, write_idx);
    }

    // traverse end node in forward direction
    {
        // get the node id to traverse from
        const int& head_id = std::get<0>(ORF_info).back();

        // parse unitig_id. Zero based, so take 1
        const auto unitig_id = abs(head_id) - 1;

        // get reference to unitig_dict object
        const auto& unitig_dict = graph_vector.at(unitig_id);

        // gather unitig information from graph_vector
        const uint8_t codon_arr = 0;
        const sdsl::bit_vector colour_arr = unitig_dict.full_colour();

        // calculate where in unitig ORF sits, to determine initial length of path
        size_t path_length = unitig_dict.size().first - std::get<1>(ORF_info).back().second;

        // check if path_length already exceeds radius
        PathVector unitig_complete_paths;
        if (path_length >= radius)
        {
            unitig_complete_paths.push_back({head_id});
        } else
        {
            // generate node tuple for iteration
            NodeTuple head_node_tuple(0, head_id, codon_arr, colour_arr, path_length);

            // recur paths
            unitig_complete_paths = iter_nodes_length(graph_vector, head_node_tuple, colour_ID, radius, repeat);
        }

        if (!unitig_complete_paths.empty())
        {
            assign_seq(graph_vector, unitig_complete_paths, kmer, is_ref, fm_idx, downstream_seq);
        }
    }

    // traverse start node in reverse direction
    {
        // get node_id to traverse from, parse unitig_id from TIS if present.
        const bool TIS_present = (!std::get<3>(ORF_info).empty()) ? true : false;
        const int head_id = ((TIS_present) ? std::get<3>(ORF_info).at(0) : std::get<0>(ORF_info).at(0)) * -1;

        // get unitig_ID Zero based, so take absolute value and take 1
        const auto unitig_id = abs(head_id) - 1;

        // get reference to unitig_dict object
        const auto& unitig_dict = graph_vector.at(unitig_id);

        // gather unitig information from graph_vector
        const uint8_t codon_arr = 0;
        const sdsl::bit_vector colour_arr = unitig_dict.full_colour();

        // calculate where in unitig ORF sits, to determine initial length of path
        size_t path_length;
        if (TIS_present)
        {
            path_length = std::get<4>(ORF_info).at(0).first + 16;
        } else
        {
            path_length = std::get<1>(ORF_info).at(0).first;
        }

        // check if path_length already exceeds radius
        PathVector unitig_complete_paths;
        if (path_length >= radius)
        {
            unitig_complete_paths.push_back({head_id});
        } else
        {
            // generate node tuple for iteration
            NodeTuple head_node_tuple(0, head_id, codon_arr, colour_arr, path_length);

            // recur paths
            unitig_complete_paths = iter_nodes_length(graph_vector, head_node_tuple, colour_ID, radius, repeat);
        }

        if (!unitig_complete_paths.empty())
        {
            assign_seq(graph_vector, unitig_complete_paths, kmer, is_ref, fm_idx, upstream_seq);
        }
    }

    return {upstream_seq, downstream_seq};
}