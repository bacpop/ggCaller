//
// Created by sth19 on 28/05/2021.
//

#include "graph.h"

GraphTuple Graph::build (const std::string& infile1,
                    const int kmer,
                    const std::vector<std::string>& stop_codons_for,
                    const std::vector<std::string>& stop_codons_rev,
                    size_t num_threads,
                    bool is_ref,
                    const bool write_graph,
                    const std::string& infile2) {
    // Set number of threads
    if (num_threads < 1)
    {
        num_threads = 1;
    }

    // read in compact coloured DBG
    cout << "Building coloured compacted DBG..." << endl;

    // set OMP number of threads
    omp_set_num_threads(num_threads);

    // initialise persistent variables
    GraphPair graph_pair;
    int overlap = kmer - 1;
    size_t nb_colours;
    std::vector<std::string> input_colours;
    NodeColourVector node_colour_vector;

    if (infile2 != "NA") {
        is_ref = 0;
    }

    // scope for ccdbg
    {
        // generate graph, writing if write_graph == true
        size_t lastindex = infile1.find_last_of(".");
        std::string outgraph = infile1.substr(0, lastindex);
        ColoredCDBG<> ccdbg = buildGraph(infile1, infile2, is_ref, kmer, num_threads, false, write_graph, outgraph);

        // get the number of colours
        nb_colours = ccdbg.getNbColors();

        // get colour names
        input_colours = ccdbg.getColorNames();

        // generate codon index for graph
        cout << "Generating graph stop codon index..." << endl;
        node_colour_vector = std::move(_index_graph(ccdbg, stop_codons_for, stop_codons_rev, kmer, nb_colours));
    }

    // make tuple containing all information needed in python back-end
    GraphTuple graph_tuple = std::make_tuple(node_colour_vector, input_colours, nb_colours, overlap);

    return graph_tuple;
}

// read existing graph and index
GraphTuple Graph::read (const std::string& graphfile,
                    const std::string& coloursfile,
                    const std::vector<std::string>& stop_codons_for,
                    const std::vector<std::string>& stop_codons_rev,
                    size_t num_threads,
                    const bool is_ref) {

    // Set number of threads
    if (num_threads < 1)
    {
        num_threads = 1;
    }

    // read in compact coloured DBG
    cout << "Reading coloured compacted DBG..." << endl;

    // set OMP number of threads
    omp_set_num_threads(num_threads);

    // initialise persistent variables
    int kmer;
    int overlap;
    size_t nb_colours;
    std::vector<std::string> input_colours;
    NodeColourVector node_colour_vector;

    // scope for ccdbg
    {
        // read in graph
        ColoredCDBG<> ccdbg;
        ccdbg.read(graphfile, coloursfile, num_threads);

        //set local variables
        kmer = ccdbg.getK();
        overlap = kmer - 1;

        // get the number of colours
        nb_colours = ccdbg.getNbColors();

        // get colour names
        input_colours = ccdbg.getColorNames();

        // generate codon index for graph
        cout << "Generating graph stop codon index..." << endl;
        node_colour_vector = std::move(_index_graph(ccdbg, stop_codons_for, stop_codons_rev, kmer, nb_colours));
    }

    // make tuple containing all information needed in python back-end
    GraphTuple graph_tuple = std::make_tuple(node_colour_vector, input_colours, nb_colours, overlap);

    return graph_tuple;
}

std::pair<ORFOverlapMap, ORFVector> Graph::findORFs (const size_t& colour_ID,
                                                     const std::vector<size_t>& node_ids,
                                                     const bool& repeat,
                                                     const size_t& overlap,
                                                     const size_t& max_path_length,
                                                     bool& is_ref,
                                                     const bool& no_filter,
                                                     const std::vector<std::string>& stop_codons_for,
                                                     const std::vector<std::string>& start_codons_for,
                                                     const size_t min_ORF_length,
                                                     const size_t max_overlap,
                                                     const bool write_idx,
                                                     const std::string& FM_fasta_file)
{
    ORFVector ORF_vector;
    // traverse graph, set scope for all_paths and fm_idx
    {
        // recursive traversal
        //cout << "Traversing graph: " << to_string(colour_ID) << endl;
        AllPaths all_paths = traverse_graph(_GraphVector, colour_ID, node_ids, repeat, max_path_length);

        // if no FM_fasta_file specified, cannot generate FM Index
        if (FM_fasta_file == "NA")
        {
            is_ref = false;
        }

        // generate FM_index if is_ref
        fm_index_coll fm_idx;
        if (is_ref)
        {
            fm_idx = index_fasta(FM_fasta_file, write_idx);
        }

        // generate ORF calls
        //cout << "Calling ORFs: " << to_string(colour_ID) << endl;
        ORF_vector = call_ORFs(all_paths, _GraphVector, stop_codons_for, start_codons_for, overlap, min_ORF_length, is_ref, fm_idx);
    }

    // if no filtering required, do not calculate overlaps
    ORFOverlapMap ORF_overlap_map;
    if (!no_filter)
    {
        // << "Determining overlaps: " << to_string(colour_ID) << endl;
        ORF_overlap_map = std::move(calculate_overlaps(_GraphVector, ORF_vector, overlap, max_overlap));
    }

    std::pair<ORFOverlapMap, ORFVector> return_pair = std::make_pair(ORF_overlap_map, ORF_vector);

    return return_pair;
}

std::unordered_set<size_t> Graph::add_ORF_info (const size_t& colour_ID,
                                                const std::vector<std::pair<size_t,size_t>>& ORF_IDs,
                                                const ORFVector& ORF_vector)
{
    // generate lists of single-node ORFs for correction of edges
    std::unordered_set<size_t> uninode_ORFs;

    for (const auto & ID_pair : ORF_IDs)
    {
        // unpack pair for source and sink nodes
        const auto & source = ID_pair.first;
        const auto & sink = ID_pair.second;

        {
            // get graph information for source node
            const auto & ORF_info = ORF_vector.at(source);

            // get start node IDs (-1 as graph is zero-based). Only add 5 prime node.
            const auto source_node_id = abs(std::get<0>(ORF_info).at(0)) - 1;
            const auto sink_node_id = abs(std::get<0>(ORF_info).back()) - 1;

            // add ORF information to graph
            _GraphVector.at(source_node_id).set_ORFs(colour_ID, source);
            // check if ORF is present only on single node
            if (source_node_id != sink_node_id)
            {
                _GraphVector.at(sink_node_id).set_ORFs(colour_ID, source);
            } else
            {
                // ORF must be uninode as source and sink nodes match
                uninode_ORFs.insert(source);
            }
        }

        // check if source and sink ORFs are the same, if not continue
        if (source != sink)
        {
            // get graph information for source node
            const auto & ORF_info = ORF_vector.at(sink);

            // get start and end node IDs (-1 as graph is zero-based)
            const auto source_node_id = abs(std::get<0>(ORF_info).at(0)) - 1;
            const auto sink_node_id = abs(std::get<0>(ORF_info).back()) - 1;

            // add ORF information to graph
            _GraphVector.at(source_node_id).set_ORFs(colour_ID, sink);
            // check if ORF is present only on single node
            if (source_node_id != sink_node_id)
            {
                _GraphVector.at(sink_node_id).set_ORFs(colour_ID, sink);
            } else
            {
                // ORF must be uninode as source and sink nodes match
                uninode_ORFs.insert(sink);
            }
        }
    }
    return uninode_ORFs;
}

std::vector<std::pair<size_t, size_t>> Graph::get_neighbouring_ORFs (const size_t& colour_ID,
                                                                     const std::vector<std::pair<size_t,size_t>>& end_ORFs,
                                                                     const ORFVector& ORF_vector,
                                                                     const std::unordered_set<size_t>& uninode_ORFs)
{
    // initialise pair of vectors (first = upstream of start_ORF, second = downstream of start_ORF)
    std::vector<std::pair<size_t, size_t>> ORF_edges;

    // iterate over each entry in end_ORFs
    for (const auto& end_ORF_pair : end_ORFs)
    {
        // scope for first item in end_ORF_pair
        {
            const auto& start_ORF = end_ORF_pair.first;

            // get ORF info
            const auto & start_ORF_info = ORF_vector.at(start_ORF);

            // get id for source_node and sink node for ORF
            auto source_node_id = std::get<0>(start_ORF_info).at(0);
            auto sink_node_id = std::get<0>(start_ORF_info).back();

            // get source_node strand and node info
            const auto& source_node_info = _GraphVector.at(abs(source_node_id) - 1);

            // check if there are any overlapping ORFs on current node and order...
            const auto& source_node_ORFs = source_node_info.get_ORFs(colour_ID);

            // get start and end nodes of the path for upstream and downstream traversal, assign to start ORF
            size_t upstream_source = start_ORF;
            size_t downstream_source = start_ORF;

            // check if there are overlapping ORFs on the start node and order
            if (source_node_ORFs.size() >= 2)
            {
                // get strand of start_ORF
                const auto& start_ORF_strand = std::get<5>(start_ORF_info);

                // if 2 or more ORFs overlapping, order and then add each connecting edge
                const auto ordered_ORFs = order_ORFs(_GraphVector, source_node_ORFs, source_node_id, ORF_vector);
                for (size_t i = 0; i < ordered_ORFs.size() - 1; i++)
                {
                    ORF_edges.push_back({ordered_ORFs.at(i), ordered_ORFs.at(i + 1)});
                }

                // assign upstream_source and downstream_source, as well as the new sink and source nodes for upstream/downstream traversal
                upstream_source = ORF_edges.at(0).first;
                const auto& upstream_ORF_info = ORF_vector.at(upstream_source);
                const auto& upstream_ORF_strand = std::get<5>(upstream_ORF_info);
                sink_node_id = std::get<0>(upstream_ORF_info).back();
                downstream_source = ORF_edges.back().second;
                const auto & downstream_ORF_info = ORF_vector.at(downstream_source);
                const auto& downstream_ORF_strand = std::get<5>(upstream_ORF_info);
                source_node_id = std::get<0>(upstream_ORF_info).at(0);

                // need to reverse sink node and source IDs if the upstream and downstream ORFs are reversed compared to original start ORF
                if (start_ORF_strand != upstream_ORF_strand)
                {
                    sink_node_id *= -1;
                }
                if (start_ORF_strand != downstream_ORF_strand)
                {
                    source_node_id *= -1;
                }
            }

            // now traverse graph using DFS, finding next ORFs upstream and downstream of sink and source node respectively
            // traverse upstream
            {
                auto next_ORFs = check_next_ORFs(_GraphVector, sink_node_id, upstream_source, colour_ID, -1, ORF_vector, uninode_ORFs);
                ORF_edges.insert(ORF_edges.end(), make_move_iterator(next_ORFs.begin()), make_move_iterator(next_ORFs.end()));
            }
            // traverse downstream
            {
                auto next_ORFs = check_next_ORFs(_GraphVector, source_node_id, downstream_source, colour_ID, 1, ORF_vector, uninode_ORFs);
                ORF_edges.insert(ORF_edges.end(), make_move_iterator(next_ORFs.begin()), make_move_iterator(next_ORFs.end()));
            }

            // scope for second item in end_ORF_pair
            if (end_ORF_pair.first != end_ORF_pair.second)
            {
                const auto& start_ORF = end_ORF_pair.second;

                // get ORF info
                const auto & start_ORF_info = ORF_vector.at(start_ORF);

                // get id for source_node and sink node for ORF
                auto source_node_id = std::get<0>(start_ORF_info).at(0);
                auto sink_node_id = std::get<0>(start_ORF_info).back();

                // get source_node strand and node info
                const auto& source_node_info = _GraphVector.at(abs(source_node_id) - 1);

                // check if there are any overlapping ORFs on current node and order...
                const auto& source_node_ORFs = source_node_info.get_ORFs(colour_ID);

                // get start and end nodes of the path for upstream and downstream traversal, assign to start ORF
                size_t upstream_source = start_ORF;
                size_t downstream_source = start_ORF;

                // check if there are overlapping ORFs on the start node and order
                if (source_node_ORFs.size() >= 2)
                {
                    // get strand of start_ORF
                    const auto& start_ORF_strand = std::get<5>(start_ORF_info);

                    // if 2 or more ORFs overlapping, order and then add each connecting edge
                    const auto ordered_ORFs = order_ORFs(_GraphVector, source_node_ORFs, source_node_id, ORF_vector);
                    for (size_t i = 0; i < ordered_ORFs.size() - 1; i++)
                    {
                        ORF_edges.push_back({ordered_ORFs.at(i), ordered_ORFs.at(i + 1)});
                    }

                    // assign upstream_source and downstream_source, as well as the new sink and source nodes for upstream/downstream traversal
                    upstream_source = ORF_edges.at(0).first;
                    const auto& upstream_ORF_info = ORF_vector.at(upstream_source);
                    const auto& upstream_ORF_strand = std::get<5>(upstream_ORF_info);
                    sink_node_id = std::get<0>(upstream_ORF_info).back();
                    downstream_source = ORF_edges.back().second;
                    const auto & downstream_ORF_info = ORF_vector.at(downstream_source);
                    const auto& downstream_ORF_strand = std::get<5>(upstream_ORF_info);
                    source_node_id = std::get<0>(upstream_ORF_info).at(0);

                    // need to reverse sink node and source IDs if the upstream and downstream ORFs are reversed compared to original start ORF
                    if (start_ORF_strand != upstream_ORF_strand)
                    {
                        sink_node_id *= -1;
                    }
                    if (start_ORF_strand != downstream_ORF_strand)
                    {
                        source_node_id *= -1;
                    }
                }

                // now traverse graph using DFS, finding next ORFs upstream and downstream of sink and source node respectively
                // traverse upstream
                {
                    auto next_ORFs = check_next_ORFs(_GraphVector, sink_node_id, upstream_source, colour_ID, -1, ORF_vector, uninode_ORFs);
                    ORF_edges.insert(ORF_edges.end(), make_move_iterator(next_ORFs.begin()), make_move_iterator(next_ORFs.end()));
                }
                // traverse downstream
                {
                    auto next_ORFs = check_next_ORFs(_GraphVector, source_node_id, downstream_source, colour_ID, 1, ORF_vector, uninode_ORFs);
                    ORF_edges.insert(ORF_edges.end(), make_move_iterator(next_ORFs.begin()), make_move_iterator(next_ORFs.end()));
                }
            }
        }
    }
    return ORF_edges;
}

std::string Graph::generate_sequence(const std::vector<int>& nodelist,
                                     const std::vector<indexPair>& node_coords,
                                     const size_t& overlap)
{
    std::string sequence;
    for (size_t i = 0; i < nodelist.size(); i++)
    {
        // initialise sequence items
        std::string unitig_seq;
        std::string substring;

        // parse information
        const auto& id = nodelist[i];
        const auto& coords = node_coords[i];
        bool strand = (id >= 0) ? true : false;

        if (strand)
        {
            unitig_seq = _GraphVector.at(abs(id) - 1).seq();
        } else {
            unitig_seq = reverse_complement(_GraphVector.at(abs(id) - 1).seq());
        }

        if (sequence.empty())
        {
            // get node_seq_len, add one as zero indexed
            int node_seq_len = (std::get<1>(coords) - std::get<0>(coords)) + 1;
            substring = unitig_seq.substr(std::get<0>(coords), node_seq_len);
        } else
        {
            // get node_seq_len, add one as zero indexed
            int node_seq_len = (std::get<1>(coords) - overlap) + 1;
            // need to account for overlap, if overlap is greater than the end of the node, sequence already accounted for
            if (node_seq_len > 0)
            {
                substring = unitig_seq.substr(overlap, node_seq_len);
            }
        }
        sequence += substring;
    }
    return sequence;
}

NodeColourVector Graph::_index_graph (const ColoredCDBG<>& ccdbg,
                         const std::vector<std::string>& stop_codons_for,
                         const std::vector<std::string>& stop_codons_rev,
                         const int& kmer,
                         const size_t& nb_colours)
{
    // get all head kmers for parrellelisation
    std::vector<Kmer> head_kmer_arr;
    for (const auto& um : ccdbg)
    {
        head_kmer_arr.push_back(um.getUnitigHead());
    }

    // structures for results
    GraphVector graph_vector(head_kmer_arr.size());
    NodeColourVector node_colour_vector(nb_colours);
    robin_hood::unordered_map<std::string, size_t> head_kmer_map;

    // run unitig indexing in parallel
    size_t unitig_id = 1;
    #pragma omp parallel
    {
        GraphVector graph_vector_private;
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
                if (unitig_dict.full_colour().at(i))
                {
                    node_colour_vector_private[i].push_back(unitig_dict.id);
                }
            }

            // add head_kmer and unitig id to map
            head_kmer_map_private[unitig_dict.head_kmer()] = unitig_dict.id;

            // add unitig to graph_vector, minus 1 as zero based
            graph_vector[unitig_dict.id - 1] = std::move(unitig_dict);
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
    update_neighbour_index(graph_vector, head_kmer_map);

    // assign the graph vector to the graph _GraphVector
    _GraphVector = std::move(graph_vector);

    // return node_colour vector
    return node_colour_vector;
}