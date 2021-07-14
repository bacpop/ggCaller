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
    std::pair<ORFVector, NodeStrandMap> ORF_pair;
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
        ORF_pair = call_ORFs(all_paths, _GraphVector, stop_codons_for, start_codons_for, overlap, min_ORF_length, is_ref, fm_idx);
    }

    // if no filtering required, do not calculate overlaps
    ORFOverlapMap ORF_overlap_map;
    if (!no_filter)
    {
        // << "Determining overlaps: " << to_string(colour_ID) << endl;
        ORF_overlap_map = std::move(calculate_overlaps(_GraphVector, ORF_pair, overlap, max_overlap));
    }

    std::pair<ORFOverlapMap, ORFVector> return_pair = std::make_pair(ORF_overlap_map, ORF_pair.first);

    return return_pair;
}

void Graph::add_ORF_info (const size_t& colour_ID,
                          const std::vector<std::pair<size_t,size_t>>& ORF_IDs,
                          const ORFVector& ORF_vector)
{
    for (const auto & ID_pair : ORF_IDs)
    {
        // unpack pair for source and sink nodes
        const auto & source = ID_pair.first;
        const auto & sink = ID_pair.second;

        {
            // get graph information for source node
            const auto & ORF_info = ORF_vector.at(source);

            // get start and end node IDs (-1 as graph is zero-based)
            const auto& source_node_id = abs(std::get<0>(ORF_info).at(0)) - 1;
            const auto& sink_node_id = abs(std::get<0>(ORF_info).back()) - 1;

            // add ORF information to graph
            _GraphVector.at(source_node_id).set_ORFs(colour_ID, source_node_id);
            _GraphVector.at(source_node_id).set_ORFs(colour_ID, sink_node_id);
        }

        // check if source and sink ORFs are the same, if not continue

        if (source != sink)
        {
            // get graph information for source node
            const auto & ORF_info = ORF_vector.at(sink);

            // get start and end node IDs (-1 as graph is zero-based)
            const auto& source_node_id = abs(std::get<0>(ORF_info).at(0)) - 1;
            const auto& sink_node_id = abs(std::get<0>(ORF_info).back()) - 1;

            // add ORF information to graph
            _GraphVector.at(source_node_id).set_ORFs(colour_ID, source_node_id);
            _GraphVector.at(source_node_id).set_ORFs(colour_ID, sink_node_id);
        }
    }
}

std::vector<std::pair<size_t, size_t> Graph::get_neighbouring_ORFs (const size_t& colour_ID,
                                                                     const std::vector<std::pair<size_t,size_t>>& ORF_IDs,
                                                                     const ORFVector& ORF_vector)
{
    std::vector<std::pair<size_t, size_t> ORF_pair_vector;

    for (const auto & ID_pair : ORF_IDs)
    {
        // unpack pair for source and sink nodes
        const auto & source = ID_pair.first;
        const auto & sink = ID_pair.second;

        {
            // get graph information for source node
            const auto & ORF_info = ORF_vector.at(source);

            // get source node IDs (as we are at source, want to look upstream)
            const auto& id = std::get<0>(ORF_info).at(0);
            // get real node id (-1 as zero indexed)
            const auto node_id = abs(id) - 1;

            // check if another ORF overlaps same node
            const auto & node_ORFs = _GraphVector.at(node_id).get_ORFs(colour_ID);

            // if there is only one ORF, then need to traverse upstream
            if (node_ORFs.size() == 1)
            {
                const auto node_vector = check_upstream_ORFs(_GraphVector, id, colour_ID);

                // order the ORFs in the node_vector, then place last and current head node together
                const auto ordered_ORFs = order_node_ends(_GraphVector, node_vector, id);

                ORF_pair_vector.push_back(std::pair<size_t, size_t>(ordered_ORFs.back().second, source));
            }
            // if there is only two ORFs, set this ORF as upstream
            else if (node_ORFs.size() == 2)
            {
                // iterate over the ORFs in the node_ORFs
                for (const auto& ORF_ID : node_ORFs)
                {
                    if (ORF_ID != source)
                    {
                        // add to return vector (first upstream node, then source node)
                        ORF_pair_vector.push_back(std::pair<size_t, size_t>(ORF_ID, source));
                    }
                }
            }
            // if there is more than one ORF, then need to check which is directly upstream
            else if (node_ORFs.size() > 2)
            {
                const auto& ordered_ORFs = order_node_ends(_GraphVector, node_ORFs, id);

                // parse the sorted vector and append to return vector (pairing this and next item, so don't traverse full vector)
                for (int i = 0; i < ordered_ORFs.size() - 1; i++)
                {
                    ORF_pair_vector.push_back(std::pair<size_t, size_t>(ordered_ORFs.at(i).second, ordered_ORFs.at(i + 1).second));
                }
            }
        }
        // do the same for the sink node
        {
            // get graph information for source node
            const auto & ORF_info = ORF_vector.at(sink);

            // get source node IDs (as we are at source, want to look upstream)
            const auto& id = std::get<0>(ORF_info).at(0);
            // get real node id (-1 as zero indexed)
            const auto node_id = abs(id) - 1;

            // check if another ORF overlaps same node
            const auto & node_ORFs = _GraphVector.at(node_id).get_ORFs(colour_ID);

            // if there is only one ORF, then need to traverse upstream
            if (node_ORFs.size() == 1)
            {
                const auto node_vector = check_upstream_ORFs(_GraphVector, id, colour_ID);

                // order the ORFs in the node_vector, then place last and current head node together
                const auto ordered_ORFs = order_node_ends(_GraphVector, node_vector, id);

                ORF_pair_vector.push_back(std::pair<size_t, size_t>(sink, ordered_ORFs.at(0).second));
            }
            // if there is only two ORFs, set this ORF as upstream
            else if (node_ORFs.size() == 2)
            {
                // iterate over the ORFs in the node_ORFs
                for (const auto& ORF_ID : node_ORFs)
                {
                    if (ORF_ID != sink)
                    {
                        // add to return vector (first upstream node, then source node)
                        ORF_pair_vector.push_back(std::pair<size_t, size_t>(sink, ORF_ID));
                    }
                }
            }
            // if there is more than one ORF, then need to check which is directly upstream
            else if (node_ORFs.size() > 2)
            {
                const auto& ordered_ORFs = order_node_ends(_GraphVector, node_ORFs, id);

                // parse the sorted vector and append to return vector (pairing this and next item, so don't traverse full vector)
                for (int i = 0; i < ordered_ORFs.size() - 1; i++)
                {
                    ORF_pair_vector.push_back(std::pair<size_t, size_t>(ordered_ORFs.at(i).second, ordered_ORFs.at(i + 1).second));
                }
            }
        }
    }
    return ORF_pair_vector;
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