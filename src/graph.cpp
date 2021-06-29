//
// Created by sth19 on 28/05/2021.
//

#include "graph.h"
// define mutex
std::mutex mtx3;

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

        // generate codon index for graph and resize _ColourGraphPaths, add to graph object
        cout << "Generating graph stop codon index..." << endl;
        _index_graph(ccdbg, stop_codons_for, stop_codons_rev, kmer, nb_colours);
        _ColourGraphPaths.resize(nb_colours);
        _NodeColourVectorTraversed.resize(nb_colours);
    }

    // make tuple containing all information needed in python back-end
    GraphTuple graph_tuple = std::make_tuple(input_colours, nb_colours, overlap);

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

        // generate codon index for graph and resize _ColourGraphPaths, add to graph object
        cout << "Generating graph stop codon index..." << endl;
        _index_graph(ccdbg, stop_codons_for, stop_codons_rev, kmer, nb_colours);
        _ColourGraphPaths.resize(nb_colours);
        _NodeColourVectorTraversed.resize(nb_colours);
    }

    // make tuple containing all information needed in python back-end
    GraphTuple graph_tuple = std::make_tuple(input_colours, nb_colours, overlap);

    return graph_tuple;
}

std::pair<ORFOverlapMap, ORFVector> Graph::findORFs (const size_t& colour_ID,
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

    // get reference to _node_ids for specific colour
    auto & node_ids = _NodeColourVector.at(colour_ID);

    // traverse graph, set scope for all_paths and fm_idx
    {
        // recursive traversal
        //cout << "Traversing graph: " << to_string(colour_ID) << endl;
        _traverse_graph(colour_ID, repeat, max_path_length);

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

        //cout << "Calling ORFs: " << to_string(colour_ID) << endl;
        ORF_pair = call_ORFs(_ColourGraphPaths.at(colour_ID), _GraphVector, stop_codons_for, start_codons_for, overlap, min_ORF_length, is_ref, fm_idx);
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

void Graph::_index_graph (const ColoredCDBG<>& ccdbg,
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
                    // check if forward and reverse stops present
                    if (unitig_dict.forward_stop())
                    {
                        node_colour_vector_private[i].insert(((int) unitig_dict.id));
                    }
                    if (unitig_dict.reverse_stop())
                    {
                        node_colour_vector_private[i].insert((((int) unitig_dict.id)) * -1);
                    }
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
                node_colour_vector[i].merge(node_colour_vector_private[i]);
            }
        }
    }
    // update neighbour index in place within graph_vector
    update_neighbour_index(graph_vector, head_kmer_map);

    // assign the graph vector to the graph _GraphVector and colour vector to _NodeColourVector
    _GraphVector = std::move(graph_vector);
    _NodeColourVector = std::move(node_colour_vector);
}

void Graph::_traverse_graph(const size_t& colour_ID,
                     const bool repeat,
                     const size_t max_path_length)
{
    // get reference to position in _ColourGraphPaths
    auto & all_paths = _ColourGraphPaths[colour_ID];

    // get reference to position in _NodeColourVector
    auto & node_ids = _NodeColourVector[colour_ID];

    // get reference to position in _NodeColourVectorTraversed
    auto & traversed_node_ids = _NodeColourVectorTraversed[colour_ID];

    // need to think about a mutex here to prevent node from interferring with traversing colour
    // traverse nodes
    for (const auto& head_id : node_ids)
    {
        // check if head has already been traversed by another thread
        mtx3.lock();
        if (traversed_node_ids.find(head_id) != traversed_node_ids.end())
        {
            mtx3.unlock();
            continue;
        }
        mtx3.unlock();

        // parse unitig_id. Zero based, so get as positive and then zero base
        int unitig_id;
        bool pos;
        if (head_id > 0)
        {
            unitig_id = head_id - 1;
            pos = true;
        } else
        {
            unitig_id = (head_id * -1) - 1;
            pos = false;
        }

        // get reference to unitig_dict object
        const auto& unitig_dict = _GraphVector.at(unitig_id);


        // gather unitig information from graph_vector
        const uint8_t codon_arr = unitig_dict.get_codon_arr(true, pos, 0);
        const size_t unitig_len = unitig_dict.size().first;
        const std::vector<bool> colour_arr = unitig_dict.full_colour();

        // generate node tuple for iteration
        NodeTuple head_node_tuple(0, head_id, codon_arr, colour_arr, unitig_len);

        // recur paths
        iter_nodes_binary(_GraphVector, _NodeColourVectorTraversed, _ColourGraphPaths, head_node_tuple, colour_ID, max_path_length, repeat);
    }
    // clear node_ids after full traversal
    node_ids.clear();
}