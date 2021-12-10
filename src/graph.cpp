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
        node_colour_vector = std::move(_index_graph(ccdbg, stop_codons_for, stop_codons_rev, kmer, nb_colours, is_ref, input_colours));
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
        node_colour_vector = std::move(_index_graph(ccdbg, stop_codons_for, stop_codons_rev, kmer, nb_colours, is_ref, input_colours));
    }

    // make tuple containing all information needed in python back-end
    GraphTuple graph_tuple = std::make_tuple(node_colour_vector, input_colours, nb_colours, overlap);

    return graph_tuple;
}

void Graph::in(std::string infile)
{
    GraphVector newg;
    std::ifstream ifs(infile);
    boost::archive::text_iarchive ia(ifs);
    ia >> newg;

    _GraphVector = newg;

    for (const auto& node : _GraphVector)
    {
        _KmerMap[node.head_kmer()] = node.id;
    }
}

void Graph::out(std::string outfile)
{
    std::ofstream ofs(outfile);
    boost::archive::text_oarchive oa(ofs);
    // write class instance to archive
    oa << _GraphVector;
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
//        cout << "Traversing graph: " << to_string(colour_ID) << endl;
        std::vector<PathVector> all_paths = traverse_graph(_GraphVector, colour_ID, node_ids, repeat, max_path_length, is_ref);

        // if no FM_fasta_file specified, cannot generate FM Index
        if (FM_fasta_file == "NA")
        {
            is_ref = false;
        }

        // generate FM_index if is_ref
        fm_index_coll fm_idx;

        if (is_ref)
        {
//            cout << "FM-indexing: " << to_string(colour_ID) << endl;
            const auto idx_file_name = FM_fasta_file + ".fm";
            if (!load_from_file(fm_idx, idx_file_name))
            {
                cout << "FM-Index not available for " << FM_fasta_file << endl;
                is_ref = false;
            }
        }

        // generate ORF calls
//        cout << "Calling ORFs: " << to_string(colour_ID) << endl;
        ORF_vector = call_ORFs(colour_ID, all_paths, _GraphVector, stop_codons_for, start_codons_for, overlap, min_ORF_length, is_ref, fm_idx);
    }

    // if no filtering required, do not calculate overlaps
    ORFOverlapMap ORF_overlap_map;
    if (!no_filter)
    {
//        cout << "Determining overlaps: " << to_string(colour_ID) << endl;
        ORF_overlap_map = std::move(calculate_overlaps(_GraphVector, ORF_vector, overlap, max_overlap));
    }

    std::pair<ORFOverlapMap, ORFVector> return_pair = std::make_pair(ORF_overlap_map, ORF_vector);

    return return_pair;
}

std::vector<std::pair<size_t, size_t>> Graph::connect_ORFs(const size_t& colour_ID,
                                                           const ORFVector& ORF_vector,
                                                           const std::vector<size_t>& target_ORFs,
                                                           const size_t& max_ORF_path_length,
                                                           const bool is_ref)
{
    std::vector<std::pair<size_t, size_t>> connected_ORFs;

    // add ORF info for colour to graph
    add_ORF_info(_GraphVector, colour_ID, target_ORFs, ORF_vector);

    // initialise prev_node_set to avoid same ORFs being traversed from again
    std::unordered_set<int> prev_node_set;

    // conduct DBG traversal for upstream...
    auto new_connections = pair_ORF_nodes(_GraphVector, colour_ID, target_ORFs, ORF_vector, max_ORF_path_length, -1, prev_node_set, is_ref);
    connected_ORFs.insert(connected_ORFs.end(), make_move_iterator(new_connections.begin()), make_move_iterator(new_connections.end()));

    // ... and downstream
    new_connections = pair_ORF_nodes(_GraphVector, colour_ID, target_ORFs, ORF_vector, max_ORF_path_length, 1, prev_node_set, is_ref);
    connected_ORFs.insert(connected_ORFs.end(), make_move_iterator(new_connections.begin()), make_move_iterator(new_connections.end()));

    return connected_ORFs;
}

std::pair<ORFMatrixVector, ORFClusterMap> Graph::generate_clusters(const ColourORFMap& colour_ORF_map,
                                                                  const size_t& overlap,
                                                                  const double& id_cutoff,
                                                                  const double& len_diff_cutoff)
{
    // group ORFs together based on single shared k-mer
    auto ORF_group_tuple = group_ORFs(colour_ORF_map, _GraphVector);

    //unpack ORF_group_tuple
    auto& ORF_mat_vector = std::get<0>(ORF_group_tuple);
    auto& ORF_group_vector = std::get<1>(ORF_group_tuple);
    auto& centroid_vector = std::get<2>(ORF_group_tuple);

    // generate clusters for ORFs based on identity
    auto cluster_map = produce_clusters(colour_ORF_map, _GraphVector, overlap, ORF_mat_vector,
                                                   ORF_group_vector, centroid_vector, id_cutoff, len_diff_cutoff);

    // generate return pair of mappings of ORF IDs and clusters
    const auto return_pair = std::make_pair(ORF_mat_vector, cluster_map);
    return return_pair;
}

RefindMap Graph::refind_gene(const size_t& colour_ID,
                             const std::unordered_map<int, std::unordered_map<std::string, ORFNodeVector>>& node_search_dict,
                             const size_t& radius,
                             bool is_ref,
                             const int kmer,
                             const std::string& FM_fasta_file,
                             const bool repeat)
{
    fm_index_coll fm_idx;
    if (is_ref)
    {
        const auto idx_file_name = FM_fasta_file + ".fm";
        if (!load_from_file(fm_idx, idx_file_name))
        {
            cout << "FM-Index not available for " << FM_fasta_file << endl;
            is_ref = false;
        }
    }

    return refind_in_nodes(_GraphVector, colour_ID, node_search_dict, radius, is_ref,
                            kmer, fm_idx, repeat);
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
            } else
            {
                break;
            }
        }
        sequence += substring;
    }
    return sequence;
}

std::tuple<std::vector<std::string>, int, std::vector<MappingCoords>> Graph::search_graph(const std::string& graphfile,
                                                                                          const std::string& coloursfile,
                                                                                          const std::vector<std::string>& query_vec,
                                                                                          const double& id_cutoff,
                                                                                          size_t num_threads)
{
    // Set number of threads
    if (num_threads < 1)
    {
        num_threads = 1;
    }

    // set OMP number of threads
    omp_set_num_threads(num_threads);

    // read in graph
    ColoredCDBG<> ccdbg;
    ccdbg.read(graphfile, coloursfile, num_threads);

    // get input colours
    std::vector<std::string> input_colours = ccdbg.getColorNames();

    //set local variables
    const int kmer = ccdbg.getK();

    std::vector<MappingCoords> query_coords(query_vec.size());

    // go through query, determine head-kmers of each node and map to _GraphVector
//    #pragma omp parallel for
    for (int i = 0; i < query_vec.size(); i++)
    {
        query_coords[i] = std::move(query_DBG(ccdbg, query_vec.at(i), kmer, _KmerMap, id_cutoff));
    }

    return {input_colours, kmer, query_coords};
}

NodeColourVector Graph::_index_graph (const ColoredCDBG<>& ccdbg,
                                     const std::vector<std::string>& stop_codons_for,
                                     const std::vector<std::string>& stop_codons_rev,
                                     const int& kmer,
                                     const size_t& nb_colours,
                                     const bool is_ref,
                                     const std::vector<std::string>& input_colours)
{
    auto node_colour_vector = index_graph(_GraphVector, _KmerMap, ccdbg, stop_codons_for, stop_codons_rev, kmer, nb_colours, is_ref, input_colours);

    // return node_colour vector
    return node_colour_vector;
}