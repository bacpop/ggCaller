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
    int overlap = kmer - 1;

    if (infile2 != "NA") {
        is_ref = 0;
    }

    // generate graph, writing if write_graph == true
    size_t lastindex = infile1.find_last_of(".");
    std::string outgraph = infile1.substr(0, lastindex);
    _ccdbg = buildGraph(infile1, infile2, is_ref, kmer, num_threads, false, write_graph, outgraph);

    // get the number of colours
    size_t nb_colours = _ccdbg.getNbColors();

    // get colour names
    std::vector<std::string> input_colours = _ccdbg.getColorNames();

    // generate codon index for graph
    cout << "Generating graph stop codon index..." << endl;
    NodeColourVector node_colour_vector = std::move(_index_graph(stop_codons_for, stop_codons_rev, kmer, nb_colours, is_ref, input_colours));

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

    // read in graph
    _ccdbg.read(graphfile, coloursfile, num_threads);

    //set local variables
    int kmer = _ccdbg.getK();
    int overlap = kmer - 1;

    // get the number of colours
    size_t nb_colours = _ccdbg.getNbColors();

    // get colour names
    std::vector<std::string> input_colours = _ccdbg.getColorNames();

    // generate codon index for graph
    cout << "Generating graph stop codon index..." << endl;
    NodeColourVector node_colour_vector = std::move(_index_graph(stop_codons_for, stop_codons_rev, kmer, nb_colours, is_ref, input_colours));

    // make tuple containing all information needed in python back-end
    GraphTuple graph_tuple = std::make_tuple(node_colour_vector, input_colours, nb_colours, overlap);

    return graph_tuple;
}

void Graph::in(const std::string& infile,
               const std::string& graphfile,
               const std::string& coloursfile,
               const size_t num_threads)
{
    std::vector<std::string> kmer_array;

    std::ifstream ifs(infile);
    boost::archive::text_iarchive ia(ifs);
    ia >> kmer_array;

    ColoredCDBG<MyUnitigMap> new_ccdbg;

    new_ccdbg.read(graphfile, coloursfile, num_threads);

    _ccdbg = std::move(new_ccdbg);

    _KmerArray.resize(kmer_array.size());

    for (int i = 0; i < kmer_array.size(); i++)
    {
        _KmerArray[i] = Kmer(kmer_array[i].c_str());

        // get a reference to the unitig map object
        auto um_pair = get_um_data(_ccdbg, _KmerArray[i]);
        auto& um = um_pair.first;
        auto& um_data = um_pair.second;

        um_data->set_id(i + 1);
    }
}

void Graph::out(const std::string& outfile)
{
    std::ofstream ofs(outfile);
    boost::archive::text_oarchive oa(ofs);
    // write class instance to archive

    std::vector<std::string> kmer_array(_KmerArray.size());

    // add all Kmers as strings into kmer_array
    for (int i = 0; i < _KmerArray.size(); i++)
    {
        kmer_array[i] = _KmerArray[i].toString();
    }

    oa << kmer_array;
}


std::pair<ORFOverlapMap, ORFVector> Graph::findORFs (const size_t colour_ID,
                                                     const std::vector<size_t>& node_ids,
                                                     const bool repeat,
                                                     const size_t overlap,
                                                     const size_t max_path_length,
                                                     bool is_ref,
                                                     const bool no_filter,
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
        std::vector<PathVector> all_paths = traverse_graph(_ccdbg, _KmerArray, colour_ID, node_ids, repeat, max_path_length, overlap, is_ref);

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
            const auto idx_file_name = FM_fasta_file + ".fmp";
            if (!load_from_file(fm_idx, idx_file_name))
            {
                cout << "FM-Index not available for " << FM_fasta_file << endl;
                is_ref = false;
            }
        }

        // generate ORF calls
//        cout << "Calling ORFs: " << to_string(colour_ID) << endl;
        ORF_vector = call_ORFs(colour_ID, all_paths, _ccdbg, _KmerArray, stop_codons_for, start_codons_for, overlap, min_ORF_length, is_ref, fm_idx);
    }

    // if no filtering required, do not calculate overlaps
    ORFOverlapMap ORF_overlap_map;
    if (!no_filter)
    {
//        cout << "Determining overlaps: " << to_string(colour_ID) << endl;
        ORF_overlap_map = std::move(calculate_overlaps(_ccdbg, _KmerArray, ORF_vector, overlap, max_overlap));
    }

    std::pair<ORFOverlapMap, ORFVector> return_pair = std::make_pair(ORF_overlap_map, ORF_vector);

    return return_pair;
}

std::vector<std::pair<size_t, size_t>> Graph::connect_ORFs(const size_t colour_ID,
                                                           const ORFVector& ORF_vector,
                                                           const std::vector<size_t>& target_ORFs,
                                                           const size_t max_ORF_path_length,
                                                           const bool is_ref,
                                                           const int overlap)
{
    std::vector<std::pair<size_t, size_t>> connected_ORFs;

    // add ORF info for colour to graph
    const auto node_to_ORFs = add_ORF_info(_KmerArray, target_ORFs, ORF_vector);

    // initialise prev_node_set to avoid same ORFs being traversed from again
    std::unordered_set<int> prev_node_set;

    // conduct DBG traversal for upstream...
    auto new_connections = pair_ORF_nodes(_ccdbg, _KmerArray, node_to_ORFs, colour_ID, target_ORFs, ORF_vector, max_ORF_path_length, -1, prev_node_set, is_ref, overlap);
    connected_ORFs.insert(connected_ORFs.end(), make_move_iterator(new_connections.begin()), make_move_iterator(new_connections.end()));

    // ... and downstream
    new_connections = pair_ORF_nodes(_ccdbg, _KmerArray, node_to_ORFs, colour_ID, target_ORFs, ORF_vector, max_ORF_path_length, 1, prev_node_set, is_ref, overlap);
    connected_ORFs.insert(connected_ORFs.end(), make_move_iterator(new_connections.begin()), make_move_iterator(new_connections.end()));

    return connected_ORFs;
}

std::pair<ORFMatrixVector, ORFClusterMap> Graph::generate_clusters(const ColourORFMap& colour_ORF_map,
                                                                   const size_t& overlap,
                                                                   const double& id_cutoff,
                                                                   const double& len_diff_cutoff)
{
    // group ORFs together based on single shared k-mer
    auto ORF_group_tuple = group_ORFs(colour_ORF_map, _KmerArray);

    //unpack ORF_group_tuple
    auto& ORF_mat_vector = std::get<0>(ORF_group_tuple);
    auto& ORF_group_vector = std::get<1>(ORF_group_tuple);
    auto& centroid_vector = std::get<2>(ORF_group_tuple);

    // generate clusters for ORFs based on identity
    auto cluster_map = produce_clusters(colour_ORF_map, _ccdbg, _KmerArray, overlap, ORF_mat_vector,
                                                   ORF_group_vector, centroid_vector, id_cutoff, len_diff_cutoff);

    // generate return pair of mappings of ORF IDs and clusters
    const auto return_pair = std::make_pair(ORF_mat_vector, cluster_map);
    return return_pair;
}

RefindMap Graph::refind_gene(const size_t& colour_ID,
                             const std::unordered_map<int, std::unordered_map<std::string, ORFNodeVector>>& node_search_dict,
                             const size_t radius,
                             bool is_ref,
                             const int kmer,
                             const std::string& FM_fasta_file,
                             const bool repeat)
{
    fm_index_coll fm_idx;
    if (is_ref)
    {
        const auto idx_file_name = FM_fasta_file + ".fmp";
        if (!load_from_file(fm_idx, idx_file_name))
        {
            cout << "FM-Index not available for " << FM_fasta_file << endl;
            is_ref = false;
        }
    }

    return refind_in_nodes(_ccdbg, _KmerArray, colour_ID, node_search_dict, radius, is_ref,
                            kmer, fm_idx, repeat);
}

std::string Graph::generate_sequence(const std::vector<int>& nodelist,
                                     const std::vector<indexPair>& node_coords,
                                     const size_t& overlap)
{
    return generate_sequence_nm(nodelist, node_coords, overlap, _ccdbg, _KmerArray);
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
    _ccdbg.read(graphfile, coloursfile, num_threads);

    // get input colours
    std::vector<std::string> input_colours = _ccdbg.getColorNames();

    //set local variables
    const int kmer = _ccdbg.getK();

    std::vector<MappingCoords> query_coords(query_vec.size());

    // go through query, determine head-kmers of each node and map to _GraphVector
    #pragma omp parallel for
    for (int i = 0; i < query_vec.size(); i++)
    {
        query_coords[i] = std::move(query_DBG(_ccdbg, query_vec.at(i), kmer, id_cutoff));
    }

    return {input_colours, kmer, query_coords};
}

std::vector<std::pair<ContigLoc, bool>> Graph::ORF_location(const std::vector<std::pair<std::vector<int>, std::vector<indexPair>>>& ORF_IDs,
                                                            const std::string& fasta_file,
                                                            const int overlap,
                                                            const bool write_idx,
                                                            size_t num_threads)
{
    // Set number of threads
    if (num_threads < 1)
    {
        num_threads = 1;
    }

    // set OMP number of threads
    omp_set_num_threads(num_threads);

    // initialise return vector
    std::vector<std::pair<ContigLoc, bool>> ORF_coords(ORF_IDs.size());

    // get the FM_index
    const auto fm_index = index_fasta(fasta_file, write_idx);

    #pragma omp parallel for
    for (int i = 0; i < ORF_IDs.size(); i++)
    {
        const auto& ORF_info = ORF_IDs.at(i);
        const auto ORF_sequence = generate_sequence_nm(ORF_info.first, ORF_info.second, overlap, _ccdbg, _KmerArray);

        // get the coordinates of the ORF
        ORF_coords[i] = get_ORF_coords(ORF_sequence, fm_index.first, fm_index.second);
    }

    return ORF_coords;
}

NodeColourVector Graph::_index_graph (const std::vector<std::string>& stop_codons_for,
                                      const std::vector<std::string>& stop_codons_rev,
                                      const int& kmer,
                                      const size_t& nb_colours,
                                      const bool is_ref,
                                      const std::vector<std::string>& input_colours)
{
    auto node_colour_vector = index_graph(_KmerArray, _ccdbg, stop_codons_for, stop_codons_rev, kmer, nb_colours, is_ref, input_colours);

    // return node_colour vector
    return node_colour_vector;
}