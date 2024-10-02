//
// Created by sth19 on 28/05/2021.
//

#include "graph.h"

GraphTuple Graph::build (const std::string& infile1,
                        const int kmer,
                        const std::vector<std::string>& stop_codons_for,
                        const std::vector<std::string>& stop_codons_rev,
                        const std::vector<std::string>& start_codons_for,
                        const std::vector<std::string>& start_codons_rev,
                        size_t num_threads,
                        bool is_ref,
                        const bool write_graph,
                        const std::string& infile2,
                        const std::unordered_set<std::string>& ref_set) {
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

    // store is_ref information in bitvector
    _RefSet.resize(nb_colours);
    // assume all colours are references
    if (is_ref && ref_set.empty())
    {
        _RefSet.set();
    } else
    {
        for (int i = 0; i < input_colours.size(); i++)
        {
            if (ref_set.find(input_colours[i]) != ref_set.end())
            {
                _RefSet[i] = 1;
            }
        }
    }

    // generate codon index for graph
    cout << "Generating graph stop codon index..." << endl;
    _index_graph(stop_codons_for, stop_codons_rev, start_codons_for, start_codons_rev, kmer, nb_colours, input_colours);

    // create vector bool for reference sequences
    std::vector<bool> ref_list(nb_colours, false);
    for (int i = 0; i < _RefSet.size(); i++)
    {
        if ((bool)_RefSet[i])
        {
            ref_list[i] = true;
        }
    }

    // make tuple containing all information needed in python back-end
    GraphTuple graph_tuple = std::make_tuple(input_colours, nb_colours, overlap, ref_list);

    return graph_tuple;
}

// read existing graph and index
GraphTuple Graph::read (const std::string& graphfile,
                    const std::string& coloursfile,
                    const std::vector<std::string>& stop_codons_for,
                    const std::vector<std::string>& stop_codons_rev,
                    const std::vector<std::string>& start_codons_for,
                    const std::vector<std::string>& start_codons_rev,
                    size_t num_threads,
                    const bool is_ref,
                    const std::unordered_set<std::string>& ref_set) {

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

    // store is_ref information in bitvector
    _RefSet.resize(nb_colours);
    // assume all colours are references
    if (is_ref && ref_set.empty())
    {
        _RefSet.set();
    } else
    {
        for (int i = 0; i < input_colours.size(); i++)
        {
            if (ref_set.find(input_colours[i]) != ref_set.end())
            {
                _RefSet[i] = 1;
            }
        }
    }

    // generate codon index for graph
    cout << "Generating graph stop codon index..." << endl;
    _index_graph(stop_codons_for, stop_codons_rev, start_codons_for, start_codons_rev, kmer, nb_colours, input_colours);

    // create vector bool for reference sequences
    std::vector<bool> ref_list(nb_colours, false);
    for (int i = 0; i < _RefSet.size(); i++)
    {
        if ((bool)_RefSet[i])
        {
            ref_list[i] = true;
        }
    }

    // make tuple containing all information needed in python back-end
    GraphTuple graph_tuple = std::make_tuple(input_colours, nb_colours, overlap, ref_list);

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


std::pair<std::map<size_t, std::string>, std::map<size_t, std::string>> Graph::findGenes (const bool repeat,
                                                                                            const size_t overlap,
                                                                                            const size_t max_path_length,
                                                                                            bool no_filter,
                                                                                            const std::vector<std::string>& stop_codons_for,
                                                                                            const std::vector<std::string>& start_codons_for,
                                                                                            const size_t min_ORF_length,
                                                                                            const size_t max_overlap,
                                                                                            const std::vector<std::string>& input_colours,
                                                                                            const std::string& ORF_model_file,
                                                                                            const std::string& TIS_model_file,
                                                                                            const float& minimum_ORF_score,
                                                                                            const float& minimum_path_score,
                                                                                            const size_t max_ORF_path_length,
                                                                                            const bool clustering,
                                                                                            const double& id_cutoff,
                                                                                            const double& len_diff_cutoff,
                                                                                            size_t num_threads,
                                                                                            const std::string& cluster_file,
                                                                                            const float& score_tolerance,
                                                                                            const std::string& tmp_dir)
{    
    // initilise map to hold file paths
    std::map<size_t, std::string> ORF_file_paths;
    std::map<size_t, std::string> Edge_file_paths;

    // initilise all colour keys
    for (int colour_ID = 0; colour_ID < input_colours.size(); colour_ID++)
    {
        ORF_file_paths[colour_ID] = "";
        Edge_file_paths[colour_ID] = "";
    }

    // Set number of threads
    if (num_threads < 1)
    {
        num_threads = 1;
    }

    // set OMP number of threads
    omp_set_num_threads(num_threads);

    // set global error for Balrog model loading
    bool error = false;

    {
        // set up progress bar
        progressbar bar(input_colours.size());
        bar.set_todo_char(" ");
        bar.set_done_char("█");
        bar.set_opening_bracket_char("|");
        bar.set_closing_bracket_char("|");

        // load Balrog model for TIS
        torch::jit::script::Module TIS_model;

        if (!no_filter)
        {
            try {
                // Deserialize the ScriptModule from a file using torch::jit::load().
                TIS_model = torch::jit::load(TIS_model_file);
            }
            catch (const c10::Error& e) {
                std::cerr << "error loading the TIS model\n";
                error = true;
            }
        }

        // initialise maps to store ORF TIS scores across threads
        tbb::concurrent_unordered_map<size_t, float> all_TIS_scores;

        // create map to hold number of times start codons chosen
        tbb::concurrent_unordered_map<size_t, tbb::concurrent_unordered_set<int>> start_chosen;

        cout << "Traversing graph to identify ORFs..." << endl;

        const int aa_kmer = std::round((float) (overlap + 1) / (float) 6) - 1;

        #pragma omp parallel for
        for (int colour_ID = 0; colour_ID < input_colours.size(); colour_ID++)
        {
            // get whether colour is reference or not
            bool is_ref = ((bool)_RefSet[colour_ID]) ? true : false;

            const auto& FM_fasta_file = input_colours.at(colour_ID);

            // if no FM_fasta_file specified, cannot generate FM Index
            if (FM_fasta_file == "NA")
            {
                is_ref = false;
            }

            // initialise ORF_vector
            ORFNodeMap ORF_map;

            // traverse graph, set scope for all_paths and fm_idx
            {
                const auto& node_ids = _NodeColourVector.at(colour_ID);

                // generate FM_index if is_ref
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

                // recursive traversal and ORF calling
                if (error)
                {
                    no_filter = true;
                }

                // convert this to map to make removal easier
                ORF_map = std::move(traverse_graph(_ccdbg, _KmerArray, _stop_freq, colour_ID, node_ids, repeat, max_path_length,
                                                   overlap, is_ref, _RefSet, fm_idx, stop_codons_for, start_codons_for, min_ORF_length,
                                                   TIS_model, minimum_ORF_score, no_filter, all_TIS_scores, _StartFreq, score_tolerance,
                                                   start_chosen, aa_kmer));

            }

            // write ORF_map file
            {
                std::string ORF_file_path = tmp_dir + "/colour_" + std::to_string(colour_ID) + "_ORFs.tmp";
                std::ofstream ofs(ORF_file_path);
                boost::archive::text_oarchive oa(ofs);
                // write class instance to archive
                oa << ORF_map;

                ORF_file_paths[colour_ID] = ORF_file_path;
            }

            // update progress bar
            #pragma omp critical
            {
                bar.update();
            }
        }
    }

    //clear objects no longer used
    _NodeColourVector.clear();
    _StartFreq.clear();

    // add new line to account for progress bar
    cout << endl;

    // generate clusters if required
    if (clustering || !no_filter)
    {
        cout << "Generating clusters of high-scoring ORFs..." << endl;
        ORFClusterMap cluster_map;

        // load Balrog gene model
        torch::jit::script::Module ORF_model;

        if (!no_filter)
        {
            try {
                // Deserialize the ScriptModule from a file using torch::jit::load().
                ORF_model = torch::jit::load(ORF_model_file);
            }
            catch (const c10::Error& e) {
                std::cerr << "error loading the ORF model\n";
                error = true;
            }
        }

        // initialise maps to store ORF scores across threads
        tbb::concurrent_unordered_map<size_t, float> all_ORF_scores;

        // scope for clustering variables
        {
            // group ORFs together based on single shared k-mer
            auto ORF_group_pair = group_ORFs(ORF_file_paths, _ccdbg, _KmerArray, overlap);

            // generate clusters for ORFs based on identity
            cluster_map = produce_clusters(ORF_file_paths, _ccdbg, _KmerArray, overlap,
                                           ORF_group_pair, id_cutoff, len_diff_cutoff);
        }

        if (!no_filter)
        {
            cout << "Scoring Centroids..." << endl;

            // keep track of clusters with low scoring centroids
            //tbb::concurrent_unordered_set<size_t> to_remove_cluster;

            // set up progress bar
            progressbar bar(ORF_file_paths.size());
            bar.set_todo_char(" ");
            bar.set_done_char("█");
            bar.set_opening_bracket_char("|");
            bar.set_closing_bracket_char("|");

            // score centroid and apply score to all ORFs in group
            std::set<std::pair<size_t, size_t>> to_remove;
            robin_hood::unordered_map<size_t, robin_hood::unordered_map<size_t, float>> ORFToScoreMap;
            
            //#pragma omp parallel for
            for (const auto& colour : ORF_file_paths)
            {
                ORFNodeMap ORF_map;
                // read in ORF_map file
                {
                    std::ifstream ifs(colour.second);
                    boost::archive::text_iarchive ia(ifs);
                    ia >> ORF_map;
                }
                
                for (auto& ORF_entry : ORF_map)
                {
                    // check if centroid, if so then score otherwise ignore
                    std::string ORF_ID_str = std::to_string(colour.first) + "_" + std::to_string(ORF_entry.first);

                    const auto& centroid_found = cluster_map.find(ORF_ID_str);
                    
                    if (centroid_found != cluster_map.end())
                    {
                        auto& centroid_info = ORF_entry.second;
                        std::string centroid_seq = translate(generate_sequence_nm(std::get<0>(centroid_info), std::get<1>(centroid_info), overlap, _ccdbg, _KmerArray));
                        
                        const auto gene_prob = score_gene(std::get<4>(centroid_info), centroid_seq, std::get<2>(centroid_info), ORF_model, all_ORF_scores);
                        
                        // if centroid score is below the min-orf score remove from cluster map
                        if (std::get<4>(centroid_info) < minimum_ORF_score)
                        {
                            // remove all ORFs (first is centroid)
                            for (const auto& homolog_ID : centroid_found->second) 
                            {
                                to_remove.insert(homolog_ID);
                            }

                            // remove from cluster map
                            cluster_map.erase(ORF_ID_str);

                        } else {
                            // add score for calculation
                            for (const auto& homolog_ID : centroid_found->second) 
                            {
                                ORFToScoreMap[homolog_ID.first][homolog_ID.second] = gene_prob;
                            }
                        }
                    }
                }
                // update progress bar
                #pragma omp critical
                {
                    bar.update();
                }
            }
            // new line for progress bar
            cout << endl;

            // now need to iterate over genes again, identify the centroid and then assign the score
            // need to hold scores for each ORF so that it can be assigned to the new ORFs from gene_prob
            
            cout << "Scoring remaining ORFs..." << endl;
            progressbar bar2(ORF_file_paths.size());
            bar2.set_todo_char(" ");
            bar2.set_done_char("█");
            bar2.set_opening_bracket_char("|");
            bar2.set_closing_bracket_char("|");
            
            for (const auto& colour : ORF_file_paths)
            {
                // set up progress bar

                ORFNodeMap ORF_map;
                // read in ORF_map file
                {
                    std::ifstream ifs(colour.second);
                    boost::archive::text_iarchive ia(ifs);
                    ia >> ORF_map;
                }

                // iterate over ORFToScoreMap for the given colour
                for (const auto& ORF_ID : ORFToScoreMap[colour.first])
                {
                    // remove if score too low and move onto next gene
                    if (to_remove.find({colour.first, ORF_ID.first}) != to_remove.end())
                    {
                        ORF_map.erase(ORF_ID.first);
                        continue;
                    }

                    std::string ORF_ID_str = std::to_string(colour.first) + "_" + std::to_string(ORF_ID.first);
                    
                    auto centroid_found = cluster_map.find(ORF_ID_str);
                    
                    // this time ignore centroids as already scored
                    if (centroid_found == cluster_map.end())
                    {
                        const auto& gene_prob = ORF_ID.second;

                        auto& ORF_info = ORF_map[ORF_ID.first];
                        
                        // get ORF sequence
                        const auto ORF_seq = generate_sequence_nm(std::get<0>(ORF_info), std::get<1>(ORF_info), overlap, _ccdbg, _KmerArray);

                        // score ORF in place
                        score_cluster(std::get<4>(ORF_info), gene_prob, ORF_seq, std::get<2>(ORF_info));

                        // remove from orf map if score too low
                        if (std::get<4>(ORF_info) < minimum_ORF_score)
                        {
                            ORF_map.erase(ORF_ID.first);
                        }
                    }
                }

                // overwrite ORF_map file
                {
                    std::ofstream ofs(colour.second);
                    boost::archive::text_oarchive oa(ofs);
                    // write class instance to archive
                    oa << ORF_map;
                }

                // update progress bar
                #pragma omp critical
                {
                    bar2.update();
                }
            }
        }

        // add line for progress bar
        cout << endl;

        // write cluster file
        {
            std::ofstream ofs(cluster_file);
            boost::archive::text_oarchive oa(ofs);
            // write class instance to archive

            oa << cluster_map;
        }
    }

    {
        // set up progress bar
        progressbar bar(input_colours.size());
        bar.set_todo_char(" ");
        bar.set_done_char("█");
        bar.set_opening_bracket_char("|");
        bar.set_closing_bracket_char("|");

        cout << "Identifying high-scoring ORFs..." << endl;
        // after clustering, determing highest scoring gene set
        #pragma omp parallel for
        for (int colour_ID = 0; colour_ID < ORF_file_paths.size(); colour_ID++)
        {
            ORFNodeMap ORF_map;
            // read in ORF_map file
            {
                std::ifstream ifs(ORF_file_paths.at(colour_ID));
                boost::archive::text_iarchive ia(ifs);
                ia >> ORF_map;
            }

            // get whether colour is reference or not
            bool is_ref = ((bool)_RefSet[colour_ID]) ? true : false;

            const auto& FM_fasta_file = input_colours.at(colour_ID);

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

            // initialise values for gene information
            std::unordered_map<size_t, std::unordered_set<size_t>> gene_edges;
            ORFNodeMap gene_map;
            std::vector<std::vector<size_t>> gene_paths;

            // if no filtering required, do not calculate overlaps, score genes or get gene_paths
            if (!no_filter)
            {
                ORFOverlapMap ORF_overlap_map;
                //        cout << "Determining overlaps: " << to_string(colour_ID) << endl;
                {
                    ORF_overlap_map = std::move(calculate_overlaps(_ccdbg, _KmerArray, ORF_map, overlap, max_overlap, is_ref, fm_idx));
                }

                if (!error)
                {
                    gene_paths = call_true_genes(ORF_map, ORF_overlap_map, minimum_path_score, _KmerArray);

                    // get high scoring genes
                    for (const auto& path : gene_paths)
                    {
                        for (const auto& ORF_ID : path)
                        {
                            if (gene_map.find(ORF_ID) == gene_map.end())
                            {
                                // simplify ORF_info
                                simplify_ORFNodeVector(ORF_map[ORF_ID], overlap);
                                gene_map[ORF_ID] = std::move(ORF_map[ORF_ID]);
                            }
                        }
                    }
                } else
                {
                    // return unfiltered genes
                    for (auto& entry : ORF_map)
                    {
                        // simplify ORF_info
                        simplify_ORFNodeVector(entry.second, overlap);
                        gene_map[entry.first] = std::move(entry.second);
                        gene_paths.push_back({entry.first});
                    }
                }
            } else
            {
                // return unfiltered genes
                for (auto& entry : ORF_map)
                {
                    // simplify ORF_info
                    simplify_ORFNodeVector(entry.second, overlap);
                    gene_map[entry.first] = std::move(entry.second);
                    gene_paths.push_back({entry.first});
                }
            }

            // deallocate ORF_map
            ORF_map.clear();

            // connect ORFs
            {
                std::set<std::pair<size_t, size_t>> connected_ORFs;

                // determine target_ORFs to connect and redundant edges
                std::set<std::pair<size_t, size_t>> redundant_edges;
                robin_hood::unordered_set<size_t> target_ORFs;
                for (const auto& path : gene_paths)
                {
                    const auto& first = path.at(0);
                    const auto& second = path.back();
                    target_ORFs.insert(first);
                    target_ORFs.insert(second);

                    // add redundant edges by ordering first and second ORF
                    if (first <= second)
                    {
                        redundant_edges.insert({first, second});
                    } else
                    {
                        redundant_edges.insert({second, first});
                    }
                }

                // add ORF info for colour to graph
                auto node_to_ORFs = add_ORF_info(_ccdbg, _KmerArray, target_ORFs, gene_map, overlap);

                // initialise prev_node_set to avoid same ORFs being traversed from again
                robin_hood::unordered_set<size_t> downstream_ORF_set;
                robin_hood::unordered_set<size_t> upstream_ORF_set;

                // conduct DBG traversal for upstream...
                auto new_connections = pair_ORF_nodes(_ccdbg, _KmerArray, node_to_ORFs, colour_ID, target_ORFs, gene_map, max_ORF_path_length, repeat, -1, downstream_ORF_set, upstream_ORF_set, overlap, is_ref, fm_idx);
                connected_ORFs.insert(std::make_move_iterator(new_connections.begin()), std::make_move_iterator(new_connections.end()));

                // ... and downstream
                new_connections = pair_ORF_nodes(_ccdbg, _KmerArray, node_to_ORFs, colour_ID, target_ORFs, gene_map, max_ORF_path_length, repeat, 1, downstream_ORF_set, upstream_ORF_set, overlap, is_ref, fm_idx);
                connected_ORFs.insert(std::make_move_iterator(new_connections.begin()), std::make_move_iterator(new_connections.end()));

                // check edges found in connected_ORFs against redundant edges
                for (const auto& edge : connected_ORFs)
                {
                    if (redundant_edges.find(edge) == redundant_edges.end())
                    {
                        std::vector<size_t> edge_vec = {edge.first, edge.second};
                        gene_paths.push_back(edge_vec);
                    }
                }
            }

            // create map to connect ORFs for panaroo and to store high scoring ORFs
            for (const auto& entry : gene_paths)
            {
                size_t last_index = entry.size() - 1;
                for (int i = 0; i < entry.size(); i++)
                {
                    const size_t& ORF = entry.at(i);
                    if (i != last_index)
                    {
                        gene_edges[ORF].insert(entry.at(i + 1));
                    }
                }
            }

            // overwrite ORF_map file with gene_map
            {
                std::ofstream ofs(ORF_file_paths.at(colour_ID));
                boost::archive::text_oarchive oa(ofs);
                // write class instance to archive
                oa << gene_map;
            }

            {
                std::string file_path = tmp_dir + "/colour_" + std::to_string(colour_ID) + "_edges.tmp";
                std::ofstream ofs(file_path);
                boost::archive::text_oarchive oa(ofs);
                // write class instance to archive
                oa << gene_edges;

                Edge_file_paths[colour_ID] = file_path;
            }

            // update colour_ORF_map and colour_edge_map
            //colour_ORF_map[colour_ID] = std::move(gene_map);
            //colour_edge_map[colour_ID] = std::move(gene_edges);

            // update progress bar
            #pragma omp critical
            {
                bar.update();
            }
        }
    }

    // add line for progress bar
    cout << endl;

    return {ORF_file_paths, Edge_file_paths};
}


std::pair<RefindMap, bool> Graph::refind_gene(const size_t& colour_ID,
                                             const NodeSearchDict& node_search_dict,
                                             const size_t radius,
                                             const int kmer,
                                             const std::string& FM_fasta_file,
                                             const bool repeat,
                                             const std::unordered_set<int>& to_avoid)
{
    // get whether colour is reference or not
    bool is_ref = ((bool)_RefSet[colour_ID]) ? true : false;

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

    return {refind_in_nodes(_ccdbg, _KmerArray, colour_ID, node_search_dict, radius, is_ref,
                            kmer, fm_idx, repeat, to_avoid), is_ref};
}

std::string Graph::generate_sequence(const std::vector<int>& nodelist,
                                     const std::vector<indexPair>& node_coords,
                                     const size_t& overlap)
{
    return generate_sequence_nm(nodelist, node_coords, overlap, _ccdbg, _KmerArray);
}

std::tuple<std::vector<std::string>, int, std::vector<std::unordered_set<int>>> Graph::search_graph(const std::vector<std::string>& query_vec,
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

    // get input colours
    std::vector<std::string> input_colours = _ccdbg.getColorNames();

    //set local variables
    const int kmer = _ccdbg.getK();

    std::vector<std::unordered_set<int>> query_nodes(query_vec.size());

    // go through query, determine head-kmers of each node and map to _GraphVector
    #pragma omp parallel for
    for (int i = 0; i < query_vec.size(); i++)
    {
        query_nodes[i] = std::move(query_DBG(_ccdbg, query_vec.at(i), kmer, id_cutoff));
    }

    return {input_colours, kmer, query_nodes};
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

void Graph::_index_graph (const std::vector<std::string>& stop_codons_for,
                          const std::vector<std::string>& stop_codons_rev,
                          const std::vector<std::string>& start_codons_for,
                          const std::vector<std::string>& start_codons_rev,
                          const int& kmer,
                          const size_t& nb_colours,
                          const std::vector<std::string>& input_colours)
{
    float stop_codon_freq = 0;
    _NodeColourVector = std::move(index_graph(_KmerArray, _ccdbg, stop_codon_freq, stop_codons_for, stop_codons_rev, start_codons_for, start_codons_rev, kmer, nb_colours, input_colours, _RefSet, _StartFreq));
    _stop_freq= stop_codon_freq;
}

ORFClusterMap read_cluster_file(const std::string& cluster_file)
{
    ORFClusterMap cluster_map;

    std::ifstream ifs(cluster_file);
    boost::archive::text_iarchive ia(ifs);
    ia >> cluster_map;

    return cluster_map;
}

void clear_graph(Graph& g)
{
    g.~Graph();
}

//void use_count(std::shared_ptr<Graph> sp)
//{
//    std::cout << "use_count() == " << sp.use_count() << endl;
//}