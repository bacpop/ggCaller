#ifndef GRAPH_H
#define GRAPH_H

#include "unitigDict.h"
#include "indexing.h"
#include "traversal.h"
#include "call_ORFs.h"
#include "match_string.h"
#include "gene_overlap.h"
#include "ORF_connection.h"
#include "ORF_clustering.h"
#include "gene_refinding.h"
#include "search_DBG.h"

class Graph {
    public:
    // build new bifrost graph and index
    GraphTuple build (const std::string& infile1,
                    const int kmer,
                    const std::vector<std::string>& stop_codons_for,
                    const std::vector<std::string>& stop_codons_rev,
                    size_t num_threads,
                    bool is_ref,
                    const bool write_graph,
                    const std::string& infile2);

    // read existing graph and index
    GraphTuple read (const std::string& graphfile,
                     const std::string& coloursfile,
                     const std::vector<std::string>& stop_codons_for,
                     const std::vector<std::string>& stop_codons_rev,
                     size_t num_threads,
                     const bool is_ref);

    // get graph object from serialised file
    void in(const std::string& infile,
            const std::string& graphfile,
            const std::string& coloursfile,
            const size_t num_threads);

    // generate serialised file from graph object
    void out(const std::string& outfile);

    // find ORFs
    std::pair<ORFOverlapMap, ORFVector> findORFs (const size_t colour_ID,
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
                                                  const std::string& FM_fasta_file);

    // search orientated paths and DBG to connect ORFs
    std::vector<std::pair<size_t, size_t>> connect_ORFs(const size_t colour_ID,
                                                        const ORFVector& ORF_vector,
                                                        const std::vector<size_t>& target_ORFs,
                                                        const size_t max_ORF_path_length,
                                                        const bool is_ref,
                                                        const int overlap);

    std::pair<ORFMatrixVector, ORFClusterMap> generate_clusters(const ColourORFMap& colour_ORF_map,
                                                                const size_t& overlap,
                                                                const double& id_cutoff,
                                                                const double& len_diff_cutoff);

    RefindMap refind_gene(const size_t& colour_ID,
                          const std::unordered_map<int, std::unordered_map<std::string, ORFNodeVector>>& node_search_dict,
                          const size_t radius,
                          bool is_ref,
                          const int kmer,
                          const std::string& FM_fasta_file,
                          const bool repeat);

    // generate sequences from ORF node_lists
    std::string generate_sequence(const std::vector<int>& nodelist,
                                  const std::vector<indexPair>& node_coords,
                                  const size_t& overlap);

    std::tuple<std::vector<std::string>, int, std::vector<MappingCoords>> search_graph(const std::string& graphfile,
                                                                                      const std::string& coloursfile,
                                                                                      const std::vector<std::string>& query_vec,
                                                                                      const double& id_cutoff,
                                                                                      size_t num_threads);

//    size_t node_size(const int& node_id)
//    {
//        return _GraphVector.at(abs(node_id) - 1).size().first;
//    };

    // clear graph object
    void clear() {_ccdbg.clear(); _KmerArray.clear();};

    private:
    // index graph
    NodeColourVector _index_graph (const std::vector<std::string>& stop_codons_for,
                                     const std::vector<std::string>& stop_codons_rev,
                                     const int& kmer,
                                     const size_t& nb_colours,
                                     const bool is_ref,
                                     const std::vector<std::string>& input_colours);


//    // stored unitigDict objects
//    GraphVector _GraphVector;
    
    // stored bifrost DBG
    ColoredCDBG<MyUnitigMap> _ccdbg;

    // mapping of head kmers to nodes
    std::vector<Kmer> _KmerArray;
};

#endif //GRAPH_H
