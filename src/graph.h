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

    // find ORFs
    std::pair<ORFOverlapMap, ORFVector> findORFs (const size_t& colour_ID,
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
                                                  const std::string& FM_fasta_file);

    // search orientated paths and DBG to connect ORFs
    std::vector<std::pair<size_t, size_t>> connect_ORFs(const size_t& colour_ID,
                                                        const ORFVector& ORF_vector,
                                                        const std::vector<size_t>& target_ORFs,
                                                        const size_t& max_ORF_path_length);

    std::pair<ORFMatrixVector, ORFClusterMap> generate_clusters(const std::unordered_map<size_t, ORFNodeMap>& colour_ORF_map,
                                                                const size_t& overlap,
                                                                const double& id_cutoff,
                                                                const double& len_diff_cutoff);

    std::pair<std::string, std::string> refind_gene(const size_t& colour_ID,
                                                     const ORFNodeVector& ORF_info,
                                                     const size_t& radius,
                                                     const bool is_ref,
                                                     const bool write_idx,
                                                     const int kmer,
                                                     const std::string& FM_fasta_file,
                                                     const bool repeat);

    // generate sequences from ORF node_lists
    std::string generate_sequence(const std::vector<int>& nodelist,
                                  const std::vector<indexPair>& node_coords,
                                  const size_t& overlap);

    // clear graph object
    void clear() {_GraphVector.clear();};

    private:
    // index graph
    NodeColourVector _index_graph (const ColoredCDBG<>& ccdbg,
                             const std::vector<std::string>& stop_codons_for,
                             const std::vector<std::string>& stop_codons_rev,
                             const int& kmer,
                             const size_t& nb_colours);

    // stored unitigDict objects
    std::vector<unitigDict> _GraphVector;
};

#endif //GRAPH_H
