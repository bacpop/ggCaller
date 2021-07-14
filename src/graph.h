#ifndef GRAPH_H
#define GRAPH_H

#include "unitigDict.h"
#include "indexing.h"
#include "traversal.h"
#include "call_ORFs.h"
#include "match_string.h"
#include "gene_overlap.h"

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

    // add ORF information to graph
    void add_ORF_info (const size_t& colour_ID,
                       const std::vector<std::pair<size_t,size_t>>& ORF_IDs,
                       const ORFVector& ORF_vector)

    // get next ORFs along from current ORF


   // generate sequences from ORF node_lists
    std::string generate_sequence(const std::vector<int>& nodelist,
                                  const std::vector<indexPair>& node_coords,
                                  const size_t& overlap);

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
