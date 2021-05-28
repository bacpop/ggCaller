#ifndef BINDINGS_H
#define BINDINGS_H

#include "ggCaller_classes.h"
#include "traversal.h"
#include "call_ORFs.h"
#include "gene_overlap.h"

// ggCaller_bindings
GraphTuple py_index_graph_exists(const std::string& graphfile,
                                 const std::string& coloursfile,
                                 const std::vector<std::string>& stop_codons_for,
                                 const std::vector<std::string>& stop_codons_rev,
                                 size_t num_threads,
                                 const bool is_ref);

GraphTuple py_index_graph_build(const std::string& infile1,
                                const int kmer,
                                const std::vector<std::string>& stop_codons_for,
                                const std::vector<std::string>& stop_codons_rev,
                                size_t num_threads,
                                bool is_ref,
                                const bool write_graph,
                                const std::string& infile2);

std::pair<ORFOverlapMap, ORFVector> py_calculate_ORFs(const UnitigVector& graph_vector,
                                                     const size_t& colour_ID,
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

#endif //BINDINGS_H
