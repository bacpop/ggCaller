#ifndef GRAPH_H
#define GRAPH_H

#include "ggCaller_classes.h"
#include "indexing.h"

class Graph {
    public:
    // build new bifrost graph and index
    GraphTuple build (const std::string& infile1,
                    const int kmer,
                    const std::vector<std::string>& stop_codons_for,
                    const std::vector<std::string>& stop_codons_rev,
                    const size_t num_threads,
                    bool is_ref,
                    const bool write_graph,
                    const std::string& infile2);

    // build from existing graph and index
    GraphTuple build (const std::string& graphfile,
                     const std::string& coloursfile,
                     const std::vector<std::string>& stop_codons_for,
                     const std::vector<std::string>& stop_codons_rev,
                     const size_t num_threads,
                     const bool is_ref);

    // add unitigDict objects to graph


    private:
    // index graph
    GraphTuple _index_graph (const ColoredCDBG<>& ccdbg,
                             const std::vector<std::string>& stop_codons_for,
                             const std::vector<std::string>& stop_codons_rev,
                             const int kmer,
                             const size_t nb_colours);

    // stored unitigDict objects
    std::vector<unitigDict> _UnitigVector;
};

#endif //GRAPH_H
