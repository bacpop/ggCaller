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
    // constructors
    // default constructor
    Graph()
    {
        cout << "Graph default constructor" << endl;
    };

    // copy
    Graph(const Graph & rhs)
    {
          std::cout << "Graph copy constructor" << std::endl;
          _GraphVector = rhs._GraphVector;
    };
    // copy operator
    Graph & operator = (const Graph & rhs)
    {
          std::cout << "Graph copy operator" << std::endl;
          if(this != &rhs)
          {
            _GraphVector = rhs._GraphVector;
          }
          return *this;
    };

    // move
    Graph(const Graph && rhs) noexcept
    {
          std::cout << "Graph move constructor" << std::endl;
          _GraphVector = std::move(rhs._GraphVector);
    };
    // move operator
    Graph & operator = (const Graph && rhs) noexcept
    {
          std::cout << "Graph move operator" << std::endl;
          if(this != &rhs)
          {
            _GraphVector = std::move(rhs._GraphVector);
          }
          return *this;
    };

    // destructor
    ~Graph()
    {
        std::cout << "Graph destructor" << std::endl;
    };

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

    std::string generate_sequence(const std::vector<int>& nodelist,
                                  const std::vector<indexPair>& node_coords,
                                  const size_t& overlap);

    private:
    // index graph
    void _index_graph (const ColoredCDBG<>& ccdbg,
                     const std::vector<std::string>& stop_codons_for,
                     const std::vector<std::string>& stop_codons_rev,
                     const int& kmer,
                     const size_t& nb_colours);

    // traverse graph
    void _traverse_graph(const size_t& colour_ID,
                         const bool repeat,
                         const size_t max_path_length);

    // stored head nodes for traversal
    NodeColourVector _NodeColourVector;

    // stored head nodes that have been traversed
    NodeColourVector _NodeColourVectorTraversed;

    // stored paths for ORF calling per colour
    std::vector<AllPaths> _ColourGraphPaths;

    // stored unitigDict objects
    std::vector<unitigDict> _GraphVector;
};

#endif //GRAPH_H
