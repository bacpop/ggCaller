#ifndef GRAPH_H
#define GRAPH_H

#include <torch/script.h>
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
#include "ORF_scoring.h"
#include "gene_graph.h"
#include "progressbar.h"

class Graph {
    public:
//    // constructors for debugging
//    // default constructor
//    Graph()
//    {
//        cout << "Graph default constructor" << endl;
//    };
//
//    // copy
//    Graph(const Graph & rhs)
//    {
//          std::cout << "Graph copy constructor" << std::endl;
//          _ccdbg = rhs._ccdbg;
//          _stop_freq = rhs._stop_freq;
//          _KmerArray = rhs._KmerArray;
//          _NodeColourVector = rhs._NodeColourVector;
//          _StartFreq = rhs._StartFreq;
//          _RefSet = rhs._RefSet;
//    };
//    // copy operator
//    Graph & operator = (const Graph & rhs)
//    {
//          std::cout << "Graph copy operator" << std::endl;
//          if(this != &rhs)
//          {
//              _ccdbg = rhs._ccdbg;
//              _stop_freq = rhs._stop_freq;
//              _KmerArray = rhs._KmerArray;
//              _NodeColourVector = rhs._NodeColourVector;
//              _StartFreq = rhs._StartFreq;
//              _RefSet = rhs._RefSet;
//          }
//          return *this;
//    };
//
//    // move
//    Graph(const Graph && rhs) noexcept
//    {
//          std::cout << "Graph move constructor" << std::endl;
//          _ccdbg = std::move(rhs._ccdbg);
//          _stop_freq = std::move(rhs._stop_freq);
//          _KmerArray = std::move(rhs._KmerArray);
//          _NodeColourVector = std::move(rhs._NodeColourVector);
//          _StartFreq = std::move(rhs._StartFreq);
//          _RefSet = std::move(rhs._RefSet);
//    };
//    // move operator
//    Graph & operator = (const Graph && rhs) noexcept
//    {
//          std::cout << "Graph move operator" << std::endl;
//          if(this != &rhs)
//          {
//              _ccdbg = std::move(rhs._ccdbg);
//              _stop_freq = std::move(rhs._stop_freq);
//              _KmerArray = std::move(rhs._KmerArray);
//              _NodeColourVector = std::move(rhs._NodeColourVector);
//              _StartFreq = std::move(rhs._StartFreq);
//              _RefSet = std::move(rhs._RefSet);
//          }
//          return *this;
//    };
//
//    // destructor
//    ~Graph()
//    {
//        std::cout << "Graph destructor" << std::endl;
//    };

    // build new bifrost graph and index
    GraphTuple build (const std::string& infile1,
                      const int kmer,
                      const std::vector<std::string>& stop_codons_for,
                      const std::vector<std::string>& stop_codons_rev,
                      const std::vector<std::string>& start_codons_for,
                      const std::vector<std::string>& start_codons_rev,
                      size_t num_threads,
                      bool is_ref,
                      const bool write_graph,
                      const std::string& infile2,
                      const std::unordered_set<std::string>& ref_set);

    // read existing graph and index
    GraphTuple read (const std::string& graphfile,
                     const std::string& coloursfile,
                     const std::vector<std::string>& stop_codons_for,
                     const std::vector<std::string>& stop_codons_rev,
                     const std::vector<std::string>& start_codons_for,
                     const std::vector<std::string>& start_codons_rev,
                     size_t num_threads,
                     const bool is_ref,
                     const std::unordered_set<std::string>& ref_set);

    // get graph object from serialised file
    void in(const std::string& infile,
            const std::string& graphfile,
            const std::string& coloursfile,
            const size_t num_threads);

    // generate serialised file from graph object
    void out(const std::string& outfile);

    // find ORFs
std::pair<std::map<size_t, std::string>, std::map<size_t, std::string>> findGenes (const bool repeat,
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
                                                                                    const std::string& tmp_dir);


std::pair<RefindMap, bool> refind_gene(const size_t& colour_ID,
                                        const NodeSearchDict& node_search_dict,
                                        const size_t radius,
                                        const int kmer,
                                        const std::string& FM_fasta_file,
                                        const bool repeat,
                                        const std::unordered_set<int>& to_avoid,
                                        const string& ORF_file_path);

    // generate sequences from ORF node_lists
    std::string generate_sequence(const std::vector<int>& nodelist,
                                  const std::vector<indexPair>& node_coords,
                                  const size_t& overlap);

    std::tuple<std::vector<std::string>, int, std::vector<std::unordered_set<int>>> search_graph(const std::vector<std::string>& query_vec,
                                                                                       const double& id_cutoff,
                                                                                       size_t num_threads);

    std::vector<std::pair<ContigLoc, bool>> ORF_location(const std::vector<std::pair<std::vector<int>, std::vector<indexPair>>>& ORF_IDs,
                                                         const std::string& fasta_file,
                                                         const int overlap,
                                                         const bool write_idx,
                                                         size_t num_threads);
    size_t rb_hash(const std::string seq)
    {
        return hasher{}(seq);
    }

    size_t node_size(const int& node_id)
    {
        // get a reference to the unitig map object
        auto um_pair = get_um_data(_ccdbg, _KmerArray, node_id);
        auto& um = um_pair.first;

        return um.size;
    };

    // clear graph object
    void clear() {_ccdbg.clear(); _KmerArray.clear();};

    private:
    // index graph
    void _index_graph (const std::vector<std::string>& stop_codons_for,
                       const std::vector<std::string>& stop_codons_rev,
                       const std::vector<std::string>& start_codons_for,
                       const std::vector<std::string>& start_codons_rev,
                       const int& kmer,
                       const size_t& nb_colours,
                       const std::vector<std::string>& input_colours);

    
    // stored bifrost DBG
    ColoredCDBG<MyUnitigMap> _ccdbg;

    // stop codon frequency
    float _stop_freq = 0;

    // mapping of head kmers to nodes
    std::vector<Kmer> _KmerArray;

    // nodes to colours
    NodeColourVector _NodeColourVector;

    // mapping of start site sequences to frequency
    robin_hood::unordered_map<std::string, size_t> _StartFreq;

    // bitset to determine if colours are refs or reads
    boost::dynamic_bitset<> _RefSet;
};

ORFClusterMap read_cluster_file(const std::string& cluster_file);

ORFNodeMap read_ORF_file(const std::string& ORF_file);

void save_ORF_file(const std::string& ORF_file_path,
                   const ORFNodeMap& ORF_map);

std::unordered_map<size_t, std::unordered_set<size_t>> read_edge_file(const std::string& cluster_file);

void clear_graph(Graph& g);

//void use_count(std::shared_ptr<Graph> sp);

#endif //GRAPH_H
