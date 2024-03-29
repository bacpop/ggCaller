//
// Created by sth19 on 07/12/2021.
//

#include "search_DBG.h"

std::unordered_set<int> query_DBG(const ColoredCDBG<MyUnitigMap>& ccdbg,
                        const std::string& query,
                        const int& kmer,
                        const double& id_cutoff)
{
    // tuple of head kmer string, strand and coordinates
    std::unordered_set<int> node_set;

    const size_t num_kmers = query.length() - kmer + 1;

    // convert query to string for search in graph
    const char *query_str = query.c_str();

    // count number of missed k-mer
    size_t missed_kmers = 0;

    for (KmerIterator it_km(query_str), it_km_end; it_km != it_km_end; ++it_km)
    {
        auto um = ccdbg.find(it_km->first);

        // if found, add to mapping coords
        if (!um.isEmpty)
        {
            auto da = um.getData();
            const MyUnitigMap* um_data = da->getData(um);

            const int node_id = um_data->get_id();

            node_set.insert(node_id);

        } else
        {
            // if not matching add to missed_kmers
            missed_kmers++;

            // check if missing this kmer means total kmer count falls below threshold
            if ((num_kmers - missed_kmers) < (num_kmers * id_cutoff))
            {
                break;
            }
        }
    }

    // check if query traversed fully
    if ((num_kmers - missed_kmers) < (num_kmers * id_cutoff))
    {
        std::unordered_set<int> empty_set;
        return empty_set;
    }
    else
    {
        return node_set;
    }
}
