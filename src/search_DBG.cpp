//
// Created by sth19 on 07/12/2021.
//

#include "search_DBG.h"

MappingCoords query_DBG(const ColoredCDBG<>& ccdbg,
                        const std::string& query,
                        const int& kmer)
{
    // tuple of head kmer string, strand and coordinates
    MappingCoords mapping_coords;

    const size_t num_kmers = query.length() - kmer + 1;

    // counter for how many k-mers into a unitig search is
    size_t kmer_counter = 0;

    // count for determining where in query sequence search is
    size_t kmer_index = 1;

    const char *query_str = query.c_str();

    for (KmerIterator it_km(query_str), it_km_end; it_km != it_km_end; ++it_km)
    {
        auto um = ccdbg.find(it_km->first);

        // if found, add to FM-index string
        if (!um.isEmpty)
        {
            const std::string head_kmer = um.getUnitigHead().toString();

            if (mapping_coords.empty())
            {
                if (um.strand)
                {
                    mapping_coords.push_back({head_kmer, um.strand, {um.dist, 0}});
                } else
                {
                    mapping_coords.push_back({head_kmer, um.strand, {0, um.dist + (kmer - 1)}});
                }
            } else
            {
                if (std::get<0>(mapping_coords.back()) != head_kmer)
                {
                    // if previous node is forward strand
                    if (std::get<1>(mapping_coords.back()))
                    {
                        std::get<2>(mapping_coords.back()).second = std::get<2>(mapping_coords.back()).first + kmer_counter + kmer - 2;
                    } else
                    {
                        std::get<2>(mapping_coords.back()).first = std::get<2>(mapping_coords.back()).second - kmer_counter - kmer + 2;
                    }

                    if (um.strand)
                    {
                        mapping_coords.push_back({head_kmer, um.strand, {um.dist, 0}});
                    } else
                    {
                        mapping_coords.push_back({head_kmer, um.strand, {0, um.dist + (kmer - 1)}});
                    }

                    kmer_counter = 0;
                }
            }

            // if at end of query, need to check final position
            // if at last kmer, need to determine position within unitig
            if (kmer_index == num_kmers)
            {
                // if previous node is forward strand, need to add on kmer length
                if (std::get<1>(mapping_coords.back()))
                {
                    std::get<2>(mapping_coords.back()).second = std::get<2>(mapping_coords.back()).first + kmer_counter + kmer - 1;
                } else
                {
                    std::get<2>(mapping_coords.back()).first = um.dist;
                }
                break;
            }

            kmer_index++;
            kmer_counter++;
        } else
        {
            break;
        }
    }

    // check if query traversed fully
    if (kmer_index != num_kmers)
    {
        return {};
    }
    else
    {
        return mapping_coords;
    }
}
