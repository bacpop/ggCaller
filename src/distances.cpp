#include "distances.h"

/*/
MIT License

Copyright (c) 2018 Gerry Tonkin-Hill

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
/*/

template <typename T>
std::vector<T> combine_vectors(const std::vector<std::vector<T>> &vec,
                               const size_t len) {
    std::vector<T> all(len);
    auto all_it = all.begin();
    for (size_t i = 0; i < vec.size(); ++i) {
        std::copy(vec[i].cbegin(), vec[i].cend(), all_it);
        all_it += vec[i].size();
    }
    return all;
}

std::vector<double> get_distances_align(const std::vector<std::string>& matrix_in,
                                        const int n_threads)
{
    omp_set_num_threads(n_threads);

    size_t n_seqs = matrix_in.size();
    size_t seq_length = matrix_in.at(0).size();

    // Shared variables for openmp loop
    std::vector<std::vector<double>> distances(n_seqs);
    uint64_t len = 0;

    // initialise bitmaps
    std::vector<boost::dynamic_bitset<>> A_snps;
    std::vector<boost::dynamic_bitset<>> C_snps;
    std::vector<boost::dynamic_bitset<>> G_snps;
    std::vector<boost::dynamic_bitset<>> T_snps;

    for (const auto& seq : matrix_in)
    {
        boost::dynamic_bitset<> As(seq_length);
        boost::dynamic_bitset<> Cs(seq_length);
        boost::dynamic_bitset<> Gs(seq_length);
        boost::dynamic_bitset<> Ts(seq_length);

        for (size_t j = 0; j < seq_length; j++) {

            char nuc = std::toupper(seq[j]);

            switch (nuc) {
                case 'A':
                    As[j] = 1;
                    break;
                case 'C':
                    Cs[j] = 1;
                    break;
                case 'G':
                    Gs[j] = 1;
                    break;
                case 'T':
                    Ts[j] = 1;
                    break;

                    // M = A or C
                case 'M':
                    As[j] = 1;
                    Cs[j] = 1;
                    break;

                    // R = A or G
                case 'R':
                    As[j] = 1;
                    Gs[j] = 1;
                    break;

                    // W = A or T
                case 'W':
                    As[j] = 1;
                    Ts[j] = 1;
                    break;

                    // S = C or G
                case 'S':
                    Cs[j] = 1;
                    Gs[j] = 1;
                    break;

                    // Y = C or T
                case 'Y':
                    Cs[j] = 1;
                    Ts[j] = 1;
                    break;

                    // K = G or T
                case 'K':
                    Gs[j] = 1;
                    Ts[j] = 1;
                    break;

                    // V = A,C or G
                case 'V':
                    As[j] = 1;
                    Cs[j] = 1;
                    Gs[j] = 1;
                    break;

                    // H = A,C or T
                case 'H':
                    As[j] = 1;
                    Cs[j] = 1;
                    Ts[j] = 1;
                    break;

                    // D = A,G or T
                case 'D':
                    As[j] = 1;
                    Gs[j] = 1;
                    Ts[j] = 1;
                    break;

                    // B = C,G or T
                case 'B':
                    Cs[j] = 1;
                    Gs[j] = 1;
                    Ts[j] = 1;
                    break;

                    // N = A,C,G or T
                default:
                    As[j] = 1;
                    Cs[j] = 1;
                    Gs[j] = 1;
                    Ts[j] = 1;
                    break;
            }
        }
        A_snps.push_back(As);
        C_snps.push_back(Cs);
        G_snps.push_back(Gs);
        T_snps.push_back(Ts);
    }

    #pragma omp parallel reduction(+:len)
    {
        #pragma omp for nowait
        for (uint64_t i = 0; i < n_seqs; i++)
        {
            // Cannot throw in an openmp block, short circuit instead
            boost::dynamic_bitset<> res(seq_length);

            for (uint64_t j = 0; j < n_seqs; j++) {

                res = A_snps[i] & A_snps[j];
                res |= C_snps[i] & C_snps[j];
                res |= G_snps[i] & G_snps[j];
                res |= T_snps[i] & T_snps[j];

                const int comp_snps = seq_length - res.count();

                distances[i].push_back(comp_snps);
            }
            len += distances[i].size();
        }
    }


    // Combine the lists from each thread
    std::vector<double> distances_all = combine_vectors(distances, len);

    return distances_all;
}

std::vector<double> get_distances_pa(const std::vector<std::vector<bool>>& matrix_in,
                                    const int n_threads)
{
    omp_set_num_threads(n_threads);

    size_t n_seqs = matrix_in.size();
    size_t seq_length = matrix_in.at(0).size();

    // Shared variables for openmp loop
    std::vector<std::vector<double>> distances(n_seqs);
    uint64_t len = 0;

    std::vector<boost::dynamic_bitset<>> present_all;
    std::vector<boost::dynamic_bitset<>> absent_all;
    for (const auto& seq : matrix_in)
    {
        boost::dynamic_bitset<> present(seq_length);
        for (size_t i = 0; i < seq_length; i++)
        {
            switch (seq[i])
            {
                case true:
                    present[i] = 1;
            }
        }
        boost::dynamic_bitset<> absent = present;
        absent.flip();
        present_all.push_back(present);
        absent_all.push_back(absent);
    }
    #pragma omp parallel reduction(+:len)
    {
        #pragma omp for nowait
        for (uint64_t i = 0; i < n_seqs; i++) {
            boost::dynamic_bitset<> res(seq_length);

            for (uint64_t j = 0; j < n_seqs; j++) {

                res = present_all[i] & present_all[j];
                res |= absent_all[i] & absent_all[j];
                const int comp_snps = seq_length - res.count();
                distances[i].push_back(comp_snps);
            }

            len += distances[i].size();
        }
    }

    // Combine the lists from each thread
    std::vector<double> distances_all = combine_vectors(distances, len);

    return distances_all;
}