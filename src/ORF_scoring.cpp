#include "ORF_scoring.h"
std::mutex mtx3;
std::mutex mtx4;

vector<size_t> sort_indexes(vector<torch::Tensor> &v) {

    // initialize original index locations
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // first sort indexes based on comparing values in v
    std::stable_sort(idx.begin(), idx.end(),
                [&v](size_t i1, size_t i2) {return v[i1].size(0) < v[i2].size(0);});

    // then sort the tensor vector in place
    std::stable_sort(v.begin(), v.end(), []
            (const torch::Tensor& i1, const torch::Tensor& i2){
        return i1.size(0) < i2.size(0);
    });

    return idx;
}

torch::Tensor tokenized_aa_seq(const std::string& aa_seq)
{
    std::vector<int> tokens = (tokenise(aa_seq).out());

    std::vector<torch::Tensor> padded_stack;

    torch::Tensor t = torch::tensor(tokens, {torch::kInt64});

    torch::Tensor zeroes = torch::zeros({t.size(0)}, torch::kInt64);

    padded_stack.push_back(std::move(t));
    padded_stack.push_back(std::move(zeroes));

    return torch::stack(padded_stack);
}

std::unordered_map<size_t, double> run_BALROG (const ColoredCDBG<MyUnitigMap>& ccdbg,
                                               const std::vector<Kmer>& head_kmer_arr,
                                               const ORFVector& ORF_vector,
                                               torch::jit::script::Module& ORF_model,
                                               torch::jit::script::Module& TIS_model,
                                               const int overlap,
                                               const float& minimum_ORF_score,
                                               const int ORF_batch_size,
                                               const int TIS_batch_size,
                                               robin_hood::unordered_map<size_t, double>& all_ORF_scores,
                                               robin_hood::unordered_map<size_t, double>& all_TIS_scores)
{
    // initialise map to return
    std::unordered_map<size_t, double> score_map;

    for (int i = 0;  i < ORF_vector.size(); i++)
    {
        const auto& ORF_info = ORF_vector.at(i);
        const auto ORF_DNA = generate_sequence_nm(std::get<0>(ORF_info), std::get<1>(ORF_info), overlap, ccdbg, head_kmer_arr);
        const auto upstream = generate_sequence_nm(std::get<3>(ORF_info), std::get<4>(ORF_info), overlap, ccdbg, head_kmer_arr);
        const auto downstream = ORF_DNA.substr (0,19);

        const auto& encode_len = (std::get<2>(ORF_info) / 3) - 2;

        double TIS_prob = 0.5;
        double gene_prob = 0;

        // score TIS
        if (upstream.size() == 16)
        {
            // reverse and encode the combined upstream + downstream sequences for scoring
            std::string combined = upstream + downstream.substr(3, downstream.size() - 3);
            std::reverse(combined.begin(), combined.end());

            size_t TIS_hash = hasher{}(combined);

            auto TIS_found = all_TIS_scores.find(TIS_hash);
            if (TIS_found != all_TIS_scores.end())
            {
                TIS_prob = TIS_found->second;
            } else
            {
                std::vector<int> encoded;

                for (const auto &c : combined)
                {
                    encoded.push_back(nuc_encode(c).out());
                }
                torch::Tensor t = torch::tensor(encoded, {torch::kInt64});

                torch::Tensor zeroes = torch::zeros({t.size(0)}, torch::kInt64);

                std::vector<torch::Tensor> padded_stack;
                padded_stack.push_back(std::move(t));
                padded_stack.push_back(std::move(zeroes));

                torch::Tensor pred = predict(TIS_model, torch::stack(padded_stack), false);
                TIS_prob = pred[0].item<double>();

                mtx3.lock();
                all_TIS_scores[TIS_hash] = TIS_prob;
                mtx4.lock();
            }
        }

        // encode start codon
        const auto start_codon = start_encode(downstream.substr(0,3)).out();

        // get gene score either from stored scores or calculate it
        {
            // get aa sequence, remove start and end codon
            const auto ORF_aa = (translate(ORF_DNA)).aa().substr(1,(ORF_DNA.size() / 3) - 2);

            size_t ORF_hash = hasher{}(ORF_aa);

            auto ORF_found = all_ORF_scores.find(ORF_hash);
            if (ORF_found != all_ORF_scores.end())
            {
                gene_prob = ORF_found->second;
            } else
            {
                torch::Tensor t = tokenized_aa_seq(ORF_aa);

                torch::Tensor pred = predict(ORF_model, t, true);

                auto sub_seq = pred.index({0, 0, torch::indexing::Slice(torch::indexing::None, encode_len)});

                gene_prob = torch::special::expit(torch::mean(torch::logit(sub_seq))).item<double>();

                mtx4.lock();
                all_ORF_scores[ORF_hash] = gene_prob;
                mtx4.unlock();
            }
        }

        // pull info from ORF_ID
        const int& ORF_len = std::get<2>(ORF_vector.at(i));

        const bool ATG = (start_codon == 0) ? 1 : 0;
        const bool GTG = (start_codon == 1) ? 1 : 0;
        const bool TTG = (start_codon == 2) ? 1 : 0;

        const double comb_prob = (gene_prob * weight_gene_prob + TIS_prob * weight_TIS_prob) +
                                 (ATG * weight_ATG) + (GTG * weight_GTG) + (TTG * weight_TTG);

        const double score = (comb_prob - probthresh) * ORF_len;
        if (score >= minimum_ORF_score)
        {
            score_map[i] = score;
        }
    }

    return score_map;
}

torch::Tensor predict(torch::jit::script::Module& module,
                      const torch::Tensor& seq_tensor,
                      const bool gene)
{
    int num_classes = 21;

    if (!gene)
    {
        num_classes = 4;
    }

    module.eval();
    torch::NoGradGuard nograd;

    auto output = torch::special::expit(module.forward({ torch::nn::functional::one_hot(seq_tensor, num_classes)
            .permute({0, 2, 1}).to(torch::kFloat32) }).toTensor());

    return output;
}
