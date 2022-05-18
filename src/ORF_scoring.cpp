#include "ORF_scoring.h"

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

std::vector<float> score_TIS (const std::vector<std::tuple<std::string, std::string, size_t>>& TIS_list,
                              torch::jit::script::Module& TIS_model,
                              const float& minimum_ORF_score,
                              tbb::concurrent_unordered_map<size_t, float>& all_TIS_scores)
{
    std::vector<float> TIS_scores(TIS_list.size());

    bool confident = true;

    // create index to hold positions of scores when added to padded_stack
    std::vector<size_t> pos_idx;
    std::vector<size_t> hash_idx;

    // hold size of vector in case padding required
    long tensor_size = 0;

    // create stack for batch processing
    std::vector<torch::Tensor> padded_stack;

    // score TIS
    for (int pos = 0; pos < TIS_list.size(); pos++)
    {
        const auto& entry = TIS_list.at(pos);
        const auto& upstream = std::get<0>(entry);
        const auto& downstream = std::get<1>(entry);

        if (upstream.size())
        {
            // reverse and encode the combined upstream + downstream sequences for scoring
            std::string combined = upstream + downstream.substr(3, downstream.size() - 3);
            std::reverse(combined.begin(), combined.end());

            size_t TIS_hash = hasher{}(combined);

            auto TIS_found = all_TIS_scores.find(TIS_hash);
            if (TIS_found != all_TIS_scores.end())
            {
                TIS_scores[pos] = TIS_found->second;
            } else
            {
                std::vector<int> encoded;

                for (const auto &c : combined)
                {
                    encoded.push_back(nuc_encode(c).out());
                }
                torch::Tensor t = torch::tensor(encoded, {torch::kInt64});

                if (!tensor_size)
                {
                    tensor_size = t.size(0);
                }

                padded_stack.push_back(std::move(t));
                pos_idx.push_back(pos);
                hash_idx.push_back(TIS_hash);
            }
        } else
        {
            // if no upstream sequence, set TIS score as 0.5
            TIS_scores[pos] = 0.5;
        }
    }

    if (!pos_idx.empty())
    {
        // pad tensor if only single sequence
        if (pos_idx.size() == 1)
        {
            torch::Tensor zeroes = torch::zeros({tensor_size}, torch::kInt64);
            padded_stack.push_back(std::move(zeroes));
        }
        torch::Tensor pred = predict(TIS_model, torch::stack(padded_stack), false);

        for (int idx = 0; idx < pos_idx.size(); idx++)
        {
            const auto& pos = pos_idx.at(idx);
            TIS_scores[pos] = pred[idx].item<float>();

            // add to shared map
            all_TIS_scores[hash_idx[idx]] = TIS_scores[pos];
        }
    }

    // check if ORF scores are high enough
    for (int pos = 0; pos < TIS_list.size(); pos++)
    {
        const auto& entry = TIS_list.at(pos);
        const auto& downstream = std::get<1>(entry);
        const auto& ORF_len = std::get<2>(entry);

        auto& TIS_prob = TIS_scores[pos];

        // encode start codon
        const auto start_codon = start_encode(downstream.substr(0,3)).out();

        const float ATG = ((start_codon == 0) ? 1 : 0) * weight_ATG;
        const float GTG = ((start_codon == 1) ? 1 : 0) * weight_GTG;
        const float TTG = ((start_codon == 2) ? 1 : 0) * weight_TTG;

        // conduct check to ensure if max score for gene prob = 1, will be higher than prev_score
        {
            // update TIS_prob to take into account start codon
            TIS_prob = (TIS_prob * weight_TIS_prob) + ATG + GTG + TTG;

            const float comb_prob = (1 * weight_gene_prob) + TIS_prob;

            float score = (comb_prob - probthresh) * ORF_len;

            // if score less than minimum, ignore
            if (score < minimum_ORF_score)
            {
                TIS_prob = 0;
            }
        }
    }

    return TIS_scores;
}

float score_gene (float& curr_prob,
                   const std::string& ORF_DNA,
                   const size_t& ORF_len,
                   torch::jit::script::Module& ORF_model,
                   tbb::concurrent_unordered_map<size_t, float>& all_ORF_scores)
{
    // pull info from sequences
    const auto downstream = ORF_DNA.substr (0,19);
    const auto& encode_len = (ORF_len / 3) - 2;

    float gene_prob = 0;

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

            gene_prob = torch::special::expit(torch::mean(torch::logit(sub_seq))).item<float>();

            all_ORF_scores.emplace(ORF_hash, gene_prob);
        }
    }

    // get weighted gene probability
    gene_prob *= weight_gene_prob;

    // curr_prob is TIS prob, use this to determine overall score
    const float comb_prob = gene_prob + curr_prob;

    // edit current probability in place
    curr_prob = (comb_prob - probthresh) * ORF_len;

    return gene_prob;
}

void score_cluster(float& curr_prob,
                   const float& gene_prob,
                   const std::string& ORF_DNA,
                   const size_t& ORF_len)
{
    // pull info from sequences
    const auto downstream = ORF_DNA.substr (0,19);
    const auto& encode_len = (ORF_len / 3) - 2;

    // curr_prob is TIS prob, use this to determine overall score
    const float comb_prob = gene_prob + curr_prob;

    // edit current probability in place
    curr_prob = (comb_prob - probthresh) * ORF_len;
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
