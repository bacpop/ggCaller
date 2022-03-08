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

std::pair<float, bool> run_BALROG (const std::string& ORF_DNA,
                                   const std::string& upstream,
                                   const size_t& ORF_len,
                                   torch::jit::script::Module& ORF_model,
                                   torch::jit::script::Module& TIS_model,
                                   const float& prev_score,
                                   tbb::concurrent_unordered_map<size_t, float>& all_ORF_scores,
                                   tbb::concurrent_unordered_map<size_t, float>& all_TIS_scores)
{
    bool confident = true;

    // pull info from sequences
    const auto downstream = ORF_DNA.substr (0,19);
    const auto& encode_len = (ORF_len / 3) - 2;

    float TIS_prob = 0.5;
    float gene_prob = 0;

    // get float for start amino acid probability
    float start_prob = 0;

    // score TIS
    if (upstream.size())
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
            TIS_prob = pred[0].item<float>();

            all_TIS_scores.emplace(TIS_hash, TIS_prob);
        }
    }

    // encode start codon
    const auto start_codon = start_encode(downstream.substr(0,3)).out();

    const float ATG = ((start_codon == 0) ? 1 : 0) * weight_ATG;
    const float GTG = ((start_codon == 1) ? 1 : 0) * weight_GTG;
    const float TTG = ((start_codon == 2) ? 1 : 0) * weight_TTG;

    // multiply by TIS_prob weight
    TIS_prob *= weight_TIS_prob;

    // conduct check to ensure if max score for gene prob = 1, will be higher than prev_score
    {
        const float comb_prob = (1 * weight_gene_prob) + TIS_prob + ATG + GTG + TTG;

        float score = (comb_prob - probthresh) * ORF_len;

        // if next
        if (score <= prev_score)
        {
            return {0, false};
        }
    }

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

            auto logit = torch::logit(sub_seq);

            start_prob = logit[-1].item<float>();

            gene_prob = torch::special::expit(torch::mean(logit)).item<float>();

            all_ORF_scores.emplace(ORF_hash, gene_prob);
        }
    }

    const float comb_prob = (gene_prob * weight_gene_prob) + TIS_prob + ATG + GTG + TTG;

    float score = (comb_prob - probthresh) * ORF_len;

    // determine if not confident about start site i.e. little support for TIS and starting residue
    if (TIS_prob <= 0.5 && start_prob <= 0)
    {
        confident = false;
    }

    return {score, confident};
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
