#include "ORF_scoring.h"

vector<size_t> sort_indexes(vector<torch::Tensor> &v) {

    // initialize original index locations
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // first sort indexes based on comparing values in v
    std::stable_sort(idx.begin(), idx.end(),
                [&v](size_t i1, size_t i2) {return v[i1].size(0) < v[i2].size(0);});

    // then sort the tensor vector in place
    std::sort(v.begin(), v.end(), []
            (const torch::Tensor& i1, const torch::Tensor& i2){
        return i1.size(0) < i2.size(0);
    });

    return idx;
}

torch::Tensor tokenized_aa_seq(const std::string& aa_seq)
{
    std::vector<int> tokens = (tokenise(aa_seq).out());

    torch::Tensor t = torch::from_blob(tokens.data(), {static_cast<long>(tokens.size())}).to(torch::kInt64);;

    return t;
}

std::pair<std::vector<torch::Tensor>, std::vector<std::pair<std::string, std::string>>> get_ORF_info (const ColoredCDBG<MyUnitigMap>& ccdbg,
                                                                                                      const std::vector<Kmer>& head_kmer_arr,
                                                                                                      const ORFVector& ORF_vector,
                                                                                                      const int overlap)
{
    std::vector<torch::Tensor> ORF_seq_enc(ORF_vector.size());
    std::vector<std::pair<std::string, std::string>> TIS_seqs(ORF_vector.size());

    for (int i = 0;  i < ORF_vector.size(); i++)
    {
        const auto& ORF_info = ORF_vector.at(i);
        const auto ORF_DNA = generate_sequence_nm(std::get<0>(ORF_info), std::get<1>(ORF_info), overlap, ccdbg, head_kmer_arr);
        const auto upstream_TIS = generate_sequence_nm(std::get<3>(ORF_info), std::get<4>(ORF_info), overlap, ccdbg, head_kmer_arr);
        const auto downstream_TIS = ORF_DNA.substr (0,18);

        const auto ORF_aa = (translate(ORF_DNA)).aa();

        TIS_seqs[i] = {upstream_TIS, downstream_TIS};
        ORF_seq_enc[i] = tokenized_aa_seq(ORF_aa);
    }

    return {ORF_seq_enc, TIS_seqs};
}

torch::jit::script::Module load_model(const std::string& model_file)
{
    torch::jit::script::Module module;
    try {
        // Deserialize the ScriptModule from a file using torch::jit::load().
        module = torch::jit::load(model_file);
    }
    catch (const c10::Error& e) {
        std::cerr << "error loading the model: " << model_file << "n";
    }

    return module;
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

    std::vector<torch::jit::IValue> inputs;
    torch::Tensor X_enc = torch::nn::functional::one_hot(seq_tensor, num_classes).permute(c10::IntArrayRef{0, 2, 1}).to(torch::kFloat32);
    inputs.push_back(std::move(X_enc));
    torch::Tensor output = torch::special::expit(module.forward(inputs).toTensor());

    return output;
}

std::unordered_map<size_t, double> run_BALROG(const ColoredCDBG<MyUnitigMap>& ccdbg,
                                               const std::vector<Kmer>& head_kmer_arr,
                                               const ORFVector& ORF_vector,
                                               torch::jit::script::Module& ORF_model,
                                               torch::jit::script::Module& TIS_model,
                                               const int overlap,
                                               const float& minimum_ORF_score,
                                               const int ORF_batch_size,
                                               const int TIS_batch_size)
{
    // initialise map to return
    std::unordered_map<size_t, double> score_map;

    const auto encoded_seqs = get_ORF_info(ccdbg, head_kmer_arr, ORF_vector, overlap);

    const auto& ORF_seq_enc = encoded_seqs.first;
    const auto& TIS_seqs = encoded_seqs.second;

    // copy and sort on size
    auto ORF_seq_enc_sorted = ORF_seq_enc;
    const auto length_idx = sort_indexes(ORF_seq_enc_sorted);

    // call TIS scores
    std::vector<double> TIS_prob_list(ORF_seq_enc.size());
    std::vector<int> start_codon_list(ORF_seq_enc.size());
    {
        std::vector<int> TIS_idx;
        std::vector<torch::Tensor> TIS_tensor;

        const int vec_size = static_cast<int>(TIS_seqs.size());
        for (int i = 0; i < TIS_seqs.size(); i++)
        {
            const auto &upstream = TIS_seqs.at(i).first;
            const auto &downstream = TIS_seqs.at(i).second;
            if (upstream.size() == 16)
            {
                // reverse and encode the combined upstream + downstream sequences for scoring
                std::string combined = upstream + downstream;
                std::reverse(combined.begin(), combined.end());
                std::vector<int> encoded;

                for (const auto &c : combined)
                {
                    encoded.push_back(nuc_encode(c).out());
                }
                torch::Tensor t = torch::from_blob(encoded.data(), {static_cast<long>(encoded.size())}).to(torch::kInt64);

                // add to TIS_tensor and keep track of indices added
                TIS_tensor.push_back(std::move(t));
                TIS_idx.push_back(i);
            } else
            {
                TIS_prob_list[i] = 0.5;
            }

            // encode start codon
            const auto start_codon = start_encode(downstream.substr(0,3)).out();
            start_codon_list[i] = start_codon;
        }

        // batch score TIS_seqs
        for (int i = 0; i < vec_size; i += TIS_batch_size)
        {
            std::vector<torch::Tensor> batch;

            // get the batch for processing
            if (i < (vec_size - TIS_batch_size))
            {
                batch = std::vector<torch::Tensor>(TIS_tensor.begin() + i, TIS_tensor.begin() + i + TIS_batch_size);
            } else
            {
                batch = std::vector<torch::Tensor>(TIS_tensor.begin() + i, TIS_tensor.end());
            }

            // initialise return tensor
            torch::Tensor seq_tensor = torch::stack(batch);

            // predict scores
            torch::Tensor pred_all = predict(TIS_model, seq_tensor, false);

            // add all to list
            double* ptr = (double*)pred_all.data_ptr();
            for (int j = 0; j < pred_all.numel(); j++)
            {
                const int& idx = TIS_idx.at(i + j);
                TIS_prob_list[idx] = *ptr++;
            }
        }
    }

    // call ORF scores
    const int vec_size = static_cast<int>(ORF_seq_enc_sorted.size());
    for (int i = 0; i < vec_size; i += ORF_batch_size)
    {
        std::vector<torch::Tensor> batch;

        // get the batch for processing
        if (i < (vec_size - ORF_batch_size))
        {
            batch = std::vector<torch::Tensor>(ORF_seq_enc_sorted.begin() + i, ORF_seq_enc_sorted.begin() + i + ORF_batch_size);
        } else
        {
            batch = std::vector<torch::Tensor>(ORF_seq_enc_sorted.begin() + i, ORF_seq_enc_sorted.end());
        }

        // determine the lengths for sequences
        std::vector<int> seq_lens;
        for (const auto& entry : batch)
        {
            seq_lens.push_back(entry.size(0));
        }

        // get largest element
        const int max_size = *std::max_element(seq_lens.begin(), seq_lens.end());

        // initialise return tensor
        torch::Tensor seq_tensor = torch::zeros({static_cast<long>(batch.size()), max_size}, torch::kInt64);
        for (int j = 0; j < batch.size(); j++)
        {
            seq_tensor.index_put_({j, torch::indexing::Slice(torch::indexing::None, seq_lens.at(j))}, batch.at(j));
        }

        // predict scores
        torch::Tensor pred_all = predict(ORF_model, seq_tensor, true);

        std::vector<double> pred(seq_lens.size());
        for (int j = 0; j < seq_lens.size(); j++)
        {
            auto sub_seq = pred_all.index({j, 0, torch::indexing::Slice(torch::indexing::None, seq_lens.at(j))});
            double gene_prob = torch::special::expit(sub_seq.mean()).item<double>();
            const int ORF_ID = i + j;
            const int& ORF_len = std::get<2>(ORF_vector.at(ORF_ID));
            const auto& TIS_prob = TIS_prob_list.at(ORF_ID);
            const auto& start_codon = start_codon_list.at(ORF_ID);

            const bool ATG = (start_codon == 0) ? 1 : 0;
            const bool GTG = (start_codon == 1) ? 1 : 0;
            const bool TTG = (start_codon == 2) ? 1 : 0;

            const double comb_prob = (gene_prob * weight_gene_prob + TIS_prob * weight_TIS_prob) +
                    (ATG * weight_ATG) + (GTG * weight_GTG) + (TTG * weight_TTG);

            const double score = (comb_prob - probthresh) * ORF_len;
            if (score >= minimum_ORF_score)
            {
                score_map[ORF_ID] = score;
            }
        }
    }

    return score_map;
}
