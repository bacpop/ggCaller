#ifndef ORF_SCORING_H
#define ORF_SCORING_H

#include <torch/torch.h>
#include <torch/script.h>
#include "definitions.h"
#include "translation.h"
#include "graph.h"
#include <numeric>

// All following are internal parameters. Change at your own risk.
const double weight_gene_prob = 0.9746869839852076;
const double weight_TIS_prob = 0.25380288790532707;
const double score_threshold = 0.47256101519707244;
const double weight_ATG = 0.84249804151264;
const double weight_GTG = 0.7083689705744909;
const double weight_TTG = 0.7512400826652517;
const double maxprob = weight_gene_prob + weight_TIS_prob + weight_ATG;
const double probthresh = score_threshold * maxprob;

vector<size_t> sort_indexes(vector<torch::Tensor> &v);

torch::Tensor tokenized_aa_seq(const std::string& aa_seq);

torch::Tensor predict(torch::jit::script::Module& module,
                      const torch::Tensor& seq_tensor,
                      const bool gene);


std::unordered_map<size_t, double> run_BALROG(const ColoredCDBG<MyUnitigMap>& ccdbg,
                                              const std::vector<Kmer>& head_kmer_arr,
                                              const ORFVector& ORF_vector,
                                              torch::jit::script::Module& ORF_model,
                                              torch::jit::script::Module& TIS_model,
                                              const int overlap,
                                              const float& minimum_ORF_score,
                                              const int ORF_batch_size,
                                              const int TIS_batch_size,
                                              robin_hood::unordered_map<size_t, double>& all_ORF_scores,
                                              robin_hood::unordered_map<size_t, double>& all_TIS_scores);

#endif //ORF_SCORING_H
