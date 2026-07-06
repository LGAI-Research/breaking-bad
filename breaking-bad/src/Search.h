#pragma once

#include "PDAG.h"
#include "ScorerInterface.h"

#include "spdlog/logger.h"

#include "EdgeQueueSet.h"

static FlatSet undirected_component_of(const PDAG& pdag, int seed);


class Search {
public:
    Search(int n_variables, ScorerInterface *scorer);
    Search(const Search &other);
    Search(const PDAG &initial_graph, ScorerInterface *scorer);

    void fit_xges(bool extended_search);
    void fit_xges0(bool extended_search);
    void fit_xges_variant2();  
    void fit_xges_variant3();
    void fit_xges_variant4();

    void fit_lges(bool score_based=false, bool prune=false);
    void fit_lges_variant1(bool score_based=false, bool prune=false);
    void fit_lges_variant2(bool score_based=false, bool prune=false);
    void fit_lges_variant3(bool score_based=false, bool prune=false);
    void fit_lges_variant4(bool score_based=false, bool prune=false);

    void heuristic_lges(std::vector<Insert> &candidate_inserts,
                           std::vector<Reverse> &candidate_reverses,
                           std::vector<Delete> &candidate_deletes,
                           UnblockedPathsMap &unblocked_paths_map,
                           bool initialize_inserts,
                            bool score_based=false,
                            bool prune=false);

    void find_inserts_to_y_for_lges(int y, std::vector<Insert> &candidate_inserts,const PDAG* DAG_ptr=nullptr ,bool score_based=false, bool prune=false,
                                    int parent_x=-1,bool positive_only=true) const;

    void update_operator_candidates_naive_for_lges(std::vector<Insert> &candidate_inserts,
                                                std::vector<Reverse> &candidate_reverses,
                                                std::vector<Delete> &candidate_deletes,
                                                const PDAG* DAG_ptr=nullptr,
                                                bool score_based=false,
                                                bool prune=false)const;

    void update_operator_candidates_efficient_for_lges(EdgeModificationsMap &edge_modifications,
                                                std::vector<Insert> &candidate_inserts,
                                                std::vector<Reverse> &candidate_reverses,
                                                std::vector<Delete> &candidate_deletes,
                                                UnblockedPathsMap &unblocked_paths_map,
                                                const PDAG* DAG_ptr=nullptr,
                                                bool score_based=false,
                                                bool prune=false);


    double recompute_total_score_from_pdag();


    enum class AlgorithmType {
        XGES,
        GES,
        OPS,
        BOSS  
    };
    void run_research(UnblockedPathsMap &unblocked_paths_map, int algo_type = 0);




    void hub_only_component_incoming_reset_and_research_edge_only(UnblockedPathsMap &unblocked_paths_map, size_t hub_threshold, double top_k_ratio = 1.0,int algo_type=0);
    void node_incoming_reset_batch(UnblockedPathsMap &unblocked_paths_map, size_t hub_threshold, double top_k_ratio = 1.0,int algo_type=0);
    void node_incoming_reset_batch_ablation_deletepa(UnblockedPathsMap &unblocked_paths_map, size_t hub_threshold, double top_k_ratio = 1.0,int algo_type=0,int max_enum_neighbors=-1);
    void single_edge_delete_clean_only(UnblockedPathsMap &unblocked_paths_map,int algo_type=0);
    void single_edge_delete_best_first(UnblockedPathsMap &unblocked_paths_map,int algo_type=0);
    


    void initialize_candidates_pure(std::vector<Insert> &candidate_inserts,
        std::vector<Reverse> &candidate_reverses,
        std::vector<Delete> &candidate_deletes) const;

    void fit_ops(bool use_reverse);
    void fit_ops_variant1();
    void fit_ops_variant2();
    void fit_ops_variant3();
    void fit_ops_variant4();

    void fit_ges(bool use_reverse);
    void fit_ges_variant1();
    void fit_ges_variant2();
    void fit_ges_variant3();
    void fit_ges_variant4();

    double set_boss_initial_graph(std::string boss_result_path,std::string input_path,double alpha);
    void heuristic_backward_only(std::vector<Delete> &candidate_deletes,UnblockedPathsMap &unblocked_paths_map);
    void fit_boss_variant0();
    void fit_boss_variant1();
    void fit_boss_variant2();
    void fit_boss_variant3();
    void fit_boss_variant4();

    double get_score() const;
    double get_initial_score() const;
    const PDAG &get_pdag() const;

    std::unique_ptr<PDAG> ground_truth_pdag;
    std::map<std::string, double> statistics;

private:
    int n_variables;
    ScorerInterface *scorer;
    PDAG pdag;
    const double initial_score = 0;
    double total_score = 0;
    std::shared_ptr<spdlog::logger> _logger;

    void heuristic_xges0(std::vector<Insert> &candidate_inserts,
                         std::vector<Reverse> &candidate_reverses,
                         std::vector<Delete> &candidate_deletes,
                         UnblockedPathsMap &unblocked_paths_map,
                         bool initialize_inserts = true);

    void update_operator_candidates_naive(std::vector<Insert> &candidate_inserts,
                                          std::vector<Reverse> &candidate_reverses,
                                          std::vector<Delete> &candidate_deletes) const;
    void update_operator_candidates_efficient(EdgeModificationsMap &edge_modifications,
                                              std::vector<Insert> &candidate_inserts,
                                              std::vector<Reverse> &candidate_reverses,
                                              std::vector<Delete> &candidate_deletes,
                                              UnblockedPathsMap &unblocked_paths_map);

    void block_each_edge_and_research(UnblockedPathsMap &unblocked_paths_map, int algo_type = 0);

  
    void find_inserts_to_y(int y, std::vector<Insert> &candidate_inserts,
                           int parent_x = -1, bool positive_only = true) const;
                           
    void find_inserts_to_y_variants(int y, std::vector<Insert> &candidate_inserts,
                                   int parent_x = -1, bool positive_only = true) const;

    void find_delete_to_y_from_x(int y, int x, std::vector<Delete> &candidate_deletes,
                                 bool positive_only = true) const;
    void find_deletes_to_y(int y, std::vector<Delete> &candidate_deletes,
                           bool positive_only = true) const;

    void find_reverse_to_y_from_x(int y, int x,
                                  std::vector<Reverse> &candidate_reverses) const;
    void find_reverse_to_y(int y, std::vector<Reverse> &candidate_reverses) const;
};