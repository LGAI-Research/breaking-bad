#include <spdlog/fmt/bundled/ostream.h> 
#include <spdlog/fmt/bundled/std.h>    

#include "Search.h"
#include "set_ops.h"
#include "utils.h"
#include "EdgeQueueSet.h"
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"
#include <stack>
#include <random>


// -------------------------------------------------------------
// [Helper] Research Dispatcher
// -------------------------------------------------------------
// algo_type: 0=XGES, 1=GES, 2=OPS, 3=BOSS
void Search::run_research(UnblockedPathsMap &unblocked_paths_map,int algo_type) {
    if (algo_type == 0) { // XGES
        std::vector<Insert> i; std::vector<Reverse> r; std::vector<Delete> d;
        heuristic_xges0(i, r, d, unblocked_paths_map, true);
    } 
    else if (algo_type == 1 || algo_type == 3) { // GES or BOSS
        fit_ges(false); 
    } 
    else if (algo_type == 2) { // OPS
        // OPS: Permutation Space Search
        fit_ops(false); 
    }
    else if(algo_type==4){ // LGES-Safe
        fit_lges(true,false);
    }
    else if(algo_type==5){ //LGES-Conservative
        fit_lges(true,true);
    }
}


// ============================================================
// Exhaustive H enumeration + global heap by approx_score
//  - For each selected hub y:
//      enumerate ALL subsets H ⊆ Ne(y)
//      keep only valid: Ne(y)\H is clique
//      compute approx_score(y,H)
//  - Push ALL candidates into one heap, pop in score order
//  - Trial: DeletePa(y,H) (forced delete ALL parents) -> repair -> research -> accept
// ============================================================

struct DeletePaCand {
    int y;
    FlatSet H;
    double approx_score; // heap key (bigger is better)

    bool operator<(const DeletePaCand& other) const {
        return approx_score < other.approx_score; // max-heap
    }
};

static inline FlatSet flatset_difference(const FlatSet& a, const FlatSet& b) {
    FlatSet out;
    out.reserve(a.size());
    for (int v : a) if (b.find(v) == b.end()) out.insert(v);
    return out;
}

static inline FlatSet flatset_union_singleton(const FlatSet& a, int s) {
    FlatSet out = a;
    out.insert(s);
    return out;
}

// Apply DeletePa(y,H) WITHOUT finalize.
//  - delete ALL directed parent edges x->y (forced)
//  - orient y-h as y->h for h in H (only if undirected exists)
//  - record edge modifications + accumulate affected edges
static inline void apply_deletepa_only_no_finalize(
    PDAG& G,
    int y,
    const FlatSet& H,
    EdgeModificationsMap& em,
    EdgeQueueSet& edges_to_check_acc
) {
    // delete all parents x -> y (forced)
    std::vector<int> parents_vec;
    {
        const auto& parents = G.get_parents(y);
        parents_vec.assign(parents.begin(), parents.end());
    }
    for (int x : parents_vec) {
        G.remove_edge_only_no_finalize(x, y, /*directed_x_to_y=*/true, em, edges_to_check_acc);
    }

    // orient y - h  to y -> h
    for (int h : H) {
        if (!G.has_undirected_edge(y, h)) continue;
        G.remove_undirected_edge(y, h);
        G.add_directed_edge(y, h);
        em.update_edge_directed(y, h, EdgeType::UNDIRECTED);
        G.add_adjacent_edges(y, h, edges_to_check_acc);
    }

    // conservative adjacency accumulation
    for (int x : parents_vec) {
        G.add_adjacent_edges(x, y, edges_to_check_acc);
    }
}

// Enumerate ALL subsets H ⊆ neighbors (including empty).
// Returns H as FlatSet via callback to avoid huge memory if you want streaming.
// Here we *do* callback-style; caller can filter/score/push to heap.
template <class F>
static inline void enumerate_all_H_subsets(
    const FlatSet& neighbors,
    int max_enum_neighbors, // -1 means no limit (DANGEROUS)
    F&& fn_take_H           // void(const FlatSet& H)
) {
    std::vector<int> nv;
    nv.reserve(neighbors.size());
    for (int v : neighbors) nv.push_back(v);

    const int k = (int)nv.size();
    if (max_enum_neighbors >= 0 && k > max_enum_neighbors) {
        // hard guard to avoid 2^k explosion
        throw std::runtime_error(
            "enumerate_all_H_subsets: |Ne(y)|=" + std::to_string(k) +
            " exceeds max_enum_neighbors=" + std::to_string(max_enum_neighbors)
        );
    }

    // NOTE: 2^k iteration: only feasible for small k
    const uint64_t total = (k >= 64) ? 0ULL : (1ULL << k);
    if (k >= 64) {
        // extremely unsafe; you can implement recursion, but it will still be astronomically large.
        throw std::runtime_error("enumerate_all_H_subsets: k>=64 is not supported.");
    }

    FlatSet H;
    for (uint64_t mask = 0; mask < total; ++mask) {
        H.clear();
        // build subset
        for (int i = 0; i < k; ++i) {
            if (mask & (1ULL << i)) H.insert(nv[i]);
        }
        fn_take_H(H);
    }
}


// ============================================================
// XGES method: Exhaustive DeletePa(y,H) ablation
// ============================================================

void Search::node_incoming_reset_batch_ablation_deletepa(
    UnblockedPathsMap &unblocked_paths_map,
    size_t hub_threshold,
    double top_k_ratio,        // hub selection ratio
    int algo_type,
    int max_enum_neighbors     // safety guard; set -1 to truly do all (danger)
) {
    constexpr double eps_gain = 1e-7;
    bool global_improved = true;

    struct HubCandidate {
        size_t indeg;
        size_t undirected_deg;
        int y;
    };

    while (global_improved) {
        global_improved = false;

        // ----------------------------------------------------
        // 1) Hub selection
        // ----------------------------------------------------
        const auto& all_nodes = pdag.get_nodes_variables();
        std::vector<HubCandidate> all_candidates;
        all_candidates.reserve(all_nodes.size());

        for (int y : all_nodes) {
            size_t indeg  = pdag.get_parents(y).size();
            size_t un_deg = pdag.get_neighbors(y).size();
            all_candidates.push_back({indeg, un_deg, y});
        }

        std::sort(all_candidates.begin(), all_candidates.end(),
                  [](const HubCandidate& a, const HubCandidate& b) {
                      size_t da = a.indeg + a.undirected_deg;
                      size_t db = b.indeg + b.undirected_deg;
                      if (da != db) return da > db;
                      return a.indeg > b.indeg;
                  });

        size_t top_k_limit = static_cast<size_t>(all_nodes.size() * top_k_ratio);
        if (!all_nodes.empty() && top_k_limit == 0) top_k_limit = 1;

        std::vector<HubCandidate> hubs;
        hubs.reserve(top_k_limit);
        for (const auto& cand : all_candidates) {
            if (hubs.size() >= top_k_limit) break;
            if (cand.indeg < hub_threshold) break;
            hubs.push_back(cand);
        }
        if (hubs.empty()) break;

        // ----------------------------------------------------
        // 2) Build GLOBAL heap of ALL valid (y,H) candidates
        //    ordered by approx_score
        // ----------------------------------------------------
        std::vector<DeletePaCand> heap;
        // reserve is hard to estimate because exhaustive can be huge

        for (const auto& hub : hubs) {
            const int y = hub.y;
            const auto& parents_y = pdag.get_parents(y);
            if (parents_y.empty()) continue;

            const auto& neighbors_y = pdag.get_neighbors(y);

            // local Δ_y = S(y, empty) - S(y, Pa(y))
            double delta_y = 0.0;
            {
                FlatSet empty;
                delta_y = scorer->local_score(y, empty) - scorer->local_score(y, parents_y);
            }

            // enumerate all H ⊆ Ne(y)
            enumerate_all_H_subsets(
                neighbors_y,
                max_enum_neighbors,
                [&](const FlatSet& H) {

                    // validity: Ne(y)\H is clique
                    FlatSet Ne_minus_H = flatset_difference(neighbors_y, H);
                    if (!pdag.is_clique(Ne_minus_H)) return;

                    // approx_score(y,H) = Δ_y + Σ_{h∈H} (S(h, Pa(h)∪{y}) - S(h, Pa(h)))
                    double approx = delta_y;
                    for (int h : H) {
                        const auto& parents_h = pdag.get_parents(h);
                        FlatSet parents_h_plus_y = flatset_union_singleton(parents_h, y);
                        approx += scorer->local_score(h, parents_h_plus_y) - scorer->local_score(h, parents_h);
                    }

                    DeletePaCand cand;
                    cand.y = y;
                    cand.H = H;
                    cand.approx_score = approx;
                    heap.push_back(std::move(cand));
                }
            );
        }

        if (heap.empty()) break;

        std::make_heap(heap.begin(), heap.end()); // max-heap by approx_score

        // ----------------------------------------------------
        // 3) Pop candidates by priority and do expensive trial
        // ----------------------------------------------------
        while (!heap.empty()) {
            std::pop_heap(heap.begin(), heap.end());
            DeletePaCand cand = std::move(heap.back());
            heap.pop_back();

            const int y = cand.y;
            const FlatSet& H = cand.H;

            // Trial copy
            statistics["deletepa_trials"]++;
            Search trial(*this);
            UnblockedPathsMap upm_copy = unblocked_paths_map;

            EdgeModificationsMap em;
            EdgeQueueSet edges_to_check;

            // Apply forced DeletePa(y,H) (no finalize)
            apply_deletepa_only_no_finalize(trial.pdag, y, H, em, edges_to_check);

            // Repair CPDAG
            trial.pdag.complete_cpdag_global(em);

            // safer scoring
            trial.total_score = trial.recompute_total_score_from_pdag();

            // Research
            if (algo_type == 0) {
                std::vector<Insert>  ci;
                std::vector<Reverse> cr;
                std::vector<Delete>  cd;

                trial.update_operator_candidates_efficient(em, ci, cr, cd, upm_copy);
                trial.heuristic_xges0(ci, cr, cd, upm_copy, false);
            } else {
                
                trial.run_research(upm_copy, algo_type);
            }

            const double gain = trial.total_score - total_score;
            const bool changed = !(pdag == trial.pdag);

            if (gain > eps_gain && changed) {
                _logger->debug(
                    "[Phase2-DeletePa-EXHAUSTIVE] Accepted y={}, |Pa(y)|={}, |H|={}, approx={}, gain={}",
                    y, pdag.get_parents(y).size(), H.size(), cand.approx_score, gain
                );

                total_score = trial.total_score;
                pdag = std::move(trial.pdag);
                unblocked_paths_map = std::move(upm_copy);

                statistics["extended_search-accepted_phase2"] += 1;
                global_improved = true;
                break; // greedy restart: recompute hubs + re-enumerate ALL
            }

            else{
                statistics["extended_search-rejected_phase2"] += 1;
            }
        }
    }
}




static FlatSet undirected_component_of(const PDAG& pdag, int seed) {
    FlatSet comp; 
    std::vector<int> stack; 

    stack.push_back(seed);
    comp.insert(seed);

    while (!stack.empty()) {
        int v = stack.back();
        stack.pop_back();

        for (int nb : pdag.get_neighbors(v)) { // undirected edges only
            if (!comp.contains(nb)) {
                comp.insert(nb); 
                stack.push_back(nb);
            }
        }
    }
    return comp;
}

using namespace std::chrono;

/** @brief Constructor of the XGES class.
 *
 * @param n_variables The number of variables in the dataset.
 * @param scorer The scorer to use to score the PDAGs.
 *
 */
Search::Search(const int n_variables, ScorerInterface *scorer)
    : n_variables(n_variables), scorer(scorer), pdag(n_variables),
      initial_score(scorer->score_pdag(pdag)) {
    total_score = initial_score;
    _logger = spdlog::get("stdout_logger");
}

Search::Search(const Search &other)
    : n_variables(other.n_variables), scorer(other.scorer), pdag(other.pdag),
      initial_score(other.initial_score), total_score(other.total_score),
      _logger(other._logger) {
    // The pointer to the ground truth PDAG is not copied. Change if needed.
}


Search::Search(const PDAG &initial_graph, ScorerInterface *scorer)
    : n_variables((int)initial_graph.get_nodes_variables().size()),
      scorer(scorer),
      pdag(initial_graph),
      initial_score(scorer->score_pdag(initial_graph)) {
    total_score = initial_score;
    _logger = spdlog::get("stdout_logger");
}


/** @brief Fit the XGES algorithm using the scorer provided in the constructor.
 *
 * @param extended_search If true, the extended search of XGES is performed. Otherwise,
 * only XGES-0 is performed.
 */
void Search::fit_xges(bool extended_search) {
    std::vector<Insert> candidate_inserts; // empty at the beginning
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map;

    candidate_inserts.reserve(100 * n_variables); // empty at the beginning
    candidate_inserts.reserve(n_variables); // empty at the beginning       
    candidate_deletes.reserve(n_variables); // empty at the beginning

    heuristic_xges0(candidate_inserts, candidate_reverses, candidate_deletes,
                    unblocked_paths_map, true);

    if (extended_search) { block_each_edge_and_research(unblocked_paths_map,0); }
}


void Search::fit_xges0(bool extended_search) {
    std::vector<Insert> candidate_inserts; // empty at the beginning
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map;

    candidate_inserts.reserve(100 * n_variables); // empty at the beginning
    candidate_inserts.reserve(n_variables); // empty at the beginning       
    candidate_deletes.reserve(n_variables); // empty at the beginning

    heuristic_xges0(candidate_inserts, candidate_reverses, candidate_deletes,
                    unblocked_paths_map, true);
}



double Search::recompute_total_score_from_pdag() {
    return scorer->score_pdag(this->pdag);
}

void Search::hub_only_component_incoming_reset_and_research_edge_only(UnblockedPathsMap &unblocked_paths_map, size_t hub_threshold, double top_k_ratio, int algo_type)
{
    // Top K% node (Example: 0.1 or 1.0)
    size_t min_hub_degree = hub_threshold; 
    constexpr double eps_gain = 1e-7;

    bool global_improved = true;

    while (global_improved) {
        global_improved = false;
        
        const auto& all_nodes = pdag.get_nodes_variables(); 
        struct HubCandidate {
            size_t indeg; 
            size_t undirected_deg; 
            int y;
        };
        std::vector<HubCandidate> all_candidates; 
        all_candidates.reserve(all_nodes.size());

        for (int y : all_nodes) { 
            size_t indeg = pdag.get_parents(y).size();     
            size_t un_deg = pdag.get_neighbors(y).size();  
            all_candidates.push_back({indeg, un_deg, y});
        }

        std::sort(all_candidates.begin(), all_candidates.end(),
                  [](const HubCandidate& a, const HubCandidate& b){ 
                      if (a.indeg != b.indeg) return a.indeg > b.indeg;
                      return (a.indeg + a.undirected_deg) > (b.indeg + b.undirected_deg); 
                  });


        // top_k_ratio
        size_t top_k_limit = (size_t)(all_nodes.size() * top_k_ratio);
        
        if (all_nodes.size() > 0 && top_k_limit == 0) top_k_limit = 1; 

        std::vector<HubCandidate> hubs;
        hubs.reserve(top_k_limit);

        for (const auto& cand : all_candidates) {
            if (hubs.size() >= top_k_limit) break; // Ratio check
            if (cand.indeg < min_hub_degree) break; // Min degree check
            hubs.push_back(cand);
        }

        if (hubs.empty()) {
            _logger->debug("EXTENDED SEARCH (hub-only): no hubs found in Top {:.1f}% with indeg >= {}", top_k_ratio * 100.0, min_hub_degree);
            break; 
        }

        _logger->debug("EXTENDED SEARCH (hub-only): {} hubs selected from Top {:.1f}% (min_deg={})",
                       hubs.size(), top_k_ratio * 100.0, min_hub_degree);

        bool local_improved = false;

        for (const auto& hub : hubs) {
            int best_y = hub.y;

            // ---- snapshot copy ----
            statistics["phase1_trials"]++; 
            Search Search_copy(*this);
            UnblockedPathsMap upm_copy; 
            const double baseline_score = total_score;

            std::vector<Insert>  cand_ins;
            std::vector<Reverse> cand_rev;
            std::vector<Delete>  cand_del;

            FlatSet comp = undirected_component_of(Search_copy.pdag, best_y);
            FlatSet external_parents;
            
            for (int v : comp) {
                for (int p : Search_copy.pdag.get_parents(v)) {
                    if (!comp.contains(p)) external_parents.insert(p);
                }
            }

            // batch edge-only delete
            EdgeModificationsMap em_batch;
            EdgeQueueSet edges_to_check_acc;
            size_t removed_cnt = 0;

            for (int x : external_parents) {
                for (int v : comp) {
                    if (!Search_copy.pdag.has_directed_edge(x, v)) continue;
                    Search_copy.pdag.remove_edge_only_no_finalize(
                        x, v, /*directed_x_to_y=*/true, em_batch, edges_to_check_acc);

                    removed_cnt++;
                    _logger->debug(
                        "  reset-step delete (hub-only): {} -> {} (hub y={})",
                        x, v, best_y
                    );
                }
            }
            if (removed_cnt == 0) {
                continue;
            }

            //Global CPDAG closure & Score & XGES-0
            Search_copy.pdag.complete_cpdag_global(em_batch);
            Search_copy.total_score = Search_copy.recompute_total_score_from_pdag();
            statistics["research_calls_phase1"]++;
            Search_copy.run_research(upm_copy,algo_type);

            //  Accept / Reject
            double gain = Search_copy.total_score - baseline_score;
            bool structure_changed = !(pdag == Search_copy.pdag);

            if (gain > eps_gain && structure_changed) {
                _logger->debug("EXTENDED SEARCH (hub-only) ACCEPTED: y={}, gain={}", best_y, gain);

                total_score = Search_copy.total_score;
                pdag = std::move(Search_copy.pdag);
                unblocked_paths_map = std::move(upm_copy);

                statistics["extended_search-accepted_phase1"] += 1;
                
                local_improved = true;
                global_improved = true; 
                break; 
            } else {
                _logger->debug("EXTENDED SEARCH (hub-only) REJECTED: y={}, gain={}", best_y, gain);
                statistics["extended_search-rejected_phase1"] += 1;
            }
        } 
        
    }
}


// Node-level Incoming Reset
void Search::node_incoming_reset_batch(UnblockedPathsMap &unblocked_paths_map, 
                                     size_t hub_threshold, double top_k_ratio, int algo_type) {
    constexpr double eps_gain = 1e-7;
    bool global_improved = true;
    struct HubCandidate {
        size_t indeg;
        size_t undirected_deg;
        int y;
    };

    while (global_improved) {
        global_improved = false;

        const auto& all_nodes = pdag.get_nodes_variables();
        std::vector<HubCandidate> all_candidates;
        all_candidates.reserve(all_nodes.size());

        for (int y : all_nodes) {
            size_t indeg = pdag.get_parents(y).size();
            size_t un_deg = pdag.get_neighbors(y).size();
            all_candidates.push_back({indeg, un_deg, y});
        }

        std::sort(all_candidates.begin(), all_candidates.end(),
                  [](const HubCandidate& a, const HubCandidate& b) {
                      size_t deg_a = a.indeg + a.undirected_deg;
                      size_t deg_b = b.indeg + b.undirected_deg;
                      if (deg_a != deg_b) return deg_a > deg_b;
                      return a.indeg > b.indeg;
                  });

        // Top K filtering
        size_t top_k_limit = (size_t)(all_nodes.size() * top_k_ratio);
        if (all_nodes.size() > 0 && top_k_limit == 0) top_k_limit = 1;

        std::vector<HubCandidate> hubs;
        hubs.reserve(top_k_limit);

        for (const auto& cand : all_candidates) {
            if (hubs.size() >= top_k_limit) break;
            if (cand.indeg < hub_threshold) break; // Threshold check
            hubs.push_back(cand);
        }

        if (hubs.empty()) {
            // _logger->debug("[Phase 2] No hubs found with threshold={}, ratio={:.2f}", hub_threshold, top_k_ratio);
            break; 
        }

        for (const auto& hub : hubs) {
            int y = hub.y;
            const auto& neighbors = pdag.get_neighbors(y);
            if(pdag.is_clique(neighbors)==false) {
                _logger->debug("[Phase 2] Rejected: Invalid CPDAG structure at y={}", y);
                continue;
            }

            const auto& parents = pdag.get_parents(y);
            if (parents.empty()) continue;

            // Snapshot Copy 
            statistics["phase2_trials"]++;
            Search Search_copy(*this);
            UnblockedPathsMap upm_copy;
            
            EdgeModificationsMap em_batch;
            EdgeQueueSet edges_to_check;
            
            // Batch Delete
            for (int x : parents) {
                Search_copy.pdag.remove_edge_only_no_finalize(x, y, true, em_batch, edges_to_check);
            }
            // Graph Repair (Meek Rule) & Validate
            Search_copy.pdag.complete_cpdag_global(em_batch);
            

            // Score Re-calculation
            Search_copy.total_score = Search_copy.recompute_total_score_from_pdag();

            // Re-search (Full Scan)
            std::vector<Insert> cand_ins;
            std::vector<Reverse> cand_rev;
            std::vector<Delete> cand_del;
            
            statistics["research_calls_phase2"]++;
            Search_copy.run_research(upm_copy,algo_type);

            // Accept Check
            double gain = Search_copy.total_score - total_score;
            bool structure_changed = !(pdag == Search_copy.pdag);

            if (gain > eps_gain && structure_changed) {
                _logger->debug("[Phase 2] Accepted: Node y={} (Ranked Hub), Cut Count={}, Gain={}", 
                               y, parents.size(), gain);
                
                total_score = Search_copy.total_score;
                pdag = std::move(Search_copy.pdag);
                unblocked_paths_map = std::move(upm_copy);
                
                statistics["extended_search-accepted_phase2"] += 1;
                global_improved = true;
                break;
            }
            else{
                statistics["extended_search-rejected_phase2"] += 1;
            }
        }
    }
}

void Search::single_edge_delete_best_first(
    UnblockedPathsMap &unblocked_paths_map,
    int algo_type
) {
    constexpr double eps_gain = 1e-7;

    struct EdgeCand {
        int x;
        int y;
        double delta;   // heap key (bigger first)
        int pa_sz;      // tie-break
        int ne_sz;      // tie-break

        // max-heap by (delta, pa_sz, ne_sz)
        bool operator<(const EdgeCand& other) const {
            if (delta != other.delta) return delta < other.delta;
            if (pa_sz != other.pa_sz) return pa_sz < other.pa_sz;
            return ne_sz < other.ne_sz;
        }
    };

    // tabu: directed edge only, so key=(x,y)
    std::set<std::pair<int,int>> rejected_edges;

    bool global_improved = true;
    while (global_improved) {
        global_improved = false;

        // Build candidate heap (DIRECTED edges only)
        std::vector<EdgeCand> heap;
        heap.reserve(pdag.get_number_of_edges() + 16);

        const auto& d_edges = pdag.get_directed_edges();

        for (auto [x, y] : d_edges) {
            if (rejected_edges.count({x, y})) continue;
            const auto& y_neighbors = pdag.get_neighbors(y);
            const auto& x_adjacent  = pdag.get_adjacent(x);

            FlatSet common;
            common.reserve(std::min(y_neighbors.size(), x_adjacent.size()));
            std::set_intersection(
                y_neighbors.begin(), y_neighbors.end(),
                x_adjacent.begin(), x_adjacent.end(),
                std::inserter(common, common.begin())
            );
            if (!pdag.is_clique(common)) {
                rejected_edges.insert({x, y});
                continue;
            }

            const auto& parents_y = pdag.get_parents(y);
            if (!parents_y.contains(x)) {
                continue;
            }

            FlatSet new_parents = parents_y;
            new_parents.erase(x);

            // delta = S(y, newPa) - S(y, Pa)
            const double delta =
                scorer->local_score(y, new_parents) - scorer->local_score(y, parents_y);

            const int pa_sz = (int)parents_y.size();
            const int ne_sz = (int)y_neighbors.size();

            heap.push_back(EdgeCand{x, y, delta, pa_sz, ne_sz});
        }

        if (heap.empty()) break;
        std::make_heap(heap.begin(), heap.end());

        // best-first trial
        while (!heap.empty()) {
            std::pop_heap(heap.begin(), heap.end());
            EdgeCand cand = heap.back();
            heap.pop_back();

            if (rejected_edges.count({cand.x, cand.y})) continue;

            // Snapshot Copy
            statistics["phase3_trials"]++;
            Search Search_copy(*this);

            EdgeModificationsMap em;
            EdgeQueueSet edges_to_check;

            // Delete (DIRECTED only, H=empty)
            Search_copy.pdag.remove_edge_only_no_finalize(
                cand.x, cand.y,
                /*directed_x_to_y=*/true,
                em, edges_to_check
            );

            // Repair (global closure)
            Search_copy.pdag.complete_cpdag_global(em);

            // Forbidden insert
            Search_copy.pdag.add_forbidden_insert(cand.x, cand.y);

            // Score recompute
            Search_copy.total_score = Search_copy.recompute_total_score_from_pdag();

            //  Research
            UnblockedPathsMap upm_copy = unblocked_paths_map;
            statistics["research_calls_phase3"]++;
            Search_copy.run_research(upm_copy, algo_type);

            //  Accept / Reject
            const double gain = Search_copy.total_score - total_score;
            const bool changed = !(pdag == Search_copy.pdag);

            if (gain > eps_gain && changed) {
                _logger->debug(
                    "[Phase 3 - Clean Heap Δ] Accepted {}->{}, Δ_local={}, gain={}",
                    cand.x, cand.y, cand.delta, gain
                );

                total_score = Search_copy.total_score;
                pdag = std::move(Search_copy.pdag);
                unblocked_paths_map = std::move(upm_copy);

                rejected_edges.clear();
                global_improved = true;
                statistics["extended_search-accepted_phase3"] += 1;
                break; // greedy restart
            } else {
                statistics["extended_search-rejected_phase3"] += 1;
                rejected_edges.insert({cand.x, cand.y});
            }
        }
    }
}

// Single Edge Deletion
void Search::single_edge_delete_clean_only(UnblockedPathsMap &unblocked_paths_map, int algo_type) {
    constexpr double eps_gain = 1e-7;
    std::set<std::pair<int, int>> rejected_edges; 

    bool global_improved = true;
    while (global_improved) {
        global_improved = false;
        
        auto edges = pdag.get_directed_edges();
        auto u_edges = pdag.get_undirected_edges();
        edges.insert(edges.end(), u_edges.begin(), u_edges.end());
        
        for (auto [x, y] : edges) {
            // Tabu
            if (rejected_edges.count({x, y})) continue;

            const auto& y_neighbors = pdag.get_neighbors(y);
            const auto& x_adjacent= pdag.get_adjacent(x);
            FlatSet common;
            common.reserve(std::min(y_neighbors.size(), x_adjacent.size()));
            std::set_intersection(
                y_neighbors.begin(), y_neighbors.end(),
                x_adjacent.begin(), x_adjacent.end(),
                std::inserter(common, common.begin())  
            );
            if(pdag.is_clique(common)==false){
                _logger->debug("[Phase 3] Rejected: Invalid Graph at {}->{}", x, y);
                rejected_edges.insert({x, y});
                continue;
            }
            
            // Edge Type Check
            bool is_directed = pdag.has_directed_edge(x, y);

            // Snapshot Copy
            statistics["phase3_trials"]++;
            Search Search_copy(*this);
            EdgeModificationsMap em;
            EdgeQueueSet edges_to_check; 

            // Delete Operation
            Search_copy.pdag.remove_edge_only_no_finalize(x, y, is_directed, em, edges_to_check);
            
            // Graph Repair (Meek Rule)
            Search_copy.pdag.complete_cpdag_global(em);

            // Forbidden
            Search_copy.pdag.add_forbidden_insert(x, y);
            if (!is_directed) {
                 Search_copy.pdag.add_forbidden_insert(y, x);
            }

            // Score Re-calculation
            Search_copy.total_score = Search_copy.recompute_total_score_from_pdag();

            // Re-search
            UnblockedPathsMap upm_copy = unblocked_paths_map;
            std::vector<Insert> cand_ins;
            std::vector<Reverse> cand_rev;
            std::vector<Delete> cand_del;
            
            statistics["research_calls_phase3"]++;
            Search_copy.run_research(upm_copy,algo_type);

            // Accept / Reject
            double gain = Search_copy.total_score - total_score;
            bool structure_changed = !(pdag == Search_copy.pdag);

            if (gain > eps_gain && structure_changed) {
                 _logger->debug("[Phase 3] Accepted: {}->{}, Gain={}", x, y, gain);
                 
                 total_score = Search_copy.total_score;
                 pdag = std::move(Search_copy.pdag);
                 unblocked_paths_map = std::move(upm_copy);
                 
                 rejected_edges.clear();
                 global_improved = true;
                 statistics["extended_search-accepted_phase3"] += 1;
                 break; // Greedy restart
            } else {
                statistics["extended_search-rejected_phase3"] += 1;
                rejected_edges.insert({x, y});
            }
        }
    }
}



void Search::fit_xges_variant2() {
    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map;

    candidate_inserts.reserve(100 * n_variables);
    candidate_reverses.reserve(n_variables);
    candidate_deletes.reserve(n_variables);

    _logger->info("=== [Variant 2] Hub-based Search with FULL Coverage (Ratio 1.0) ===");

    heuristic_xges0(candidate_inserts, candidate_reverses, candidate_deletes, unblocked_paths_map, true);

    _logger->info("[Variant 2] Running Hub-based Component Reset checking ALL nodes (Ratio 1.0)...");
    
    hub_only_component_incoming_reset_and_research_edge_only(unblocked_paths_map, 1,1.0,0);

    _logger->info("=== [Variant 2] Finished. Final Score: {} ===", total_score);
}

void Search::fit_xges_variant3() {
    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map;

    candidate_inserts.reserve(100 * n_variables);
    candidate_reverses.reserve(n_variables);
    candidate_deletes.reserve(n_variables);

    _logger->info("=== [Variant 3] Scale-down Extended Search (Hub -> Node -> Edge) ===");

    // -------------------------------------------------------------------
    // Phase 0: Base XGES-0
    // -------------------------------------------------------------------
    _logger->info("[Variant 3] Phase 0: Base XGES-0");
    heuristic_xges0(candidate_inserts, candidate_reverses, candidate_deletes, unblocked_paths_map, true);


    // -------------------------------------------------------------------
    // Phase 1: Hub Component Reset (Coarse-grained)
    // -------------------------------------------------------------------
    _logger->info("[Variant 3] Phase 1: Hub Component Reset");
    double score_p1 = total_score;
    // Parameters: threshold=3, ratio=0.1 (Top 10% hubs)
    hub_only_component_incoming_reset_and_research_edge_only(unblocked_paths_map,2,1.0,0);
    
    _logger->info(" > Phase 1 Gain: {}", total_score - score_p1);


    // -------------------------------------------------------------------
    // Phase 2: Node Incoming Reset (Medium-grained)
    // Strategy: Batch delete incoming edges -> Repair(Meek) -> Check Valid(CPDAG) -> Rescore
    // -------------------------------------------------------------------
    _logger->info("[Variant 3] Phase 2: Node Incoming Reset Batch (Try & Validate)");
    double score_p2 = total_score;
    
    node_incoming_reset_batch(unblocked_paths_map,2,1.0,0);
    
    _logger->info(" > Phase 2 Gain: {}", total_score - score_p2);


    // -------------------------------------------------------------------
    // Phase 3: Single Edge Deletion (Fine-grained)
    // Strategy: Single delete -> Repair(Meek) -> Check Valid(CPDAG) -> Rescore
    // -------------------------------------------------------------------
    _logger->info("[Variant 3] Phase 3: Single Edge Deletion (Try & Validate)");
    double score_p3 = total_score;
    
    single_edge_delete_clean_only(unblocked_paths_map,0);
    
    _logger->info(" > Phase 3 Gain: {}", total_score - score_p3);
}


void Search::fit_xges_variant4() {
    _logger->info("=== [XGES Variant 2] Base XGES-0 + per-node DeletePa(y, H) (single-node baseline) ===");

    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map;

    candidate_inserts.reserve(100 * n_variables);
    candidate_reverses.reserve(n_variables);
    candidate_deletes.reserve(n_variables);

    // 1. Run base XGES-0
    _logger->info("[XGES Variant 2] Phase 0: Base XGES-0");
    heuristic_xges0(candidate_inserts, candidate_reverses, candidate_deletes, unblocked_paths_map, true);

    // 2. Single-node DeletePa(y, H) loop (NO component grouping).
    //    algo_type=0 matches the post-DeletePa research path used for XGES.
    _logger->info("[XGES Variant 2] Phase 1: Per-node DeletePa(y, H) with no component grouping");
    double score_before = total_score;
    node_incoming_reset_batch_ablation_deletepa(
        unblocked_paths_map,
        /*hub_threshold=*/1,
        /*top_k_ratio=*/1.0,
        /*algo_type=*/0,            // XGES dispatcher
        /*max_enum_neighbors=*/30   // safety cap; throws for hubs with |Ne(y)|>10
    );
    _logger->info(" > Variant 2 Gain: {}", total_score - score_before);

    // Logging equivalent to GES variant4
    _logger->info("DeletePa trials: {}", statistics["deletepa_trials"]);
    _logger->info("research_calls_phase1: {}", statistics["research_calls_phase1"]);
    _logger->info("research_calls_phase2: {}", statistics["research_calls_phase2"]);
    _logger->info("research_calls_phase3: {}", statistics["research_calls_phase3"]);

    _logger->info("extended_search-accepted_phase1: {}", statistics["extended_search-accepted_phase1"]);
    _logger->info("extended_search-accepted_phase2: {}", statistics["extended_search-accepted_phase2"]);
    _logger->info("extended_search-accepted_phase3: {}", statistics["extended_search-accepted_phase3"]);

    _logger->info("extended_search-rejected_phase1: {}", statistics["extended_search-rejected_phase1"]);
    _logger->info("extended_search-rejected_phase2: {}", statistics["extended_search-rejected_phase2"]);
    _logger->info("extended_search-rejected_phase3: {}", statistics["extended_search-rejected_phase3"]);

    _logger->info("=== [XGES Variant 2] Finished. Final Score: {} ===", total_score);
}




std::string edge_type_to_string(EdgeType type) {
    switch (type) {
        case EdgeType::DIRECTED_TO_Y: return "->";
        case EdgeType::DIRECTED_TO_X: return "<-";
        case EdgeType::UNDIRECTED:    return "-";
        case EdgeType::NONE:          return "x"; 
        default:                      return "?";
    }
}



/** @brief Extended search of XGES.
 *
 * Extended search is performed after XGES-0 has exhausted all possible operators.
 * It successively delete some edges, and resume the search with XGES-0 (preventing
 * the deleted edges from being inserted again). If the score increases, the new PDAG
 * is kept.
 *
 * At the beggining of the search, `all_edge_deletes` contains all the possible deletes
 * for the current PDAG. If the PDAG is updated, we keep iterating through `all_edge_deletes`
 * until all the deletes have been tried. Only then we recompute `all_edge_deletes` on the
 * last updated PDAG. Until no new deletes are found. This is more efficient than recomputing
 * all the deletes at each PDAG update, and re-testing all of them (as most of them will be
 * the same).
 *
 *
 * @param unblocked_paths_map The unblocked paths map used by XGES-0.
 */

// [Modified] Standard Extended Search with Algorithm Dispatcher
// algo_type: 0=XGES, 1=GES, 2=OPS, 3=BOSS, 4=LGES-Safe, 5=LGES-Conservative
void Search::block_each_edge_and_research(UnblockedPathsMap &unblocked_paths_map, int algo_type) {
    std::vector<Delete> all_edge_deletes;
    bool deletes_of_pdag_are_updated = false;

    // (x, y, directed_bool, frozenset(C))
    std::set<std::tuple<int, int, bool, FlatSet>> rejected_deletes; 

    while (!all_edge_deletes.empty() || !deletes_of_pdag_are_updated) {
        if (all_edge_deletes.empty()) {
            
            for (const int y: pdag.get_nodes_variables()) {
                find_deletes_to_y(y, all_edge_deletes, false);
            }

            std::vector<Delete> filtered;
            filtered.reserve(all_edge_deletes.size());
            for (const auto &del : all_edge_deletes) {
                std::tuple<int, int, bool, FlatSet> key =
                    std::make_tuple(del.x, del.y, del.directed, del.C);
                if (rejected_deletes.find(key) == rejected_deletes.end()) { 
                    filtered.push_back(del);
                }
            }
            all_edge_deletes = std::move(filtered);

            if (all_edge_deletes.empty()) { break; }
            deletes_of_pdag_are_updated = true;
        }

        std::pop_heap(all_edge_deletes.begin(), all_edge_deletes.end());
        auto delete_ = std::move(all_edge_deletes.back());
        all_edge_deletes.pop_back();

        if (!pdag.is_delete_valid(delete_)) { continue; } 

        // Apply the delete
        statistics["phase3_trials"]++;
        Search Search_copy(*this);
        EdgeModificationsMap edge_modifications;
        
        Search_copy.pdag.apply_delete(delete_, edge_modifications);
        Search_copy.total_score += delete_.score;
        Search_copy.pdag.add_forbidden_insert(delete_.x, delete_.y);

        if (!delete_.directed) {
            Search_copy.pdag.add_forbidden_insert(delete_.y, delete_.x);
        }

        UnblockedPathsMap blocked_paths_map_copy = unblocked_paths_map;
        _logger->debug("EXTENDED SEARCH: {}", fmt::streamed(delete_));
        for (auto &[fst, snd]: edge_modifications) { _logger->trace("\tEdge {}", fmt::streamed(snd)); }

        
        if (algo_type == 0) { // XGES
             std::vector<Insert> candidate_inserts;
             std::vector<Reverse> candidate_reverses;
             std::vector<Delete> candidate_deletes;
             
             // Efficient Update (XGES only optimization)
             Search_copy.update_operator_candidates_efficient(
                    edge_modifications, candidate_inserts, candidate_reverses,
                    candidate_deletes, blocked_paths_map_copy);
             Search_copy.heuristic_xges0(candidate_inserts, candidate_reverses,
                                      candidate_deletes, blocked_paths_map_copy, false);
        } 
        else { // GES, OPS, BOSS, LGES-Safe, LGES-Conservative -> General Research
             Search_copy.run_research(blocked_paths_map_copy,algo_type);
        }

        // =========================================================================

        if (pdag == Search_copy.pdag) { 
            _logger->debug("No changes by EXTENDED SEARCH");
            statistics["extended_search-rejected_phase3"] += 1;
            continue; 
        }
        if (Search_copy.total_score - total_score > 1e-7) {
            _logger->debug("EXTENDED SEARCH ACCEPTED: with increase {} and {}",
                           Search_copy.total_score - total_score, fmt::streamed(delete_));
            total_score = Search_copy.total_score;
            pdag = Search_copy.pdag;
            unblocked_paths_map = std::move(blocked_paths_map_copy);
            deletes_of_pdag_are_updated = false;

            rejected_deletes.clear(); 
            statistics["extended_search-accepted_phase3"] += 1;
        } else {
            _logger->debug("EXTENDED SEARCH REJECTED: {} {}", fmt::streamed(delete_),
                           Search_copy.total_score);
            statistics["extended_search-rejected_phase3"] += 1;

            std::tuple<int, int, bool, FlatSet> key =
                std::make_tuple(delete_.x, delete_.y, delete_.directed, delete_.C);
            rejected_deletes.insert(key);
        }
    }
}

/** @brief Run the XGES-0 heuristic.
 *
 * XGES-0 is the main heuristic of the XGES algorithm. It maintains lists of all possible
 * operators and apply them in order: delete, reverse, and insert. XGES-0 stops when no
 * more operators can be applied.
 * XGES-0 uses update_operator_candidates_efficient to update the operators after each
 * PDAG update.
 *
 * @param candidate_inserts The vector of candidate inserts.
 * @param candidate_reverses The vector of candidate reverses.
 * @param candidate_deletes The vector of candidate deletes.
 * @param unblocked_paths_map The unblocked paths map used by XGES-0.
 * @param initialize_inserts If true, the candidate inserts are initialized at the beginning
 * of the heuristic.
 */
void Search::heuristic_xges0(std::vector<Insert> &candidate_inserts,
                           std::vector<Reverse> &candidate_reverses,
                           std::vector<Delete> &candidate_deletes,
                           UnblockedPathsMap &unblocked_paths_map,
                           bool initialize_inserts) {
    if (initialize_inserts) {
        // find all possible inserts
        auto start_init_inserts = high_resolution_clock::now();
        for (int y = 0; y < n_variables; ++y) {
            find_inserts_to_y(y, candidate_inserts, -1, true);
            find_deletes_to_y(y, candidate_deletes, true);
            find_reverse_to_y(y, candidate_reverses);
        }
        statistics["initialize_inserts-time"] += measure_time(start_init_inserts);
    }
    EdgeModificationsMap edge_modifications;
    int i_operations = 1;

    Insert last_insert(-1, -1, FlatSet{}, -1, FlatSet{});

    // XGES-0 main loop, in order: delete, reverse, insert; one operator per iteration
    while (!candidate_inserts.empty() || !candidate_reverses.empty() || !candidate_deletes.empty()) {
        edge_modifications.clear();

        if (!candidate_deletes.empty()) {
            // apply the best delete if possible
            std::pop_heap(candidate_deletes.begin(), candidate_deletes.end());
            auto best_delete = std::move(candidate_deletes.back());
            candidate_deletes.pop_back();
            if (pdag.is_delete_valid(best_delete)) {
                pdag.apply_delete(best_delete, edge_modifications);
                total_score += best_delete.score;
                _logger->debug("{}: {}, current_score: {}", i_operations, fmt::streamed(best_delete), total_score);
            } else {
                continue;
            }

        } else if (!candidate_reverses.empty()) {
            // apply the best reverse if possible (no delete available)
            std::pop_heap(candidate_reverses.begin(), candidate_reverses.end());
            auto best_reverse = std::move(candidate_reverses.back());
            candidate_reverses.pop_back();
            if (pdag.is_reverse_valid(best_reverse, unblocked_paths_map)) {
                pdag.apply_reverse(best_reverse, edge_modifications);
                total_score += best_reverse.score;
                _logger->debug("{}: {}, current_score: {}", i_operations, fmt::streamed(best_reverse), total_score);
            } else {
                continue;
            }

        } else if (!candidate_inserts.empty()) {
            // apply the best insert if possible (no delete or reverse available)
            std::pop_heap(candidate_inserts.begin(), candidate_inserts.end());
            auto best_insert = std::move(candidate_inserts.back());
            candidate_inserts.pop_back();
            if (best_insert.y == last_insert.y &&
                abs(best_insert.score - last_insert.score) < 1e-10 &&
                best_insert.x == last_insert.x && best_insert.T == last_insert.T) {
                statistics["probable_insert_duplicates"] += 1;
                continue;
            }
            last_insert = std::move(best_insert);
            if (pdag.is_insert_valid(last_insert, unblocked_paths_map)) {
                pdag.apply_insert(last_insert, edge_modifications);
                total_score += last_insert.score;
                _logger->debug("{}: {}, current_score: {}", i_operations, fmt::streamed(last_insert), total_score);
            } else {
                continue;
            }
        }
        // if we reach this point, we have applied an operator
        i_operations++;
        for (auto &edge_modification: edge_modifications) {
            _logger->trace("\tEdge {}", fmt::streamed(edge_modification.second));
        }

        // update the new possible operators
        // update_operator_candidates_naive(candidate_inserts, candidate_reverses,
        //                                  candidate_deletes);
        update_operator_candidates_efficient(edge_modifications, candidate_inserts,
                                             candidate_reverses, candidate_deletes,
                                             unblocked_paths_map);
    }
}



void Search::fit_lges(bool score_based, bool prune) {
    std::vector<Insert> candidate_inserts; // empty at the beginning
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map;

    candidate_inserts.reserve(100 * n_variables); // empty at the beginning
    candidate_inserts.reserve(n_variables); // empty at the beginning       
    candidate_deletes.reserve(n_variables); // empty at the beginning
    //_logger->info("Safe Insert: {}",score_based);
    heuristic_lges(candidate_inserts, candidate_reverses, candidate_deletes,
                    unblocked_paths_map, true, score_based, prune);
}

void Search::fit_lges_variant1(bool score_based, bool prune) {
    fit_lges(score_based, prune);
    // Priming Step (Candidate Initialization ONLY)
    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map; // Start empty

    candidate_inserts.reserve(100 * n_variables);
    candidate_reverses.reserve(n_variables);
    candidate_deletes.reserve(n_variables);

    _logger->debug("[LGES Variant 1] Initializing candidates (Pure Mode)...");
    
    initialize_candidates_pure(candidate_inserts, candidate_reverses, candidate_deletes);

    // Standard Extended Search
    _logger->info("[LGES Variant 1] Starting Standard Extended Search...");
    double score_before = total_score;
    if(score_based&&!prune){
        block_each_edge_and_research(unblocked_paths_map,4);
    }
    else if(score_based&&prune){
        block_each_edge_and_research(unblocked_paths_map,5);
    }
    _logger->info("=== [LGES Variant 1] Finished. Gain: {} ===", total_score - score_before);
}



void Search::fit_lges_variant2(bool score_based,bool prune) {
    _logger->info("=== [LGES Variant 2] Running LGES + hub node Extended Search===");

    // Run GES
    fit_lges(score_based,prune);


    // Priming Step (Candidate Initialization ONLY)
    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map; // Start empty

    candidate_inserts.reserve(100 * n_variables);
    candidate_reverses.reserve(n_variables);
    candidate_deletes.reserve(n_variables);

    _logger->debug("[LGES Variant 2] Initializing candidates (Pure Mode)...");
    
    initialize_candidates_pure(candidate_inserts, candidate_reverses, candidate_deletes);

    // Phase 1: Hub-based Extended Search
    _logger->info("[LGES Variant 2] Phase 1: Running Hub-based Component Incoming Reset...");
    double score_before_hub = total_score;
    
    if(score_based&&!prune){
        hub_only_component_incoming_reset_and_research_edge_only(unblocked_paths_map,1,1.0,4);
    }
    else if(score_based&&prune){
        hub_only_component_incoming_reset_and_research_edge_only(unblocked_paths_map,1,1.0,5);
    }
    
    _logger->info("[LGES Variant 2] Phase 1 Finished. Score Gain: {}", total_score - score_before_hub);
}

void Search::fit_lges_variant3(bool score_based, bool prune) {
    // -------------------------------------------------------------------
    // Phase 0: Base OPS
    // -------------------------------------------------------------------
    _logger->info("[Variant 3] Phase 0: Base LGES");
    // 1. Run GES
    fit_lges(score_based,prune);

    // 2. Priming Step (Candidate Initialization ONLY)
    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map; // Start empty

    candidate_inserts.reserve(100 * n_variables);
    candidate_reverses.reserve(n_variables);
    candidate_deletes.reserve(n_variables);

    // -------------------------------------------------------------------
    // Phase 1: Hub Component Reset (Coarse-grained)
    // -------------------------------------------------------------------
    _logger->info("[Variant 3] Phase 1: Hub Component Reset");
    double score_p1 = total_score;
    // Parameters: threshold=3, ratio=0.1 (Top 10% hubs)
    if(score_based&&!prune){
        hub_only_component_incoming_reset_and_research_edge_only(unblocked_paths_map,2,1.0,4);
    }
    else if(score_based&&prune){
        hub_only_component_incoming_reset_and_research_edge_only(unblocked_paths_map,2,1.0,5);
    }
    
    _logger->info(" > Phase 1 Gain: {}", total_score - score_p1);


    // -------------------------------------------------------------------
    // Phase 2: Node Incoming Reset (Medium-grained)
    // Strategy: Batch delete incoming edges -> Repair(Meek) -> Check Valid(CPDAG) -> Rescore
    // -------------------------------------------------------------------
    _logger->info("[Variant 3] Phase 2: Node Incoming Reset Batch (Try & Validate)");
    double score_p2 = total_score;
    
    if(score_based&&!prune){
        node_incoming_reset_batch(unblocked_paths_map,2,1.0,4);
    }
    else if(score_based&&prune){
        node_incoming_reset_batch(unblocked_paths_map,2,1.0,5);
    }
    
    _logger->info(" > Phase 2 Gain: {}", total_score - score_p2);


    // -------------------------------------------------------------------
    // Phase 3: Single Edge Deletion (Fine-grained)
    // Strategy: Single delete -> Repair(Meek) -> Check Valid(CPDAG) -> Rescore
    // -------------------------------------------------------------------
    _logger->info("[Variant 3] Phase 3: Single Edge Deletion (Try & Validate)");
    double score_p3 = total_score;
    
    if(score_based&&!prune){
        single_edge_delete_clean_only(unblocked_paths_map,4);
    }
    else if(score_based&&prune){
        single_edge_delete_clean_only(unblocked_paths_map,5);
    }
    
    _logger->info(" > Phase 3 Gain: {}", total_score - score_p3);
    _logger->info("=== [Variant 3] Finished. Final Score: {} ===", total_score);
}


void Search::fit_lges_variant4(bool score_based, bool prune) {
    _logger->info("=== [LGES Variant 2] Base LGES + per-node DeletePa(y, H) (single-node baseline) ===");

    // 1. Run base LGES
    _logger->info("[LGES Variant 2] Phase 0: Base LGES");
    fit_lges(score_based, prune);

    // 2. Priming Step (Candidate Initialization ONLY)
    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map; // Start empty

    candidate_inserts.reserve(100 * n_variables);
    candidate_reverses.reserve(n_variables);
    candidate_deletes.reserve(n_variables);

    _logger->debug("[LGES Variant 2] Initializing candidates (Pure Mode)...");
    initialize_candidates_pure(candidate_inserts, candidate_reverses, candidate_deletes);

    // 3. Single-node DeletePa(y, H) loop (NO component grouping).
    //    Determine dispatcher ID based on LGES mode
    int algo_type = 1; // fallback to GES if not score_based
    if (score_based && !prune) {
        algo_type = 4; // LGES-Safe
    } else if (score_based && prune) {
        algo_type = 5; // LGES-Conservative
    }

    _logger->info("[LGES Variant 2] Phase 1: Per-node DeletePa(y, H) with no component grouping");
    double score_before = total_score;
    node_incoming_reset_batch_ablation_deletepa(
        unblocked_paths_map,
        /*hub_threshold=*/1,
        /*top_k_ratio=*/1.0,
        /*algo_type=*/algo_type,
        /*max_enum_neighbors=*/30   // safety cap; throws for hubs with |Ne(y)|>10
    );
    _logger->info(" > Variant 2 Gain: {}", total_score - score_before);

    // Logging equivalent to other variant4s
    _logger->info("DeletePa trials: {}", statistics["deletepa_trials"]);
    _logger->info("research_calls_phase1: {}", statistics["research_calls_phase1"]);
    _logger->info("research_calls_phase2: {}", statistics["research_calls_phase2"]);
    _logger->info("research_calls_phase3: {}", statistics["research_calls_phase3"]);

    _logger->info("extended_search-accepted_phase1: {}", statistics["extended_search-accepted_phase1"]);
    _logger->info("extended_search-accepted_phase2: {}", statistics["extended_search-accepted_phase2"]);
    _logger->info("extended_search-accepted_phase3: {}", statistics["extended_search-accepted_phase3"]);

    _logger->info("extended_search-rejected_phase1: {}", statistics["extended_search-rejected_phase1"]);
    _logger->info("extended_search-rejected_phase2: {}", statistics["extended_search-rejected_phase2"]);
    _logger->info("extended_search-rejected_phase3: {}", statistics["extended_search-rejected_phase3"]);

    _logger->info("=== [LGES Variant 2] Finished. Final Score: {} ===", total_score);
}


void Search::heuristic_lges(std::vector<Insert> &candidate_inserts,
                           std::vector<Reverse> &candidate_reverses,
                           std::vector<Delete> &candidate_deletes,
                           UnblockedPathsMap &unblocked_paths_map,
                           bool initialize_inserts,
                            bool score_based,
                            bool prune) {
    if (initialize_inserts) {
        // find all possible inserts
        auto start_init_inserts = high_resolution_clock::now();
        for (int y = 0; y < n_variables; ++y) {
            find_inserts_to_y(y, candidate_inserts, -1, true);
            //also check delete, reversal candiates
            find_deletes_to_y(y, candidate_deletes, true);
            find_reverse_to_y(y, candidate_reverses);
        }
        statistics["initialize_inserts-time"] += measure_time(start_init_inserts);
    }
    EdgeModificationsMap edge_modifications;
    int i_operations = 1;

    Insert last_insert(-1, -1, FlatSet{}, -1, FlatSet{});

    // XGES-0 main loop, in order: delete, reverse, insert; one operator per iteration
    while (!candidate_inserts.empty() || !candidate_reverses.empty() || !candidate_deletes.empty()) {
        edge_modifications.clear();

        if (!candidate_deletes.empty()) {
            // apply the best delete if possible
            std::pop_heap(candidate_deletes.begin(), candidate_deletes.end());
            auto best_delete = std::move(candidate_deletes.back());
            candidate_deletes.pop_back();
            if (pdag.is_delete_valid(best_delete)) {
                pdag.apply_delete(best_delete, edge_modifications);
                total_score += best_delete.score;
                _logger->debug("{}: {}, current_score: {}", i_operations, fmt::streamed(best_delete), total_score);
            } else {
                continue;
            }

        } else if (!candidate_reverses.empty()) {
            // apply the best reverse if possible (no delete available)
            std::pop_heap(candidate_reverses.begin(), candidate_reverses.end());
            auto best_reverse = std::move(candidate_reverses.back());
            candidate_reverses.pop_back();
            if (pdag.is_reverse_valid(best_reverse, unblocked_paths_map)) {
                pdag.apply_reverse(best_reverse, edge_modifications);
                total_score += best_reverse.score;
                _logger->debug("{}: {}, current_score: {}", i_operations, fmt::streamed(best_reverse), total_score);
            } else {
                continue;
            }

        } else if (!candidate_inserts.empty()) {
            // apply the best insert if possible (no delete or reverse available)
            std::pop_heap(candidate_inserts.begin(), candidate_inserts.end());
            auto best_insert = std::move(candidate_inserts.back());
            candidate_inserts.pop_back();
            if (best_insert.y == last_insert.y &&
                abs(best_insert.score - last_insert.score) < 1e-10 &&
                best_insert.x == last_insert.x && best_insert.T == last_insert.T) {
                statistics["probable_insert_duplicates"] += 1;
                continue;
            }
            last_insert = std::move(best_insert);
            if (pdag.is_insert_valid(last_insert, unblocked_paths_map)) {
                pdag.apply_insert(last_insert, edge_modifications);
                total_score += last_insert.score;
                _logger->debug("{}: {}, current_score: {}", i_operations, fmt::streamed(last_insert), total_score);
            } else {
                continue;
            }
        }
        // if we reach this point, we have applied an operator
        i_operations++;
        for (auto &edge_modification: edge_modifications) {
            _logger->trace("\tEdge {}", fmt::streamed(edge_modification.second));
        }

        std::unique_ptr<PDAG> sampled_DAG;
        PDAG* DAG_ptr=nullptr;
        //SafeInsert: DAG Sampling
        if(score_based){
            sampled_DAG=std::make_unique<PDAG>(pdag.get_dag_extension());
            DAG_ptr=sampled_DAG.get();
        }
        //_logger->info("DAG sampling: {}",DAG_ptr->check_is_cpdag());

        // update the new possible operators
         update_operator_candidates_naive_for_lges(candidate_inserts, candidate_reverses,
                                      candidate_deletes,DAG_ptr,score_based,prune);
        //update_operator_candidates_efficient_for_lges(edge_modifications, candidate_inserts,
          //                                   candidate_reverses, candidate_deletes,
            //                                 unblocked_paths_map,DAG_ptr,score_based,prune);
    }
}

void Search::fit_ops(bool use_reverse) {
    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map;

    candidate_inserts.reserve(100 * n_variables);
    candidate_inserts.reserve(n_variables);
    candidate_deletes.reserve(n_variables);

    for (int y = 0; y < n_variables; ++y) {
        find_inserts_to_y(y, candidate_inserts, -1, true);
    }
    EdgeModificationsMap edge_modifications;
    int i_operations = 1;

    Insert last_insert(-1, -1, FlatSet{}, -1, FlatSet{});

    while (!candidate_inserts.empty() || !candidate_reverses.empty() ||
           !candidate_deletes.empty()) {

        // Look for the operator with the highest score
        // Each of candidate_inserts, candidate_reverses, candidate_deletes are heaps
        // with the highest score at the front.
        double best_score = 0;
        int best_type = -1;// 0: insert, 1: reverse, 2: delete

        if (!candidate_inserts.empty()) {
            if (candidate_inserts.front().score > best_score) {
                best_type = 0;
                best_score = candidate_inserts.front().score;
            }
        }
        if (!candidate_reverses.empty() && use_reverse) {
            if (candidate_reverses.front().score > best_score) {
                best_type = 1;
                best_score = candidate_reverses.front().score;
            }
        }
        if (!candidate_deletes.empty()) {
            if (candidate_deletes.front().score > best_score) { best_type = 2; }
        }

        edge_modifications.clear();
        if (best_type == 0) {
            std::pop_heap(candidate_inserts.begin(), candidate_inserts.end());
            auto best_insert = std::move(candidate_inserts.back());
            candidate_inserts.pop_back();
            if (best_insert.y == last_insert.y &&
                abs(best_insert.score - last_insert.score) < 1e-10 &&
                best_insert.x == last_insert.x && best_insert.T == last_insert.T) {
                statistics["probable_insert_duplicates"] += 1;
                continue;
            }
            last_insert = std::move(best_insert);
            if (pdag.is_insert_valid(last_insert, unblocked_paths_map)) {
                pdag.apply_insert(last_insert, edge_modifications);
                total_score += last_insert.score;
                _logger->debug("{}: {}", i_operations, fmt::streamed(last_insert));
            } else {
                continue;
            }
        } else if (best_type == 1) {
            // apply the best reverse if possible (no delete available)
            std::pop_heap(candidate_reverses.begin(), candidate_reverses.end());
            auto best_reverse = std::move(candidate_reverses.back());
            candidate_reverses.pop_back();
            if (pdag.is_reverse_valid(best_reverse, unblocked_paths_map)) {
                pdag.apply_reverse(best_reverse, edge_modifications);
                total_score += best_reverse.score;
                _logger->debug("{}: {}", i_operations, fmt::streamed(best_reverse));
            } else {
                continue;
            }
        } else if (best_type == 2) {
            std::pop_heap(candidate_deletes.begin(), candidate_deletes.end());
            auto best_delete = std::move(candidate_deletes.back());
            candidate_deletes.pop_back();
            if (pdag.is_delete_valid(best_delete)) {
                pdag.apply_delete(best_delete, edge_modifications);
                total_score += best_delete.score;
                _logger->debug("{}: {}", i_operations, fmt::streamed(best_delete));
            } else {
                continue;
            }
        } else {
            continue;
        }


        i_operations++;
        for (auto &edge_modification: edge_modifications) {
            _logger->trace("\tEdge {}", fmt::streamed(edge_modification.second));
        }

        update_operator_candidates_naive(candidate_inserts, candidate_reverses,
                                         candidate_deletes);

        if (!use_reverse) { candidate_reverses.clear(); }
    }
}


void Search::fit_ops_variant1() {
    _logger->info("=== [OPS Variant 1] Running OPS + Standard Extended Search ===");

    fit_ops(false);
    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map; // Start empty

    candidate_inserts.reserve(100 * n_variables);
    candidate_reverses.reserve(n_variables);
    candidate_deletes.reserve(n_variables);

    _logger->debug("[OPS Variant 1] Initializing candidates (Pure Mode)...");
    
    initialize_candidates_pure(candidate_inserts, candidate_reverses, candidate_deletes);

    _logger->info("[OPS Variant 1] Starting Standard Extended Search...");
    double score_before = total_score;
    block_each_edge_and_research(unblocked_paths_map,2);
    _logger->info("=== [OPS Variant 1] Finished. Gain: {} ===", total_score - score_before);
}


void Search::fit_ops_variant2() {
    _logger->info("=== [OPS Variant 2] Running OPS + Hybrid Extended Search (Hub -> Std) ===");

    fit_ops(false);

    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map; // Start empty

    candidate_inserts.reserve(100 * n_variables);
    candidate_reverses.reserve(n_variables);
    candidate_deletes.reserve(n_variables);

    _logger->debug("[OPS Variant 2] Initializing candidates (Pure Mode)...");
    
    initialize_candidates_pure(candidate_inserts, candidate_reverses, candidate_deletes);

    _logger->info("[OPS Variant 2] Phase 1: Running Hub-based Component Incoming Reset...");
    double score_before_hub = total_score;
    
    hub_only_component_incoming_reset_and_research_edge_only(unblocked_paths_map,1,1.0,2);
}

void Search::fit_ops_variant3() {
    // -------------------------------------------------------------------
    // Phase 0: Base OPS
    // -------------------------------------------------------------------
    _logger->info("[Variant 3] Phase 0: Base XGES-0");
    // 1. Run OPS
    fit_ops(false);

    // 2. Priming Step (Candidate Initialization ONLY)
    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map; // Start empty

    candidate_inserts.reserve(100 * n_variables);
    candidate_reverses.reserve(n_variables);
    candidate_deletes.reserve(n_variables);

    // -------------------------------------------------------------------
    // Phase 1: Hub Component Reset (Coarse-grained)
    // -------------------------------------------------------------------
    _logger->info("[Variant 3] Phase 1: Hub Component Reset");
    double score_p1 = total_score;
    // Parameters: threshold=3, ratio=0.1 (Top 10% hubs)
    hub_only_component_incoming_reset_and_research_edge_only(unblocked_paths_map,2,1.0,2);
    
    _logger->info(" > Phase 1 Gain: {}", total_score - score_p1);


    // -------------------------------------------------------------------
    // Phase 2: Node Incoming Reset (Medium-grained)
    // Strategy: Batch delete incoming edges -> Repair(Meek) -> Check Valid(CPDAG) -> Rescore
    // -------------------------------------------------------------------
    _logger->info("[Variant 3] Phase 2: Node Incoming Reset Batch (Try & Validate)");
    double score_p2 = total_score;
    
    node_incoming_reset_batch(unblocked_paths_map,2,1.0,2);
    
    _logger->info(" > Phase 2 Gain: {}", total_score - score_p2);


    // -------------------------------------------------------------------
    // Phase 3: Single Edge Deletion (Fine-grained)
    // Strategy: Single delete -> Repair(Meek) -> Check Valid(CPDAG) -> Rescore
    // -------------------------------------------------------------------
    _logger->info("[Variant 3] Phase 3: Single Edge Deletion (Try & Validate)");
    double score_p3 = total_score;
    
    single_edge_delete_clean_only(unblocked_paths_map,2);
    
    _logger->info(" > Phase 3 Gain: {}", total_score - score_p3);


    _logger->info("=== [Variant 3] Finished. Final Score: {} ===", total_score);
}

void Search::fit_ops_variant4() {
    _logger->info("=== [OPS Variant 2] Base OPS + per-node DeletePa(y, H) (single-node baseline) ===");

    // 1. Run base OPS
    _logger->info("[OPS Variant 2] Phase 0: Base OPS");
    fit_ops(false);

    // 2. Priming Step (Candidate Initialization ONLY)
    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map; // Start empty

    candidate_inserts.reserve(100 * n_variables);
    candidate_reverses.reserve(n_variables);
    candidate_deletes.reserve(n_variables);

    _logger->debug("[OPS Variant 2] Initializing candidates (Pure Mode)...");
    initialize_candidates_pure(candidate_inserts, candidate_reverses, candidate_deletes);

    // 3. Single-node DeletePa(y, H) loop (NO component grouping).
    //    algo_type=2 matches the post-DeletePa research path used for OPS.
    _logger->info("[OPS Variant 2] Phase 1: Per-node DeletePa(y, H) with no component grouping");
    double score_before = total_score;
    node_incoming_reset_batch_ablation_deletepa(
        unblocked_paths_map,
        /*hub_threshold=*/1,
        /*top_k_ratio=*/1.0,
        /*algo_type=*/2,            // OPS dispatcher
        /*max_enum_neighbors=*/30   // safety cap; throws for hubs with |Ne(y)|>10
    );
    _logger->info(" > Variant 2 Gain: {}", total_score - score_before);

    // Logging equivalent to GES/XGES variant4
    _logger->info("DeletePa trials: {}", statistics["deletepa_trials"]);
    _logger->info("research_calls_phase1: {}", statistics["research_calls_phase1"]);
    _logger->info("research_calls_phase2: {}", statistics["research_calls_phase2"]);
    _logger->info("research_calls_phase3: {}", statistics["research_calls_phase3"]);

    _logger->info("extended_search-accepted_phase1: {}", statistics["extended_search-accepted_phase1"]);
    _logger->info("extended_search-accepted_phase2: {}", statistics["extended_search-accepted_phase2"]);
    _logger->info("extended_search-accepted_phase3: {}", statistics["extended_search-accepted_phase3"]);

    _logger->info("extended_search-rejected_phase1: {}", statistics["extended_search-rejected_phase1"]);
    _logger->info("extended_search-rejected_phase2: {}", statistics["extended_search-rejected_phase2"]);
    _logger->info("extended_search-rejected_phase3: {}", statistics["extended_search-rejected_phase3"]);

    _logger->info("=== [OPS Variant 2] Finished. Final Score: {} ===", total_score);
}



void Search::fit_ges(bool use_reverse) {
    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map;

    candidate_inserts.reserve(100 * n_variables);
    candidate_inserts.reserve(n_variables);
    candidate_deletes.reserve(n_variables);


    EdgeModificationsMap edge_modifications;
    int i_operations = 1;
    int last_operation = -1;

    while (i_operations > last_operation) {
        last_operation = i_operations;

        // only apply inserts
        for (int y = 0; y < n_variables; ++y) { find_inserts_to_y(y, candidate_inserts); }
        while (!candidate_inserts.empty()) {
            edge_modifications.clear();
            std::pop_heap(candidate_inserts.begin(), candidate_inserts.end());
            auto best_insert = std::move(candidate_inserts.back());
            candidate_inserts.pop_back();
            if (pdag.is_insert_valid(best_insert, unblocked_paths_map)) {
                pdag.apply_insert(best_insert, edge_modifications);
                total_score += best_insert.score;
                _logger->debug("{}: {}, current_score: {}", i_operations, fmt::streamed(best_insert), total_score);
                i_operations++;
                
                
                for (auto &edge_modification: edge_modifications) {
                    _logger->trace("\tEdge {}", fmt::streamed(edge_modification.second));
                }


                candidate_inserts.clear();
                for (int y = 0; y < n_variables; ++y) {
                    find_inserts_to_y(y, candidate_inserts);
                }
            }
        }
        // only apply deletes
        for (int y = 0; y < n_variables; ++y) { find_deletes_to_y(y, candidate_deletes); }
        while (!candidate_deletes.empty()) {
            edge_modifications.clear();
            std::pop_heap(candidate_deletes.begin(), candidate_deletes.end());
            auto best_delete = std::move(candidate_deletes.back());
            candidate_deletes.pop_back();
            if (pdag.is_delete_valid(best_delete)) {
                pdag.apply_delete(best_delete, edge_modifications);
                total_score += best_delete.score;
                _logger->debug("{}: {}, current_score: {}", i_operations, fmt::streamed(best_delete), total_score);
                i_operations++;
                
                for (auto &edge_modification: edge_modifications) {
                    _logger->trace("\tEdge {}", fmt::streamed(edge_modification.second));
                }


                candidate_deletes.clear();
                for (int y = 0; y < n_variables; ++y) {
                    find_deletes_to_y(y, candidate_deletes);
                }
            }
        }
        // only apply reverses, if use_reverse is true
        // note that without reverse, there is no multiple
        // iterations of forward and backward
        if (!use_reverse) { return; }
        for (int y = 0; y < n_variables; ++y) {
            find_reverse_to_y(y, candidate_reverses);
        }
        while (!candidate_reverses.empty()) {
            edge_modifications.clear();
            std::pop_heap(candidate_reverses.begin(), candidate_reverses.end());
            auto best_reverse = std::move(candidate_reverses.back());
            candidate_reverses.pop_back();
            if (pdag.is_reverse_valid(best_reverse, unblocked_paths_map)) {
                pdag.apply_reverse(best_reverse, edge_modifications);
                total_score += best_reverse.score;
                _logger->debug("{}: {}, current_score: {}", i_operations, fmt::streamed(best_reverse), total_score);
                i_operations++;

                for (auto &edge_modification: edge_modifications) {
                    _logger->trace("\tEdge {}", fmt::streamed(edge_modification.second));
                }

                candidate_reverses.clear();
                for (int y = 0; y < n_variables; ++y) {
                    find_reverse_to_y(y, candidate_reverses);
                }
            }
        }
    }
}

void Search::initialize_candidates_pure(std::vector<Insert> &candidate_inserts,
    std::vector<Reverse> &candidate_reverses,
    std::vector<Delete> &candidate_deletes) const {
// find all possible inserts, deletes, reverses
for (int y = 0; y < n_variables; ++y) {
find_inserts_to_y(y, candidate_inserts, -1, true);
find_deletes_to_y(y, candidate_deletes, true);
find_reverse_to_y(y, candidate_reverses);
}
}


void Search::fit_ges_variant1() {
    _logger->info("=== [GES Variant 1] Running GES + Standard Extended Search ===");

    // 1. Run GES (Enabling reverse phase by default for extended search variants)
    fit_ges(false);

    // 2. Priming Step (Candidate Initialization ONLY)
    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map; // Start empty

    candidate_inserts.reserve(100 * n_variables);
    candidate_reverses.reserve(n_variables);
    candidate_deletes.reserve(n_variables);

    _logger->debug("[GES Variant 1] Initializing candidates (Pure Mode)...");
    
    initialize_candidates_pure(candidate_inserts, candidate_reverses, candidate_deletes);

    // 4. Standard Extended Search
    _logger->info("[GES Variant 1] Starting Standard Extended Search...");
    double score_before = total_score;
    block_each_edge_and_research(unblocked_paths_map,1);
    _logger->info("=== [GES Variant 1] Finished. Gain: {} ===", total_score - score_before);

    _logger->info("extended_search-accepted_phase1: {}", statistics["extended_search-accepted_phase1"]);
    _logger->info("extended_search-accepted_phase2: {}", statistics["extended_search-accepted_phase2"]);
    _logger->info("extended_search-accepted_phase3: {}", statistics["extended_search-accepted_phase3"]);

    _logger->info("extended_search-rejected_phase1: {}", statistics["extended_search-rejected_phase1"]);
    _logger->info("extended_search-rejected_phase2: {}", statistics["extended_search-rejected_phase2"]);
    _logger->info("extended_search-rejected_phase3: {}", statistics["extended_search-rejected_phase3"]);
}

void Search::fit_ges_variant2() {
    _logger->info("=== [GES Variant 2] Running GES + hub node Extended Search===");

    // 1. Run GES
    fit_ges(false);

    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map; // Start empty

    candidate_inserts.reserve(100 * n_variables);
    candidate_reverses.reserve(n_variables);
    candidate_deletes.reserve(n_variables);

    _logger->debug("[GES Variant 2] Initializing candidates (Pure Mode)...");
    
    initialize_candidates_pure(candidate_inserts, candidate_reverses, candidate_deletes);

    _logger->info("[GES Variant 2] Phase 1: Running Hub-based Component Incoming Reset...");
    double score_before_hub = total_score;
    
    hub_only_component_incoming_reset_and_research_edge_only(unblocked_paths_map,1,1.0,1);


    _logger->info("Phase2 trials: {}", statistics["phase2_trials"]);
    _logger->info("Phase3 trials: {}", statistics["phase3_trials"]);
    _logger->info("DeletePa trials: {}", statistics["deletepa_trials"]);
    _logger->info("research_calls_phase1: {}", statistics["research_calls_phase1"]);
    _logger->info("research_calls_phase2: {}", statistics["research_calls_phase2"]);
    _logger->info("research_calls_phase2: {}", statistics["research_calls_phase3"]);

    _logger->info("extended_search-accepted_phase1: {}", statistics["extended_search-accepted_phase1"]);
    _logger->info("extended_search-accepted_phase2: {}", statistics["extended_search-accepted_phase2"]);
    _logger->info("extended_search-accepted_phase3: {}", statistics["extended_search-accepted_phase3"]);

    _logger->info("extended_search-rejected_phase1: {}", statistics["extended_search-rejected_phase1"]);
    _logger->info("extended_search-rejected_phase2: {}", statistics["extended_search-rejected_phase2"]);
    _logger->info("extended_search-rejected_phase3: {}", statistics["extended_search-rejected_phase3"]);
}

void Search::fit_ges_variant3() {
    // -------------------------------------------------------------------
    // Phase 0: Base GES
    // -------------------------------------------------------------------
    _logger->info("[Variant 3] Phase 0: Base GES");
    // 1. Run GES
    fit_ges(false);

    // 2. Priming Step (Candidate Initialization ONLY)
    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map; // Start empty

    candidate_inserts.reserve(100 * n_variables);
    candidate_reverses.reserve(n_variables);
    candidate_deletes.reserve(n_variables);

    // -------------------------------------------------------------------
    // Phase 1: Hub Component Reset (Coarse-grained)
    // -------------------------------------------------------------------
    _logger->info("[Variant 3] Phase 1: Hub Component Reset");
    double score_p1 = total_score;
    hub_only_component_incoming_reset_and_research_edge_only(unblocked_paths_map,2,1.0,1);
    
    _logger->info(" > Phase 1 Gain: {}", total_score - score_p1);


    // -------------------------------------------------------------------
    // Phase 2: Node Incoming Reset (Medium-grained)
    // Strategy: Batch delete incoming edges -> Repair(Meek) -> Check Valid(CPDAG) -> Rescore
    // -------------------------------------------------------------------
    _logger->info("[Variant 3] Phase 2: Node Incoming Reset Batch (Try & Validate)");
    double score_p2 = total_score;
    
    node_incoming_reset_batch(unblocked_paths_map,2,1.0,1);
    
    _logger->info(" > Phase 2 Gain: {}", total_score - score_p2);


    // -------------------------------------------------------------------
    // Phase 3: Single Edge Deletion (Fine-grained)
    // Strategy: Single delete -> Repair(Meek) -> Check Valid(CPDAG) -> Rescore
    // -------------------------------------------------------------------
    _logger->info("[Variant 3] Phase 3: Single Edge Deletion (Try & Validate)");
    double score_p3 = total_score;
    
    single_edge_delete_clean_only(unblocked_paths_map,1);
    
    _logger->info(" > Phase 3 Gain: {}", total_score - score_p3);


    _logger->info("=== [Variant 3] Finished. Final Score: {} ===", total_score);
    _logger->info("Phase2 trials: {}", statistics["phase2_trials"]);
    _logger->info("Phase3 trials: {}", statistics["phase3_trials"]);
    _logger->info("DeletePa trials: {}", statistics["deletepa_trials"]);
    _logger->info("research_calls_phase1: {}", statistics["research_calls_phase1"]);
    _logger->info("research_calls_phase2: {}", statistics["research_calls_phase2"]);
    _logger->info("research_calls_phase2: {}", statistics["research_calls_phase3"]);

    _logger->info("extended_search-accepted_phase1: {}", statistics["extended_search-accepted_phase1"]);
    _logger->info("extended_search-accepted_phase2: {}", statistics["extended_search-accepted_phase2"]);
    _logger->info("extended_search-accepted_phase3: {}", statistics["extended_search-accepted_phase3"]);

    _logger->info("extended_search-rejected_phase1: {}", statistics["extended_search-rejected_phase1"]);
    _logger->info("extended_search-rejected_phase2: {}", statistics["extended_search-rejected_phase2"]);
    _logger->info("extended_search-rejected_phase3: {}", statistics["extended_search-rejected_phase3"]);
}

// =====================================================================================
// fit_ges_variant4: Base GES + per-node DeletePa(y, H) baseline (NO component grouping)
// =====================================================================================
void Search::fit_ges_variant4() {
    _logger->info("=== [GES Variant 2] Base GES + per-node DeletePa(y, H) (single-node baseline) ===");

    // 1. Run base GES (matches variants 1/2/3)
    fit_ges(false);

    // 2. Priming Step (Candidate Initialization ONLY — mirrors variant 2's "Pure Mode")
    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map; // Start empty

    candidate_inserts.reserve(100 * n_variables);
    candidate_reverses.reserve(n_variables);
    candidate_deletes.reserve(n_variables);

    _logger->debug("[GES Variant 2] Initializing candidates (Pure Mode)...");
    initialize_candidates_pure(candidate_inserts, candidate_reverses, candidate_deletes);

    // 3. Single-node DeletePa(y, H) loop (NO component grouping).
    //    algo_type=1 matches the post-DeletePa research path used by variants 2/3 so that
    //    the only material difference between variant 2 and Variant 2 is the
    //    component-wise vs. single-node operator choice.
    _logger->info("[GES Variant 2] Phase 1: Per-node DeletePa(y, H) with no component grouping");
    double score_before = total_score;
    node_incoming_reset_batch_ablation_deletepa(
        unblocked_paths_map,
        /*hub_threshold=*/1,
        /*top_k_ratio=*/1.0,
        /*algo_type=*/1,
        /*max_enum_neighbors=*/30   // safety cap; throws for hubs with |Ne(y)|>10
    );
    _logger->info(" > Variant 2 Gain: {}", total_score - score_before);

    _logger->info("DeletePa trials: {}", statistics["deletepa_trials"]);
    _logger->info("research_calls_phase1: {}", statistics["research_calls_phase1"]);
    _logger->info("research_calls_phase2: {}", statistics["research_calls_phase2"]);
    _logger->info("research_calls_phase3: {}", statistics["research_calls_phase3"]);

    _logger->info("extended_search-accepted_phase1: {}", statistics["extended_search-accepted_phase1"]);
    _logger->info("extended_search-accepted_phase2: {}", statistics["extended_search-accepted_phase2"]);
    _logger->info("extended_search-accepted_phase3: {}", statistics["extended_search-accepted_phase3"]);

    _logger->info("extended_search-rejected_phase1: {}", statistics["extended_search-rejected_phase1"]);
    _logger->info("extended_search-rejected_phase2: {}", statistics["extended_search-rejected_phase2"]);
    _logger->info("extended_search-rejected_phase3: {}", statistics["extended_search-rejected_phase3"]);

    _logger->info("=== [GES Variant 2] Finished. Final Score: {} ===", total_score);
}

std::string to_string_with_precision(double value, int n = 1) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(n) << value;
    return out.str();
}

double Search::set_boss_initial_graph(std::string boss_result_path,
                                    std::string input_path,
                                    double alpha) {

    _logger->debug("[BOSS] Initializing PDAG from BOSS results...");

    std::filesystem::path full_path_obj(boss_result_path);
    std::filesystem::path parent_dir = full_path_obj.parent_path();
    std::string parent_name = parent_dir.filename().string();

    std::string boss_root_dir;
    bool is_realworld = false;

    if (boss_result_path.find("Realworld") != std::string::npos) {
        is_realworld = true;

        // .../outputs/Realworld/ARTH150/boss/run-...
        boss_root_dir = parent_dir.parent_path().parent_path().string();
    }
    else {
        // ER or SF
        // .../outputs/boss/SF/m4/var50/deg2.0/run-...
        std::filesystem::path deg_dir = parent_dir;
        std::filesystem::path var_dir = deg_dir.parent_path();
        std::filesystem::path maybe_m = var_dir.parent_path();

        std::string maybe_m_name = maybe_m.filename().string();

        if (maybe_m_name.size() > 1 && maybe_m_name[0] == 'm') {
            boss_root_dir = maybe_m.parent_path().string();  // SF
        } else {
            boss_root_dir = maybe_m.string();                // ER
        }
    }

    _logger->debug("[BOSS] Root Dir: {}", boss_root_dir);
    std::filesystem::path p_path(input_path);
    std::string filename = p_path.stem().string();

    int p_val = 0;
    double deg_val = 0.0;
    int n_val = 0;
    int id_val = 0;
    int m_val = -1;

    std::string graph_type_str = "ER";
    std::string dataset_name;

    // Real-world
    if (is_realworld) {

        graph_type_str = "REAL";

        char dataset[128];
        sscanf(filename.c_str(),
               "Data_%[^_]_%d_%d",
               dataset, &n_val, &id_val);

        dataset_name = std::string(dataset);

        std::string alpha_str = to_string_with_precision(alpha, 1);
        std::string n_str     = std::to_string(n_val);
        std::string graph_str = std::to_string(id_val);

        std::string matrix_filename =
            "run-" + dataset_name + "_" + n_str +
            "-seed" + graph_str +
            "-alpha" + alpha_str +
            "-boss-variant0-search_result.csv";

        std::string matrix_path =
            boss_root_dir + "/" + dataset_name +
            "/boss/" + matrix_filename;

        _logger->debug("[BOSS] Matrix Path: {}", matrix_path);

        if (!std::filesystem::exists(matrix_path)) {
            _logger->error("[BOSS] Matrix file not found: {}", matrix_path);
            return -1.0;
        }

        // -------- Time Log --------
        std::string log_filename =
            "execution_log_v" + dataset_name +
            "_n" + n_str +
            "_alpha" + alpha_str + ".csv";

        std::string log_path =
            boss_root_dir + "/" + dataset_name +
            "/boss/time_log/" + log_filename;

        double elapsed_time = 0.0;

        std::ifstream log_file(log_path);
        if (log_file.is_open()) {

            std::string line;
            std::getline(log_file, line);

            while (std::getline(log_file, line)) {

                std::stringstream ss(line);
                std::string segment;
                std::vector<std::string> row;

                while (std::getline(ss, segment, ',')) {
                    row.push_back(segment);
                }

                if (row.size() >= 7) {
                    try {
                        if (std::stoi(row[4]) == id_val) {
                            elapsed_time = std::stod(row[6]);
                            break;
                        }
                    } catch (...) { continue; }
                }
            }

            log_file.close();
        }

        // -------- PDAG Load --------
        try {
            this->pdag = PDAG::from_file_pdag(matrix_path);
            EdgeModificationsMap em;
            this->pdag.complete_cpdag_global(em);
            total_score = scorer->score_pdag(pdag);
        }
        catch (const std::exception& e) {
            _logger->error("[BOSS] PDAG load failed: {}", e.what());
            return -1.0;
        }

        _logger->info("[BOSS] Real-world initial graph loaded. Time: {}s",
                      elapsed_time);

        return elapsed_time;
    }

    // Scale-Free
    if (filename.find("ScaleFree") != std::string::npos) {

        graph_type_str = "SF";

        sscanf(filename.c_str(),
               "ScaleFree_m%d_Data_var%d_n%d_seed%d",
               &m_val, &p_val, &n_val, &id_val);

        std::filesystem::path deg_dir = parent_dir;
        std::string deg_dir_name = deg_dir.filename().string();

        if (deg_dir_name.rfind("deg", 0) == 0)
            deg_val = std::stod(deg_dir_name.substr(3));
    }
    else if (filename.find("nonlinear") != std::string::npos) {
        graph_type_str = "NL_ER";
        sscanf(filename.c_str(),
               "Data_nonlinear_var%d_avg_deg%lf_n_samples%d_graph_num%d",
               &p_val, &deg_val, &n_val, &id_val);
    }
    else if (filename.find("gumbel") != std::string::npos) {
        graph_type_str = "gumbel";
        sscanf(filename.c_str(),
               "Data_gumbel_var%d_avg_deg%lf_n_samples%d_graph_num%d",
               &p_val, &deg_val, &n_val, &id_val);
    }
    else {
        // ER
        sscanf(filename.c_str(),
               "Data_var%d_avg_deg%lf_n_samples%d_graph_num%d",
               &p_val, &deg_val, &n_val, &id_val);
    }

    std::string var_str   = std::to_string(p_val);
    std::string deg_str   = to_string_with_precision(deg_val, 1);
    std::string n_str     = std::to_string(n_val);
    std::string alpha_str = to_string_with_precision(alpha, 1);
    std::string graph_str = std::to_string(id_val);

    std::string base_path;

    if (graph_type_str == "SF") {
        base_path = boss_root_dir + "/m" + std::to_string(m_val) +
                    "/var" + var_str + "/deg" + deg_str;
    } else {
        base_path = boss_root_dir +
                    "/var" + var_str + "/deg" + deg_str;
    }

    std::string matrix_filename;

    if (graph_type_str == "SF") {
        matrix_filename =
            "run-SF-m" + std::to_string(m_val) +
            "-" + var_str + "-" + deg_str +
            "-" + n_str + "-alpha" + alpha_str +
            "-boss-variant0-graph" + graph_str +
            "-search_result.csv";
    } 
    else if (graph_type_str == "NL_ER") {
        matrix_filename =
            "run-NL_ER-" + var_str + "-" + deg_str +
            "-" + n_str + "-alpha" + alpha_str +
            "-boss-variant0-graph" + graph_str +
            "-search_result.csv";
    }
    else if (graph_type_str == "gumbel") {
        matrix_filename =
            "run-gumbel-" + var_str + "-" + deg_str +
            "-" + n_str + "-alpha" + alpha_str +
            "-boss-variant0-graph" + graph_str +
            "-search_result.csv";
    }
    
    else {
        matrix_filename =
            "run-ER-" + var_str + "-" + deg_str +
            "-" + n_str + "-alpha" + alpha_str +
            "-boss-variant0-graph" + graph_str +
            "-search_result.csv";
    }

    std::string matrix_path = base_path + "/" + matrix_filename;

    _logger->debug("[BOSS] Matrix Path: {}", matrix_path);

    if (!std::filesystem::exists(matrix_path)) {
        _logger->error("[BOSS] Matrix file not found: {}", matrix_path);
        return -1.0;
    }

    // -------- Time Log --------
    double elapsed_time = 0.0;

    std::string log_filename =
        "execution_log_v" + var_str +
        "_deg" + deg_str +
        "_n" + n_str +
        "_alpha" + alpha_str + ".csv";

    std::string log_path =
        base_path + "/time_log/" + log_filename;

    std::ifstream log_file(log_path);
    if (log_file.is_open()) {

        std::string line;
        std::getline(log_file, line);

        while (std::getline(log_file, line)) {

            std::stringstream ss(line);
            std::string segment;
            std::vector<std::string> row;

            while (std::getline(ss, segment, ',')) {
                row.push_back(segment);
            }

            if (row.size() >= 7) {
                try {
                    if (std::stoi(row[4]) == id_val) {
                        elapsed_time = std::stod(row[6]);
                        break;
                    }
                } catch (...) { continue; }
            }
        }

        log_file.close();
    }

    // -------- PDAG Load --------
    try {
        this->pdag = PDAG::from_file_pdag(matrix_path);
        EdgeModificationsMap em;
        this->pdag.complete_cpdag_global(em);
        total_score = scorer->score_pdag(pdag);
    }
    catch (const std::exception& e) {
        _logger->error("[BOSS] PDAG load failed: {}", e.what());
        return -1.0;
    }

    _logger->info("[BOSS] Initial graph loaded. Time: {}s",
                  elapsed_time);

    return elapsed_time;
}

void Search::heuristic_backward_only(std::vector<Delete> &candidate_deletes,
                                   UnblockedPathsMap &unblocked_paths_map) {
    EdgeModificationsMap edge_modifications;
    int i_operations = 1;
    bool improved = true;

    // GES style: Keep deleting while improvements are possible
    while (improved) {
        improved = false;
        candidate_deletes.clear();

        // 1. Find all possible deletes for current graph
        // (Note: finding deletes is relatively cheap compared to inserts)
        for (int y = 0; y < n_variables; ++y) {
            find_deletes_to_y(y, candidate_deletes, true); // true = positive score only
        }

        if (candidate_deletes.empty()) break;

        // 2. Apply deletes greedily
        // GES usually applies one best move and re-evaluates.
        // Or we can try to apply multiple non-conflicting moves (like XGES heuristic).
        // Here we follow the safer "apply best and update" approach to match GES logic.
        
        // candidate_deletes is a heap (max score at top)
        while (!candidate_deletes.empty()) {
            std::pop_heap(candidate_deletes.begin(), candidate_deletes.end());
            auto best_delete = std::move(candidate_deletes.back());
            candidate_deletes.pop_back();

            if (pdag.is_delete_valid(best_delete)) {
                edge_modifications.clear();
                pdag.apply_delete(best_delete, edge_modifications);
                total_score += best_delete.score;
                _logger->debug("[Backward Only] {}: {}, current_score: {}", 
                               i_operations, fmt::streamed(best_delete), total_score);
                i_operations++;
                improved = true;
                
                // After applying one delete, the scores of other candidates might change.
                // In exact GES, we should re-evaluate. 
                // But efficient GES implementations often update locally.
                // Here, to be simple and consistent with XGES structure, we break and re-find deletes.
                // This corresponds to one edge removal per 'outer' loop iteration.
                // (Optimization: use update_operator_candidates_efficient if desired, 
                //  but for backward phase restarting find is acceptable overhead)
                break; 
            }
        }
    }
}

void Search::fit_boss_variant0(){
    _logger->info("=== [BOSS Variant 1] Running BOSS + Backward Pruning + Standard Extended Search ===");

    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    candidate_inserts.reserve(100 * n_variables);
    candidate_reverses.reserve(n_variables);
    candidate_deletes.reserve(n_variables);
    UnblockedPathsMap unblocked_paths_map;

    // 1. Backward Pruning (Deletion Only Phase)
    _logger->info("[BOSS Variant 0] Step 1: Running Backward Phase (Deletion Only)...");
    double score_before_pruning = total_score;
    heuristic_backward_only(candidate_deletes, unblocked_paths_map);

}



void Search::fit_boss_variant1(){
    _logger->info("=== [BOSS Variant 1] Running BOSS + Backward Pruning + Standard Extended Search ===");

    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;

    candidate_inserts.reserve(100 * n_variables);
    candidate_reverses.reserve(n_variables);
    candidate_deletes.reserve(n_variables);
    UnblockedPathsMap unblocked_paths_map;

    // 1. Backward Pruning (Deletion Only Phase)
    _logger->info("[BOSS Variant 1] Step 1: Running Backward Phase (Deletion Only)...");
    double score_before_pruning = total_score;

    heuristic_backward_only(candidate_deletes, unblocked_paths_map);
    
    _logger->info("[BOSS Variant 1] Backward Pruning Finished. Gain: {}", total_score - score_before_pruning);
    _logger->debug("[BOSS Variant 1] Step 2: Re-initializing candidates for Extended Search...");
    
    candidate_inserts.clear(); candidate_reverses.clear(); candidate_deletes.clear();
    initialize_candidates_pure(candidate_inserts, candidate_reverses, candidate_deletes);

    _logger->info("[BOSS Variant 1] Step 3: Starting Standard Extended Search...");
    double score_before_ext = total_score;
    block_each_edge_and_research(unblocked_paths_map,3);
    _logger->info("=== [BOSS Variant 1] Finished. Final Gain after Ext: {} ===", total_score - score_before_ext);
}


void Search::fit_boss_variant2() {
    _logger->info("=== [BOSS Variant 2] Running BOSS + hub node Extended Search===");

    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map; // Start empty

    candidate_inserts.reserve(100 * n_variables);
    candidate_reverses.reserve(n_variables);
    candidate_deletes.reserve(n_variables);

    // 1. Backward Pruning (Deletion Only Phase)
    _logger->info("[BOSS Variant 1] Step 1: Running Backward Phase (Deletion Only)...");
    double score_before_pruning = total_score;

    heuristic_backward_only(candidate_deletes, unblocked_paths_map);
    
    _logger->info("[BOSS Variant 1] Backward Pruning Finished. Gain: {}", total_score - score_before_pruning);
    _logger->debug("[BOSS Variant 1] Step 2: Re-initializing candidates for Extended Search...");
    
    candidate_inserts.clear(); candidate_reverses.clear(); candidate_deletes.clear();
    initialize_candidates_pure(candidate_inserts, candidate_reverses, candidate_deletes);

    _logger->info("[BOSS Variant 2] Phase 1: Running Hub-based Component Incoming Reset...");
    double score_before_hub = total_score;
    
    hub_only_component_incoming_reset_and_research_edge_only(unblocked_paths_map,1,1.0,3);
}


void Search::fit_boss_variant3() {
    // -------------------------------------------------------------------
    // Phase 0: Base BOSS
    // -------------------------------------------------------------------
    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map; // Start empty

    candidate_inserts.reserve(100 * n_variables);
    candidate_reverses.reserve(n_variables);
    candidate_deletes.reserve(n_variables);

    // 1. Backward Pruning (Deletion Only Phase)
    _logger->info("[BOSS Variant 1] Step 1: Running Backward Phase (Deletion Only)...");
    double score_before_pruning = total_score;
    heuristic_backward_only(candidate_deletes, unblocked_paths_map);
    
    _logger->info("[BOSS Variant 1] Backward Pruning Finished. Gain: {}", total_score - score_before_pruning);

    _logger->debug("[BOSS Variant 1] Step 2: Re-initializing candidates for Extended Search...");
    
    candidate_inserts.clear(); candidate_reverses.clear(); candidate_deletes.clear();
    initialize_candidates_pure(candidate_inserts, candidate_reverses, candidate_deletes);


    // -------------------------------------------------------------------
    // Phase 1: Hub Component Reset (Coarse-grained)
    // -------------------------------------------------------------------
    _logger->info("[Variant 3] Phase 1: Hub Component Reset");
    double score_p1 = total_score;
    hub_only_component_incoming_reset_and_research_edge_only(unblocked_paths_map,2,1.0,3);
    
    _logger->info(" > Phase 1 Gain: {}", total_score - score_p1);


    // -------------------------------------------------------------------
    // Phase 2: Node Incoming Reset (Medium-grained)
    // Strategy: Batch delete incoming edges -> Repair(Meek) -> Check Valid(CPDAG) -> Rescore
    // -------------------------------------------------------------------
    _logger->info("[Variant 3] Phase 2: Node Incoming Reset Batch (Try & Validate)");
    double score_p2 = total_score;
    
    node_incoming_reset_batch(unblocked_paths_map,2,1.0,3);
    
    _logger->info(" > Phase 2 Gain: {}", total_score - score_p2);


    // -------------------------------------------------------------------
    // Phase 3: Single Edge Deletion (Fine-grained)
    // Strategy: Single delete -> Repair(Meek) -> Check Valid(CPDAG) -> Rescore
    // -------------------------------------------------------------------
    _logger->info("[Variant 3] Phase 3: Single Edge Deletion (Try & Validate)");
    double score_p3 = total_score;
    
    single_edge_delete_clean_only(unblocked_paths_map,3);
    
    _logger->info(" > Phase 3 Gain: {}", total_score - score_p3);


    _logger->info("=== [Variant 3] Finished. Final Score: {} ===", total_score);
}


void Search::fit_boss_variant4() {
    _logger->info("=== [BOSS Variant 2] Base BOSS + per-node DeletePa(y, H) (single-node baseline) ===");

    std::vector<Insert> candidate_inserts;
    std::vector<Reverse> candidate_reverses;
    std::vector<Delete> candidate_deletes;
    UnblockedPathsMap unblocked_paths_map; // Start empty

    candidate_inserts.reserve(100 * n_variables);
    candidate_reverses.reserve(n_variables);
    candidate_deletes.reserve(n_variables);

    // 1. Backward Pruning (Deletion Only Phase)
    _logger->info("[BOSS Variant 2] Phase 0: Base BOSS (Backward Phase)...");
    double score_before_pruning = total_score;
    heuristic_backward_only(candidate_deletes, unblocked_paths_map);
    _logger->info("[BOSS Variant 2] Backward Pruning Finished. Gain: {}", total_score - score_before_pruning);

    // 2. Priming Step (Candidate Initialization ONLY)
    _logger->debug("[BOSS Variant 2] Initializing candidates (Pure Mode)...");
    candidate_inserts.clear(); candidate_reverses.clear(); candidate_deletes.clear();
    initialize_candidates_pure(candidate_inserts, candidate_reverses, candidate_deletes);

    // 3. Single-node DeletePa(y, H) loop (NO component grouping).
    //    algo_type=3 matches the post-DeletePa research path used for BOSS.
    _logger->info("[BOSS Variant 2] Phase 1: Per-node DeletePa(y, H) with no component grouping");
    double score_before = total_score;
    node_incoming_reset_batch_ablation_deletepa(
        unblocked_paths_map,
        /*hub_threshold=*/1,
        /*top_k_ratio=*/1.0,
        /*algo_type=*/3,            // BOSS dispatcher
        /*max_enum_neighbors=*/30   // safety cap; throws for hubs with |Ne(y)|>10
    );
    _logger->info(" > Variant 2 Gain: {}", total_score - score_before);

    // Logging equivalent to other variant4s
    _logger->info("DeletePa trials: {}", statistics["deletepa_trials"]);
    _logger->info("research_calls_phase1: {}", statistics["research_calls_phase1"]);
    _logger->info("research_calls_phase2: {}", statistics["research_calls_phase2"]);
    _logger->info("research_calls_phase3: {}", statistics["research_calls_phase3"]);

    _logger->info("extended_search-accepted_phase1: {}", statistics["extended_search-accepted_phase1"]);
    _logger->info("extended_search-accepted_phase2: {}", statistics["extended_search-accepted_phase2"]);
    _logger->info("extended_search-accepted_phase3: {}", statistics["extended_search-accepted_phase3"]);

    _logger->info("extended_search-rejected_phase1: {}", statistics["extended_search-rejected_phase1"]);
    _logger->info("extended_search-rejected_phase2: {}", statistics["extended_search-rejected_phase2"]);
    _logger->info("extended_search-rejected_phase3: {}", statistics["extended_search-rejected_phase3"]);

    _logger->info("=== [BOSS Variant 2] Finished. Final Score: {} ===", total_score);
}

void Search::update_operator_candidates_naive(
        std::vector<Insert> &candidate_inserts, std::vector<Reverse> &candidate_reverses,
        std::vector<Delete> &candidate_deletes) const {
    candidate_inserts.clear();
    candidate_reverses.clear();
    candidate_deletes.clear();
    for (int y = 0; y < n_variables; ++y) {
        find_inserts_to_y(y, candidate_inserts);
        find_reverse_to_y(y, candidate_reverses);
        find_deletes_to_y(y, candidate_deletes);
    }
}

void Search::update_operator_candidates_naive_for_lges(
        std::vector<Insert> &candidate_inserts, std::vector<Reverse> &candidate_reverses,
        std::vector<Delete> &candidate_deletes,
        const PDAG* DAG_ptr,
        bool score_based,
        bool prune) const {
    candidate_inserts.clear();
    candidate_reverses.clear();
    candidate_deletes.clear();
    for (int y = 0; y < n_variables; ++y) {
        find_inserts_to_y_for_lges(y, candidate_inserts,DAG_ptr,score_based,prune);
        find_reverse_to_y(y, candidate_reverses);
        find_deletes_to_y(y, candidate_deletes);
    }
}


void add_pairs(std::set<std::pair<int, int>> &pairs, const FlatSet &xs, int y) {
    for (auto x: xs) { pairs.emplace(x, y); }
}
void add_pairs(std::set<std::pair<int, int>> &pairs, int x, const FlatSet &ys) {
    for (auto y: ys) { pairs.emplace(x, y); }
}
void add_pairs(std::set<std::pair<int, int>> &pairs, const FlatSet &xs,
               const FlatSet &ys) {
    for (auto x: xs) {
        for (auto y: ys) { pairs.emplace(x, y); }
    }
}

/** @brief Update the candidate operators after a PDAG update.
 *
 * After a PDAG update with a set of edge modifications, the candidate operators
 * are updated following the conditions detailed in the XGES paper. The insert, reverse,
 * and delete operators are updated according to the type of each edge modification.
 *
 * @param edge_modifications The edge modifications that have been applied to the PDAG.
 * @param candidate_inserts The vector of candidate inserts.
 * @param candidate_reverses The vector of candidate reverses.
 * @param candidate_deletes The vector of candidate deletes.
 * @param unblocked_paths_map The unblocked paths map used by XGES-0.
 */

void Search::update_operator_candidates_efficient(EdgeModificationsMap &edge_modifications,
                                                std::vector<Insert> &candidate_inserts,
                                                std::vector<Reverse> &candidate_reverses,
                                                std::vector<Delete> &candidate_deletes,
                                                UnblockedPathsMap &unblocked_paths_map) {

    auto start_time = high_resolution_clock::now();
    // First, undo all the edge modifications
    for (auto &[fst, edge_modification]: edge_modifications) {
        pdag.apply_edge_modification(edge_modification, true);
    }

    std::set<int> full_insert_to_y;
    std::map<int, std::set<int>> partial_insert_to_y;
    std::set<int> full_delete_to_y;
    std::set<int> full_delete_from_x;
    std::set<std::pair<int, int>> delete_x_y;
    std::set<int> full_reverse_to_y;
    std::set<int> full_reverse_from_x;
    std::set<std::pair<int, int>> reverse_x_y;

    // Re-apply the edge modifications one by one and update the operators
    for (auto &[fst, edge_modification]: edge_modifications) {
        int a;
        int b;
        if (edge_modification.is_old_directed()) {
            a = edge_modification.get_old_source();
            b = edge_modification.get_old_target();
        } else if (edge_modification.is_new_directed()) {
            a = edge_modification.get_new_source();
            b = edge_modification.get_new_target();
        } else {
            a = edge_modification.x;
            b = edge_modification.y;
        }
        // Track inserts
        switch (edge_modification.get_modification_id()) {
            case 1:// a  b becomes a -- b
                // y = a
                full_insert_to_y.insert(a);
            case 2:// a  b becomes a → b
                // y = b
                full_insert_to_y.insert(b);
                // y \in Ne(a) ∩ Ne(b)
                std::ranges::set_intersection(
                        pdag.get_neighbors(a), pdag.get_neighbors(b),
                        std::inserter(full_insert_to_y, full_insert_to_y.begin()));
                // x=a and y \in Ne(b)
                for (auto target: pdag.get_neighbors(b)) {
                    partial_insert_to_y[target].insert(a);
                }
                // x=b and y \in Ne(a)
                for (auto target: pdag.get_neighbors(a)) {
                    partial_insert_to_y[target].insert(b);
                }
                break;
            case 3:// a -- b becomes a  b
                // x=a and y \in Ne(b) u {b}
                for (auto target: pdag.get_neighbors(b)) {
                    if (target == a) { continue; }
                    partial_insert_to_y[target].insert(a);
                }
                partial_insert_to_y[b].insert(a);
                // y=a and x \in Ad(b)
                partial_insert_to_y[a].insert(pdag.get_adjacent(b).begin(),
                                              pdag.get_adjacent(b).end());
                // x=b and y \in Ne(a) u {a}
                for (auto target: pdag.get_neighbors(a)) {
                    if (target == b) { continue; }
                    partial_insert_to_y[target].insert(b);
                }
                partial_insert_to_y[a].insert(b);
                // y=b and x \in Ad(a)
                partial_insert_to_y[b].insert(pdag.get_adjacent(a).begin(),
                                              pdag.get_adjacent(a).end());
                //SD(x,y,a,b)
                if (unblocked_paths_map.contains({a, b})) {
                    for (auto [x, y]: unblocked_paths_map[{a, b}]) {
                        partial_insert_to_y[y].insert(x);
                    }
                }
                // unblocked_paths_map.erase({a, b}); // Leave them for the reverse
                // SD(x,y,b,a)
                if (unblocked_paths_map.contains({b, a})) {
                    for (auto [x, y]: unblocked_paths_map[{b, a}]) {
                        partial_insert_to_y[y].insert(x);
                    }
                }
                // unblocked_paths_map.erase({b, a}); // Leave them for the reverse
                break;
            case 4:// a -- b becomes a → b
                // y = a and x \in Ad(b)
                partial_insert_to_y[a].insert(pdag.get_adjacent(b).begin(),
                                              pdag.get_adjacent(b).end());
                // y = b
                full_insert_to_y.insert(b);
                // SD(x,y,b,a)
                if (unblocked_paths_map.contains({b, a})) {
                    for (auto [x, y]: unblocked_paths_map[{b, a}]) {
                        partial_insert_to_y[y].insert(x);
                    }
                }
                // unblocked_paths_map.erase({b, a}); // Leave them for the reverse
                break;
            case 5:// a → b becomes a  b
                // x=a and y \in Ne(b) u {b}
                for (auto target: pdag.get_neighbors(b)) {
                    partial_insert_to_y[target].insert(a);
                }
                partial_insert_to_y[b].insert(a);
                // x=b and y \in Ne(a) u {a}
                for (auto target: pdag.get_neighbors(a)) {
                    partial_insert_to_y[target].insert(b);
                }
                partial_insert_to_y[a].insert(b);
                full_insert_to_y.insert(b);
                // SD(x,y,a,b)
                if (unblocked_paths_map.contains({a, b})) {
                    for (auto [x, y]: unblocked_paths_map[{a, b}]) {
                        partial_insert_to_y[y].insert(x);
                    }
                }
                // unblocked_paths_map.erase({a, b}); // Leave them for the reverse
                break;
            case 6:// a → b becomes a -- b
                full_insert_to_y.insert(a);
                full_insert_to_y.insert(b);
                break;
            case 7:// a → b becomes a ← b
                full_insert_to_y.insert(a);
                full_insert_to_y.insert(b);
                // SD(x,y,a,b)
                if (unblocked_paths_map.contains({a, b})) {
                    for (auto [x, y]: unblocked_paths_map[{a, b}]) {
                        partial_insert_to_y[y].insert(x);
                    }
                }
                // unblocked_paths_map.erase({a, b}); // Leave them for the reverse
                break;
        }
        // Track deletes
        FlatSet x_intersection;
        switch (edge_modification.get_modification_id()) {
            case 1:// a  b becomes a -- b
                full_delete_to_y.insert(a);
            case 2:// a  b becomes a → b
                full_delete_to_y.insert(b);
                full_delete_from_x.insert(a);
                full_delete_from_x.insert(b);
                // x \in  Ad(a) ∩ Ad(b) and y \in Ne(a) ∩ Ne(b) [almost never happens]
                std::ranges::set_intersection(
                        pdag.get_adjacent(a), pdag.get_adjacent(b),
                        std::inserter(x_intersection, x_intersection.begin()));
                if (!x_intersection.empty()) {
                    FlatSet y_intersection;
                    std::ranges::set_intersection(
                            pdag.get_neighbors(a), pdag.get_neighbors(b),
                            std::inserter(y_intersection, y_intersection.begin()));
                    add_pairs(delete_x_y, x_intersection, y_intersection);
                }
                break;
            case 3:// a -- b becomes a  b
                break;
            case 4:// a -- b becomes a → b
            case 5:// a → b becomes a  b
                full_delete_to_y.insert(b);
                break;
            case 6:// a → b becomes a -- b
            case 7:// a → b becomes a ← b
                full_delete_to_y.insert(a);
                full_delete_to_y.insert(b);
                break;
        }

        // Track reverse
        switch (edge_modification.get_modification_id()) {
            case 1:// a  b becomes a -- b
                // y \in {a, b}
                full_reverse_to_y.insert(a);
            case 2:// a  b becomes a → b
                // y = b
                full_reverse_to_y.insert(b);
                // y \in Ne(a) ∩ Ne(b)
                std::ranges::set_intersection(
                        pdag.get_neighbors(a), pdag.get_neighbors(b),
                        std::inserter(full_reverse_to_y, full_reverse_to_y.begin()));
                // x \in {a, b}
                full_reverse_from_x.insert(a);
                full_reverse_from_x.insert(b);
                break;
            case 3:// a -- b becomes a  b
                // y \in {a, b} or x \in {a, b}
                full_reverse_to_y.insert(a);
                full_reverse_to_y.insert(b);
                full_reverse_from_x.insert(a);
                full_reverse_from_x.insert(b);
                // SD(x,y,a,b)
                if (unblocked_paths_map.contains({a, b})) {
                    for (auto [x, y]: unblocked_paths_map[{a, b}]) {
                        reverse_x_y.emplace(x, y);
                    }
                    unblocked_paths_map.erase({a, b});
                }
                // SD(x,y,b,a)
                if (unblocked_paths_map.contains({b, a})) {
                    for (auto [x, y]: unblocked_paths_map[{b, a}]) {
                        reverse_x_y.emplace(x, y);
                    }
                    unblocked_paths_map.erase({b, a});
                }
                break;
            case 4:// a -- b becomes a → b
                // y \in {a, b} or x = b
                full_reverse_to_y.insert(a);
                full_reverse_to_y.insert(b);
                full_reverse_from_x.insert(b);
                // SD(x,y,b,a)
                if (unblocked_paths_map.contains({b, a})) {
                    for (auto [x, y]: unblocked_paths_map[{b, a}]) {
                        reverse_x_y.emplace(x, y);
                    }
                    unblocked_paths_map.erase({b, a});
                }
                break;
            case 5:// a → b becomes a  b
                // y = b or x \in {a, b}
                full_reverse_to_y.insert(b);
                full_reverse_from_x.insert(a);
                full_reverse_from_x.insert(b);
                // SD(x,y,a,b)
                if (unblocked_paths_map.contains({a, b})) {
                    for (auto [x, y]: unblocked_paths_map[{a, b}]) {
                        reverse_x_y.emplace(x, y);
                    }
                    unblocked_paths_map.erase({a, b});
                }
                break;
            case 6:// a → b becomes a -- b
                // y \in {a, b} or x = b
                full_reverse_to_y.insert(a);
                full_reverse_to_y.insert(b);
                full_reverse_from_x.insert(b);
                break;
            case 7:// a → b becomes a ← b
                // y \in {a, b} or x \in {a, b}
                full_reverse_to_y.insert(a);
                full_reverse_to_y.insert(b);
                full_reverse_from_x.insert(a);
                full_reverse_from_x.insert(b);
                // SD(x,y,a,b)
                if (unblocked_paths_map.contains({a, b})) {
                    for (auto [x, y]: unblocked_paths_map[{a, b}]) {
                        reverse_x_y.emplace(x, y);
                    }
                    unblocked_paths_map.erase({a, b});
                }
                break;
        }
        pdag.apply_edge_modification(edge_modification);
    }

    // Find the inserts
    // step 1: remove partial inserts that are now full
    std::vector<int> keys_to_erase;
    for (const auto &[y, xs]: partial_insert_to_y) {
        if (full_insert_to_y.contains(y)) { keys_to_erase.push_back(y); }
    }
    for (auto key: keys_to_erase) { partial_insert_to_y.erase(key); }
    // step 2: find the partial inserts
    for (const auto &[y, xs]: partial_insert_to_y) {
        for (auto x: xs) {
            // check that x is not adjacent to y
            if (!pdag.get_adjacent(y).contains(x) && x != y) {
                find_inserts_to_y(y, candidate_inserts, x);
            }
        }
    }
    // step 3: find the full inserts
    for (auto y: full_insert_to_y) { find_inserts_to_y(y, candidate_inserts); }

    // Find the deletes
    // step 1: form the edges to delete
    for (auto x: full_delete_from_x) { add_pairs(delete_x_y, x, pdag.get_neighbors(x)); }
    for (auto x: full_delete_from_x) { add_pairs(delete_x_y, x, pdag.get_children(x)); }
    for (auto y: full_delete_to_y) { add_pairs(delete_x_y, pdag.get_parents(y), y); }
    for (auto y: full_delete_to_y) { add_pairs(delete_x_y, pdag.get_neighbors(y), y); }
    // step 2: find the deletes
    for (auto [x, y]: delete_x_y) {
        if (x != y) { find_delete_to_y_from_x(y, x, candidate_deletes); }
    }

    // Find the reverses
    // step 1: form the edges to reverse
    for (auto x: full_reverse_from_x) { add_pairs(reverse_x_y, x, pdag.get_parents(x)); }
    for (auto y: full_reverse_to_y) { add_pairs(reverse_x_y, pdag.get_children(y), y); }
    // step 2: find the reverses
    for (auto [x, y]: reverse_x_y) {
        // check that x ← y
        if (pdag.has_directed_edge(y, x) && x != y) {
            find_reverse_to_y_from_x(y, x, candidate_reverses);
        }
    }
    statistics["update_operators-time"] += measure_time(start_time);
}

void Search::update_operator_candidates_efficient_for_lges(EdgeModificationsMap &edge_modifications,
                                                std::vector<Insert> &candidate_inserts,
                                                std::vector<Reverse> &candidate_reverses,
                                                std::vector<Delete> &candidate_deletes,
                                                UnblockedPathsMap &unblocked_paths_map,
                                                const PDAG* DAG_ptr,
                                                bool score_based,
                                                bool prune) {

    auto start_time = high_resolution_clock::now();
    // First, undo all the edge modifications
    for (auto &[fst, edge_modification]: edge_modifications) {
        pdag.apply_edge_modification(edge_modification, true);
    }

    std::set<int> full_insert_to_y;
    std::map<int, std::set<int>> partial_insert_to_y;
    std::set<int> full_delete_to_y;
    std::set<int> full_delete_from_x;
    std::set<std::pair<int, int>> delete_x_y;
    std::set<int> full_reverse_to_y;
    std::set<int> full_reverse_from_x;
    std::set<std::pair<int, int>> reverse_x_y;

    // Re-apply the edge modifications one by one and update the operators
    for (auto &[fst, edge_modification]: edge_modifications) {
        int a;
        int b;
        if (edge_modification.is_old_directed()) {
            a = edge_modification.get_old_source();
            b = edge_modification.get_old_target();
        } else if (edge_modification.is_new_directed()) {
            a = edge_modification.get_new_source();
            b = edge_modification.get_new_target();
        } else {
            a = edge_modification.x;
            b = edge_modification.y;
        }
        // Track inserts
        switch (edge_modification.get_modification_id()) {
            case 1:// a  b becomes a -- b
                // y = a
                full_insert_to_y.insert(a);
            case 2:// a  b becomes a → b
                // y = b
                full_insert_to_y.insert(b);
                // y \in Ne(a) ∩ Ne(b)
                std::ranges::set_intersection(
                        pdag.get_neighbors(a), pdag.get_neighbors(b),
                        std::inserter(full_insert_to_y, full_insert_to_y.begin()));
                // x=a and y \in Ne(b)
                for (auto target: pdag.get_neighbors(b)) {
                    partial_insert_to_y[target].insert(a);
                }
                // x=b and y \in Ne(a)
                for (auto target: pdag.get_neighbors(a)) {
                    partial_insert_to_y[target].insert(b);
                }
                break;
            case 3:// a -- b becomes a  b
                // x=a and y \in Ne(b) u {b}
                for (auto target: pdag.get_neighbors(b)) {
                    if (target == a) { continue; }
                    partial_insert_to_y[target].insert(a);
                }
                partial_insert_to_y[b].insert(a);
                // y=a and x \in Ad(b)
                partial_insert_to_y[a].insert(pdag.get_adjacent(b).begin(),
                                              pdag.get_adjacent(b).end());
                // x=b and y \in Ne(a) u {a}
                for (auto target: pdag.get_neighbors(a)) {
                    if (target == b) { continue; }
                    partial_insert_to_y[target].insert(b);
                }
                partial_insert_to_y[a].insert(b);
                // y=b and x \in Ad(a)
                partial_insert_to_y[b].insert(pdag.get_adjacent(a).begin(),
                                              pdag.get_adjacent(a).end());
                //SD(x,y,a,b)
                if (unblocked_paths_map.contains({a, b})) {
                    for (auto [x, y]: unblocked_paths_map[{a, b}]) {
                        partial_insert_to_y[y].insert(x);
                    }
                }
                // unblocked_paths_map.erase({a, b}); // Leave them for the reverse
                // SD(x,y,b,a)
                if (unblocked_paths_map.contains({b, a})) {
                    for (auto [x, y]: unblocked_paths_map[{b, a}]) {
                        partial_insert_to_y[y].insert(x);
                    }
                }
                // unblocked_paths_map.erase({b, a}); // Leave them for the reverse
                break;
            case 4:// a -- b becomes a → b
                // y = a and x \in Ad(b)
                partial_insert_to_y[a].insert(pdag.get_adjacent(b).begin(),
                                              pdag.get_adjacent(b).end());
                // y = b
                full_insert_to_y.insert(b);
                // SD(x,y,b,a)
                if (unblocked_paths_map.contains({b, a})) {
                    for (auto [x, y]: unblocked_paths_map[{b, a}]) {
                        partial_insert_to_y[y].insert(x);
                    }
                }
                // unblocked_paths_map.erase({b, a}); // Leave them for the reverse
                break;
            case 5:// a → b becomes a  b
                // x=a and y \in Ne(b) u {b}
                for (auto target: pdag.get_neighbors(b)) {
                    partial_insert_to_y[target].insert(a);
                }
                partial_insert_to_y[b].insert(a);
                // x=b and y \in Ne(a) u {a}
                for (auto target: pdag.get_neighbors(a)) {
                    partial_insert_to_y[target].insert(b);
                }
                partial_insert_to_y[a].insert(b);
                full_insert_to_y.insert(b);
                // SD(x,y,a,b)
                if (unblocked_paths_map.contains({a, b})) {
                    for (auto [x, y]: unblocked_paths_map[{a, b}]) {
                        partial_insert_to_y[y].insert(x);
                    }
                }
                // unblocked_paths_map.erase({a, b}); // Leave them for the reverse
                break;
            case 6:// a → b becomes a -- b
                full_insert_to_y.insert(a);
                full_insert_to_y.insert(b);
                break;
            case 7:// a → b becomes a ← b
                full_insert_to_y.insert(a);
                full_insert_to_y.insert(b);
                // SD(x,y,a,b)
                if (unblocked_paths_map.contains({a, b})) {
                    for (auto [x, y]: unblocked_paths_map[{a, b}]) {
                        partial_insert_to_y[y].insert(x);
                    }
                }
                // unblocked_paths_map.erase({a, b}); // Leave them for the reverse
                break;
        }
        // Track deletes
        FlatSet x_intersection;
        switch (edge_modification.get_modification_id()) {
            case 1:// a  b becomes a -- b
                full_delete_to_y.insert(a);
            case 2:// a  b becomes a → b
                full_delete_to_y.insert(b);
                full_delete_from_x.insert(a);
                full_delete_from_x.insert(b);
                // x \in  Ad(a) ∩ Ad(b) and y \in Ne(a) ∩ Ne(b) [almost never happens]
                std::ranges::set_intersection(
                        pdag.get_adjacent(a), pdag.get_adjacent(b),
                        std::inserter(x_intersection, x_intersection.begin()));
                if (!x_intersection.empty()) {
                    FlatSet y_intersection;
                    std::ranges::set_intersection(
                            pdag.get_neighbors(a), pdag.get_neighbors(b),
                            std::inserter(y_intersection, y_intersection.begin()));
                    add_pairs(delete_x_y, x_intersection, y_intersection);
                }
                break;
            case 3:// a -- b becomes a  b
                break;
            case 4:// a -- b becomes a → b
            case 5:// a → b becomes a  b
                full_delete_to_y.insert(b);
                break;
            case 6:// a → b becomes a -- b
            case 7:// a → b becomes a ← b
                full_delete_to_y.insert(a);
                full_delete_to_y.insert(b);
                break;
        }

        // Track reverse
        switch (edge_modification.get_modification_id()) {
            case 1:// a  b becomes a -- b
                // y \in {a, b}
                full_reverse_to_y.insert(a);
            case 2:// a  b becomes a → b
                // y = b
                full_reverse_to_y.insert(b);
                // y \in Ne(a) ∩ Ne(b)
                std::ranges::set_intersection(
                        pdag.get_neighbors(a), pdag.get_neighbors(b),
                        std::inserter(full_reverse_to_y, full_reverse_to_y.begin()));
                // x \in {a, b}
                full_reverse_from_x.insert(a);
                full_reverse_from_x.insert(b);
                break;
            case 3:// a -- b becomes a  b
                // y \in {a, b} or x \in {a, b}
                full_reverse_to_y.insert(a);
                full_reverse_to_y.insert(b);
                full_reverse_from_x.insert(a);
                full_reverse_from_x.insert(b);
                // SD(x,y,a,b)
                if (unblocked_paths_map.contains({a, b})) {
                    for (auto [x, y]: unblocked_paths_map[{a, b}]) {
                        reverse_x_y.emplace(x, y);
                    }
                    unblocked_paths_map.erase({a, b});
                }
                // SD(x,y,b,a)
                if (unblocked_paths_map.contains({b, a})) {
                    for (auto [x, y]: unblocked_paths_map[{b, a}]) {
                        reverse_x_y.emplace(x, y);
                    }
                    unblocked_paths_map.erase({b, a});
                }
                break;
            case 4:// a -- b becomes a → b
                // y \in {a, b} or x = b
                full_reverse_to_y.insert(a);
                full_reverse_to_y.insert(b);
                full_reverse_from_x.insert(b);
                // SD(x,y,b,a)
                if (unblocked_paths_map.contains({b, a})) {
                    for (auto [x, y]: unblocked_paths_map[{b, a}]) {
                        reverse_x_y.emplace(x, y);
                    }
                    unblocked_paths_map.erase({b, a});
                }
                break;
            case 5:// a → b becomes a  b
                // y = b or x \in {a, b}
                full_reverse_to_y.insert(b);
                full_reverse_from_x.insert(a);
                full_reverse_from_x.insert(b);
                // SD(x,y,a,b)
                if (unblocked_paths_map.contains({a, b})) {
                    for (auto [x, y]: unblocked_paths_map[{a, b}]) {
                        reverse_x_y.emplace(x, y);
                    }
                    unblocked_paths_map.erase({a, b});
                }
                break;
            case 6:// a → b becomes a -- b
                // y \in {a, b} or x = b
                full_reverse_to_y.insert(a);
                full_reverse_to_y.insert(b);
                full_reverse_from_x.insert(b);
                break;
            case 7:// a → b becomes a ← b
                // y \in {a, b} or x \in {a, b}
                full_reverse_to_y.insert(a);
                full_reverse_to_y.insert(b);
                full_reverse_from_x.insert(a);
                full_reverse_from_x.insert(b);
                // SD(x,y,a,b)
                if (unblocked_paths_map.contains({a, b})) {
                    for (auto [x, y]: unblocked_paths_map[{a, b}]) {
                        reverse_x_y.emplace(x, y);
                    }
                    unblocked_paths_map.erase({a, b});
                }
                break;
        }
        pdag.apply_edge_modification(edge_modification);
    }

    // Find the inserts
    // step 1: remove partial inserts that are now full
    std::vector<int> keys_to_erase;
    for (const auto &[y, xs]: partial_insert_to_y) {
        if (full_insert_to_y.contains(y)) { keys_to_erase.push_back(y); }
    }
    for (auto key: keys_to_erase) { partial_insert_to_y.erase(key); }
    // step 2: find the partial inserts
    for (const auto &[y, xs]: partial_insert_to_y) {
        for (auto x: xs) {
            // check that x is not adjacent to y
            if (!pdag.get_adjacent(y).contains(x) && x != y) {
                find_inserts_to_y_for_lges(y, candidate_inserts,DAG_ptr,score_based,prune,x);
            }
        }
    }
    // step 3: find the full inserts
    for (auto y: full_insert_to_y) { find_inserts_to_y_for_lges(y, candidate_inserts,DAG_ptr,score_based,prune); }

    // Find the deletes
    // step 1: form the edges to delete
    for (auto x: full_delete_from_x) { add_pairs(delete_x_y, x, pdag.get_neighbors(x)); }
    for (auto x: full_delete_from_x) { add_pairs(delete_x_y, x, pdag.get_children(x)); }
    for (auto y: full_delete_to_y) { add_pairs(delete_x_y, pdag.get_parents(y), y); }
    for (auto y: full_delete_to_y) { add_pairs(delete_x_y, pdag.get_neighbors(y), y); }
    // step 2: find the deletes
    for (auto [x, y]: delete_x_y) {
        if (x != y) { find_delete_to_y_from_x(y, x, candidate_deletes); }
    }

    // Find the reverses
    // step 1: form the edges to reverse
    for (auto x: full_reverse_from_x) { add_pairs(reverse_x_y, x, pdag.get_parents(x)); }
    for (auto y: full_reverse_to_y) { add_pairs(reverse_x_y, pdag.get_children(y), y); }
    // step 2: find the reverses
    for (auto [x, y]: reverse_x_y) {
        // check that x ← y
        if (pdag.has_directed_edge(y, x) && x != y) {
            find_reverse_to_y_from_x(y, x, candidate_reverses);
        }
    }
    statistics["update_operators-time"] += measure_time(start_time);
}


/**
 * @brief Find and score all possible inserts to y.
 *
 * The candidate  Insert(x, y, T, E) are such that:
 *  1. x is not adjacent to y (x ∉ Ad(y))
 *  2. T ⊆ Ne(y) \ Ad(x)
 *  3. [Ne(y) ∩ Ad(x)] ∪ T is a clique
 *  Not enforced at that stage: [Ne(y) ∩ Ad(x)] ∪ T blocks all semi-directed paths from y to x
 *  5. E = [Ne(y) ∩ Ad(x)] ∪ T ∪ Pa(y)
 *
 * @param y The target node.
 * @param candidate_inserts The vector of candidate inserts to append inserts to.
 * @param parent_x The parent of y, if known. If -1, no pre-selection is done.
 * @param positive_only If true, only inserts with a positive score are considered.
 */
void Search::find_inserts_to_y(int y, std::vector<Insert> &candidate_inserts, int parent_x,
                             bool positive_only) const {
    auto &adjacent_y = pdag.get_adjacent(y);
    auto &parents_y = pdag.get_parents(y);

    std::set<int> possible_parents;

    if (parent_x != -1) {
        possible_parents.insert(parent_x);
    } else {
        // for now: no pre-selection
        auto &nodes = pdag.get_nodes_variables();

        // 1. x is not adjacent to y (x ∉ Ad(y))
        std::ranges::set_difference(
                nodes, adjacent_y,
                std::inserter(possible_parents, possible_parents.begin()));
        possible_parents.erase(y);// only needed because we don't have pre-selection
    }


    for (int x: possible_parents) {
        // 3. [Ne(y) ∩ Ad(x)] ∪ T is a clique
        // So in particular, [Ne(y) ∩ Ad(x)] must be a clique.
        auto neighbors_y_adjacent_x = pdag.get_neighbors_adjacent(y, x);
        if (!pdag.is_clique(neighbors_y_adjacent_x)) { continue; }

        // 2. T ⊆ Ne(y) \ Ad(x)
        auto neighbors_y_not_adjacent_x = pdag.get_neighbors_not_adjacent(y, x);

        // We enumerate all T ⊆ Ne(y) \ Ad(x) such that [Ne(y) ∩ Ad(x)] ∪ T is a clique (noted C(x,y,T)).
        // If the condition fails for some T, it will fail for all its supersets.
        // Hence, we enumerate the T in inclusion order, to minimize the number of T we need to check.
        // We simulate a recursive search using a stack.
        // The stack contains valid T, and an iterator for the next entries in neighbors_y_not_adjacent_x to consider.

        // The effective parents_y are [Ne(y) ∩ Ad(x)] ∪ T ∪ Pa(y).
        FlatSet &effective_parents_y =
                neighbors_y_adjacent_x;// just renaming it, no copy necessary
        effective_parents_y.insert(parents_y.begin(), parents_y.end());
        // <T: set of nodes, iterator over neighbors_y_not_adjacent_x, effective_parents: set of nodes>
        std::stack<std::tuple<FlatSet, FlatSet::iterator, FlatSet>> stack;
        // we know that T = {} is valid
        stack.emplace(FlatSet{}, neighbors_y_not_adjacent_x.begin(), effective_parents_y);

        while (!stack.empty()) {
            auto top = std::move(stack.top());
            stack.pop();
            auto &T = std::get<0>(top);
            auto it = std::get<1>(top);
            auto &effective_parents = std::get<2>(top);

            // change if we parallelize
            double score = scorer->score_insert(y, effective_parents, x);
            if (score > 0 || !positive_only) {
                // using move(T)/move(effective_parents) should also work even though we look them up
                // later. but we play it safe for now.
                candidate_inserts.emplace_back(x, y, T, score, effective_parents);
                std::push_heap(candidate_inserts.begin(), candidate_inserts.end());
            }

            // Look for other candidate T using the iterator, which gives us the next elements z to consider.
            while (it != neighbors_y_not_adjacent_x.end()) {
                // We define T' = T ∪ {z} and we check C(x,y,T') is a clique.
                // Since C(x,y,T) was a clique, we only need to check that z is adjacent to all nodes in C(x,y,T).
                auto z = *it;
                ++it;
                auto &adjacent_z = pdag.get_adjacent(z);
                // We check that C(x,y,T) ⊆ Ad(z); i.e. T ⊆ Ad(z) and [Ne(y) ∩ Ad(x)] ⊆ Ad(z).
                if (std::ranges::includes(adjacent_z, T) &&
                    std::ranges::includes(adjacent_z, neighbors_y_adjacent_x)) {
                    // T' is a candidate
                    auto T_prime = T;
                    T_prime.insert(z);
                    auto effective_parents_prime = effective_parents;
                    effective_parents_prime.insert(z);
                    stack.emplace(std::move(T_prime), it,
                                  std::move(effective_parents_prime));
                }
            }
        }
    }
}

void Search::find_inserts_to_y_for_lges(int y, std::vector<Insert> &candidate_inserts,const PDAG* DAG_ptr,bool score_based, bool prune, int parent_x,
                             bool positive_only) const {
    auto &adjacent_y = pdag.get_adjacent(y);
    auto &parents_y = pdag.get_parents(y);

    std::set<int> possible_parents;

    if (parent_x != -1) {
        possible_parents.insert(parent_x);
    } else {
        // for now: no pre-selection
        auto &nodes = pdag.get_nodes_variables();

        // 1. x is not adjacent to y (x ∉ Ad(y))
        std::ranges::set_difference(
                nodes, adjacent_y,
                std::inserter(possible_parents, possible_parents.begin()));
        possible_parents.erase(y);// only needed because we don't have pre-selection
    }

    //LGES Safe Insert Logic
    FlatSet origin_parents;
    double current_score=0.0;
    if(score_based&&DAG_ptr!=nullptr){
        origin_parents=DAG_ptr->get_parents(y);
        current_score=scorer->local_score(y,origin_parents);
    }

    for (int x: possible_parents) {
        //Safe Insert Logic
        if(score_based&&DAG_ptr!=nullptr){
            auto new_parents=origin_parents;
            new_parents.insert(x);
            double new_score=scorer->local_score(y,new_parents);
            if(new_score-current_score<=1e-9){
                //if(x==0&&y==7)
                //_logger->debug("Safe Insert: {}->{}",x,y);
                continue;
            }
        }
        // 3. [Ne(y) ∩ Ad(x)] ∪ T is a clique
        // So in particular, [Ne(y) ∩ Ad(x)] must be a clique.
        auto neighbors_y_adjacent_x = pdag.get_neighbors_adjacent(y, x);
        if (!pdag.is_clique(neighbors_y_adjacent_x)) { continue; }

        // 2. T ⊆ Ne(y) \ Ad(x)
        auto neighbors_y_not_adjacent_x = pdag.get_neighbors_not_adjacent(y, x);

        // We enumerate all T ⊆ Ne(y) \ Ad(x) such that [Ne(y) ∩ Ad(x)] ∪ T is a clique (noted C(x,y,T)).
        // If the condition fails for some T, it will fail for all its supersets.
        // Hence, we enumerate the T in inclusion order, to minimize the number of T we need to check.
        // We simulate a recursive search using a stack.
        // The stack contains valid T, and an iterator for the next entries in neighbors_y_not_adjacent_x to consider.

        // The effective parents_y are [Ne(y) ∩ Ad(x)] ∪ T ∪ Pa(y).
        FlatSet &effective_parents_y =
                neighbors_y_adjacent_x;// just renaming it, no copy necessary
        effective_parents_y.insert(parents_y.begin(), parents_y.end());
        // <T: set of nodes, iterator over neighbors_y_not_adjacent_x, effective_parents: set of nodes>
        std::stack<std::tuple<FlatSet, FlatSet::iterator, FlatSet>> stack;
        // we know that T = {} is valid
        stack.emplace(FlatSet{}, neighbors_y_not_adjacent_x.begin(), effective_parents_y);

        std::vector<Insert> conservative_candidates;
        bool conservative=true;

        while (!stack.empty()) {
            auto top = std::move(stack.top());
            stack.pop();
            auto &T = std::get<0>(top);
            auto it = std::get<1>(top);
            auto &effective_parents = std::get<2>(top);

            // change if we parallelize
            double score = scorer->score_insert(y, effective_parents, x);
            //LGES Conservative Insert Logic
            if(prune){
                if(score< -1e-9){
                    conservative=false;
                    break;
                }
                if (score > 0 || !positive_only) {
                    conservative_candidates.emplace_back(x, y, T, score, effective_parents);
                }
            }
            else {
                if (score > 0 || !positive_only) {
                    // using move(T)/move(effective_parents) should also work even though we look them up
                    // later. but we play it safe for now.
                    candidate_inserts.emplace_back(x, y, T, score, effective_parents);
                    std::push_heap(candidate_inserts.begin(), candidate_inserts.end());
                }
            }

            // Look for other candidate T using the iterator, which gives us the next elements z to consider.
            while (it != neighbors_y_not_adjacent_x.end()) {
                // We define T' = T ∪ {z} and we check C(x,y,T') is a clique.
                // Since C(x,y,T) was a clique, we only need to check that z is adjacent to all nodes in C(x,y,T).
                auto z = *it;
                ++it;
                auto &adjacent_z = pdag.get_adjacent(z);
                // We check that C(x,y,T) ⊆ Ad(z); i.e. T ⊆ Ad(z) and [Ne(y) ∩ Ad(x)] ⊆ Ad(z).
                if (std::ranges::includes(adjacent_z, T) &&
                    std::ranges::includes(adjacent_z, neighbors_y_adjacent_x)) {
                    // T' is a candidate
                    auto T_prime = T;
                    T_prime.insert(z);
                    auto effective_parents_prime = effective_parents;
                    effective_parents_prime.insert(z);
                    stack.emplace(std::move(T_prime), it,
                                  std::move(effective_parents_prime));
                }
            }
        }
        if(prune&&conservative){
            for(auto &iter: conservative_candidates){
                candidate_inserts.push_back(std::move(iter));
                std::push_heap(candidate_inserts.begin(), candidate_inserts.end());
            }
        }
    }
}

void Search::find_inserts_to_y_variants(
    int y,
    std::vector<Insert> &candidate_inserts,
    int parent_x,
    bool positive_only
) const {
    auto &adjacent_y = pdag.get_adjacent(y);
    auto &parents_y  = pdag.get_parents(y);

    std::set<int> possible_parents;

    if (parent_x != -1) {
        possible_parents.insert(parent_x);
    } else {
        // for now: no pre-selection
        auto &nodes = pdag.get_nodes_variables();

        // 1. x is not adjacent to y (x ∉ Ad(y))
        std::ranges::set_difference(
            nodes, adjacent_y,
            std::inserter(possible_parents, possible_parents.begin())
        );
        possible_parents.erase(y); // only needed because we don't have pre-selection
    }

    // --- case 1: positive_only == false
    if (!positive_only) {
        for (int x : possible_parents) {
            // 3. [Ne(y) ∩ Ad(x)] ∪ T is a clique
            auto neighbors_y_adjacent_x = pdag.get_neighbors_adjacent(y, x);
            if (!pdag.is_clique(neighbors_y_adjacent_x)) { continue; }

            // 2. T ⊆ Ne(y) \ Ad(x)
            auto neighbors_y_not_adjacent_x = pdag.get_neighbors_not_adjacent(y, x);

            FlatSet &effective_parents_y = neighbors_y_adjacent_x; // just renaming it, no copy necessary
            effective_parents_y.insert(parents_y.begin(), parents_y.end());

            std::stack<std::tuple<FlatSet, FlatSet::iterator, FlatSet>> stack;
            stack.emplace(FlatSet{}, neighbors_y_not_adjacent_x.begin(), effective_parents_y);

            while (!stack.empty()) {
                auto top = std::move(stack.top());
                stack.pop();
                auto &T                = std::get<0>(top);
                auto it                = std::get<1>(top);
                auto &effective_parents = std::get<2>(top);

                double score = scorer->score_insert(y, effective_parents, x);
                candidate_inserts.emplace_back(x, y, T, score, effective_parents);
                std::push_heap(candidate_inserts.begin(), candidate_inserts.end());

                while (it != neighbors_y_not_adjacent_x.end()) {
                    auto z = *it;
                    ++it;
                    auto &adjacent_z = pdag.get_adjacent(z);
                    if (std::ranges::includes(adjacent_z, T) &&
                        std::ranges::includes(adjacent_z, neighbors_y_adjacent_x)) {
                        auto T_prime = T;
                        T_prime.insert(z);
                        auto effective_parents_prime = effective_parents;
                        effective_parents_prime.insert(z);
                        stack.emplace(std::move(T_prime), it, std::move(effective_parents_prime));
                    }
                }
            }
        }
        return;
    }

    // --- case 2: positive_only == true -> edge-level filter
    constexpr double TAU_POS_FRAC = 0.5;

    for (int x : possible_parents) {
        // 3. [Ne(y) ∩ Ad(x)] ∪ T is a clique
        auto neighbors_y_adjacent_x = pdag.get_neighbors_adjacent(y, x);
        if (!pdag.is_clique(neighbors_y_adjacent_x)) { continue; }

        // 2. T ⊆ Ne(y) \ Ad(x)
        auto neighbors_y_not_adjacent_x = pdag.get_neighbors_not_adjacent(y, x);

        FlatSet &effective_parents_y = neighbors_y_adjacent_x; // just renaming it, no copy necessary
        effective_parents_y.insert(parents_y.begin(), parents_y.end());

        std::stack<std::tuple<FlatSet, FlatSet::iterator, FlatSet>> stack;
        stack.emplace(FlatSet{}, neighbors_y_not_adjacent_x.begin(), effective_parents_y);

        std::vector<Insert> local_inserts;
        local_inserts.reserve(32);

        while (!stack.empty()) {
            auto top = std::move(stack.top());
            stack.pop();
            auto &T                 = std::get<0>(top);
            auto it                 = std::get<1>(top);
            auto &effective_parents = std::get<2>(top);

            double score = scorer->score_insert(y, effective_parents, x);

            local_inserts.emplace_back(x, y, T, score, effective_parents);

            while (it != neighbors_y_not_adjacent_x.end()) {
                auto z = *it;
                ++it;
                auto &adjacent_z = pdag.get_adjacent(z);
                if (std::ranges::includes(adjacent_z, T) &&
                    std::ranges::includes(adjacent_z, neighbors_y_adjacent_x)) {
                    auto T_prime = T;
                    T_prime.insert(z);
                    auto effective_parents_prime = effective_parents;
                    effective_parents_prime.insert(z);
                    stack.emplace(std::move(T_prime), it, std::move(effective_parents_prime));
                }
            }
        }
        int pos = 0;
        int neg = 0;
        double best_pos_score = -std::numeric_limits<double>::infinity();
        std::size_t best_pos_idx = static_cast<std::size_t>(-1);

        for (std::size_t i = 0; i < local_inserts.size(); ++i) {
            double s = local_inserts[i].score;
            if (s > 0) {
                pos++;
                if (s > best_pos_score) {
                    best_pos_score = s;
                    best_pos_idx   = i;
                }
            } else if (s < 0) {
                neg++;
            }
        }
        if (pos == 0) { continue; }
        double pos_frac = (pos + neg > 0) ? static_cast<double>(pos) / (pos + neg) : 1.0;
        if (neg > 0 && pos_frac < TAU_POS_FRAC) {
            continue;
        }
        if (best_pos_idx != static_cast<std::size_t>(-1)) {
            const auto &ins = local_inserts[best_pos_idx];
            candidate_inserts.push_back(ins);
            std::push_heap(candidate_inserts.begin(), candidate_inserts.end());
        }
    }
}

/** Find all possible deletes to y from x.
 *
 * The candidate deletes (x, y, C, E) are such that:
 *  1. x is a parent or a neighbor of y
 *  2. C ⊆ [Ne(y) ∩ Ad(x)]
 *  3. C is a clique
 *  4. E = Pa(y) ∪ [Ne(y) ∩ Ad(x)] \ C ∪ {x}
 *
 *
 * @param y The target node.
 * @param x The source node.
 * @param candidate_deletes The vector of candidate deletes to append deletes to.
 * @param positive_only If true, only deletes with a positive score are considered.
 */


void Search::find_delete_to_y_from_x(int y, int x, std::vector<Delete> &candidate_deletes,
                                   bool positive_only) const {
    const FlatSet &parents_y = pdag.get_parents(y);
    auto neighbors_y_adjacent_x = pdag.get_neighbors_adjacent(y, x);
    bool directed_xy = pdag.has_directed_edge(x, y);

    //duplicated undirected delete is not allowed, only when x<y ㅣ \text{If edge is undirected and } x>y,\ \text{skip it.}
    if (!directed_xy && x > y) {return;}


    // find all possible C ⊆ [Ne(y) ∩ Ad(x)] that are cliques
    // <C set of nodes, iterator over neighbors_y_adjacent_x, set of effective_parents>
    std::stack<std::tuple<FlatSet, FlatSet::iterator, FlatSet>> stack;
    FlatSet effective_parents_init;
    effective_parents_init.reserve(parents_y.size() + neighbors_y_adjacent_x.size() + 1);
    // note: Chickering Corollary 18 is incorrect. Pa(y) might not contain x, it has to be added.
    union_with_single_element(parents_y, x, effective_parents_init);
    // note that C = {} is valid
    stack.emplace(FlatSet{}, neighbors_y_adjacent_x.begin(), effective_parents_init);

    while (!stack.empty()) {
        auto top = std::move(stack.top());
        stack.pop();
        auto C = std::get<0>(top);
        auto it = std::get<1>(top);
        auto effective_parents = std::get<2>(top);

        double score = scorer->score_delete(y, effective_parents, x);
        if (score > 0 || !positive_only) {
            candidate_deletes.emplace_back(x, y, C, score, effective_parents,
                                           directed_xy);
            std::push_heap(candidate_deletes.begin(), candidate_deletes.end());
        }

        // Look for other candidate C using the iterator, which gives us the next elements
        // z to consider.
        while (it != neighbors_y_adjacent_x.end()) {
            // We define C' = C ∪ {z} and we check if C' is a clique.
            // Since C is a clique, we only need to check that z is adjacent to all nodes
            // in C.
            auto z = *it;
            ++it;
            auto &adjacent_z = pdag.get_adjacent(z);
            // We check C ⊆ Ad(z)
            if (std::ranges::includes(adjacent_z, C)) {
                // C' is a candidate
                auto C_prime = C;
                C_prime.insert(z);
                auto effective_parents_prime = effective_parents;
                effective_parents_prime.insert(z);
                stack.emplace(std::move(C_prime), it, std::move(effective_parents_prime));
            }
        }
    }
}


/** @brief Find all possible deletes to y.
 *
 * @see find_delete_to_y_from_x
 *
 * @param y The target node.
 * @param candidate_deletes The vector of candidate deletes to append deletes to.
 * @param positive_only If true, only deletes with a positive score are considered.
 */
void Search::find_deletes_to_y(const int y, std::vector<Delete> &candidate_deletes,
                             bool positive_only) const {
    auto &neighbors_y = pdag.get_neighbors(y);
    auto &parents_y = pdag.get_parents(y);

    for (int x: parents_y) {
        find_delete_to_y_from_x(y, x, candidate_deletes, positive_only);
    }
    for (int x: neighbors_y) {
        find_delete_to_y_from_x(y, x, candidate_deletes, positive_only);
    }
}

/**
 * Only used for the naive update of the operators.
 */
void Search::find_reverse_to_y(int y, std::vector<Reverse> &candidate_reverses) const {
    // look for all possible x ← y
    auto &children_y = pdag.get_children(y);

    for (int x: children_y) {
        auto &parents_x = pdag.get_parents(x);
        std::vector<Insert> candidate_inserts;
        find_inserts_to_y(y, candidate_inserts, x, false);

        for (auto &insert: candidate_inserts) {
            // change if we parallelize
            double score = insert.score + scorer->score_delete(x, parents_x, y);

            if (score > 0) {
                candidate_reverses.emplace_back(insert, score, parents_x);
                std::push_heap(candidate_reverses.begin(), candidate_reverses.end());
            }
        }
    }
}

/** @brief Find all possible reverses to y from x.
 *
 * @see find_inserts_to_y
 *
 * The candidate Reverse(x, y, T, E, F) are such that:
 * 1. x ← y
 * 2. T ⊆ Ne(y) \ Ad(x)
 * 3. [Ne(y) ∩ Ad(x)] ∪ T is a clique
 * 4. All semi-directed paths from y to x are blocked by [Ne(y) ∩ Ad(x)] ∪ T ∪ Ne(x)
 * 5. E = [Ne(y) ∩ Ad(x)] ∪ T ∪ Pa(y)
 * 6. F = Pa(x)
 *
 * A reverse is composed of a delete and an insert. The function first finds all possible
 * inserts to y from x using `find_inserts_to_y` and then evaluate the delete part.
 *
 */
void Search::find_reverse_to_y_from_x(int y, int x,
                                    std::vector<Reverse> &candidate_reverses) const {
    if (!pdag.has_directed_edge(y, x)) { return; }
    std::vector<Insert> candidate_inserts;
    find_inserts_to_y(y, candidate_inserts, x, false);
    auto &parents_x = pdag.get_parents(x);
    for (auto &insert: candidate_inserts) {
        double score = insert.score + scorer->score_delete(x, parents_x, y);
        if (score > 0) {
            candidate_reverses.emplace_back(insert, score, parents_x);
            std::push_heap(candidate_reverses.begin(), candidate_reverses.end());
        }
    }
}


const PDAG &Search::get_pdag() const { return pdag; }
double Search::get_score() const { return total_score; }
double Search::get_initial_score() const{ return initial_score; }