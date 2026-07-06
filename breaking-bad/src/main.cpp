#include <spdlog/fmt/bundled/ostream.h> 
#include <spdlog/fmt/bundled/std.h> 


#include "PDAG.h"
#include "Search.h"
#include "utils.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>

#include "cxxopts.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include <spdlog/sinks/basic_file_sink.h>



namespace fs = std::filesystem;
using namespace std::chrono;

// Builds the Gaussian BIC scorer from an input .npy file. Implemented in
// BICScoreloader.cpp — complete that file before building
// (see the README section "BIC scorer").
std::unique_ptr<ScorerInterface>
make_bic_scorer_from_npy(const std::string &npy_path, double alpha,
                         int &n_variables, int &n_samples);


int main(int argc, char *argv[]) {
    srand(0);
    std::cout << std::setprecision(16);
    // Set up the logger
    auto logger = spdlog::stdout_color_mt("stdout_logger");
    logger->set_pattern("[%^%l%$] %v");

    // Define the command line options
    cxxopts::Options options("xges", "Run on given data");
    auto option_adder = options.add_options();
    option_adder("input", "Input data numpy file (must be a contiguous C array)",
                 cxxopts::value<std::string>());
    option_adder("output", "Output file (default `search-graph.txt`)",
                 cxxopts::value<std::string>()->default_value("search-graph.txt"));
    option_adder("alpha,a", "Alpha parameter for the BIC score (default 2.)",
                 cxxopts::value<double>()->default_value("2."));
    option_adder("stats", "File to save statistics (default `search-stats.txt`)",
                 cxxopts::value<std::string>()->default_value("search-stats.txt"));
    option_adder("0,xges0", "Do not perform the extended search of XGES, just XGES-0.",
                 cxxopts::value<bool>()->default_value("false"));
    option_adder("baseline,b", "Run a baseline instead of XGES",
                 cxxopts::value<std::string>()->default_value(""));
    option_adder("graph_truth,g", "Graph truth file", cxxopts::value<std::string>());


    option_adder("initial_pdag,i", "Initial PDAG file", cxxopts::value<std::string>()->default_value(""));

    option_adder("verbose,v", "Level of verbosity (0-3)",
                 cxxopts::value<int>()->default_value("1"));

    option_adder("variant,x", "Variant (1-3)",
                 cxxopts::value<int>()->default_value("0"));
    option_adder("k_extended,k", "K value for K-extended search",
                 cxxopts::value<int>()->default_value("2"));
    option_adder("version,e", "Version (1-4)",
                 cxxopts::value<int>()->default_value("1"));
    option_adder("delete_op,d", "Hub node deletion options",
                 cxxopts::value<int>()->default_value("1"));

    auto args = options.parse(argc, argv);

    if (int verbose_level = args["verbose"].as<int>(); verbose_level == 0) {
        logger->set_level(spdlog::level::off);
    } else if (verbose_level == 1) {
        logger->set_level(spdlog::level::info);
    } else if (verbose_level == 2) {
        logger->set_level(spdlog::level::debug);
    } else if (verbose_level == 3) {
        logger->set_level(spdlog::level::trace);
    } else {
        throw std::runtime_error("Invalid verbose level");
    }


    // Parse the command line options
    fs::path data_path = args["input"].as<std::string>();
    fs::path output_path = args["output"].as<std::string>();
    double alpha = args["alpha"].as<double>();

    logger->info("Loading input: {}", data_path.string());
    if (data_path.extension() != ".npy") {
        throw std::runtime_error("Input file must be a .npy file");
    }

    logger->info("Computing covariance matrix");
    auto start_time = high_resolution_clock::now();
    int n_variables = 0, n_samples = 0;
    std::unique_ptr<ScorerInterface> scorer_ptr =
            make_bic_scorer_from_npy(data_path.string(), alpha, n_variables, n_samples);
    ScorerInterface &scorer = *scorer_ptr;
    double elapsed_secs = measure_time(start_time);
    logger->info("Input loaded. Shape: {} x {}", n_samples, n_variables);
    logger->info("Covariance computed in {} seconds", elapsed_secs);


    std::unique_ptr<Search> Search_ptr;

    std::string init_path = args["initial_pdag"].as<std::string>();
    if (!init_path.empty()) {
        logger->info("Loading initial graph from file: {}", init_path);
        PDAG init_pdag = PDAG::from_file_pdag(init_path);
    

        if ((int)init_pdag.get_nodes_variables().size() != n_variables) {
            throw std::runtime_error("Initial PDAG node count mismatch with data columns.");
        }
    
        Search_ptr = std::make_unique<Search>(init_pdag, &scorer);
    } else {
        Search_ptr = std::make_unique<Search>(n_variables, &scorer);
    }

    Search& Search = *Search_ptr;
    logger->info("initialized with initial score: {}", Search.get_initial_score());


    if (args.count("graph_truth") > 0) {
        auto ground_truth_pdag = std::make_unique<PDAG>(
                PDAG::from_file(args["graph_truth"].as<std::string>()));
        Search.ground_truth_pdag = std::move(ground_truth_pdag);

        //scoring ground_truth_pdag
        double ground_truth_score = scorer.score_pdag(*Search.ground_truth_pdag);
        logger->info("Score of ground truth PDAG: {}", ground_truth_score);
        int n_edges_truth=Search.ground_truth_pdag->get_number_of_edges();
        logger->info("Number of edges in ground truth PDAG: {}", n_edges_truth);
    }

    if (args.count("baseline") > 0) {
        int variant = args["variant"].as<int>();
        int k_extended = args["k_extended"].as<int>();
        int version = args["version"].as<int>();
        bool extended_search = !args["0"].as<bool>();
        int deletion_options = args["delete_op"].as<int>();
        std::string baseline = args["baseline"].as<std::string>();
        start_time = high_resolution_clock::now();
        if  (baseline == "xges"){
            //set variant
            int variant = args["variant"].as<int>();
            int k_extended = args["k_extended"].as<int>();
            int version = args["version"].as<int>();
            bool extended_search = !args["0"].as<bool>();
            int deletion_options = args["delete_op"].as<int>();
            start_time = high_resolution_clock::now();
            
            if (variant == 0) {
                logger->info("Fitting XGES0 without extended search: {}", extended_search);
                Search.fit_xges0(extended_search);
            }
            else if (variant == 1) {
                logger->info("Fitting XGES with extended search: {}", extended_search);
                Search.fit_xges(extended_search);
            }
            
            else if (variant == 2) {
                logger->info("Fitting XGES with Full Hub Scan (Ratio 1.0): {}", extended_search);
                Search.fit_xges_variant2();
            }
            else if (variant == 3) {
                logger->info("Fitting XGES with Full Hub Scan (Ratio 1.0): {}", extended_search);
                Search.fit_xges_variant3();
            }
            else if (variant == 4) {
                logger->info("Fitting XGES Variant 4 (Single-node DeletePa)");
                Search.fit_xges_variant4();
            }
            else {
                logger->info("Fitting with extended search: {}", extended_search);
                
                Search.fit_xges(extended_search);
    
            }
            elapsed_secs = measure_time(start_time);
            logger->info("search completed in {} seconds", elapsed_secs);
        }
        else if (baseline == "ops") {
            if (variant == 1) {
                logger->info("Baseline GES with Variant 1 (Standard Extended Search)");
                Search.fit_ops_variant1();
            }else if (variant == 2) {
                logger->info("Baseline GES with Variant 2 (DP)");
                Search.fit_ops_variant2();
            }else if (variant == 3) {
                logger->info("Baseline GES with Variant 3 (DP+)");
                Search.fit_ops_variant3();
            } else if (variant == 4) {
                logger->info("Baseline OPS with Variant 4 (Single-node)");
                Search.fit_ops_variant4(); 
            }
            else {
                // Default GES (No Reverse)
                Search.fit_ops(false);
            }
        } else if (baseline == "ops-r") {
            Search.fit_ops(true);
        } else if (baseline == "ges") {
            
            if (variant == 1) {
                logger->info("Baseline GES with Variant 1 (Standard Extended Search)");
                Search.fit_ges_variant1();
            }else if (variant == 2) {
                logger->info("Baseline GES with Variant 2 (DP)");
                Search.fit_ges_variant2();
            }else if (variant == 3) {
                logger->info("Baseline GES with Variant 3 (DP+)");
                Search.fit_ges_variant3();
            }else if (variant == 4) {
                logger->info("Baseline GES with Variant 4 (Single-node)");
                Search.fit_ges_variant4();
            } else {
                // Default GES (No Reverse)
                Search.fit_ges(false);
            }
        } else if (baseline == "lges-safe") {
            if(variant==1){
                Search.fit_lges_variant1(true,false);
            } else if(variant==2){
                Search.fit_lges_variant2(true,false);
            } else if(variant==3){
                Search.fit_lges_variant3(true,false);
            } else if(variant==4){
                Search.fit_lges_variant4(true,false);
            }
            else {
                Search.fit_lges(true,false);
            } 
        } else if (baseline == "lges-cons") {
            if(variant==1){
                Search.fit_lges_variant1(true,true);
            } else if(variant==2){
                Search.fit_lges_variant2(true,true);
            } else if(variant==3){
                Search.fit_lges_variant3(true,true);
            } else if(variant==4){
                Search.fit_lges_variant4(true,true);
            }
            else {
                Search.fit_lges(true,true);
            } 
        } else if (baseline == "ges-r") {
            Search.fit_ges(true);
        } else if (baseline == "boss"){
            double boss_time=Search.set_boss_initial_graph(output_path,data_path,alpha);
            start_time = high_resolution_clock::now();
            if(variant==0){
                logger->info("Baseline BOSS");
                Search.fit_boss_variant0();
            }
            else if(variant == 1){
                logger->info("Baseline BOSS with Variant 1 (Standard Extended Search)");
                Search.fit_boss_variant1();
            }
            else if(variant == 2){
                logger->info("Baseline BOSS with Variant 2 (DP)");
                Search.fit_boss_variant2();
            }
            else if(variant == 3){
                logger->info("Baseline BOSS with Variant 3 (DP+)");
                Search.fit_boss_variant3();
            }
            else if(variant == 4){
                logger->info("Baseline BOSS with Variant 4 (Single-node)");
                Search.fit_boss_variant4();
            }
            else{
                throw std::runtime_error("Invalid boss variant");
            }
            elapsed_secs=measure_time(start_time);
            elapsed_secs+=boss_time;
        } else {
            throw std::runtime_error("Invalid baseline");
        }
        if(baseline!="boss") {
            elapsed_secs = measure_time(start_time);
        }
        logger->info("Baseline mode: {}", baseline);
        logger->info("Baseline completed in {} seconds", elapsed_secs);
        
    } 

    logger->info("Score: {}", Search.get_score());
    int n_edges_pred=Search.get_pdag().get_number_of_edges();
    logger->info("Number of edges in predicted PDAG: {}", n_edges_pred);

    // Save the output
    std::ofstream out_file(output_path);
    out_file << Search.get_pdag().get_adj_string();
    out_file.close();

    // Save the statistics
    std::ofstream stats_file(args["stats"].as<std::string>());
    stats_file << std::setprecision(16);
    stats_file << "time, " << elapsed_secs << std::endl;
    stats_file << "score, " << Search.get_score() << std::endl;
    stats_file << "score check, " << scorer.score_pdag(Search.get_pdag()) << std::endl;
    stats_file << "score_empty, " << Search.get_initial_score() << std::endl;
    stats_file << "score_increase, " << Search.get_score() - Search.get_initial_score()
               << std::endl;
    for (auto &kv: Search.statistics) {
        stats_file << kv.first << ", " << kv.second << std::endl;
    }
    for (auto &kv: Search.get_pdag().statistics) {
        stats_file << kv.first << ", " << kv.second << std::endl;
    }
    for (auto &kv: scorer.statistics) {
        stats_file << kv.first << ", " << kv.second << std::endl;
    }
    return 0;
}


void test_pdag() {
    PDAG pdag(10);
    pdag.add_undirected_edge(0, 1);
    pdag.add_undirected_edge(1, 2);

    PDAG dag_extension_true1(10);
    dag_extension_true1.add_directed_edge(0, 1);
    dag_extension_true1.add_directed_edge(1, 2);
    PDAG dag_extension_true2(10);
    dag_extension_true2.add_directed_edge(2, 1);
    dag_extension_true2.add_directed_edge(1, 0);
    assert(pdag.get_dag_extension() == dag_extension_true1 ||
           pdag.get_dag_extension() == dag_extension_true2);
}