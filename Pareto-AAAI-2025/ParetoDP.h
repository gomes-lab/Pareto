//
// Created by Jonathan Gomes Selman on 7/6/17.
//

#ifndef AMAZON_PROJECT_DP_ALGORITHM_H
#define AMAZON_PROJECT_DP_ALGORITHM_H

#include <cstdlib>
#include "Network.h"
#include <stack>
#include <cmath>
#include <thread>
#include <algorithm>
#include <mutex>
#include <iomanip>
#include <unordered_set>
#include <cmath>
#include "Solution.h"
#include "Consts.h"
#include "IOHandler.h"
#include "Transform.h"
#include "DPBase.h"
#include "DPTransforms.h"
#include "DPOriginal.h"
#include "DPBruteForce.h"
#include "Config.h"
#include "SolutionCompressor.h"
#include <bitset>

#define COMPLETION_FLAG 1

struct RunInfo {
    unsigned long num_policies_generated = 0;
    unsigned long num_policies_pruned = 0;
    unsigned long num_policies_final = 0;
    unsigned long num_comparisons = 0;

    unsigned long transforms_considered;
    unsigned long transforms_pruned;
    unsigned long policies_avoided;

    unsigned long max_node_policies = 0;
    float wall_time = 0;
    float cpu_time = 0;

    RunInfo(unsigned long num_generated, unsigned long num_pruned, unsigned long num_final, unsigned long max_num,
            unsigned long num_comparisons, float wall_time, float cpu_time, unsigned long transforms_considered,
            unsigned long transforms_pruned, unsigned long policies_avoided) :
        num_policies_generated(num_generated), num_policies_pruned(num_pruned), num_policies_final(num_final),
        num_comparisons(num_comparisons), max_node_policies(max_num), wall_time(wall_time), cpu_time(cpu_time),
        transforms_considered(transforms_considered), transforms_pruned(transforms_pruned), policies_avoided(policies_avoided) {

    }

    RunInfo(RunInfo &src) = default;

    RunInfo() {
        num_policies_generated = 0;
        num_policies_pruned = 0;
        num_policies_final = 0;
        num_comparisons = 0;

        transforms_considered = 0;
        transforms_pruned = 0;
        policies_avoided = 0;

        max_node_policies = 0;
        wall_time = 0;
        cpu_time = 0;
    }
};

class ParetoDP {
public:
    ParetoDP(Network *net, const IOHandler& io, const Config& config, bool brute_force); // Constructor
    virtual ~ParetoDP();


    /*
     * Calculate the k values for each node.
     * -------------------------------------
     * The k values are calculated by taking the local value for each criteria and
     * multiplying by the constant K_EPSILON.
     *      k_energy = K_EPSILON * associated_energy
     */
    void calc_theoretical_ks_tree(TreeNode* node);

    /*
     * Implements the dynamic programming algorithm with a divide-and-conquer method
     * for determining non-dominated solutions as well as batched pruning.
     */
    void computeDP(TreeNode *node, int dynamic_depth, int static_depth);

    void writeFrontier(TreeNode *node);

    void writeMeta(TreeNode *node, float duration);

    void build_DP();
    void brute_force();

    /*
     * Recursively traverse the binary-tree with a post-order traversal calculating
     * the Pareto_Opt_List for each parent node.
     * -----------------------------------------
     * 1) Recurse on left/right children
     * 2) Compute DP for parent
     */
    void build_DP_table_recursive_tree(TreeNode* node, int dynamic_depth, int static_depth);

    // Below are verious methods to print out information based on each experimental run
    // ---------------------------------------------------------------------------------
    void print_DP_Output(std::ostream &stream = std::cout); // If no output stream is provided then prints  to console

    void print_dams_generic_tree_vec(TreeNode* node, std::ostream &out = std::cout);

    void print_dams_helper(Solution *result, std::ostream &out, TreeNode *node);
    // --------------------------------------------------------------------------------

    bool readPolicies(TreeNode* node);

    /*
     * Runs an experiment for a given epsilon and river basin.
     * Outputs all results to a file with a standardized format.
     * Keeps track of several experiment metrics such as the execution
     * time and time spent doing particular tasks.
     */
    void run_experiment();
    void run_brute_force();

    Network *net; // Hyper node network with the all the data
    RunInfo run_info;
    DPBase *dp;

private:
    // Used to track how much progress
    float node_processed_counter;

    Config config;

    // Used to calculate the theoretical k values
    int num_criteria;
    int true_num_criteria;

    // Used for experimentation on algorithm complexity
    float time_sorting;
    float time_generating;
    float time_nlogn;
    float time_copying;
    float time_generating_transforms;
    float time_sorting_transforms;
    float time_transforming;
    float time_merging;

    IOHandler io;

    SolutionCompressor* compressor;

    std::atomic<unsigned long long> clocks_sorting;
    std::atomic<unsigned long long> clocks_generating;
    std::atomic<unsigned long long> clocks_nlogn;
    std::atomic<unsigned long long> clocks_copying;
    std::atomic<unsigned long long> clocks_generating_transforms;
    std::atomic<unsigned long long> clocks_sorting_transforms;
    std::atomic<unsigned long long> clocks_transforming;
    std::atomic<unsigned long long> clocks_merging;

    // Used for comparing binary vs. non-binary
    unsigned long max_policies_considered;
    unsigned long total_policies_considered;
    unsigned long total_pruned_policies;
    unsigned long total_policies_compared;
};

#endif //AMAZON_PROJECT_DP_ALGORITHM_H
