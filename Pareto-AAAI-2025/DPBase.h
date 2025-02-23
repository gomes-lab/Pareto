//
// Created by marc on 10/24/21.
//

#ifndef DP_PARALLEL_GENERATION_DPBASE_H
#define DP_PARALLEL_GENERATION_DPBASE_H

#include "Network.h"
#include "Solution.h"
#include "Transform.h"
#include "DivideAndConquer.h"
#include "ThreadPool.h"
#include "Config.h"
#include <ctime>
#include <unordered_set>
#include <cstdlib>
#include <bitset>

class DPBase {
public:
    DPBase(Network *net, int num_criteria, const Config &config);
    ~DPBase();
    virtual void computeDP(TreeNode* node);

    Solution* generateSolution(int node_id, int tree_node_id, const double *criteria_all,
                               const std::pair<double, double> *criteria, std::vector<double> dam_decisions,
                               std::vector<Solution *> pareto_decisions);
    bool checkSolution(Solution* soln, TreeNode* node);

    unsigned long total_policies_considered;
    unsigned long total_policies_compared;
    unsigned long total_pruned_policies;

    unsigned long num_comparisons;

    unsigned long max_policies_considered;

    unsigned long clocks_generating;
    unsigned long clocks_generating_transforms;
    unsigned long clocks_transforming;
    unsigned long clocks_merging;
    unsigned long clocks_copying;
    unsigned long clocks_nlogn;
    unsigned long clocks_sorting;
    unsigned long clocks_sorting_transforms;

    unsigned long transforms_considered;
    unsigned long transforms_pruned;
    unsigned long policies_avoided;

    std::vector<std::vector<double>> ks;
    DivideAndConquer *dnc;
protected:
    std::atomic<unsigned int> batch_number;
    int num_criteria;

    Network *net;

    ThreadPool pool;
    std::vector< std::future<int> > results;
    std::vector<std::unordered_set<Solution *, SolutionHash, EqualSolutions> >  solution_sets;

    Config config;
};


#endif //DP_PARALLEL_GENERATION_DPBASE_H
