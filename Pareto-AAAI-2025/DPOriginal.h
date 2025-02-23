//
// Created by marc on 10/24/21.
//

#ifndef DP_PARALLEL_GENERATION_DPORIGINAL_H
#define DP_PARALLEL_GENERATION_DPORIGINAL_H

#include "Network.h"
#include "Solution.h"
#include "Transform.h"
#include "DivideAndConquer.h"
#include "ThreadPool.h"
#include <ctime>
#include <unordered_set>
#include <cstdlib>
#include "DPBase.h"
#include "Config.h"

class DPOriginal : public DPBase {
public:
    DPOriginal(Network *net, int num_criteria, const Config &config);
    virtual void computeDP(TreeNode* node);
protected:
    void computeDPTwoChild(TreeNode* node, unsigned long &left_considered,
                           unsigned long &right_considered, unsigned long &num_considered,
                           std::unordered_set<Solution*, SolutionHash, EqualSolutions> &all_solution_set,
                           clock_t &t);

    void computeDPOneChild(TreeNode* node, unsigned long &left_considered, unsigned long &num_considered,
                           std::unordered_set<Solution*, SolutionHash, EqualSolutions> &all_solution_set,
                           clock_t &t);
    void computeLeaf(TreeNode* node);

    void generate_criteria(std::pair<double, double> *criteria, double *criteria_all, std::vector<Solution *> &grouping,
                           std::vector<double> &decisionVec, TreeNode *node);
    bool add_solution(Solution* adding, std::unordered_set<Solution *, SolutionHash, EqualSolutions> &solution_set,
                      int &num_generated, clock_t &t, std::atomic<unsigned int> &batch_num,
                      int thread_id, unsigned long &num_compared, int &local_batches,
                      std::chrono::duration<float> &local_duration);
    void initialize_criteria(TreeNode *node, std::pair<double, double> *criteria, double *criteria_all);
    void update_criteria(TreeNode *node, Solution *solution_decision, std::pair<double, double> *criteria,
                         double *criteria_all, double decision, int child_index, std::array<double, MAX_CRITERIA> &criteria_built);
    void assign_node_k(TreeNode *node, std::array<double, MAX_CRITERIA> &built_criteria, double* k_vals);
    void round_criteria(TreeNode *node, std::pair<double, double> *criteria, double *k_vals);
    double round_with_k(double value, double k);
};


#endif //DP_PARALLEL_GENERATION_DPORIGINAL_H
