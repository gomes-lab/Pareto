//
// Created by marc on 10/24/21.
//

#ifndef DP_PARALLEL_GENERATION_DPBRUTEFORCE_H
#define DP_PARALLEL_GENERATION_DPBRUTEFORCE_H

#include "Network.h"
#include "Solution.h"
#include "Transform.h"
#include "DivideAndConquer.h"
#include <ctime>
#include <unordered_set>
#include <cstdlib>
#include "DPOriginal.h"
#include "Config.h"

class DPBruteForce : public DPOriginal {
public:
    DPBruteForce(Network *net, int num_criteria, const Config &config);
    virtual void computeDP(TreeNode* node);
    void divide_and_conquer(TreeNode* node);
private:
    void computeLeaf(TreeNode* node);
    void computeDPOneChild(TreeNode* node, unsigned long &left_considered, unsigned long &num_considered,
                           std::unordered_set<Solution*, SolutionHash, EqualSolutions> &all_solution_set,
                           clock_t &t);
    void computeDPTwoChild(TreeNode* node, unsigned long &left_considered,
                           unsigned long &right_considered, unsigned long &num_considered,
                           std::unordered_set<Solution*, SolutionHash, EqualSolutions> &all_solution_set,
                           clock_t &t);
    bool add_solution(Solution* adding, std::unordered_set<Solution *, SolutionHash, EqualSolutions> &solution_set,
                      int &num_generated, clock_t &t, std::atomic<unsigned int> &batch_num,
                      int thread_id, unsigned long &num_compared, int &local_batches,
                      std::chrono::duration<float> &local_duration);
};


#endif //DP_PARALLEL_GENERATION_DPORIGINAL_H
