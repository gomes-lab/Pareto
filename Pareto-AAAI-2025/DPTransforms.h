//
// Created by marc on 10/24/21.
//

#ifndef DP_PARALLEL_GENERATION_DPTransforms_H
#define DP_PARALLEL_GENERATION_DPTransforms_H

#include "Network.h"
#include "Solution.h"
#include "Transform.h"
#include "DivideAndConquer.h"
#include "ThreadPool.h"
#include <ctime>
#include <unordered_set>
#include <cstdlib>
#include "DPBase.h"
#include "Transform.h"
#include "Config.h"

class DPTransforms : public DPBase {
public:
    DPTransforms(Network *net, int num_criteria, const Config &config);
    virtual void computeDP(TreeNode* node);
private:
    void computeDPTwoChild(TreeNode* node, unsigned long &left_considered,
                           unsigned long &right_considered, unsigned long &num_considered,
                           std::unordered_set<Solution*, SolutionHash, EqualSolutions> &all_solution_set,
                           clock_t &t);

    void computeDPOneChild(TreeNode* node, unsigned long &left_considered, unsigned long &num_considered,
                           std::unordered_set<Solution*, SolutionHash, EqualSolutions> &all_solution_set,
                           clock_t &t);
    void computeLeaf(TreeNode* node);

    std::vector<Transform*> divide_and_conquer_transforms(std::vector<Transform*> &curr_set, int dimension);
    std::vector<Transform*> divide_and_conquer_transforms_helper(std::vector<Transform*> &curr_set, unsigned long low,
                                                                 unsigned long high, int dimension);
    std::vector<Transform*> L2DTransforms(std::vector<Transform*>& all_possible);

    std::vector<Transform*> marry_transform(std::vector<Transform*> &curr_set, unsigned long low,
                                            unsigned long high, int dimension);
    std::vector<Transform*> marry_2d_transform(std::vector<Transform*> &curr_set);
};


#endif //DP_PARALLEL_GENERATION_DPTransforms_H
