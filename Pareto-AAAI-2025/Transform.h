//
// Created by marc on 11/16/21.
//

#ifndef DP_PARALLEL_GENERATION_TRANSFORM_H
#define DP_PARALLEL_GENERATION_TRANSFORM_H

#include "Network.h"
#include "Solution.h"
#include <unordered_set>
#include <cstdlib>
#include <atomic>
#include "DivideAndConquer.h"
#include "Config.h"

class Transform {
public:
    Transform(TreeNode* parentNode, TreeNode* inner_node, Solution* outer_solution, const Dam &outer_dam, const Dam &inner_dam,
              double outer_dec, double inner_dec, int num_criteria, double k_epsilon,
              const std::vector<std::vector<double> > &w, int batch_size, Network *net, DivideAndConquer *dnc, bool verbose,
              Config config);

    Transform(TreeNode* parentNode, TreeNode* inner_node, const Dam &inner_dam,
              double inner_dec, int num_criteria, double k_epsilon, const std::vector<std::vector<double> > &w,
              int batch_size, Network *net, DivideAndConquer *dnc, bool verbose, Config config);

    bool add_solution(Solution* adding, std::unordered_set<Solution *, SolutionHash, EqualSolutions> &solution_set, int &num_generated, std::atomic<unsigned int> &batch_num);

    std::vector<Solution *> get_grouping(int j);
    void generate_criteria(std::pair<double, double> *criteria, double *criteria_all, int j);
    void initialize_criteria(std::pair<double, double> *criteria, double* criteria_all);
    void update_criteria(Solution *inner_solution, std::pair<double, double> *criteria, double* criteria_all);
    void round_criteria(std::pair<double, double> *criteria);
    double round_with_k(double value, double k);

    double base_criteria_minmax[MAX_CRITERIA];
    bool inferior;

    double k_vals[MAX_CRITERIA]{};

    double outer_dec;
    double inner_dec;

    std::vector<double> decisionVec;

    TreeNode* parentNode;
    TreeNode* inner_node;
    Solution* outer_solution;
private:
    std::pair<double, double> base_criteria[MAX_CRITERIA];
    double base_criteria_all[MAX_CRITERIA];
    double built_criteria[MAX_CRITERIA]{};
    Dam outer_dam;
    Dam inner_dam;
    int num_criteria;
    double k_epsilon;
    Network *net;
    DivideAndConquer *dnc;
    std::vector<std::vector<double> > w;
    bool verbose;
    int batch_size;
    Config config;
};

struct Compare_Transforms_Ordering {
    int dimension;
    int total_dimensions;
//    Compare_Dimensions_Reverse cd;
    Compare_Dimensions cd;

    explicit Compare_Transforms_Ordering(int dimension, int total_dimensions, bool use_linear_preferences) : cd(dimension, total_dimensions, use_linear_preferences) {
        this->dimension = dimension;
        this->total_dimensions = total_dimensions;
    }

    bool operator () (Transform* t1, Transform* t2) const {
        if (t1->outer_dec != t2->outer_dec || t1->inner_dec != t2->inner_dec) {
            // Calculate an "ordering" - put outer first
            if (t1->outer_dec == -1 || t2->outer_dec == -1) {
                return t1->inner_dec < t2->inner_dec;
            }
            else {
                return (2 * t1->inner_dec + t1->outer_dec) < (2 * t2->inner_dec + t2->outer_dec);
            }
        }
        else {
            return this->cd(t1->outer_solution, t2->outer_solution);
        }
    }
};

struct Compare_Transforms_New {
    explicit Compare_Transforms_New(int dimension, int total_dimensions) {
        this->dimension = dimension;
        this->total_dimensions = total_dimensions;
    }
    bool operator () (Transform* t1, Transform* t2) const {
        // Because arrays are 0 based index we subtract from the dimension
        for (int i = dimension - 1; i >= 0; i--) {
            // If the solutions are equal in the given dimension, try next dimension
            if (t1->base_criteria_minmax[i] == t2->base_criteria_minmax[i]) {
                continue;
            }
            return t1->base_criteria_minmax[i] > t2->base_criteria_minmax[i]; // Sort in descending order
        }
        // We will never have solutions equal in all dimensions so now we compare across all dimensions
        // Resolve tie in all lower dimensions by comparing upper dimensions
        int tie = dimension;
        while (tie <= total_dimensions - 1) {
            if (t1->base_criteria_minmax[tie] == t2->base_criteria_minmax[tie]) {
                tie++;
            } else {
                return t1->base_criteria_minmax[tie] > t2->base_criteria_minmax[tie];
            }
        }
        return false;
    }

    int dimension;
    int total_dimensions;
};

#endif //DP_PARALLEL_GENERATION_TRANSFORM_H
