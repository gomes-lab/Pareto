//
// Created by marc on 10/24/21.
//

#include "DPBase.h"

DPBase::DPBase(Network *net, int num_criteria, const Config &config) :
        net(net), num_criteria(num_criteria), config(config),
        total_policies_compared(0), total_policies_considered(0), total_pruned_policies(0), max_policies_considered(0),
        clocks_generating_transforms(0), clocks_generating(0), clocks_copying(0), clocks_nlogn(0), clocks_sorting(0),
        clocks_sorting_transforms(0), clocks_transforming(0), clocks_merging(0), num_comparisons(0), batch_number(1),
        transforms_considered(0), transforms_pruned(0), policies_avoided(0), pool(config.num_threads)  {
    dnc = new DivideAndConquer(num_criteria, config.use_linear_preferences);

    for (int i = 0; i < config.num_threads; ++i){
        solution_sets.emplace_back(std::unordered_set<Solution *, SolutionHash, EqualSolutions>());
    }
}

void DPBase::computeDP(TreeNode *node) {
    return;
}

DPBase::~DPBase() {
    delete(dnc);
}

Solution* DPBase::generateSolution(int node_id, int tree_node_id, const double *criteria_all,
                                   const std::pair<double, double> *criteria, std::vector<double> dam_decisions,
                                   std::vector<Solution *> pareto_decisions) {
    auto soln = new Solution(node_id, tree_node_id, criteria_all, criteria, dam_decisions, pareto_decisions,
                             (config.dcip_index >= 0 ? net->num_criteria + 1 : net->num_criteria),
                             num_criteria, config.use_linear_preferences, config.w, net->minmax, net->normalizing_factor);
    if (config.compress_distance_mode == COMPRESS_DISTANCE_DECISION_SPACE) {
        soln->generate_representation();
    }

    return soln;
}

bool DPBase::checkSolution(Solution* soln, TreeNode* node) {
    bool valid_solution = true;

    // Check upper bounds
    double min_forced;
    if (config.upper_bounds.size() > 0) {
        for (int i = 0; i < net->num_criteria; i++) {
            if (config.upper_bounds[i] < 0) {
                continue;
            }

            min_forced = net->min_forced[i] - net->hyper_nodes[soln->tree_node_id]->min_forced[i];
            if (soln->criteria[i] > 0 && soln->criteria[i] + min_forced > config.upper_bounds[i]) {
                node->num_ub_broken++;
                valid_solution = false;
            } else if (soln->criteria[i] < 0 && soln->criteria[i] - min_forced < -config.upper_bounds[i]) {
                node->num_ub_broken++;
                valid_solution = false;
            }
        }
    }

    // Check lower bounds
    double max_possible;
    if (config.lower_bounds.size() > 0) {
        for (int i = 0; i < net->num_criteria; i++) {
            if (config.lower_bounds[i] < 0) {
                continue;
            }

            max_possible = net->max_possible[i] - net->hyper_nodes[soln->tree_node_id]->max_possible[i];
            if (soln->criteria[i] > 0 && soln->criteria[i] + max_possible < config.lower_bounds[i]) {
                node->num_lb_broken++;
                valid_solution = false;
            }
            else if (soln->criteria[i] < 0 && soln->criteria[i] - max_possible > -config.lower_bounds[i]) {
                node->num_lb_broken++;
                valid_solution = false;
            }
        }
    }

    // Check mutually exclusive dams

    if (!valid_solution) {
        delete(soln);
    }
    
    return valid_solution;
}