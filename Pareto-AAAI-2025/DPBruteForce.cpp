//
// Created by marc on 10/24/21.
//

#include "DPBruteForce.h"

DPBruteForce::DPBruteForce(Network *net, int num_criteria, const Config &config) : DPOriginal(net, num_criteria, config) {
}

void DPBruteForce::divide_and_conquer(TreeNode* node) {
    std::vector<Solution *> solution_set(node->vec_frontier.begin(), node->vec_frontier.end());
    node->vec_frontier = dnc->divide_and_conquer(solution_set, num_criteria, false, num_comparisons);
}

/**
 * Main DP step that determines the node type (leaf, one child, two children, three+ children, etc.) and executes the
 * appropriate sub-function for that node type.
 *
 * @param node
 */
void DPBruteForce::computeDP(TreeNode* node) {
    clock_t total_t = clock();
    unsigned long num_considered = 0;
    unsigned long left_considered = 0;
    unsigned long right_considered = 0;
    unsigned long num_compared = 0;
    this->batch_number = 1;
    // Check if we are a leaf node
    if (node->is_leaf) {
        num_considered = 1;
        computeLeaf(node);
    } else { // Not leaf node
        clock_t t = clock();

        std::unordered_set<Solution *, SolutionHash, EqualSolutions> all_solution_set;

        std::pair<TreeNode *, TreeNode *> child_pair = node->get_children();
        if (child_pair.first != nullptr && child_pair.second != nullptr) { // Has both children
            // Use transformations to reduce policies to search by finding dominated transformations
            computeDPTwoChild(node, left_considered, right_considered, num_considered, all_solution_set, t);
        } else if (child_pair.first != nullptr) { // Just left
            computeDPOneChild(node, left_considered, num_considered, all_solution_set, t);
        }
        t = clock() - t;
        clocks_generating += (unsigned long long) t;
    }

    std::sort(node->vec_frontier.begin(), node->vec_frontier.end(), Compare_Dimensions(num_criteria, num_criteria, config.use_linear_preferences));

    total_policies_considered += num_considered;
    total_policies_compared += num_compared;
    max_policies_considered = std::max(max_policies_considered, num_considered);
    total_pruned_policies += num_considered - node->vec_frontier.size();
    node->max_frontier_size = num_considered;

    total_t = clock() - total_t;
    node->clocks_spent += (unsigned long long)total_t;
    node->time_spent = (double)node->clocks_spent/CLOCKS_PER_SEC;
}

/**
 * Compute the solution for a leaf - only one solution containing the criteria values
 *
 * @param node
 */
void DPBruteForce::computeLeaf(TreeNode* node) {
    // List of criteria values
    double criteria_all[MAX_CRITERIA];
    std::pair<double, double> criteria[MAX_CRITERIA];

    initialize_criteria(node, criteria, criteria_all);

    node->vec_frontier.push_back(new Solution(node->node_data.id, node->node_id, criteria_all, criteria, {}, {}, net->num_criteria,
                                              num_criteria, config.use_linear_preferences, config.w, net->minmax, net->normalizing_factor));

    node->max_frontier_size = 0;
    if (this->config.verbose) {
        std::cout << "Leaf node - one solution generated" << std::endl << std::endl;
    }
}

void DPBruteForce::computeDPOneChild(TreeNode* node, unsigned long &left_considered, unsigned long &num_considered,
                                     std::unordered_set<Solution*, SolutionHash, EqualSolutions> &all_solution_set,
                                     clock_t &t) {
    std::pair<TreeNode*, TreeNode*> child_pair = node->get_children();
    // Track policies
    int num_generated = 0;
    TreeNode *left = child_pair.first;
    left_considered = left->vec_frontier.size() * left->parentDam.decision.size();
    num_considered = left_considered;
    if (this->config.verbose) {
        std::cout << "Generating solutions for one child: " << num_considered << std::endl;
    }
    unsigned long num_compared;
    int local_batches = 0;
    unsigned long total_generated = 0;
    std::chrono::duration<float> local_duration(0);
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    // Loop through left dam
    for (double left_dam_plan: left->parentDam.decision) {
        std::vector<double> decision = {left_dam_plan, -1};

        for (Solution *currL: left->vec_frontier) {
            std::vector<Solution *> pair = {currL, nullptr};

            double criteria_all[MAX_CRITERIA];
            std::pair<double, double> criteria[MAX_CRITERIA];
            generate_criteria(criteria, criteria_all, pair, decision, node);
            Solution* soln = generateSolution(node->node_data.id, node->node_id, criteria_all, criteria, decision, pair);

            total_generated++;

            if (total_generated % 1000000 == 0) {
                std::chrono::high_resolution_clock::time_point stop = std::chrono::high_resolution_clock::now();

                std::chrono::duration<float> duration = stop - start;
                std::cout << "Generated " << total_generated << " after " << duration.count() << std::endl;
            }

            if (!checkSolution(soln, node)) {
                continue;
            }

            // Call method to generate new solution
            add_solution(soln, all_solution_set, num_generated, t,
                         batch_number, -1, num_compared, local_batches, local_duration);
        }
    }

    // Copying
    t = clock();
    // Add the solutions from the hash set to a vector.
    // No ties exist
    std::vector<Solution *> solution_set(all_solution_set.begin(), all_solution_set.end());
    t = clock() - t;
    clocks_copying += (unsigned long long) t;

    node->vec_frontier = solution_set;
}

/**
 * Original code for executing the Dynamic Programming step - does not take advantage of the transform objects
 * or the fact that child sub-frontiers are already sorted.
 *
 * @param node
 * @param left_considered
 * @param right_considered
 * @param num_considered
 * @param all_solution_set
 * @param t
 */
void DPBruteForce::computeDPTwoChild(TreeNode* node, unsigned long &left_considered,
                                     unsigned long &right_considered, unsigned long &num_considered,
                                     std::unordered_set<Solution*, SolutionHash, EqualSolutions> &all_solution_set,
                                     clock_t &t)
{
    // Calculate the number that will be considered
    std::pair<TreeNode*, TreeNode*> child_pair = node->get_children();
    TreeNode *left = child_pair.first;
    TreeNode *right = child_pair.second;
    left_considered = left->vec_frontier.size() * left->parentDam.decision.size();
    right_considered = right->vec_frontier.size() * right->parentDam.decision.size();
    num_considered = left_considered * right_considered;
    std::chrono::high_resolution_clock::time_point start, stop;
    std::chrono::duration<float> duration;

    if (this->config.verbose) {
        std::cout << "Considering two children: " << num_considered << std::endl;
    }

    start = std::chrono::high_resolution_clock::now();
    unsigned long total_compared = 0;

    // Take all the combinations of the dam placements
    // Note we go *backwards* because when doing the cartesian
    // product it does it backwards

    int num_generated = 0;
    int local_batches = 0;
    std::chrono::duration<float> local_duration(0);

    unsigned long ties = 0;
    unsigned long total_generated = 0;

    for (double right_dam_plan: right->parentDam.decision) {
        for (double left_dam_plan: left->parentDam.decision) {
            // Create decision vector
            std::vector<double> decision = {left_dam_plan, right_dam_plan};
            // Take all combinations of the left and right tables
            for (Solution *currR: right->vec_frontier) {
                for (Solution *currL: left->vec_frontier) {
                    // Create pair
                    std::vector<Solution *> pair = {currL, currR};

                    double criteria_all[MAX_CRITERIA];
                    std::pair<double, double> criteria[MAX_CRITERIA];
                    generate_criteria(criteria, criteria_all, pair, decision, node);
                    Solution* soln = generateSolution(node->node_data.id, node->node_id, criteria_all, criteria, decision, pair);

                    total_generated++;

                    if (total_generated % 1000000 == 0) {
                        stop = std::chrono::high_resolution_clock::now();

                        duration = stop - start;
                        std::cout << "Generated " << total_generated << " after " << duration.count() << std::endl;
                    }

                    if (!checkSolution(soln, node)) {
                        continue;
                    }

                    // Call method to generate new solution
                    add_solution(soln, all_solution_set, num_generated, t,
                                 batch_number, -1, total_compared, local_batches, local_duration);
                }
            }
        }
    }

    stop = std::chrono::high_resolution_clock::now();
    duration = stop - start;
    if (this->config.verbose) {
        std::cout << "Generating " << duration.count() << std::endl;
    }

    // Copying
    t = clock();
    // Add the solutions from the hash set to a vector.
    // No ties exist
    if (this->config.verbose) {
        std::cout << "Copying set" << std::endl;
    }
    std::vector<Solution *> solution_set(all_solution_set.begin(), all_solution_set.end());
    t = clock() - t;
    clocks_copying += (unsigned long long) t;

    node->vec_frontier = solution_set;
}

/**
 * Add a solution to the set of Policies, if it should be added
 *
 * @param grouping
 * @param decisionVec
 * @param node
 * @param solution_set
 * @param num_generated
 * @param t
 * @param batch_num
 */
bool
DPBruteForce::add_solution(Solution* adding, std::unordered_set<Solution *, SolutionHash, EqualSolutions> &solution_set,
                       int &num_generated, clock_t &t, std::atomic<unsigned int> &batch_num,
                       int thread_id, unsigned long &num_compared, int &local_batches,
                       std::chrono::duration<float> &local_duration) {
    bool added_new;
    auto find = solution_set.find(adding);
    if (find != solution_set.end()) { // Solution tie exists!
        int rand_num = rand() % 2; // Only want to do if epsilon != 0 i.e. we are rounding!!!!

        if (rand_num == 0) { // Remove the current and replace it with the new node
            Solution* old = *find;

            solution_set.erase(find);
            delete (old);
            solution_set.insert(adding);
        } else {
            delete(adding);
        }

        added_new = false;
    } else {
        solution_set.insert(adding);
        num_generated++;

        added_new = true;
    }

    return added_new;
}