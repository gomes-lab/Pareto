//
// Created by marc on 10/24/21.
//

#include "DPOriginal.h"

DPOriginal::DPOriginal(Network *net, int num_criteria, const Config &config) : DPBase(net, num_criteria, config) {
}

/**
 * Main DP step that determines the node type (leaf, one child, two children, three+ children, etc.) and executes the
 * appropriate sub-function for that node type.
 *
 * @param node
 */
void DPOriginal::computeDP(TreeNode* node) {
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

    if (this->config.verbose) {
        for (int i = 0; i < num_criteria; i++) {
            std::cout << "Level " << i << " sorts " << dnc->sorts[i] << " divides " << dnc->divides[i] << " marries " << dnc->marries[i] << std::endl;
            dnc->sorts[i] = 0;
            dnc->divides[i] = 0;
            dnc->marries[i] = 0;
        }

        std::cout << "End - kept: " << node->vec_frontier.size() << "\tdropped: " <<
                  (num_considered - node->vec_frontier.size()) << std::endl << std::endl;
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
void DPOriginal::computeLeaf(TreeNode* node) {
    // List of criteria values
    double criteria_all[MAX_CRITERIA];
    std::pair<double, double> criteria[MAX_CRITERIA];

    initialize_criteria(node, criteria, criteria_all);

    node->vec_frontier.push_back(new Solution(node->node_data.id, node->node_id, criteria_all, criteria, {}, {},
                                              (config.dcip_index >= 0 ? net->num_criteria + 1 : net->num_criteria),
                                              num_criteria, config.use_linear_preferences, config.w,
                                              net->minmax, net->normalizing_factor));

    node->max_frontier_size = 0;
    if (this->config.verbose) {
        std::cout << "Leaf node - one solution generated" << std::endl << std::endl;
    }
}

void DPOriginal::computeDPOneChild(TreeNode* node, unsigned long &left_considered, unsigned long &num_considered,
                                           std::unordered_set<Solution*, SolutionHash, EqualSolutions> &all_solution_set,
                                           clock_t &t) {
    std::pair<TreeNode*, TreeNode*> child_pair = node->get_children();
    // Track policies
    int num_generated = 0;
    TreeNode *left = child_pair.first;
    left_considered = left->vec_frontier.size() * left->parentDam.decision.size();
    num_considered = left_considered;
    if (this->config.verbose) {
        std::cout << "Considering one child: " << num_considered << std::endl;
    }
    unsigned long num_compared;
    int local_batches = 0;
    std::chrono::duration<float> local_duration(0);
    // Loop through left dam
    for (double left_dam_plan: left->parentDam.decision) {
        std::vector<double> decision = {left_dam_plan, -1};

        for (Solution *currL: left->vec_frontier) {
            std::vector<Solution *> pair = {currL, nullptr};

            double criteria_all[MAX_CRITERIA];
            std::pair<double, double> criteria[MAX_CRITERIA];
            generate_criteria(criteria, criteria_all, pair, decision, node);
            Solution* soln = generateSolution(node->node_data.id, node->node_id, criteria_all, criteria, decision, pair);

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

    // Time nlogn
    t = clock();
    // Perform the divide and conquer
    node->vec_frontier = dnc->divide_and_conquer(solution_set, num_criteria, false, num_comparisons);
    solution_set.clear();
    t = clock() - t;
    clocks_nlogn += (unsigned long long) t;
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
void DPOriginal::computeDPTwoChild(TreeNode* node, unsigned long &left_considered,
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
    if (right->vec_frontier.size() > 2 * config.num_threads) {
        int n = (int) right->vec_frontier.size() / config.num_threads;

        unsigned long total_ties = 0;

        for (int i = 0; i < config.num_threads; ++i) {
            results.emplace_back(pool.enqueue([i, this, node, right, left, all_solution_set, n, &total_compared, &total_ties] {
                                                  // mtx.lock();
                                                  auto t = clock();
                                                  // printf("[%d] Start\n",i );
                                                  auto it1 = right->vec_frontier.begin();

                                                  unsigned long num_compared = 0;
                                                  unsigned long local_ties = 0;
                                                  std::advance(it1, i * n);

                                                  auto it2 = right->vec_frontier.begin();
                                                  if (i == config.num_threads - 1) {
                                                      it2 = right->vec_frontier.end();
                                                  } else {
                                                      std::advance(it2, (i + 1) * n);
                                                  }
                                                  std::vector<Solution *> partial_right(it1, it2);

                                                  int num_generated = 0;
                                                  int local_batches = 0;
                                                  std::chrono::duration<float> local_duration(0);

                                                  for (double right_dam_plan: right->parentDam.decision) {
                                                      for (double left_dam_plan: left->parentDam.decision) {
                                                          // Create decision vector
                                                          std::vector<double> decision = {left_dam_plan, right_dam_plan};
                                                          // Take all combinations of the left and right tables

                                                          for (Solution *currL: left->vec_frontier) {
                                                              for (Solution *currR: partial_right) {
                                                                  // Create pair

                                                                  std::vector<Solution *> pair = {currL, currR};

                                                                  double criteria_all[MAX_CRITERIA];
                                                                  std::pair<double, double> criteria[MAX_CRITERIA];
                                                                  generate_criteria(criteria, criteria_all, pair, decision, node);
                                                                  Solution* soln = generateSolution(node->node_data.id, node->node_id, criteria_all, criteria, decision, pair);

                                                                  if (!checkSolution(soln, node)) {
                                                                      continue;
                                                                  }

                                                                  // Call method to generate new solution
                                                                  if (!add_solution(soln, solution_sets[i], num_generated, t,
                                                                               this->batch_number, i, num_compared, local_batches, local_duration)) {
                                                                      local_ties++;
                                                                  }

                                                              }
                                                          }

                                                      }
                                                  }

                                                  num_compared += solution_sets[i].size();

                                                  std::chrono::high_resolution_clock::time_point startb = std::chrono::high_resolution_clock::now();
                                                  unsigned long size = solution_sets[i].size();

                                                  t = clock();
                                                  std::vector<Solution* > batched_solutions(solution_sets[i].begin(), solution_sets[i].end());
                                                  t = clock() - t;
                                                  clocks_copying += (unsigned long long)t;
                                                  t = clock();
                                                  std::vector<Solution* > temp_answers = dnc->divide_and_conquer(batched_solutions, num_criteria, false, num_comparisons);
                                                  batched_solutions.clear();
                                                  t = clock() - t;
                                                  clocks_nlogn += (unsigned long long)t;
                                                  t = clock();
                                                  solution_sets[i] = std::unordered_set<Solution*, SolutionHash, EqualSolutions> (temp_answers.begin(), temp_answers.end());
                                                  temp_answers.clear();
                                                  t = clock() - t;
                                                  std::chrono::high_resolution_clock::time_point stopb = std::chrono::high_resolution_clock::now();


                                                  std::chrono::duration<float> durationb = stopb - startb;
                                                  local_batches++;
                                                  local_duration = local_duration + durationb;

                                                  total_compared += num_compared;
                                                  if (config.verbose) {
                                                      std::cout << "(" << i << ") Num compared " << num_compared
                                                                << std::endl;

                                                      std::cout << "(" << i << ") Batch " << this->batch_number++ << " (" << size << " to " << solution_sets[i].size() << ") took "
                                                                << durationb.count() << std::endl;

                                                      std::cout << "(" << i << ") Ran " << local_batches << " batches, taking "
                                                                << local_duration.count() << " (" << (local_duration.count() / (float) local_batches) << " per batch)" << std::endl;

                                                      std::cout << "(" << i << ") Ties: " << local_ties << std::endl;
                                                  }

                                                  total_ties += local_ties;

                                                  return i;
                                              }
            ));
        }
        for (auto &&result: results) {
            int i = result.get();

        }
        results.clear();

        if (this->config.verbose) {
            std::cout << "Total compared " << total_compared << std::endl;
            std::cout << "Total ties " << total_ties << std::endl;
        }

        unsigned long merge_ties = 0;

        for (int i = 0; i < config.num_threads; ++i) {
            for (auto &soln : solution_sets[i]) {
                auto find = all_solution_set.find(soln);
                if (find != all_solution_set.end()) { // Solution tie exists!
                    int rand_num = rand() % 2; // Only want to do if epsilon != 0 i.e. we are rounding!!!!

                    if (rand_num == 0) { // Remove the current and replace it with the new node
                        Solution *old = *find;

                        all_solution_set.erase(find);
                        delete (old);
                        all_solution_set.insert(soln);
                    } else {
                        delete (soln);
                    }
                    merge_ties++;
                } else {
                    all_solution_set.insert(soln);
                }
            }
            solution_sets[i].clear();
        }

        std::cout << "Merge ties " << merge_ties << std::endl;
    } else {
        int num_generated = 0;
        int local_batches = 0;
        std::chrono::duration<float> local_duration(0);

        unsigned long ties = 0;

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

                        if (!checkSolution(soln, node)) {
                            continue;
                        }

                        // Call method to generate new solution
                        if (!add_solution(soln, all_solution_set, num_generated, t,
                                     batch_number, -1, total_compared, local_batches, local_duration)) {
                            ties++;
                        }
                    }
                }
            }
        }

        std::cout << "Total ties " << ties << std::endl;
    }

    stop = std::chrono::high_resolution_clock::now();
    duration = stop - start;
    if (this->config.verbose) {
        std::cout << "Batching and generating " << duration.count() << std::endl;
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

    start = std::chrono::high_resolution_clock::now();
    // Time nlogn
    t = clock();
    // Perform the divide and conquer
    node->vec_frontier = dnc->divide_and_conquer(solution_set, num_criteria, false, num_comparisons);
    solution_set.clear();
    t = clock() - t;
    clocks_nlogn += (unsigned long long) t;
    stop = std::chrono::high_resolution_clock::now();
    duration = stop - start;
    if (this->config.verbose) {
        std::cout << "Final D&C " << duration.count() << std::endl;
    }
}

void DPOriginal::generate_criteria(std::pair<double, double> *criteria, double *criteria_all, std::vector<Solution *> &grouping,
                       std::vector<double> &decisionVec, TreeNode *node) {
    double k_vals[MAX_CRITERIA];
    initialize_criteria(node, criteria, criteria_all);

    std::array<double, MAX_CRITERIA> criteria_built{};
    for (int i = 0; i < MAX_CRITERIA; i++) {
        criteria_built[i] = 0;
    }

    // Add the corresponding values for each of the Pareto pairs from the children
    for (int i = 0; i < decisionVec.size(); i++) {
        if (decisionVec[i] != -1) {
            update_criteria(node, grouping[i], criteria, criteria_all, decisionVec[i], i, criteria_built);

            // If we have the DCIP, update it specially, and add it to the end of the criteria_all
            if (config.dcip_index >= 0) {
                criteria[config.dcip_index].first +=
                        decisionVec[i] * pow(grouping[i]->criteria_2[config.connectivity_index].first, 2.0);

                criteria[config.dcip_index].second += decisionVec[i] * pow(grouping[i]->criteria_2[config.connectivity_index].second, 2.0);
            }
        }
    }

    if (config.dcip_index >= 0) {
        criteria_all[net->num_criteria] = criteria[config.dcip_index].first;
    }

    // Later we should not update this everytime!
    assign_node_k(node, criteria_built, k_vals);

    if (config.dcip_index >= 0) {
        k_vals[config.dcip_index] = net->dcip_k;
    }

    round_criteria(node, criteria, k_vals); // Round the criteria values based on local ks
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
DPOriginal::add_solution(Solution* adding, std::unordered_set<Solution *, SolutionHash, EqualSolutions> &solution_set,
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

    // Lets do a trick to dynamically prune
    // Process a batch
    if (config.batch_size != 0 && num_generated >= config.batch_size) {
        std::chrono::high_resolution_clock::time_point startb = std::chrono::high_resolution_clock::now();
        unsigned long size = solution_set.size();

        t = clock() - t;
        clocks_generating += (unsigned long long)t;
        t = clock();
        std::vector<Solution* > batched_solutions(solution_set.begin(), solution_set.end());
        t = clock() - t;
        clocks_copying += (unsigned long long)t;
        t = clock();
        std::vector<Solution* > temp_answers = dnc->divide_and_conquer(batched_solutions, num_criteria, false, num_comparisons);
        batched_solutions.clear();
        t = clock() - t;
        clocks_nlogn += (unsigned long long)t;
        t = clock();
        solution_set = std::unordered_set<Solution*, SolutionHash, EqualSolutions> (temp_answers.begin(), temp_answers.end());
        temp_answers.clear();
        t = clock() - t;
        clocks_copying += (unsigned long long)t;
        t = clock();
        num_generated = 0;

        std::chrono::high_resolution_clock::time_point stopb = std::chrono::high_resolution_clock::now();
        std::chrono::duration<float> durationb = stopb - startb;

        if (this->config.verbose) {
            if (thread_id != -1) {
                std::cout << "(" << thread_id << ") ";
            }
            std::cout << "Batch " << batch_num++ << " (" << size << " to " << solution_set.size() << ") took "
                      << durationb.count() << std::endl;
        }
        num_compared += size;
        local_batches++;
        local_duration = local_duration + durationb;
    }

    return added_new;
}


/**
 * Initialize the criteria array from a node
 *
 * @param node
 * @param criteria
 */
void DPOriginal::initialize_criteria(TreeNode *node, std::pair<double, double> *criteria, double *criteria_all) {
    for (int i = 0; i < num_criteria; i++) {
        criteria[i] = node->node_data.r_vals[i];
    }

    for (int i = 0; i < net->num_criteria; i++) {
        criteria_all[i] = node->node_data.r_all[i];
    }
}


/**
 * Update the criteria array for a node given a decision made on a child node and solution
 *
 * @param node
 * @param solution_decision
 * @param criteria
 * @param decision
 * @param child_index
 * @param criteria_built
 */
void DPOriginal::update_criteria(TreeNode *node, Solution *solution_decision, std::pair<double, double> *criteria,
                                 double *criteria_all, double decision, int child_index, std::array<double, MAX_CRITERIA> &criteria_built) {
    Dam dam;

    dam = node->children[child_index]->parentDam;

    for (int i = 0; i < num_criteria; i++) {
        criteria[i].first += decision * (dam.s_vals[i].first + solution_decision->criteria_2[i].first * dam.p_vals[i]) +
                             (1 - decision) * solution_decision->criteria_2[i].first * dam.q_vals[i];

        double solution_crit = std::abs(solution_decision->criteria_2[i].second);

        criteria[i].second += decision * (dam.s_vals[i].second + solution_crit * dam.p_vals[i]) +
                              (1 - decision) * solution_crit * dam.q_vals[i];

        criteria_built[i] += decision * dam.s_vals[i].first;
    }

    for (int i = 0; i < net->num_criteria; i++) {
        criteria_all[i] += decision * (dam.s_all[i] + solution_decision->criteria[i] * dam.p_all[i]) +
                (1 - decision) * solution_decision->criteria[i] * dam.q_all[i];
    }
}


/**
 * Assign a node's k values based on the node's values and the built dams
 *
 * @param node
 * @param built_criteria
 */
void DPOriginal::assign_node_k(TreeNode *node, std::array<double, MAX_CRITERIA> &built_criteria, double *k_vals) {
    for (int i = 0; i < num_criteria; ++i) {
        k_vals[i] = node->node_data.r_vals[i].first * config.epsilon + std::floor((built_criteria[i] + std::numeric_limits<double>::epsilon()) / net->min_vals[i]) * net->min_vals[i] * config.epsilon / 2.0;
    }
}


/**
 * Round the criteria based on the node's k values
 *
 * @param node
 * @param criteria
 */
void DPOriginal::round_criteria(TreeNode *node, std::pair<double, double> *criteria, double *k_vals) {
    for (int i = 0; i < num_criteria; i++) {
        criteria[i].second = net->minmax[i] * round_with_k(criteria[i].second, k_vals[i]);
    }
}


/**
 * Round a value given a k value
 *
 * @param value
 * @param k
 * @return
 */
double DPOriginal::round_with_k(double value, double k) {
    if (k == 0) {
        return value;
    }

    return std::floor((value + std::numeric_limits<double>::epsilon()) / k) * k;
}