//
// Created by marc on 10/24/21.
//

#include "DPTransforms.h"

DPTransforms::DPTransforms(Network *net, int num_criteria, const Config &config) : DPBase(net, num_criteria, config) {
}

/**
 * Main DP step that determines the node type (leaf, one child, two children, three+ children, etc.) and executes the
 * appropriate sub-function for that node type.
 *
 * @param node
 */
void DPTransforms::computeDP(TreeNode *node) {
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
    node->clocks_spent += (unsigned long long) total_t;
    node->time_spent = (double) node->clocks_spent / CLOCKS_PER_SEC;
}

/**
 * Compute the solution for a leaf - only one solution containing the criteria values
 *
 * @param node
 */
void DPTransforms::computeLeaf(TreeNode *node) {
    // List of criteria values
    std::pair<double, double> criteria[MAX_CRITERIA];
    double criteria_all[MAX_CRITERIA];

    for (int i = 0; i < num_criteria; i++) {
        criteria[i] = node->node_data.r_vals[i];
    }

    for (int i = 0; i < net->num_criteria; i++) {
        criteria_all[i] = node->node_data.r_all[i];
    }

    ks.clear();
    std::vector<double> ik;
    for (int i = 0; i < num_criteria; i++) {
        ik.emplace_back(node->node_data.k_vals[i]);
    }
    ks.emplace_back(ik);

    node->vec_frontier.push_back(
            new Solution(node->node_data.id, node->node_id, criteria_all, criteria, {}, {},
                         (config.dcip_index >= 0 ? net->num_criteria + 1 : net->num_criteria),
                         num_criteria, config.use_linear_preferences, config.w, net->minmax,
                         net->normalizing_factor));

    node->max_frontier_size = 0;
    if (this->config.verbose) {
        std::cout << "Leaf node - one solution generated" << std::endl << std::endl;
    }
}

void DPTransforms::computeDPOneChild(TreeNode *node, unsigned long &left_considered, unsigned long &num_considered,
                                     std::unordered_set<Solution *, SolutionHash, EqualSolutions> &all_solution_set,
                                     clock_t &t) {
    std::pair<TreeNode *, TreeNode *> child_pair = node->get_children();
    // Track policies
    int num_generated = 0;
    TreeNode *left = child_pair.first;
    left_considered = left->vec_frontier.size() * left->parentDam.decision.size();
    num_considered = left_considered;
    if (this->config.verbose) {
        std::cout << "Considering one child: " << num_considered << std::endl;
    }
    // Loop through left dam
    ks.clear();
    for (double left_dam_plan: left->parentDam.decision) {
        std::vector<double> decision = {left_dam_plan, -1};

        Transform *transform = new Transform(node, left, left->parentDam, left_dam_plan, num_criteria,
                                             config.epsilon, config.w, config.batch_size, net, dnc, config.verbose, config);
        transforms_considered++;

        std::vector<double> ik;
        for (int i = 0; i < num_criteria; i++) {
            ik.emplace_back(transform->k_vals[i]);
        }
        ks.emplace_back(ik);

        for (int j = 0; j < left->vec_frontier.size(); j++) {
            std::pair<double, double> criteria[MAX_CRITERIA];
            double criteria_all[MAX_CRITERIA];
            transform->generate_criteria(criteria, criteria_all, j);
            Solution* soln = generateSolution(node->node_data.id, node->node_id, criteria_all,
                                              criteria, transform->decisionVec,
                                              transform->get_grouping(j));

            if (!checkSolution(soln, node)) {
                continue;
            }

            transform->add_solution(soln, all_solution_set, num_generated, batch_number);
        }

        delete(transform);
    }

    // Copying
    t = clock();
    // Add the solutions from the hash set to a vector.
    // No ties exist
    std::vector<Solution *> solution_set(all_solution_set.begin(), all_solution_set.end());
    all_solution_set.clear();
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
void DPTransforms::computeDPTwoChild(TreeNode *node, unsigned long &left_considered,
                                     unsigned long &right_considered, unsigned long &num_considered,
                                     std::unordered_set<Solution *, SolutionHash, EqualSolutions> &all_solution_set,
                                     clock_t &t) {
    // Calculate the number that will be considered
    std::pair<TreeNode *, TreeNode *> child_pair = node->get_children();
    TreeNode *left = child_pair.first;
    TreeNode *right = child_pair.second;
    left_considered = left->vec_frontier.size() * left->parentDam.decision.size();
    right_considered = right->vec_frontier.size() * right->parentDam.decision.size();
    num_considered = left_considered * right_considered;
    std::chrono::high_resolution_clock::time_point start, stop;
    std::chrono::duration<float> duration;

    t = clock();

    ks.clear();

    start = std::chrono::high_resolution_clock::now();
    std::vector<Transform*> transforms;
    for (double right_dam_plan: right->parentDam.decision) {
        std::vector<Transform*> sub_transforms;

        for (double left_dam_plan: left->parentDam.decision) {
            bool first = true;
            for (Solution *currL: left->vec_frontier) {
                Transform* transform =
                        new Transform(node, right, currL, left->parentDam,
                                         right->parentDam, left_dam_plan, right_dam_plan, num_criteria,
                                      config.epsilon, config.w, config.batch_size, net, dnc, config.verbose, config);

                sub_transforms.emplace_back(transform);

                if (first) {
                    first = false;
                    std::vector<double> ik;
                    for (int i = 0; i < num_criteria; i++) {
                        ik.emplace_back(transform->k_vals[i]);
                    }
                    ks.emplace_back(ik);
                }
            }
        }

        unsigned long sub_transforms_considered = sub_transforms.size();
        transforms_considered += sub_transforms_considered;
        sub_transforms = divide_and_conquer_transforms(sub_transforms, num_criteria);
        transforms_pruned += sub_transforms_considered - sub_transforms.size();

        transforms.insert(transforms.end(), sub_transforms.begin(), sub_transforms.end());
    }
    t = clock() - t;
    clocks_generating_transforms += (unsigned long long) t;
    stop = std::chrono::high_resolution_clock::now();
    duration = stop - start;

    t = clock();
    if (this->config.verbose) {
        std::cout << "Transform domination check " << duration.count() << std::endl;
        std::cout << "Considering two children: " << transforms.size() * right->vec_frontier.size() << " (" << transforms.size() << " transforms)" << std::endl;
    }

    std::sort(transforms.begin(), transforms.end(), Compare_Transforms_Ordering(num_criteria, num_criteria, config.use_linear_preferences));

    unsigned long total_compared = 0;
    start = std::chrono::high_resolution_clock::now();
    if (right->vec_frontier.size() > 2 * config.num_threads) {
        int n = (int) right->vec_frontier.size() / config.num_threads;
        unsigned long total_ties = 0;

        for (int i = 0; i < config.num_threads; ++i) {
            results.emplace_back(pool.enqueue([i, this, node, right, left, all_solution_set, n, transforms, &total_compared, &total_ties] {
                // mtx.lock();
                auto t = clock();
                // printf("[%d] Start\n",i );
                int start = i * n;

                int end = (i + 1) * n;
                if (i == config.num_threads - 1) {
                    end = right->vec_frontier.size();
                }

                int num_generated = 0;
                unsigned long num_compared = 0;
                int local_batches = 0;
                unsigned long local_ties = 0;

                std::chrono::high_resolution_clock::time_point startb, stopb;
                std::chrono::duration<float> duration;
                std::chrono::duration<float> local_duration(0);
                for (auto &transform: transforms) {
                    for (int j = start; j < end; j++) {
                        std::pair<double, double> criteria[MAX_CRITERIA];
                        double criteria_all[MAX_CRITERIA];
                        transform->generate_criteria(criteria, criteria_all, j);
                        Solution* soln = generateSolution(node->node_data.id, node->node_id, criteria_all,
                                                          criteria, transform->decisionVec,
                                                          transform->get_grouping(j));

                        if (!checkSolution(soln, node)) {
                            continue;
                        }

                        if (!transform->add_solution(soln, solution_sets[i], num_generated, this->batch_number)) {
                            local_ties++;
                        }

                        // Lets do a trick to dynamically prune
                        // Process a batch
                        if (config.batch_size != 0 && num_generated >= config.batch_size) {
                            startb = std::chrono::high_resolution_clock::now();

                            unsigned long num_comparisons;
                            unsigned long size = solution_sets[i].size();
                            num_compared += size;

                            t = clock() - t;
                            clocks_generating += (unsigned long long) t;
                            t = clock();
                            std::vector<Solution *> batched_solutions(solution_sets[i].begin(), solution_sets[i].end());
                            t = clock() - t;
                            clocks_copying += (unsigned long long) t;
                            t = clock();
                            std::vector<Solution *> temp_answers = dnc->divide_and_conquer(batched_solutions,
                                                                                         num_criteria, false,
                                                                                         num_comparisons);
                            batched_solutions.clear();
                            t = clock() - t;
                            clocks_nlogn += (unsigned long long) t;
                            t = clock();
                            solution_sets[i] = std::unordered_set<Solution *, SolutionHash, EqualSolutions>(
                                    temp_answers.begin(),
                                    temp_answers.end());
                            temp_answers.clear();
                            t = clock() - t;
                            clocks_copying += (unsigned long long) t;
                            t = clock();
                            num_generated = 0;
                            stopb = std::chrono::high_resolution_clock::now();
                            duration = stopb - startb;
                            local_duration = local_duration + duration;
                            local_batches++;

                            if (this->config.verbose) {
                                std::cout << "(" << i << ") Batch " << this->batch_number++ << " (" << size << " to "
                                          << solution_sets[i].size() << ") took " << duration.count() << std::endl;
                            }
                        }
                    }
                }
                t = clock() - t;
                clocks_generating += (unsigned long long)t;

                num_compared += solution_sets[i].size();

                unsigned long size = solution_sets[i].size();

                startb = std::chrono::high_resolution_clock::now();
                t = clock();
                std::vector<Solution *> batched_solutions(solution_sets[i].begin(), solution_sets[i].end());
                t = clock() - t;

                clocks_copying += (unsigned long long) t;
                t = clock();
                std::vector<Solution *> temp_answers = dnc->divide_and_conquer(batched_solutions, num_criteria, false,
                                                                           num_comparisons);
                batched_solutions.clear();
                t = clock() - t;
                clocks_nlogn += (unsigned long long) t;
                t = clock();
                solution_sets[i] = std::unordered_set<Solution *, SolutionHash, EqualSolutions>(temp_answers.begin(),
                                                                                            temp_answers.end());
                temp_answers.clear();
                t = clock() - t;
                clocks_copying += (unsigned long long) t;

                stopb = std::chrono::high_resolution_clock::now();

                total_compared += num_compared;
                duration = stopb - startb;
                local_batches++;
                local_duration = local_duration + duration;
                if (this->config.verbose) {
                    std::cout << "(" << i << ") Batch " << this->batch_number++ << " (" << size << " to "
                              << solution_sets[i].size() << ") took " << duration.count() << std::endl;

                    std::cout << "(" << i << ") Num compared " << num_compared << "(" << (float)num_compared / (float) local_batches
                              << "per batch)" << std::endl;

                    std::cout << "(" << i << ") Ran " << local_batches << " batches, taking " << local_duration.count()
                              << " (" << (local_duration.count() / (float) local_batches) << " per batch)" << std::endl;

                    std::cout << "(" << i << ") Ties: " << local_ties << std::endl;
                }

                total_ties += local_ties;

                return i;
            }));
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
        t = clock();
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
        t = clock() - t;
        clocks_copying += (unsigned long long)t;

        if (this->config.verbose) {
            std::cout << "Merge ties " << merge_ties << std::endl;
        }
    } else {
        int num_generated = 0;
        unsigned long ties = 0;

        for (auto &transform : transforms) {
            for (int j = 0; j < right->vec_frontier.size(); j++) {
                std::pair<double, double> criteria[MAX_CRITERIA];
                double criteria_all[MAX_CRITERIA];
                transform->generate_criteria(criteria, criteria_all, j);
                Solution* soln = generateSolution(node->node_data.id, node->node_id, criteria_all,
                                                  criteria, transform->decisionVec,
                                                  transform->get_grouping(j));

                if (!checkSolution(soln, node)) {
                    continue;
                }

                if (!transform->add_solution(soln, all_solution_set, num_generated, this->batch_number)) {
                    ties++;
                }

                // Lets do a trick to dynamically prune
                // Process a batch
                if (config.batch_size != 0 && num_generated >= config.batch_size) {
                    std::chrono::high_resolution_clock::time_point startb = std::chrono::high_resolution_clock::now();
                    unsigned long num_comparisons;
                    unsigned long size = all_solution_set.size();

                    t = clock() - t;
                    clocks_generating += (unsigned long long) t;
                    t = clock();
                    std::vector<Solution *> batched_solutions(all_solution_set.begin(), all_solution_set.end());
                    t = clock() - t;
                    clocks_copying += (unsigned long long) t;
                    t = clock();
                    std::vector<Solution *> temp_answers = dnc->divide_and_conquer(batched_solutions,
                                                                                 num_criteria, false,
                                                                                 num_comparisons);
                    batched_solutions.clear();
                    t = clock() - t;
                    clocks_nlogn += (unsigned long long) t;
                    t = clock();
                    all_solution_set = std::unordered_set<Solution *, SolutionHash, EqualSolutions>(
                            temp_answers.begin(),
                            temp_answers.end());
                    temp_answers.clear();
                    t = clock() - t;
                    clocks_copying += (unsigned long long) t;
                    t = clock();
                    num_generated = 0;
                    std::chrono::high_resolution_clock::time_point stopb = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<float> durationb = stopb - startb;

                    if (this->config.verbose) {
                        std::cout << "Batch " << batch_number++ << " (" << size << " to " << all_solution_set.size()
                                  << ") took " << durationb.count() << std::endl;
                    }
                }
            }
        }
        t = clock() - t;
        clocks_generating += (unsigned long long)t;

        if (this->config.verbose) {
            std::cout << "Total ties " << ties << std::endl;
        }
    }

    for (auto &transform : transforms) {
        delete(transform);
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
    all_solution_set.clear();
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

/**
 * Two dimensional sorting of Transforms that removes dominated transformations.
 *
 * @param all_possible
 * @return
 */
std::vector<Transform* > DPTransforms::L2DTransforms(std::vector<Transform* >& all_possible) {
    std::sort(all_possible.begin(), all_possible.end(), Compare_Transforms_New(2, num_criteria));

    std::vector<Transform*> transforms;
    double second_criteria_max;

    transforms.push_back(all_possible[0]);
    second_criteria_max = all_possible[0]->base_criteria_minmax[0];
    int deleted = 0;

    for (unsigned long curr_index = 1; curr_index < all_possible.size(); curr_index++) {
        // If the last criteria is larger than the max so far update and include solution
        if (all_possible[curr_index]->base_criteria_minmax[0] > second_criteria_max) {
            {
                transforms.push_back(all_possible[curr_index]);
            }
            second_criteria_max = all_possible[curr_index]->base_criteria_minmax[0];
        } else { // Free the memory for the dominated solution
            delete(all_possible[curr_index]);

            deleted++;
        }
    }

    if (this->config.verbose) {
        std::cout << "Transforms deleted: " << deleted << std::endl;
    }
    return transforms;
}

/**
 * Divide and conquer algorithm for determining Pareto Dominance of the Matrix Transforms
 *
 * @param curr_set
 * @param dimension
 * @return
 */
std::vector<Transform*> DPTransforms::divide_and_conquer_transforms(std::vector<Transform*> &curr_set, int dimension) {
    if (dimension == 2) {
        return L2DTransforms(curr_set);
    }

    // Sort the vector by the first criteria and then call helper function
    std::sort(curr_set.begin(), curr_set.end(), Compare_Transforms_New(dimension, num_criteria));

    unsigned long num = curr_set.size();
    std::vector<Transform*> r = divide_and_conquer_transforms_helper(curr_set, 0, curr_set.size() - 1, dimension);

    if (this->config.verbose) {
        std::cout << "Transforms deleted: " << (num - r.size()) << std::endl;
    }

    return r;
}


/**
 * Helper function for the divide and conquer algorithm, splits problems into inferior and superior sets and
 * determines dominance
 *
 * @param curr_set
 * @param low
 * @param high
 * @param dimension
 * @return
 */
std::vector<Transform*> DPTransforms::divide_and_conquer_transforms_helper(std::vector<Transform*> &curr_set, unsigned long low,
                                                                 unsigned long high, int dimension) {
    if (high - low == 0) {
        // We want to return the single value
        return std::vector<Transform* > (curr_set.begin() + low, curr_set.begin() + high + 1);
    }

    unsigned long midpoint = (high + low) / 2;

    std::vector<Transform* > superior = divide_and_conquer_transforms_helper(curr_set, low, midpoint, dimension);
    std::vector<Transform* > inferior = divide_and_conquer_transforms_helper(curr_set, midpoint + 1, high, dimension);

    std::vector<Transform* > sup_inf;
    sup_inf.reserve(superior.size() + inferior.size());

    // Copy those from superior
    for (Transform* sol: superior) {
        sol->inferior = false;
        sup_inf.push_back(sol);
    }
    // Copy from inferior
    for (Transform* sol: inferior) {
        sol->inferior = true;
        sup_inf.push_back(sol);
    }
    // Sort by d - 1 before calling marry
    std::sort(sup_inf.begin(), sup_inf.end(), Compare_Transforms_New(dimension - 1, num_criteria));

    std::vector<Transform* > keep_inferior = marry_transform(sup_inf, 0, sup_inf.size() - 1, dimension - 1);
    // Take the union of superior and the inferior points kept
    superior.insert(superior.end(), keep_inferior.begin(), keep_inferior.end());

    return superior;
}


/**
 * Marry sets of sorted Transforms together to identify dominated Transforms
 *
 * @param curr_set
 * @param low
 * @param high
 * @param dimension
 * @return
 */
std::vector<Transform*> DPTransforms::marry_transform(std::vector<Transform*> &curr_set, unsigned long low, unsigned long high, int dimension) {
    if (dimension == 2) { // We can use marry_2D to quickly find the optimal pairs
        return marry_2d_transform(curr_set);
    }
    if (curr_set.empty()){
        return curr_set;
    }
    // Create divide plane and count the inferior and superior in each section
    unsigned long median = (low + high) / 2;
    // We need to keep the superior points from the first half for future use
    std::vector<Transform*> sup_first;
    // Lets count the number in each as a base case!
    unsigned long sup1 = 0, inf1 = 0, sup2 = 0, inf2 = 0;
    // First half
    for (unsigned long i = low; i <= median; i++) {

        if (curr_set[i]->inferior) {
            inf1++;
        } else {
            sup1++;
            sup_first.push_back(curr_set[i]);
        }

    }
    // Second half
    for (unsigned long i = median + 1; i <= high; i++) {
        if (curr_set[i]->inferior) {
            inf2++;
        }
        else {
            sup2++;
        }
    }
    // Check base cases!!
    if (inf1 + inf2 == high - low + 1) { // All inferior so we return all!
        return std::vector<Transform* > (curr_set.begin() + low, curr_set.begin() + high + 1);
    }
    else if (sup1 + sup2 == high - low + 1) { // All superior so return none
        return {};
    }
    else if (inf1 == 0 && sup2 == 0) { // All superior points remain superior and inf remain inferior -- drop dimension
        std::vector<Transform* > immediate_drop(curr_set.begin() + low, curr_set.begin() + high + 1);
        // Sort and drop dimension on curr set
        std::sort(immediate_drop.begin(), immediate_drop.end(), Compare_Transforms_New(dimension - 1, num_criteria));

        return marry_transform(immediate_drop, 0, immediate_drop.size() - 1, dimension - 1);
    }
    else if (sup1 == 0 && inf2 == 0){ // Inferior points are in first half and superior in second
        return std::vector<Transform* > (curr_set.begin() + low, curr_set.begin() + median + 1);
    }

    std::vector<Transform* > inf_first = marry_transform(curr_set, low, median, dimension);
    std::vector<Transform* > inf_second = marry_transform(curr_set, median + 1, high, dimension);

    // Now we must merge the sup_first with inf_second by looking at d - 1
    // Take the union of the two and sort
    sup_first.insert(sup_first.end(), inf_second.begin(), inf_second.end());
    // Sort by d - 1 before calling marry
    std::sort(sup_first.begin(), sup_first.end(), Compare_Transforms_New(dimension - 1, num_criteria));

    std::vector<Transform* > keep_inf_second = marry_transform(sup_first, 0, sup_first.size() - 1, dimension - 1);
    // Take union of inf_first and keep_int_second
    inf_first.insert(inf_first.end(), keep_inf_second.begin(), keep_inf_second.end());
    return inf_first;
}

/**
 * Two dimensional marry function for transforms
 * @param curr_set
 * @return
 */
std::vector<Transform*> DPTransforms::marry_2d_transform(std::vector<Transform*> &curr_set) {
    // Solutions to be saved
    std::vector<Transform* > saved_solutions;
    // Used to track the max first criteria seen so far
    double max_first = std::numeric_limits<double>::lowest();
    for (Transform* sol: curr_set) {
        // compare max_first or the first criteria
        if (sol->base_criteria_minmax[0] > max_first) {
            // Keep solution if is inferior
            if (sol->inferior) {
                saved_solutions.push_back(sol);
            } else {
                max_first = sol->base_criteria_minmax[0];
            }
        } else if (sol->inferior) { // Free the memory pointed to by sol because it is dominated
            delete(sol);
        }
    }

    return saved_solutions;
}