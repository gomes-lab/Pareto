//
// Created by Jonathan Gomes Selman on 7/6/17.
//

#include "ParetoDP.h"

#include <utility>

/**
 * Constructor for the ParetoDP class. This class is used to execute the full Dynamic Programming algorithm
 * to calculate the Pareto frontier for the problem. After constructing, the user should execute run_expertiment()
 * to calculate the frontier.
 *
 * @param Network net - Network Object
 * @param int root - Root node nodeID
 * @param optimized_criteria - vector of criteria being used in determining the Pareto frontier
 * @param epsilon - epsilon value for rounding solutions
 * @param seed - seed to use for the random number generator - only consistent if running with one thread
 * @param batch_size - batch size to reach before pruning
 * @param num_threads - number of threads to use
 * @param io - IO object for reading and writing data
 * @param use_transforms - use the Transforms objects when merging child frontiers
 * @param merge_on_transform - run the merge step on transformations, instead of computing an entire divide-and-conquer step
 */
ParetoDP::ParetoDP(Network *net, const IOHandler& io, const Config& config, bool brute_force)
        : num_criteria(config.optimized_criteria.size()), net(net), io(io), config(config) {

    this->true_num_criteria = num_criteria;

    if (config.use_linear_preferences) {
        this->num_criteria = (int) config.w.size() + 1;
    }

    // Initialize counters
    this->total_pruned_policies = 0;
    this->total_policies_considered = 0;
    this->max_policies_considered = 0;
    this->total_policies_compared = 0;

    // Initialize timers
    this->time_copying = 0;
    this->time_generating = 0;
    this->time_sorting = 0;
    this->time_nlogn = 0;
    this->time_generating_transforms = 0;
    this->time_transforming = 0;
    this->time_sorting_transforms = 0;
    this->time_merging = 0;

    this->clocks_copying = 0;
    this->clocks_generating = 0;
    this->clocks_sorting = 0;
    this->clocks_nlogn = 0;
    this->clocks_generating_transforms = 0;
    this->clocks_transforming = 0;
    this->clocks_sorting_transforms = 0;
    this->clocks_merging = 0;

    this->node_processed_counter = 0;

    if (config.compress_mode != COMPRESS_MODE_NONE) {
        compressor = new SolutionCompressor(net, (config.dcip_index >= 0 ? net->num_criteria + 1 : net->num_criteria), num_criteria,
                                            config.compress_gamma, config.compress_distance_mode, config.compress_linkage,
                                            config.num_threads, config.compression_chunk_size, config.compress_objectives);
    }

    if (brute_force) {
        this->dp = new DPBruteForce(this->net, this->num_criteria, config);
    }
    else if (config.use_transforms) {
        this->dp = new DPTransforms(this->net, this->num_criteria, config);
    }
    else {
        this->dp = new DPOriginal(this->net, this->num_criteria, config);
    }

    // Calculate rounding constants
    calc_theoretical_ks_tree(this->net->root_node);
}

ParetoDP::~ParetoDP() {
    delete dp;
}

/**
 * Calculates the theoretical k values for the tree recursively
 *
 * @param node
 */
void ParetoDP::calc_theoretical_ks_tree(TreeNode* node) {
    if (node != nullptr) {
        for (int i = 0; i < num_criteria; ++i) {
            node->node_data.k_vals[i] = node->node_data.r_vals[i].first * config.epsilon;
        }

        // Recurse
        for (auto &child : node->children) {
            calc_theoretical_ks_tree(child);
        }
    }
}

/**
 * Read the policies for a given node - used on a restart to re-read the data.
 *
 * @param node
 * @return
 */
bool ParetoDP::readPolicies(TreeNode* node) {
    int node_id = node->node_id;
    std::stringstream reader;
    std::string line;
    if (this->io.check_file(META_FILE_NAME, node_id)) {
        int flag;
        int max_solution_size;
        std::fstream meta_file = this->io.get_file(false, META_FILE_NAME, node_id);
        this->io.get_line(meta_file, reader, line);
        reader >> flag;
        this->io.get_line(meta_file, reader, line);
        reader >> max_solution_size;
        // This node has been completed, get the policies
        if (flag) {
            node->max_frontier_size = max_solution_size;
            std::fstream solution_file = this->io.get_file(false, SOLUTION_FILE_NAME, node_id);

            int vec_size;

            this->io.get_line(solution_file, reader, line);
            reader >> vec_size;

            for (int i = 0; i < vec_size; i++) {
                int solution_id, solution_node_id, solution_tree_node_id, crit_size, dam_dec_size, pareto_dec_size;
                bool inferior;
                this->io.get_line(solution_file, reader, line);
                reader >> solution_node_id >> solution_tree_node_id >> solution_id >> inferior >> crit_size;

                double solution_crit_all[MAX_CRITERIA];
                std::pair<double, double> solution_crit[crit_size];

                for (int j = 0; j < crit_size; j++) {
                    std::string crit_vals_tmp;
                    reader >> crit_vals_tmp;

                    size_t pos;
                    pos = crit_vals_tmp.find('|');

                    solution_crit[j] = {std::stof(crit_vals_tmp.substr(0, pos)), std::stof(crit_vals_tmp.substr(pos+1))};
                }

                std::vector<double> dam_decisions;

                reader >> dam_dec_size;
                for (int j = 0; j < dam_dec_size; j++) {
                    double dec;
                    reader >> dec;
                    dam_decisions.emplace_back(dec);
                }

                int num_children;
                reader >> num_children;
                std::vector<Solution*> children;
                for (int j = 0; j < num_children; j++) {
                    std::string child_tmp;
                    reader >> child_tmp;

                    size_t pos = child_tmp.find('|');

                    int child_id, child_solution_id;
                    child_id = std::stoi(child_tmp.substr(0, pos));
                    child_solution_id = std::stoi(child_tmp.substr(pos+1));
                    Solution* child = node->children[j]->vec_frontier[child_solution_id];
                    children.emplace_back(child);
                }

                auto p = new Solution(solution_node_id, solution_tree_node_id, solution_crit_all, solution_crit,
                                      dam_decisions, children, net->num_criteria, num_criteria,
                                      config.use_linear_preferences, config.w, net->minmax, net->normalizing_factor);
                p->solution_id = solution_id;

                node->vec_frontier.emplace_back(p);
            }
        }

        return true;
    }
    else {
        return false;
    }
}


/**
 * Main DP step that determines the node type (leaf, one child, two children, three+ children, etc.) and executes the
 * appropriate sub-function for that node type.
 *
 * @param node
 */
void ParetoDP::computeDP(TreeNode *node, int dynamic_depth, int static_depth) {
    std::cout << "Working on node " << node->node_id << " at depth (dynamic) " << dynamic_depth;
    std::cout << " (static) " << static_depth << " - ";
    if (node->intermediate) {
        std::cout << "intermediate of " << node->node_data.id << " - ";
    }
    std::cout << node->children.size() << " children";
    if (node->children.size() == 1) {
        std::cout << " - child " << node->children[0]->node_id << ": "
                  << node->children[0]->vec_frontier.size() << " (" << node->children[0]->parentDam.decision.size() << ")";
    }
    else if (node->children.size() == 2) {
        std::cout << " - left " << node->children[0]->node_id << ": "
                  << node->children[0]->vec_frontier.size() << " (" << node->children[0]->parentDam.decision.size() << ")";

        std::cout << " - right " << node->children[1]->node_id << ": "
                  << node->children[1]->vec_frontier.size() << " (" << node->children[1]->parentDam.decision.size() << ")";
    }
    std::cout << std::endl;

    std::chrono::high_resolution_clock::time_point start, stop;
    std::chrono::duration<float> duration;

    start = std::chrono::high_resolution_clock::now();
    this->dp->computeDP(node);
    stop = std::chrono::high_resolution_clock::now();
    duration = stop - start;

    // If static depth less than depth compression, compress - need at least 2 solutions
    if (this->config.compress_mode == COMPRESS_MODE_DEPTH_STATIC && static_depth <= this->config.compression_depth
        && node->vec_frontier.size() >= 2) {
        std::cout << "Compressing solutions due to static depth" << std::endl;
        long orig_size = node->vec_frontier.size();
        compressor->compress_solutions(node);
        std::cout << "Compressed " << orig_size << " solutions down to " << node->vec_frontier.size();
        std::cout << std::endl << std::endl;
    }

    // If dynamic depth less than depth compression, compress - need at least 2 solutions
    if (this->config.compress_mode == COMPRESS_MODE_DEPTH_DYNAMIC && dynamic_depth <= this->config.compression_depth
        && node->vec_frontier.size() >= 2) {
        std::cout << "Compressing solutions due to dynamic depth" << std::endl;
        long orig_size = node->vec_frontier.size();
        compressor->compress_solutions(node);
        std::cout << "Compressed " << orig_size << " solutions down to " << node->vec_frontier.size();
        std::cout << std::endl << std::endl;
    }

    // If frontier size is greater than threshold, compress
    if (this->config.compress_mode == COMPRESS_MODE_THRESHOLD
        && node->vec_frontier.size() >= this->config.compression_threshold) {
        std::cout << "Compressing solutions due to frontier size" << std::endl;
        long orig_size = node->vec_frontier.size();
        compressor->compress_solutions(node);
        std::cout << "Compressed " << orig_size << " solutions down to " << node->vec_frontier.size();
        std::cout << std::endl << std::endl;
    }

    // Re-write the node data and tree meta in case it has changed or this is an intermediate node
    net->writeMeta();
    net->writeTreeNode(node, false);

    if (this->io.config.save_progress) {
        writeFrontier(node);
    }

    writeMeta(node, duration.count());

    std::cout << "Completed node " << node->node_id << " with " << node->vec_frontier.size() << " solutions" << std::endl;
}

/**
 * Write out meta data for a node
 *
 * @param node
 */
void ParetoDP::writeMeta(TreeNode *node, float duration) {
    std::fstream meta_file = io.get_file(true, META_FILE_NAME, node->node_id);

    meta_file << COMPLETION_FLAG << std::endl;
    meta_file << node->node_id << " " << (node->intermediate ? 1 : 0);

    meta_file << " " << node->children.size();
    for (auto &child : node->children) {
        meta_file << " " << child->node_id;
    }

    meta_file << std::endl;
    meta_file << node->max_frontier_size << " " << node->vec_frontier.size();
    meta_file << " " << duration << std::endl;

    if (!dp->ks.empty()) {
        meta_file << dp->ks.size() << " " << num_criteria;
        for (auto &k : dp->ks) {
            for (float v : k) {
                meta_file << " " << v;
            }
        }

        meta_file << std::endl;
    }
}

/**
 * Save the frontier for a given node
 *
 * @param node
 */
void ParetoDP::writeFrontier(TreeNode *node) {
    std::fstream solution_file = io.get_file(true, SOLUTION_FILE_NAME, node->node_id);

    int i = 0;
    solution_file << node->vec_frontier.size() << std::endl;
    for (Solution *p : node->vec_frontier) {
        p->solution_id = i; // Set the solution ID here so we can retrieve the info on a restart
        solution_file << p->node_id << " " << p->tree_node_id << " " << p->solution_id << " " << (p->inferior ? "1" : "0") << " ";

        solution_file << MAX_CRITERIA << " ";
        for (auto const &crit : p->criteria_2) {
            solution_file << crit.first << "|" << crit.second << " ";
        }

        solution_file << p->dam_decisions.size() << " ";
        for (auto const &dec : p->dam_decisions) {
            solution_file << dec << " ";
        }

        solution_file << node->children.size() << " ";
        int j = 0;
        for (auto &child : node->children) {
            Solution* child_solution = p->pareto_decisions[j];
            solution_file << child_solution->node_id << "|" << child_solution->solution_id << " ";
            j++;
        }

        solution_file << std::endl;

        i++;
    }

    solution_file.close();
}

/**
 * Build the DP table
 */
void ParetoDP::build_DP() {
    node_processed_counter = 0;

    build_DP_table_recursive_tree(net->root_node, 0, 0);

    // Compress if we are either doing posthoc compression or the compression depth is set to 0 (equivalent of
    // final node)
    if (this->config.compress_mode == COMPRESS_MODE_FINAL ||
            ((this->config.compress_mode == COMPRESS_MODE_DEPTH_DYNAMIC
            || this->config.compress_mode == COMPRESS_MODE_DEPTH_STATIC)
            && this->config.compression_depth == 0)) {
        std::cout << "Compressing solutions at final node" << std::endl;
        long orig_size = net->root_node->vec_frontier.size();
        compressor->compress_solutions(net->root_node);
        std::cout << "Compressed " << orig_size << " solutions down to " << net->root_node->vec_frontier.size();
        std::cout << std::endl << std::endl;
    }
}

void ParetoDP::brute_force() {
    node_processed_counter = 0;

    build_DP_table_recursive_tree(net->root_node, 0, 0);

    std::vector<Solution *> solution_set(net->root_node->vec_frontier.begin(), net->root_node->vec_frontier.end());
    net->root_node->vec_frontier = this->dp->dnc->divide_and_conquer(solution_set, num_criteria, false,
                                                                     this->dp->num_comparisons);
}


/**
 * Recursively build the DP table - also dynamically generates intermediate nodes for  nodes with more than two children
 *
 * @param node
 */
void ParetoDP::build_DP_table_recursive_tree(TreeNode *node, int dynamic_depth, int static_depth) {
    if (node != nullptr) {
        for (auto &child : node->children) {
            build_DP_table_recursive_tree(child, dynamic_depth + node->children.size(), static_depth + 1);
        }

        // Helpful tracker
        node_processed_counter++;

        // Subtract 2 - one for extra count at size, and an extra because the final two children
        // is a single join
        int extra_depth = node->children.size() - 2;
        while (node->has_next()) {
            TreeNode* next_node = node->get_next(net->num_nodes_total, net->compare_operator);
            net->hyper_nodes[next_node->node_id] = next_node;
            // Children are at this node's depth + 1
            computeDP(next_node, dynamic_depth + extra_depth + 1, static_depth + 1);
            extra_depth--;
        }
    }
}


/**
 * Print the DP output to a file
 *
 * @param stream
 */
void ParetoDP::print_DP_Output(std::ostream &stream) {
    // Output num_solutions
    stream << "num_solutions: ";
    stream << net->root_node->vec_frontier.size() << std::endl;
    // Give data about pruning
    stream << "# pruning steps (# nodes): " << net->num_nodes_total << std::endl;
    stream << "Max policies considered: " << max_policies_considered << std::endl;
    stream << "Policies considered: " << total_policies_considered << std::endl;
    stream << "Pruned policies: " << total_pruned_policies << std::endl;

    // Output the epsilon
    stream << "epsilon: " << this->config.epsilon << std::endl;
    // ThreadData
    stream << "batch size: " << this->config.batch_size << std::endl;
    // Output the relevant criteria
    stream << "criteria considered: ";
    for (int i = 0; i < true_num_criteria; i++) {
        if (i != 0) {
            stream << ", ";
        }
        stream << config.optimized_criteria[i];
    }
    stream << std::endl;

    stream << "criteria all: ";
    for (int i = 0; i < net->num_criteria; i++) {
        if (i != 0) {
            stream << ", ";
        }
        stream << net->criteria[i];
    }
    if (config.dcip_index >= 0) {
        stream << ", dcip";
    }
    stream << std::endl;

    // Output the info for each solution
    print_dams_generic_tree_vec(net->root_node, stream);
}


/**
 * Print the dams to an output buffer
 *
 * @param node
 * @param out
 */
void ParetoDP::print_dams_generic_tree_vec(TreeNode *node, std::ostream &out) {
    int i = 0;
    for (Solution *result: node->vec_frontier) {
        print_dams_helper(result, out, node);
    }
}


/**
 * Helper function for printing dams to an output buffer
 *
 * @param result
 * @param out
 * @param node
 */
void ParetoDP::print_dams_helper(Solution *result, std::ostream &out, TreeNode *node) {
    std::stack<std::pair<TreeNode* , Solution* > > paths;
    paths.emplace(node, result);

    // Explore the "solution" to see what dams are built for the given solution
    std::vector<int> dams_built;

    while (!paths.empty()) {
        TreeNode* decision_node = paths.top().first;
        Solution* decision = paths.top().second;
        paths.pop();

        int i = 0;
        for (auto &child : decision_node->children) {
            if (decision->dam_decisions[i] == 1)
            {
                dams_built.push_back(child->parentDam.id);
            }
            paths.emplace(child, decision->pareto_decisions[i]);

            i++;
        }
    }

    // Print the solutions as described by the output file format -- ***.word

    // Print the criteria values
    for (int i = 0; i < net->num_criteria; i++) {
        out << result->criteria[i] << ", ";
    }
    if (config.dcip_index >= 0) {
        out << pow(result->criteria[net->num_criteria], 0.5) << ", ";
    }
    // Print number of dam
    out << dams_built.size() << ",";
    // Sort dam ids
    std::sort(dams_built.begin(), dams_built.end());
    // Print the dam ids
    for (int dam_id : dams_built) {
        out << " " << dam_id;
    }

    out << std::endl;
}


void ParetoDP::run_brute_force() {
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    clock_t t = clock();

    brute_force();

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    t = clock() - t;

    this->time_generating = (double)this->dp->clocks_generating/CLOCKS_PER_SEC;
    this->time_copying = (double)this->dp->clocks_copying/CLOCKS_PER_SEC;
    this->time_nlogn = (double)this->dp->clocks_nlogn/CLOCKS_PER_SEC;
    this->time_sorting = (double)this->dp->clocks_sorting/CLOCKS_PER_SEC;
    this->time_generating_transforms = (double)this->dp->clocks_generating_transforms/CLOCKS_PER_SEC;
    this->time_sorting_transforms = (double)this->dp->clocks_sorting_transforms/CLOCKS_PER_SEC;
    this->time_transforming = (double)this->dp->clocks_transforming/CLOCKS_PER_SEC;
    this->time_merging = (double)this->dp->clocks_merging/CLOCKS_PER_SEC;

    std::chrono::duration<float> duration = t2 - t1;

    run_info = RunInfo(this->dp->total_policies_considered, this->dp->total_pruned_policies, net->root_node->vec_frontier.size(),
                       this->dp->max_policies_considered, this->dp->total_policies_compared, duration.count(), ((float)t)/CLOCKS_PER_SEC,
                       this->dp->transforms_considered, this->dp->transforms_pruned, this->dp->policies_avoided);

    // Output times!
    if (!config.use_transforms) {
        std::cout << "Time generating: " << time_generating << std::endl;
        std::cout << "Time copying: " << time_copying << std::endl;
        std::cout << "Time sorting: " << time_sorting << std::endl;
        std::cout << "Time nlogn: " << time_nlogn << std::endl << std::endl;

    }
    else {
        std::cout << "Time generating transforms: " << time_generating_transforms << std::endl;
        std::cout << "Time sorting transforms: " << time_sorting_transforms << std::endl;
        std::cout << "Time transforming: " << time_transforming << std::endl;
        std::cout << "Time copying: " << time_copying << std::endl;
        std::cout << "Time sorting: " << time_sorting << std::endl;
    }

    std::cout << "Total policies considered: " << this->dp->total_policies_considered << std::endl;
    std::cout << "Total policies pruned: " << this->dp->total_pruned_policies << std::endl;
    std::cout << "Total comparisons made (approximately): " << this->dp->total_policies_compared << std::endl;
    std::cout << "Max policies considered: " << this->dp->max_policies_considered << std::endl << std::endl;

    // Output to exp file time
    std::cout << "Wall time experiment: " << duration.count() << std::endl;
    std::cout << "CPU time experiment: " << ((float) t) / CLOCKS_PER_SEC << std::endl;

    // Get the time
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [80];

    time (&rawtime);
    timeinfo = localtime (&rawtime);

    strftime (buffer,80,"%c",timeinfo);
    std::string time_stamp = buffer;
    // Get rid of the spaces
    std::replace(time_stamp.begin(), time_stamp.end(), ' ', '_');
    std::replace(time_stamp.begin(), time_stamp.end(), ':', '_');

    if (io.config.save_results) {
        // Print the result of the trail
        std::fstream output = io.get_file(true, SOLUTION_FILE_NAME);
        if (output.is_open()) {
            // Print date and time for trial
            output << "Date/time: " << time_stamp << "\n";

            // File used
            output << "Data file: " << io.config.experiment_name << "\n";

            // Set the precision
            std::ios::fmtflags old_settings = output.flags(); //save previous format flags
            long old_precision = output.precision();

            // NEW
            output << std::fixed << std::setprecision(4) << "Wall time: " << duration.count() << " seconds.\n";
            output << std::fixed << std::setprecision(4) << "CPU time: " << ((float) t) / CLOCKS_PER_SEC
                   << " seconds.\n";


            //Restore precision
            output.flags(old_settings);
            output.precision(old_precision);
            // Print the seed
            output << "seed: " << this->config.seed << "\n";

            print_DP_Output(output);
            output.close();
        } else {
            std::cerr << "COULD NOT OPEN SOLUTION FILE" << std::endl;
            exit(1);
        }
    }
}


/**
 * Runs the experiment to determine the Pareto Frontier for a given problem
 */
void ParetoDP::run_experiment() {
    // Run the experiment
    // Monitor wall and CPU time
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    clock_t t = clock();

    build_DP();

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    t = clock() - t;

    this->time_generating = (double)this->dp->clocks_generating/CLOCKS_PER_SEC;
    this->time_copying = (double)this->dp->clocks_copying/CLOCKS_PER_SEC;
    this->time_nlogn = (double)this->dp->clocks_nlogn/CLOCKS_PER_SEC;
    this->time_sorting = (double)this->dp->clocks_sorting/CLOCKS_PER_SEC;
    this->time_generating_transforms = (double)this->dp->clocks_generating_transforms/CLOCKS_PER_SEC;
    this->time_sorting_transforms = (double)this->dp->clocks_sorting_transforms/CLOCKS_PER_SEC;
    this->time_transforming = (double)this->dp->clocks_transforming/CLOCKS_PER_SEC;
    this->time_merging = (double)this->dp->clocks_merging/CLOCKS_PER_SEC;

    std::chrono::duration<float> duration = t2 - t1;

    run_info = RunInfo(this->dp->total_policies_considered, this->dp->total_pruned_policies, net->root_node->vec_frontier.size(),
                       this->dp->max_policies_considered, this->dp->total_policies_compared, duration.count(), ((float)t)/CLOCKS_PER_SEC,
                       this->dp->transforms_considered, this->dp->transforms_pruned, this->dp->policies_avoided);

    // Output times!
    if (!config.use_transforms) {
        std::cout << "Time generating: " << time_generating << std::endl;
        std::cout << "Time copying: " << time_copying << std::endl;
        std::cout << "Time sorting: " << time_sorting << std::endl;
        std::cout << "Time nlogn: " << time_nlogn << std::endl << std::endl;

    }
    else {
        std::cout << "Time generating transforms: " << time_generating_transforms << std::endl;
        std::cout << "Time sorting transforms: " << time_sorting_transforms << std::endl;
        std::cout << "Time transforming: " << time_transforming << std::endl;
        std::cout << "Time copying: " << time_copying << std::endl;
        std::cout << "Time sorting: " << time_sorting << std::endl;
    }

    std::cout << "Total policies considered: " << this->dp->total_policies_considered << std::endl;
    std::cout << "Total policies pruned: " << this->dp->total_pruned_policies << std::endl;
    std::cout << "Total comparisons made (approximately): " << this->dp->total_policies_compared << std::endl;
    std::cout << "Max policies considered: " << this->dp->max_policies_considered << std::endl << std::endl;

    // Output to exp file time
    std::cout << "Wall time experiment: " << duration.count() << std::endl;
    std::cout << "CPU time experiment: " << ((float) t) / CLOCKS_PER_SEC << std::endl;

    // Get the time
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [80];

    time (&rawtime);
    timeinfo = localtime (&rawtime);

    strftime (buffer,80,"%c",timeinfo);
    std::string time_stamp = buffer;
    // Get rid of the spaces
    std::replace(time_stamp.begin(), time_stamp.end(), ' ', '_');
    std::replace(time_stamp.begin(), time_stamp.end(), ':', '_');

    if (io.config.save_results) {
        // Print the result of the trail
        std::fstream output = io.get_file(true, SOLUTION_FILE_NAME);
        if (output.is_open()) {
            // Print date and time for trial
            output << "Date/time: " << time_stamp << "\n";

            // File used
            output << "Data file: " << io.config.experiment_name << "\n";

            // Set the precision
            std::ios::fmtflags old_settings = output.flags(); //save previous format flags
            long old_precision = output.precision();

            // NEW
            output << std::fixed << std::setprecision(4) << "Wall time: " << duration.count() << " seconds.\n";
            output << std::fixed << std::setprecision(4) << "CPU time: " << ((float) t) / CLOCKS_PER_SEC
                   << " seconds.\n";


            //Restore precision
            output.flags(old_settings);
            output.precision(old_precision);
            // Print the seed
            output << "seed: " << this->config.seed << "\n";

            print_DP_Output(output);
            output.close();
        } else {
            std::cerr << "COULD NOT OPEN SOLUTION FILE" << std::endl;
            exit(1);
        }
    }
}
