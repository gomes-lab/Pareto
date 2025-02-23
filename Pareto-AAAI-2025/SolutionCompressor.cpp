//
// Created by Marc Grimson on 11/27/23.
//

#include "SolutionCompressor.h"

/**
 * Solution Compressor for selecting represenative solutions to reduce the solution set size
 *
 * @param Network* net - the tree network
 * @param int num_criteria - number of criteria
 * @param float gamma - gamma value for compression strength
 * @param int distance_mode - distance mode (decision space distance or objective space distance)
 * @param int compress_linkage - linkage mode (single, average, or complete)
 */
SolutionCompressor::SolutionCompressor(Network* net, int num_criteria_all, int num_criteria_optim,
                                       float gamma, const int distance_mode,
                                       const int compress_linkage, const int threads, const long chunk_size,
                                       const int compress_objectives)
    : num_criteria_all(num_criteria_all), num_criteria_optim(num_criteria_optim), gamma(gamma), net(net), distance_mode(distance_mode),
        compress_linkage(compress_linkage), threads(threads), chunk_size(chunk_size),
        compress_objectives(compress_objectives) {

}

/**
 * Add a cluster to the clustering algorithm
 *
 * @param node
 * @param cluster
 * @param leaves
 * @param leaves_by_cluster
 */
void SolutionCompressor::add_cluster(std::vector<Solution*> &solutions, const int& cluster, std::vector<int>* leaves,
                                     std::vector<std::vector<int>> &leaves_by_cluster) {
    if (cluster < solutions.size()) {
        leaves->push_back(cluster);
    } else {
        std::vector<int> to_add = leaves_by_cluster[cluster - solutions.size()];
        for (const auto& i : to_add) {
            leaves->push_back(i);
        }
    }
}

/**
 * Find the leaves of the hierarchical cluster
 *
 * @param node
 * @param linkage_matrix
 * @param leaves_by_cluster
 */
void SolutionCompressor::find_leaves(std::vector<Solution*> &solutions, std::vector<std::vector<int>> &linkage_matrix,
                                     std::vector<std::vector<int>> &leaves_by_cluster) {
    for (int i = 0; i < linkage_matrix.size(); i++) {
        std::vector<int> leaves;
        int cluster_1 = linkage_matrix[i][0];
        int cluster_2 = linkage_matrix[i][1];
        add_cluster(solutions, cluster_1, &leaves, leaves_by_cluster);
        add_cluster(solutions, cluster_2, &leaves, leaves_by_cluster);
        leaves_by_cluster.push_back(leaves);
    }
}

/**
 * Run linkage for clustering using the hclust_fast algorithm
 *
 * @param node
 * @param distances
 * @param linkage_matrix
 */
void SolutionCompressor::run_linkage(std::vector<Solution*> &solutions, double* distances, std::vector<std::vector<int>> &linkage_matrix) {
    double* height = new double[solutions.size()-1];
    int* merge = new int[2*(solutions.size()-1)];

    switch(compress_linkage) {
        case COMPRESS_LINKAGE_SINGLE:
            hclust_fast(solutions.size(), distances, HCLUST_METHOD_SINGLE, merge, height);
            break;
        case COMPRESS_LINKAGE_AVERAGE:
            hclust_fast(solutions.size(), distances, HCLUST_METHOD_AVERAGE, merge, height);
            break;
        case COMPRESS_LINKAGE_COMPLETE:
            hclust_fast(solutions.size(), distances, HCLUST_METHOD_COMPLETE, merge, height);
            break;
        default:
            throw "Unknown linkage type";
    }

    for (int i = 0; i < solutions.size()-1; ++i) {
        int j = i + solutions.size()-1;

        int first, second;
        if (*(merge+i) < 0) {
            first = abs(*(merge+i)) - 1;
        } else {
            first = *(merge+i) - 1 + solutions.size();
        }
        if (*(merge+j) < 0) {
            second = abs(*(merge+j)) - 1;
        } else {
            second = *(merge+j) - 1 + solutions.size();
        }

        std::vector<int> linkage = {first, second};

        linkage_matrix.push_back(linkage);
    }

    delete[] distances;
    delete[] height;
    delete[] merge;

    distances = NULL;
    height = NULL;
    merge = NULL;
}

/**
 * Find coverage based on the gamma value to identify representative solutions
 *
 * @param node
 * @param linkage_matrix
 * @param leaves_by_cluster
 * @param covered
 */
void SolutionCompressor::cover(
        std::vector<Solution*> &solutions,
        const std::vector<std::vector<int>>& linkage_matrix,
        const std::vector<std::vector<int>>& leaves_by_cluster,
        int* covered, int start)
{
    const int n = linkage_matrix.size();
    assert(n == leaves_by_cluster.size());
    std::vector<int> leaves = leaves_by_cluster[n-1];
    if (leaves.size() == 1) {
        int leave = leaves[0];
        assert(covered[leave + start] == -1);
        covered[leave + start] = leave;
        return;
    }

    int potential_center = -1;

    potential_center = leaves[rand() % leaves.size()];

    if (potential_center == -1) {
        throw "Center not found";
    }

    Solution* obj_cen = solutions[potential_center];
    bool is_center = true;
    int nc = 0;
    if (compress_objectives == COMPRESS_OBJECTIVES_ALL) {
        nc = num_criteria_all;
    }
    else {
        nc = num_criteria_optim;
    }
    for (const auto& leave : leaves) {
        Solution* obj_leave = solutions[leave];
        for (int i = 0; i < nc; i++) {
            if (net->minmax[i] < 0) {
                if ((compress_objectives == COMPRESS_OBJECTIVES_OPTIMIZED
                        && obj_cen->criteria_2[i].first*(1-gamma) > obj_leave->criteria_2[i].first)
                    || (compress_objectives == COMPRESS_OBJECTIVES_ALL
                        && obj_cen->criteria[i]*(1-gamma) > obj_leave->criteria[i])) {
                    is_center = false;
                    break;
                }
            } else {
                if ((compress_objectives == COMPRESS_OBJECTIVES_OPTIMIZED
                        && obj_cen->criteria_2[i].first < obj_leave->criteria_2[i].first*(1-gamma))
                    || (compress_objectives == COMPRESS_OBJECTIVES_ALL
                        && obj_cen->criteria[i] < obj_leave->criteria[i]*(1-gamma))) {
                    is_center = false;
                    break;
                }
            }
        }
        if (!is_center) {
            break;
        }
    }

    if (!is_center) {
        int left, right;
        left = linkage_matrix[n-1][0];
        right = linkage_matrix[n-1][1];
        if (left < solutions.size()) {
            if (covered[left + start] != -1) {
                std::cerr << "Double center for: " << left << "\t" << covered[left + start] << "\t" << left << endl;
                abort();
            }
            covered[left + start] = left;
        } else {
            left = left - solutions.size();
            std::vector<std::vector<int>> linkage_matrix_left;
            std::vector<std::vector<int>> leaves_by_cluster_left;
            for (size_t i = 0; i <= left; ++i) {
                linkage_matrix_left.push_back(linkage_matrix[i]);
                leaves_by_cluster_left.push_back(leaves_by_cluster[i]);
            }
            cover(solutions, linkage_matrix_left, leaves_by_cluster_left, covered, start);
        }
        if (right < solutions.size()) {
            if (covered[right + start] != -1) {
                std::cerr << "Double center for: " << right << "\t" << covered[right + start] << "\t" << potential_center << endl;
                abort();
            }
            covered[right + start] = right;
        } else {
            right = right - solutions.size();
            std::vector<std::vector<int>> linkage_matrix_right;
            std::vector<std::vector<int>> leaves_by_cluster_right;
            for (size_t i = 0; i <= right; ++i) {
                linkage_matrix_right.push_back(linkage_matrix[i]);
                leaves_by_cluster_right.push_back(leaves_by_cluster[i]);
            }
            cover(solutions, linkage_matrix_right, leaves_by_cluster_right, covered, start);
        }
    } else {
        for (const auto& leave : leaves) {
            if (covered[leave + start] != -1) {
                std::cerr << "Double center for: " << leave << "\t" << covered[leave + start] << "\t" << potential_center << endl;
                abort();
            }
            covered[leave + start] = potential_center;
        }
    }
}

/**
 * Compress the solutions using the hierarchical clustering algorithm
 *
 * @param node
 */
void SolutionCompressor::compress_solutions(TreeNode *node) {
    std::cout << "Normalizing solutions" << std::endl;
    // First we normalize the solution values
    normalize_solutions(node);

    long frontier_size = node->vec_frontier.size();

    int chunks = ceil((double) frontier_size / (double) chunk_size);
    int local_threads = threads;

    if (chunks < local_threads) {
        local_threads = chunks;
    }

    std::vector<std::vector<std::vector<int>>> linkages(chunks);
    std::vector<std::vector<std::vector<int>>> linkage_clusters(chunks);
    // Parallel process
    #pragma omp parallel for num_threads(local_threads)
    for (int i = chunks - 1; i >= 0; i--) {
        int tid = omp_get_thread_num();

        int start = i * chunk_size;
        int end = (i + 1) * chunk_size;
        if (end > frontier_size) {
            end = frontier_size;
        }

        std::cout << "From " << start << " to " << end << std::endl;

        std::vector<Solution *> solns_subset(node->vec_frontier.begin() + start, node->vec_frontier.begin() + end);

        std::cout << tid << ": Calculating distances" << std::endl;
        long long size = (solns_subset.size() - 1) * solns_subset.size() / 2;
        double *distances = new double[size];
        calculate_distances(solns_subset, distances);

        std::cout << tid << ": Running linkage" << std::endl;
        std::vector<std::vector<int>> linkage_matrix;

        run_linkage(solns_subset, distances, linkage_matrix);
        linkages[i] = linkage_matrix;

        std::cout << tid << ": Finding leaves" << std::endl;
        std::vector<std::vector<int>> leaves_by_cluster;
        find_leaves(solns_subset, linkage_matrix, leaves_by_cluster);

        linkage_clusters[i] = leaves_by_cluster;
    }

    int* covered = new int[frontier_size];
    for (int i = 0; i < frontier_size; i++) {
        covered[i] = -1;
    }

    std::vector<std::set<int>> centers(chunks);

    #pragma omp parallel for num_threads(local_threads)
    for (int i = chunks - 1; i >= 0; i--) {
        int tid = omp_get_thread_num();

        int start = i * chunk_size;
        int end = (i + 1) * chunk_size;
        if (end > frontier_size) {
            end = frontier_size;
        }

        std::vector<Solution *> solns_subset(node->vec_frontier.begin() + start, node->vec_frontier.begin() + end);

        std::cout << "Running coverage" << std::endl;
        cover(solns_subset, linkages[i], linkage_clusters[i], covered, start);
        for (int j = 0; j < solns_subset.size(); j++) {
            int cover = covered[j + start] + start;

            centers[i].insert(cover);
        }
    }


    // Construct a set to allow for deleting of unused solutions more easily
    std::unordered_set<Solution*, SolutionHash, EqualSolutions> kept_solutions;
    for (auto thread_centers : centers) {
        for (auto center : thread_centers) {
            kept_solutions.insert(node->vec_frontier[center]);
        }
    }

    // Delete all solutions that aren't in the new set
    for (auto *soln : node->vec_frontier) {
        // Solution not in the kept solutions set
        if (kept_solutions.find(soln) == kept_solutions.end()) {
            delete(soln);
        }
    }

    // Copy to a vector
    std::vector<Solution* > new_solutions(kept_solutions.begin(), kept_solutions.end());

    // Now rerun a single-threaded clustering
    std::cout << "Calculating distances" << std::endl;
    long long size = (new_solutions.size() - 1) * new_solutions.size() / 2;
    double *distances = new double[size];
    calculate_distances(new_solutions, distances);

    std::cout << "Running linkage" << std::endl;
    std::vector<std::vector<int>> linkage_matrix;

    run_linkage(new_solutions, distances, linkage_matrix);

    std::cout << "Finding leaves" << std::endl;
    std::vector<std::vector<int>> leaves_by_cluster;
    find_leaves(new_solutions, linkage_matrix, leaves_by_cluster);

    int* final_covered = new int[new_solutions.size()];
    for (int i = 0; i < new_solutions.size(); i++) {
        final_covered[i] = -1;
    }

    std::set<int> final_centers;

    std::cout << "Running coverage" << std::endl;
    cover(new_solutions, linkage_matrix, leaves_by_cluster, final_covered, 0);
    for (int j = 0; j < new_solutions.size(); j++) {
        int cover = final_covered[j];

        final_centers.insert(cover);
    }

    kept_solutions.clear();

    for (auto center : final_centers) {
        kept_solutions.insert(new_solutions[center]);
    }

    // Delete all solutions that aren't in the new set
    for (auto *soln : new_solutions) {
        // Solution not in the kept solutions set
        if (kept_solutions.find(soln) == kept_solutions.end()) {
            delete(soln);
        }
    }

    // Copy to a vector
    std::vector<Solution* > final(kept_solutions.begin(), kept_solutions.end());

    node->vec_frontier = final;

    // Clean up
    for (auto linkage_matrix : linkages) {
        linkage_matrix.clear();
    }
    for (auto leaves_by_cluster : linkage_clusters) {
        leaves_by_cluster.clear();
    }
    delete[] (covered);
    delete[] (final_covered);
}

/**
 * Calculate euclidean distances between two solutions based on the objectives
 *
 * @param soln1
 * @param soln2
 * @return
 */
float SolutionCompressor::calculate_euclidean_distance(Solution *soln1, Solution *soln2) {
    float sum_of_difference = 0;
    int nc = 0;
    if (compress_objectives == COMPRESS_OBJECTIVES_ALL) {
        nc = num_criteria_all;
    }
    else {
        nc = num_criteria_optim;
    }
    for (int i = 0; i < nc; i++) {
        sum_of_difference += pow((soln1->criteria_normalized[i]-soln2->criteria_normalized[i]), 2);
    }
    return sqrt(sum_of_difference);
}

/**
 * Calculate hamming distance between two solutions based on the decisions
 *
 * @param soln1
 * @param soln2
 * @return
 */
float SolutionCompressor::calculate_hamming_distance(Solution *soln1, Solution *soln2) {
    float distance = 0.0;

    boost::dynamic_bitset<> temp;
    temp.operator=(soln1->representation);
    temp.operator^=(soln2->representation);

    distance = temp.count();

    return distance;
}

/**
 * Calculate distance matrix between all solutions
 *
 * @param node
 * @param distances
 */
void SolutionCompressor::calculate_distances(std::vector<Solution*> &solutions, double* distances) {
    int k = 0;
    for (int i = 0; i < solutions.size(); i++) {
        for (int j = i + 1; j < solutions.size(); j++) {
            auto* soln1 = solutions[i];
            auto* soln2 = solutions[j];

            switch(distance_mode) {
                case COMPRESS_DISTANCE_OBJECTIVE_SPACE:
                    distances[k] = static_cast<double>(calculate_euclidean_distance(soln1, soln2));
                    break;
                case COMPRESS_DISTANCE_DECISION_SPACE:
                    distances[k] = static_cast<double>(calculate_hamming_distance(soln1, soln2));
                    break;
            }

            k++;
        }
    }
}

/**
 * Normalize solution values to be between 0 and 1
 *
 * @param node
 */
void SolutionCompressor::normalize_solutions(TreeNode *node) {
    if (compress_objectives == COMPRESS_OBJECTIVES_ALL) {
        std::vector<float> mins(num_criteria_all, std::numeric_limits<float>::max());
        std::vector<float> maxs(num_criteria_all, std::numeric_limits<float>::lowest());

        for (auto *soln : node->vec_frontier) {
            for (int i = 0; i < num_criteria_all; i++) {
                if (mins[i] > soln->criteria[i]) {
                    mins[i] = soln->criteria[i];
                }
                if (maxs[i] < soln->criteria[i]) {
                    maxs[i] = soln->criteria[i];
                }
            }
        }

        for (auto *soln : node->vec_frontier) {
            for (int i = 0; i < num_criteria_all; i++) {
                if (maxs[i] - mins[i] == 0) {
                    soln->criteria_normalized[i] = 0.0;
                } else {
                    soln->criteria_normalized[i] = (soln->criteria[i] - mins[i]) / (maxs[i] - mins[i]);
                }
            }
        }
    }
    else {
        std::vector<float> mins(num_criteria_optim, std::numeric_limits<float>::max());
        std::vector<float> maxs(num_criteria_optim, std::numeric_limits<float>::lowest());

        for (auto *soln : node->vec_frontier) {
            for (int i = 0; i < num_criteria_optim; i++) {
                if (mins[i] > soln->criteria_2[i].first) {
                    mins[i] = soln->criteria_2[i].first;
                }
                if (maxs[i] < soln->criteria_2[i].first) {
                    maxs[i] = soln->criteria_2[i].first;
                }
            }
        }

        for (auto *soln : node->vec_frontier) {
            for (int i = 0; i < num_criteria_optim; i++) {
                if (maxs[i] - mins[i] == 0) {
                    soln->criteria_normalized[i] = 0.0;
                } else {
                    soln->criteria_normalized[i] = (soln->criteria_2[i].first - mins[i]) / (maxs[i] - mins[i]);
                }
            }
        }
    }
}