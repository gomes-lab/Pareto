//
// Created by Marc Grimson on 11/27/23.
//

#ifndef HYDRODAMDP_SOLUTIONCOMPRESSOR_H
#define HYDRODAMDP_SOLUTIONCOMPRESSOR_H

#include <fenv.h>

#include "Solution.h"
#include "Network.h"
#include "Consts.h"
#include "Config.h"
#include "fastcluster.h"

#include <set>
#include <unordered_set>
#include <vector>
#include <omp.h>

class SolutionCompressor {
private:
    int num_criteria_all;
    int num_criteria_optim;
    float gamma;
    int distance_mode;
    int compress_linkage;
    int compress_objectives;
    int threads;
    long chunk_size;
    Network* net;

    float calculate_euclidean_distance(Solution* soln1, Solution* soln2);
    float calculate_hamming_distance(Solution* soln1, Solution* soln2);
    void run_linkage(std::vector<Solution*> &solutions, double* distances, std::vector<std::vector<int>> &linkage_matrix);
    void find_leaves(std::vector<Solution*> &solutions, std::vector<std::vector<int>> &linkage_matrix,
                     std::vector<std::vector<int>> &leaves_by_cluster);
    void add_cluster(std::vector<Solution*> &solutions, const int& cluster, std::vector<int>* leaves,
                     std::vector<std::vector<int>> &leaves_by_cluster);
    void cover(std::vector<Solution*> &solutions, const std::vector<std::vector<int>> &linkage_matrix,
               const std::vector<std::vector<int>> &leaves_by_cluster,
               int* covered, int start);

    void calculate_distances(std::vector<Solution*> &solutions, double* distances);
    void normalize_solutions(TreeNode* node);

public:
    SolutionCompressor(Network* net, int num_criteria_all, int num_criteria_optim, float gamma, const int distance_mode, const int compress_linkage,
                       const int threads, const long chunk_size, const int compress_objectives);
    void compress_solutions(TreeNode* node);
};


#endif //HYDRODAMDP_SOLUTIONCOMPRESSOR_H
