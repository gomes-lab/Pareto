//
// Created by marc on 7/22/21.
//

#ifndef DP_PARALLEL_GENERATION_DIVIDEANDCONQUER_H
#define DP_PARALLEL_GENERATION_DIVIDEANDCONQUER_H

#include "Solution.h"
#include <iostream>

class DivideAndConquer {
public:
    int num_criteria;
    bool use_linear_preferences;

    std::array<unsigned long, MAX_CRITERIA> sorts{};
    std::array<unsigned long, MAX_CRITERIA> divides{};
    std::array<unsigned long, MAX_CRITERIA> marries{};
    std::array<std::array<unsigned long, MAX_CRITERIA>, 4> marray_cases{};

    DivideAndConquer(int num_criteria, bool use_linear_preferences);

    std::vector<Solution *> divide_and_conquer(std::vector<Solution *> &curr_set, int dimension, bool presorted, unsigned long &num_comparisons);
    std::vector<Solution *> divide_and_conquer_helper(std::vector<Solution *> &curr_set, unsigned long low, unsigned long high,
                                                                      int dimension);
    std::vector<Solution *> marry(std::vector<Solution *> &curr_set, unsigned long low, unsigned long high, int dimension);
    std::vector<Solution *> marry_2d(std::vector<Solution *> &curr_set);
    std::vector<Solution *> L2D(std::vector<Solution *> &all_possible);
};


#endif //DP_PARALLEL_GENERATION_DIVIDEANDCONQUER_H
