//
// Created by marc on 12/3/20.
//

#ifndef DP_PARALLEL_GENERATION_NODE_H
#define DP_PARALLEL_GENERATION_NODE_H

#include <vector>
#include <array>
#include "Consts.h"

// Represents an edge in the tree
// Used by the adjacency list implementation
struct Edge {
    int endNode;
    int damId;
};

// Represents a Hyper Node in the tree
struct Node {
    std::array<std::pair<double, double>, MAX_CRITERIA> r_vals;
    std::array<double, MAX_CRITERIA> r_all{};
    std::array<double, MAX_CRITERIA> k_vals{};

    int id{};
};

// Represents a Dam in the Tree
struct Dam {
    std::array<std::pair<double, double>, MAX_CRITERIA> s_vals;
    std::array<double, MAX_CRITERIA> p_vals{};
    std::array<double, MAX_CRITERIA> q_vals{};

    std::array<double, MAX_CRITERIA> s_all{};
    std::array<double, MAX_CRITERIA> p_all{};
    std::array<double, MAX_CRITERIA> q_all{};

    int status{};
    int id{};

    // Holds the possible decisions for the dam
    //      {0,1} = planned dam that we are analyzing
    //      {1} = dam is already built
    //      {0} = this is a dummy dam that exists in the binary-tree representation of the graph and is never built!
    std::vector<double> decision;
};

#endif //DP_PARALLEL_GENERATION_NODE_H
