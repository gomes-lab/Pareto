//
// Created by Jonathan Gomes Selman on 7/6/17.
//

#ifndef AMAZON_PROJECT_HYPERNET_H
#define AMAZON_PROJECT_HYPERNET_H

#include <utility>
#include <vector>
#include <unordered_map>
#include "string"
#include <fstream>
#include <iostream>
#include <sstream>
#include "Solution.h"
#include <algorithm>
#include "Node.h"
#include "Consts.h"
#include "IOHandler.h"
#include "Config.h"

// Represents a Node in the tree representation of the network
struct TreeNode {
    int node_id;
    int rand_id; // Used only for random sort
    Node node_data; // The actual data for the tree node

    bool downstream_build;

    std::vector<TreeNode*> children;

    std::array<double, MAX_CRITERIA> max_possible{};
    std::array<double, MAX_CRITERIA> min_forced{};

    int size;
    int height;

    // Incoming Dam, if it exists
    Dam parentDam;

    bool intermediate; // Indicates whether the node is an intermediate node --- created to make the tree binary
    bool policies_complete; // Policies have been completed
    bool node_complete;
    bool sorted; // Flag to determine if this has been sorted once before

    unsigned long max_frontier_size;
    unsigned long num_ub_broken;
    unsigned long num_lb_broken;
    unsigned long long clocks_spent;
    double time_spent;

    std::vector<Solution*> vec_frontier;

    bool is_leaf;

    // Compares Node size as defined by the size of the subtree rooted at that node
    static bool Compare_Node_Size(TreeNode* t1, TreeNode* t2) {
        return t1->size < t2->size;
    }

    // Compares the Frontier size of the nodes
    static bool Compare_Frontier_Size(TreeNode* t1, TreeNode* t2) {
        return t1->vec_frontier.size() * t1->parentDam.decision.size() < t2->vec_frontier.size() * t2->parentDam.decision.size();
    }

    // Compares Node IDs and keeps it in the defined node order
    static bool Compare_Node_IDs(TreeNode* t1, TreeNode* t2) {
        return t1->node_id < t2->node_id;
    }

    static bool Compare_Random(TreeNode* t1, TreeNode* t2) {
        if (t1->rand_id == -1) {
            t1->rand_id = rand();
        }
        if (t2->rand_id == -1) {
            t2->rand_id = rand();
        }

        return t1->rand_id < t2->rand_id;
    }

    static bool Compare_Node_Size_Reverse(TreeNode* t1, TreeNode* t2) {
        return Compare_Node_Size(t2, t1);
    }

    static bool Compare_Frontier_Size_Reverse(TreeNode* t1, TreeNode* t2) {
        return Compare_Frontier_Size(t2, t1);
    }

    bool has_next() {
        bool hn = !node_complete;

        if (children.size() <= 2) {
            node_complete = true;
        }

        return hn;
    }

    struct Post_Sort_Operator {
        virtual void post_sort(std::vector<TreeNode*> &sorted_children) = 0;
    };

    struct Shift_Index_Post_Sort : Post_Sort_Operator {
        int index;
        int shift;

        void post_sort(std::vector<TreeNode*> &sorted_children) override {
            if (index < sorted_children.size()) {
                if (shift > 0) {
                    unsigned long shift_end = std::max(sorted_children.size() - (shift + index + 1), (unsigned long) 0);
                    for (unsigned long i = sorted_children.size() - index - 1; i > shift_end; i--) {
                        TreeNode* tmp = sorted_children[i];
                        sorted_children[i] = sorted_children[i-1];
                        sorted_children[i - 1] = tmp;
                    }
                }
                else if (shift < 0) {
                    unsigned long shift_end = std::min(sorted_children.size() - (shift + index + 1), sorted_children.size() - 1);
                    for (unsigned long i = sorted_children.size() - index - 1; i < shift_end; i++) {
                        TreeNode* tmp = sorted_children[i];
                        sorted_children[i] = sorted_children[i+1];
                        sorted_children[i + 1] = tmp;
                    }
                }
            }
        }

        Shift_Index_Post_Sort(int index, int shift) : index(index), shift(shift) {

        }
    };

    struct LKO_Post_Sort : Post_Sort_Operator {
        std::vector<int> indices; // List of indices in the original sorted list to "leave out" by placing at the end

        void post_sort(std::vector<TreeNode*> &sorted_children) override {
            for (int i = 0; i < indices.size(); i++) {
                int idx = indices[indices.size() - 1 - i];
                if (idx < sorted_children.size()) {
                    for (unsigned long k = sorted_children.size() - idx - 1; k > 0; k--) {
                        TreeNode* tmp = sorted_children[k];
                        sorted_children[k] = sorted_children[k-1];
                        sorted_children[k-1] = tmp;
                    }
                }
            }
        }

        explicit LKO_Post_Sort(std::vector<int> indices) : indices(std::move(indices)) {

        }
    };

    struct Compare_Operator {
        std::function<bool(TreeNode*, TreeNode*)> compare_func;
        Post_Sort_Operator *post_sort = nullptr;
        bool resort; // Re-sort the children for intermediate nodes
    };

    TreeNode* get_next(int &next_id, Compare_Operator compare) {
        // More than 2 children, need to make an intermediate node out of the first two
        if (children.size() > 2) {
            if (!sorted || (sorted && compare.resort)) {
                std::sort(children.begin(), children.end(), std::move(compare.compare_func));

                if (nullptr != compare.post_sort) {
                    compare.post_sort->post_sort(children);
                }

                sorted = true;
            }

            TreeNode* right = children.back();
            children.pop_back();
            TreeNode* left = children.back();
            children.pop_back();
            auto* intermediate_node = new TreeNode(node_data, next_id++);
            intermediate_node->parentDam.status = 2;
            intermediate_node->parentDam.decision = {0};
            for (double & q_val : intermediate_node->parentDam.q_vals) {
                q_val = 1;
            }
            for (double & q_all : intermediate_node->parentDam.q_all) {
                q_all = 1;
            }

            intermediate_node->downstream_build = downstream_build;
            intermediate_node->add_child(left);
            intermediate_node->add_child(right);
            intermediate_node->height = height;
            this->height++;

            for (int i = 0; i < MAX_CRITERIA; i++) {
                intermediate_node->max_possible[i] = left->max_possible[i] + right->max_possible[i]
                         + left->parentDam.s_all[i] + right->parentDam.s_all[i]
                         + intermediate_node->node_data.r_all[i];
                intermediate_node->min_forced[i] = left->min_forced[i] + right->min_forced[i]
                         + left->parentDam.s_all[i] + right->parentDam.s_all[i]
                         + intermediate_node->node_data.r_all[i];
            }

            intermediate_node->size = left->size + right->size + 1;

            add_child(intermediate_node);

            return intermediate_node;
        }

        return this;
    }

    std::pair<TreeNode*, TreeNode*> get_children() {
        std::pair<TreeNode*, TreeNode*> child_pair = {nullptr, nullptr};

        // With at most 2 children left, we will be done after this iteration
        if (children.size() > 2) {
            std::cerr << "Error getting children of node " << node_id << ". Attempted to get children while there are more than two." << std::endl;
            exit(1);
        }

        if (!children.empty()) {
            child_pair.first = children.front();
        }
        if (children.size() >= 2) {
            child_pair.second = children.back();
        }

        return child_pair;
    }

    void add_child(TreeNode* child) {
        children.emplace_back(child);
        is_leaf = false;
    }

    explicit TreeNode(const Node& node) {
        node_id = node.id;
        node_data = node;
        intermediate = false;
        is_leaf = true;
        size = 1;
        rand_id = -1;
        clocks_spent = 0;
        time_spent = 0;
        num_ub_broken = 0;
        num_lb_broken = 0;

        max_frontier_size = 0;
        policies_complete = false;
        node_complete = false;
        sorted = false;
        height = 0;
        downstream_build = false;

        for (int i = 0; i < MAX_CRITERIA; i++) {
            max_possible[i] = 0.0;
            min_forced[i] = 0.0;
        }
    }

    TreeNode(const Node& node, int id) {
        node_id = id;
        node_data = node;
        intermediate = true;
        is_leaf = true;
        size = 1;
        rand_id = -1;
        clocks_spent = 0;
        time_spent = 0;
        num_ub_broken = 0;
        num_lb_broken = 0;

        max_frontier_size = 0;
        policies_complete = false;
        node_complete = false;
        sorted = false;
        height = 0;
        downstream_build = false;

        for (int i = 0; i < MAX_CRITERIA; i++) {
            max_possible[i] = 0.0;
            min_forced[i] = 0.0;
        }
    }

    ~TreeNode() {
        for (auto & p : vec_frontier) {
            delete(p);
        }
    }

    TreeNode(const TreeNode& src) {
        node_id = src.node_id;
        node_data = src.node_data;
        children = src.children;
        parentDam = src.parentDam;
        intermediate = src.intermediate;
        vec_frontier = src.vec_frontier;
        max_frontier_size = src.max_frontier_size;
        policies_complete = src.policies_complete;
        node_complete = src.node_complete;
        is_leaf = src.is_leaf;
        size = src.size;
        rand_id = src.rand_id;
        clocks_spent = src.clocks_spent;
        time_spent = src.time_spent;
        sorted = src.sorted;
        height = src.height;
        num_ub_broken = src.num_ub_broken;
        num_lb_broken = src.num_lb_broken;
        downstream_build = src.downstream_build;

        for (int i = 0; i < MAX_CRITERIA; i++) {
            max_possible[i] = src.max_possible[i];
            min_forced[i] = src.min_forced[i];
        }
    }
};

class Network {
public:

    /*
     * Constructor used by the most current version
     */
    Network(const Config &config, IOHandler &io);

    virtual ~Network(); // Destructor

    /*
    * Follow the rule of 3:
    *      Destructor
    *      Copy constructor
    *      Assignment operator
    */
    void destructor_helper(TreeNode* node);
    Network(const Network& src); // Copy constructor
    Network& operator=(const Network& src); // Assignment operator

    /*
     * Read the adjacent list from the input file
     */
    void readAdjList(std::fstream& input);

    /*
     * Read the data for each hyper_node
     */
    void readNodeData(std::fstream& input);


    /*
     * Read info for each dam
     */
    void readDamInfo(std::fstream& input);

    void generateStatistics(TreeNode* node);

    void generateTree();

    void writeMeta();
    void writeTree();
    void writeTreeNode(TreeNode* node, bool writeChildren);

    void writeDam(const Dam& d, int node_id);

    /*
     * Print the Hyper-Tree in-order
     */
    void print_adjacent_tree(int node);

    /*
     * Print the Hyper-Tree in-order
     */
    void print_tree(TreeNode* root);

    /*
     * Calculates the rounding constants for the
     * preliminary rounding of energy and seismic risk.
     * ------------------------------------------------
     * These rounding constants are determined by the minimal
     * criteria value for energy and seismic risk.
     * -------------------------------------------
     * The formula is:
     *      epsilon * minimum_criteria_val / 2
     */
    void calculate_k();


    /*
     * Round all of the dam energy and seismic risk value
     * based on the respective rounding constants. This
     * is the preliminary rounding phase before the
     * dynamic programming algorithm is run.
     */
    void round_criteria(int node);

    /*
     * Helper method to round specific values:
     *      floor ( energy / k_factor) * k_factor
     * Note: If the k_factor is zero, the value of
     * 'energy' is returned.
     */
    double round_energy(double energy, double k_factor);

    // Data structures for different attributes of Hyper Tree
    std::unordered_map<int, std::vector<Edge> > adj_lists; // Adjacency list of edges for each node
    std::vector<Node> nodes; // Info about each node
    std::unordered_map<int, Dam> dams; // Info for dams
    std::unordered_map<int, TreeNode*> hyper_nodes;

    std::array<double, MAX_CRITERIA> max_possible{};
    std::array<double, MAX_CRITERIA> min_forced{};

    // Basic data about the tree
    int num_nodes;
    int num_edges;
    int num_dams;
    int num_criteria;
    std::vector<std::string> criteria;
    int root; // Defines the root of the tree

    double dcip_k;

    TreeNode* root_node{};

    int num_nodes_total;
    bool normalized = false;

    TreeNode::Compare_Operator compare_operator;

    IOHandler io;
    Config config;

    std::array<double, MAX_CRITERIA> minmax{};

    std::array<double, MAX_CRITERIA> rounding_vals{};
    std::array<double, MAX_CRITERIA> node_criteria_indicator{};
    std::array<double, MAX_CRITERIA> min_vals{};

    std::array<double, MAX_CRITERIA> max_normalize_r{};
    std::array<double, MAX_CRITERIA> max_normalize_s{};
    std::array<double, MAX_CRITERIA> normalizing_factor{};

private:

    /*
     * Calculates the number of intermediate nodes that will be needed to represent
     * a parent node with a given number of children in a binary tree. The number
     * of intermediate nodes is num_children - 2.
     */
    static int calculate_number_intermediate_nodes(int num_children);

    // Helper method to perform the operations of the default constructor
    void construct_new();

    void generate_tree_helper(TreeNode* node);

    /*
     * Used to copy the tree into a new Network object
     */
    void copy_constructor_helper(TreeNode*& newTree, TreeNode* oldTree);

    /*
     * Updates the criteria of the Node to reflect the split
     * that is being made when introducing intermediary nodes.
     * For example, if the number of intermediates is 5, we
     * will update a connectivity value of 12 as follows:
     *      connectivity = 12 / (5 + 1)
     * NOTE: we add one to the number of intermediates to include the
     * original parent.
     */
    void split_node_data(Node& data, int num_intermediates);

    /*
     * This method tests to see if a line is a comment line. If the line is a comment line then it is
     * skipped until a non-comment line is reached. This function then consumes the flag characters
     * of that line and returns with the stringstream holding the contents of that line.
     */
    void read_comments(std::fstream& input, std::stringstream& reader);

};

#endif //AMAZON_PROJECT_HYPERNET_H
