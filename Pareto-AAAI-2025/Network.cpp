//
// Created by Jonathan Gomes Selman on 7/6/17.
//
// Refer to the file describing the input file for hyper_net.
// File: format-dams-v1_xiaojian.docx, found in cpp_inputs file

#include "Network.h"

#include <utility>

Network::Network(const Config &config, IOHandler &io) : config(config),  io(io) {
    switch (config.sort_type) {
        case SORT_TYPE_NONE:
            compare_operator.compare_func = TreeNode::Compare_Node_IDs;
            break;
        case SORT_TYPE_TREE:
            compare_operator.compare_func = TreeNode::Compare_Node_Size;
            break;
        case SORT_TYPE_TREE_REVERSE:
            compare_operator.compare_func = TreeNode::Compare_Node_Size_Reverse;
            break;
        case SORT_TYPE_FRONTIER:
            compare_operator.compare_func = TreeNode::Compare_Frontier_Size;
            break;
        case SORT_TYPE_FRONTIER_REVERSE:
            compare_operator.compare_func = TreeNode::Compare_Frontier_Size_Reverse;
            break;
        case SORT_TYPE_RANDOM:
            compare_operator.compare_func = TreeNode::Compare_Random;
            break;
    }
    compare_operator.resort = false;

    construct_new();
}

void Network::construct_new() {
    std::fstream input = this->io.get_file(false, DATA_FILE_NAME);

    if (input.fail()) { // Check to make sure the file was properly opened
        std::cerr << "error opening file" << std::endl;
        exit(1);
    }

    for (int i = 0; i < MAX_CRITERIA; i++) {
        max_possible[i] = 0.0;
        min_forced[i] = 0.0;
    }

    std::stringstream reader; // Used for processing lines

    // First line
    read_comments(input, reader); // Check for comments and process the flag character
    reader >> this->num_nodes >> this->num_edges >> this->num_criteria; // Assign values from first line
    this->num_dams = this->num_edges;
    this->num_nodes_total = num_nodes;

    // Criterion line
    read_comments(input, reader); // Check for comments and process the flag characters
    criteria.resize(this->num_criteria);
    for (int i = 0; i < this->num_criteria; i++) {
        std::string crit;
        reader >> crit;

        double mm = 1;
        if (crit.at(0) == '-') {
            mm = -1;
            crit = crit.substr(1);
        } else if (crit.at(0) == '+') {
            mm = 1;
            crit = crit.substr(1);
        }

        this->criteria[i] = crit;

        for (int j = 0; j < config.optimized_criteria.size(); j++) {
            if (config.optimized_criteria[j] == crit) {
                minmax[j] = mm;
            }
        }
    }

    for (int j = 0; j < config.optimized_criteria.size(); j++) {
        if (config.optimized_criteria[j] == "dcip") {
            minmax[j] = 1.0;
        }
    }

    for (int j = 0; j < config.optimized_criteria.size(); ++j) {
        min_vals[j] = std::numeric_limits<double>::max();
    }

    dcip_k = 0.0;

    // Dam Data
    readDamInfo(input);
    // Node Data
    readNodeData(input);
    // Edge List
    readAdjList(input);

    dcip_k = dcip_k * config.epsilon / (2.0 * num_nodes);

    // After reading in the data we want to update the energies for the dams
    // Using the rounding scheme proposed for AAAI paper. This is
    // the pre-rounding step.
    calculate_k();
    round_criteria(root);

    root_node = nullptr;

    generateTree();

    for (int i = 0; i < MAX_CRITERIA; i++) {
        max_normalize_r[i] = 0;
        max_normalize_s[i] = 0;
    }

    generateStatistics(root_node);

    for (int i = 0; i < MAX_CRITERIA; i++) {
        normalizing_factor[i] = std::fmax(max_normalize_r[i], max_normalize_s[i]);
    }

    if (io.config.save_progress) {
        writeTree();
    }
}


Network::~Network() {
    destructor_helper(root_node);
    root_node = nullptr;
}

Network::Network(const Network &src): num_nodes(src.num_nodes), num_edges(src.num_edges), num_dams(src.num_dams),
                                         num_criteria(src.num_criteria), root(src.root), adj_lists(src.adj_lists),
                                         nodes(src.nodes), dams(src.dams), config(src.config),
                                         minmax(src.minmax), node_criteria_indicator(src.node_criteria_indicator),
                                         min_vals(src.min_vals), normalizing_factor(src.normalizing_factor),
                                         max_normalize_s(src.max_normalize_s), max_normalize_r(src.max_normalize_r),
                                         num_nodes_total(src.num_nodes_total), dcip_k(src.dcip_k), io(src.io), compare_operator(src.compare_operator) {

    copy_constructor_helper(root_node, src.root_node);
}

void Network::copy_constructor_helper(TreeNode *&newTree, TreeNode *oldTree) {
    // Copy the oldTree into newTree
    if (oldTree != nullptr) {
        newTree = new TreeNode(*oldTree);
        // Recurse
        for (int i = 0; i < oldTree->children.size(); i++) {
            copy_constructor_helper(newTree->children[i], oldTree->children[i]);
        }
    }
}

Network &Network::operator=(const Network &src) {
    // Generate new tree
    TreeNode* tmpRoot;
    copy_constructor_helper(tmpRoot, src.root_node);

    // Copy instance variables
    num_nodes = src.num_nodes;
    num_edges = src.num_edges;
    num_dams = src.num_dams;
    num_criteria = src.num_criteria;
    root = src.root;
    adj_lists = src.adj_lists;
    nodes = src.nodes;
    dams = src.dams;
    config = src.config;
    node_criteria_indicator = src.node_criteria_indicator;
    io = src.io;
    compare_operator = src.compare_operator;

    dcip_k = src.dcip_k;

    // Free old tree
    destructor_helper(root_node);

    this->root_node = tmpRoot;
    return *this;
}

void Network::destructor_helper(TreeNode* node) {
    // Free tree with post-order traversal
    if (node != nullptr) {
        for (auto &child : node->children) {
            destructor_helper(child);
        }
        delete node;
    }
}

void Network::generateStatistics(TreeNode* node) {
    if (node == nullptr) return;

    // If we have reached a leaf stop
    node->size = 1;

    for (auto &child : node->children) {
        generateStatistics(child);

        node->size += child->size;

        for (int i = 0; i < num_criteria; i++) {
            node->max_possible[i] += child->max_possible[i] + child->parentDam.s_all[i];
            node->min_forced[i] += child->min_forced[i];
            if (child->parentDam.status == 1)
            {
                node->min_forced[i] += child->parentDam.s_all[i];
            }
        }
    }

    for (int i = 0; i < config.optimized_criteria.size(); i++) {
        max_normalize_r[i] += node->node_data.r_vals[i].first;
        max_normalize_s[i] += node->parentDam.s_vals[i].first;
    }

    for (int i = 0; i < num_criteria; i++) {
        max_possible[i] += node->node_data.r_all[i] + node->parentDam.s_all[i];
        node->max_possible[i] += node->node_data.r_all[i];

        min_forced[i] += node->node_data.r_all[i];
        if (node->parentDam.status == 1) {
            min_forced[i] += node->parentDam.s_all[i];
        }
        node->min_forced[i] += node->node_data.r_all[i];
    }

    // If more than 2 children, will eventually grow by #children - 2 intermediate nodes
    if (node->children.size() > 2) {
        node->size += (int)node->children.size() - 2;
    }
}

void Network::readAdjList(std::fstream& input) {
    // Read in the root line
    std::stringstream reader;
    read_comments(input, reader); // Check for comments and process the root flag
    reader >> root; // Get root id

    // Read through each edge
    for (int i = 0; i < num_edges; i++) {
        // Get each line to be processed
        read_comments(input, reader); // Check for comments and process the root flag

        int parent;
        Edge edge_to_child{};
        reader >> parent >> edge_to_child.endNode >> edge_to_child.damId;

        // Put the new edge into the adjacency list map
        adj_lists[parent].push_back(edge_to_child);
    }
}


void Network::readNodeData(std::fstream& input) {
    // Initialize the size of the node vector
    nodes.resize(num_nodes);
    for (int i = 0; i < num_nodes; i++) {
        Node newNode;

        // Get each line to be processed
        std::stringstream reader;
        read_comments(input, reader); // Check for comments and process the root flag

        reader >> newNode.id; // Read in the node id
        int k = 0;
        for (const std::string& criterion: this->criteria) { // Read the criterion in order dictated by file
            double r_val;
            reader >> r_val;

            newNode.r_all[k] = r_val;
            k += 1;

            for (int j = 0; j < config.optimized_criteria.size(); j++) {
                if (config.optimized_criteria[j] == criterion) {
                    newNode.r_vals[j] = {r_val, r_val};
                }
                if (config.optimized_criteria[j] == "connectivity") {
                    dcip_k += pow(r_val, 2.0);
                }
            }
        }

        // Special handling for DCIP
        for (int j = 0; j < config.optimized_criteria.size(); j++) {
            if (config.optimized_criteria[j] == "dcip") {
                newNode.r_vals[j] = {0.0, 0.0};
            }
        }

        this->nodes[newNode.id] = newNode;
    }
}

void Network::readDamInfo(std::fstream& input) {
    dams.reserve(num_dams);
    for (int i = 0; i < num_dams; i++){
        Dam newDam;

        // Get each line to be processed
        std::stringstream reader;
        read_comments(input, reader); // Check for comments and process the root flag

        reader >> newDam.id;
        reader >> newDam.status;

        int k = 0;
        for (const std::string& criterion: this->criteria) {
            double s_val, p_val, q_val;
            std::string val;
            reader >> val;

            std::string token;
            size_t pos;
            // Find the s val
            pos = val.find('|');
            token = val.substr(0, pos);
            s_val = std::stod(token);
            val.erase(0, pos + 1);

            // Find the p val
            pos = val.find('|');
            token = val.substr(0, pos);
            p_val = std::stod(token);
            val.erase(0, pos + 1);

            // The rest is the q val
            q_val = std::stod(val);

            newDam.s_all[k] = s_val;
            newDam.p_all[k] = p_val;
            newDam.q_all[k] = q_val;
            k += 1;

            for (int j = 0; j < config.optimized_criteria.size(); j++) {
                if (config.optimized_criteria[j] == criterion) {
                    newDam.s_vals[j] = {s_val, s_val};
                    newDam.p_vals[j] = p_val;
                    newDam.q_vals[j] = q_val;

                    if (newDam.s_vals[j].first != 0) {
                        min_vals[j] = std::fmin(min_vals[j], newDam.s_vals[j].first);
                    }
                }
            }
        }

        // Special handling for DCIP
        for (int j = 0; j < config.optimized_criteria.size(); j++) {
            if (config.optimized_criteria[j] == "dcip") {
                newDam.s_vals[j] = {0.0, 0.0};
                newDam.p_vals[j] = 1.0;
                newDam.q_vals[j] = 1.0;
            }
        }

        // Record the dam status (i.e. planned or already built)
        if (newDam.status == 0) {
            newDam.decision = {0, 1};
        } else if (newDam.status == 1) {
            newDam.decision = {1};
        } else {
            newDam.decision = {0};
        }
        this->dams[newDam.id] = newDam;
    }
}

void Network::calculate_k() {
    for (int j = 0; j < config.optimized_criteria.size(); ++j) {
        rounding_vals[j] = config.epsilon * min_vals[j] / 2.0;
    }
}

void Network::round_criteria(int node) {
    for (Edge& e: adj_lists[node]){
        // Round the edge value
        for (int j = 0; j < config.optimized_criteria.size(); ++j) {
            dams[e.damId].s_vals[j].second = round_energy(dams[e.damId].s_vals[j].second, rounding_vals[j]);
        }
        round_criteria(e.endNode);
    }
}

double Network::round_energy(double energy, double k_factor) {
    // Avoid the case of divide by zero
    if (k_factor == 0) {
        return energy;
    }
    return std::floor(energy / k_factor) * k_factor;
}

void Network::read_comments(std::fstream &input, std::stringstream &reader) {
    // Get the line to be checked
    std::string line;
    getline(input, line);

    // Clear and re-assign the contents of the stringstream
    reader.str(line);
    reader.clear();
    // Test the first characters to see if the line starts with '%' and is thus a comment.
    // Keep getting new lines until all comments are processed.
    // In the end the flag character(s) that appear at the beginning of the line
    // containing the data to be parsed is/are consumed.
    std::string dummy;
    reader >> dummy;
    while (dummy.compare(0, 1, "%") == 0) {
        getline(input, line);
        reader.str(line);
        reader.clear();
        reader >> dummy;
    }
}

void Network::generateTree() {
    Node root_data = nodes[root];
    int num_intermediates_root = calculate_number_intermediate_nodes(adj_lists[root].size());
    split_node_data(root_data, num_intermediates_root);

    root_node = new TreeNode(root_data);
    root_node->parentDam.status = -1;
    for (int i = 0; i < MAX_CRITERIA; i++) {
        root_node->parentDam.p_vals[i] = 0.0;
        root_node->parentDam.q_vals[i] = 0.0;
        std::pair<double, double> p = std::make_pair(0.0, 0.0);
        root_node->parentDam.s_vals[i] = p;
    }
    generate_tree_helper(root_node);
}

void Network::generate_tree_helper(TreeNode* node) {
    // If we have reached a leaf stop
    if (adj_lists[node->node_data.id].empty()) {
        return;
    }

    // Make the left child be the child at current children_index
    for (auto edge : adj_lists[node->node_data.id]) {
        Node node_data = nodes[edge.endNode];
        int num_intermediates = calculate_number_intermediate_nodes(adj_lists[edge.endNode].size());
        split_node_data(node_data, num_intermediates);

        TreeNode* new_node = new TreeNode(node_data);
        new_node->parentDam = dams[edge.damId];
        node->add_child(new_node);

        if (node->downstream_build || node->parentDam.status == 1) {
            new_node->downstream_build = true;
        }

        generate_tree_helper(new_node);

        node->height = std::max(node->height, new_node->height + 1);
    }
}

void Network::writeMeta() {
    std::fstream tree_file = this->io.get_file(true, TREE_FILE_NAME);

    tree_file << config.river_network << std::endl;
    tree_file << this->num_nodes << " " << " " << this->num_nodes_total << " "
              << num_edges << " " << num_dams
              << std::endl;
    tree_file << config.epsilon << std::endl;

    tree_file << this->num_criteria;
    for (const auto &crit : criteria) {
        tree_file << " " << crit;
    }
    tree_file << std::endl;

    tree_file << this->config.optimized_criteria.size();
    for (const auto &crit : config.optimized_criteria) {
        tree_file << " " << crit;
    }
    tree_file << std::endl;

    for (const auto &mm : minmax) {
        tree_file << mm << " ";
    }
    tree_file << std::endl;

    for (const auto &rv : rounding_vals) {
        tree_file << rv << " ";
    }
    tree_file << std::endl;

    for (const auto &nci : node_criteria_indicator) {
        tree_file << nci << " ";
    }
    tree_file << std::endl;

    for (const auto &mv : min_vals) {
        tree_file << mv << " ";
    }
    tree_file << std::endl;

    tree_file << root << std::endl;
    tree_file.close();
}

void Network::writeTree() {
    writeMeta();
    writeTreeNode(root_node, true);
}

void Network::writeTreeNode(TreeNode *node, bool writeChildren) {
    if (node == nullptr) return;

    this->io.make_node_dir(node->node_id);
    std::fstream node_file = io.get_file(true, NODE_FILE_NAME, node->node_id);
    node_file << node->node_data.id << " " << node->node_id << " " << node->intermediate << " " << node->size << " ";
    node_file << node->clocks_spent << " " << node->time_spent << " " << node->height << " " << node->num_ub_broken;
    node_file << " " << node->num_lb_broken << " ";

    node_file << node->children.size() << " ";
    for (auto &child : node->children) {
        node_file << child->node_id << " ";
    }

    node_file << node->node_data.r_vals.size() << " ";
    for (auto const &r_val : node->node_data.r_vals) {
        node_file << r_val.first << "|" << r_val.second << " ";
    }

    node_file << node->node_data.k_vals.size() << " ";
    for (auto const &k_val : node->node_data.k_vals) {
        node_file << k_val << " ";
    }
    node_file << std::endl;

    if (node->parentDam.status != -1) {
        writeDam(node->parentDam, node->node_id);
    }

    node_file.close();

    if (writeChildren) {
        for (auto &child : node->children) {
            writeTreeNode(child, writeChildren);
        }
    }
}

void Network::writeDam(const Dam& d, int node_id) {
    std::fstream dam_file = this->io.get_file(true, DAM_FILE_NAME, node_id);

    dam_file << d.id << " " << d.status << " ";
    dam_file << d.decision.size() << " ";
    for (const auto &dec : d.decision) {
        dam_file << dec << " ";
    }

    dam_file << d.s_vals.size() << " ";
    for (const auto &s_val : d.s_vals) {
        dam_file << s_val.first << "|" << s_val.second << " ";
    }

    dam_file << d.p_vals.size() << " ";
    for (const auto &p_val : d.p_vals) {
        dam_file << p_val << " ";
    }

    dam_file << d.q_vals.size() << " ";
    for (const auto &q_val : d.q_vals) {
        dam_file << q_val << " ";
    }

    dam_file << std::endl;
}

int Network::calculate_number_intermediate_nodes(int num_children) {
    return num_children - 2;
}

void Network::print_adjacent_tree(int node) {
    for (Edge& e: adj_lists[node]) {
        print_adjacent_tree(e.endNode);
    }
    std::cout << nodes[node].id << std::endl;
}

void Network::print_tree(TreeNode *node) {
    if (node == nullptr) {
        return;
    }

    for (auto &child : node->children) {
        print_tree(child);
    }

    std::cout << node->node_data.id << " is intermediate: " << node->intermediate << std::endl;
}


void Network::split_node_data(Node &data, int num_intermediates) {
    if (num_intermediates > 0) {
        for (auto & r_val : data.r_vals) {
            r_val.first = r_val.first / (double)(num_intermediates + 1);
            r_val.second = r_val.second / (double)(num_intermediates + 1);
        }

        for (auto & r_all : data.r_all) {
            r_all = r_all / (double)(num_intermediates + 1);
        }
    }
}
