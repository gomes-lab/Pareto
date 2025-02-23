//
// Created by marc on 11/16/21.
//

#include "Transform.h"

Transform::Transform(TreeNode* parentNode, TreeNode* inner_node, Solution* outer_solution, const Dam &outer_dam, const Dam &inner_dam,
                           double outer_dec, double inner_dec, int num_criteria, double k_epsilon,
                           const std::vector<std::vector<double> > &w, int batch_size, Network *net, DivideAndConquer *dnc, bool verbose,
                           Config config)
                           : parentNode(parentNode), inner_node(inner_node), outer_solution(outer_solution),
                           outer_dam(outer_dam), inner_dam(inner_dam), outer_dec(outer_dec), inner_dec(inner_dec), num_criteria(num_criteria),
                           k_epsilon(k_epsilon), batch_size(batch_size), net(net),
                           dnc(dnc), w(w), verbose(verbose), inferior(false), config(config) {

    decisionVec = {outer_dec, inner_dec};

    for (int i = 0; i < num_criteria; i++) {
        built_criteria[i] = this->outer_dec * this->outer_dam.s_vals[i].first + this->inner_dec * this->inner_dam.s_vals[i].first;

        k_vals[i] = this->parentNode->node_data.r_vals[i].first * k_epsilon + std::floor(
                (built_criteria[i] + std::numeric_limits<double>::epsilon()) / net->min_vals[i])
                        * net->min_vals[i] * k_epsilon / 2.0;

        base_criteria[i].first = parentNode->node_data.r_vals[i].first
                                 + this->outer_dec * (this->outer_dam.s_vals[i].first + this->outer_solution->criteria_2[i].first * this->outer_dam.p_vals[i])
                                 + (1 - this->outer_dec) * this->outer_solution->criteria_2[i].first * this->outer_dam.q_vals[i]
                                 + this->inner_dec * inner_dam.s_vals[i].first;
        base_criteria[i].second = parentNode->node_data.r_vals[i].second
                                  + this->outer_dec * (this->outer_dam.s_vals[i].second + std::abs(this->outer_solution->criteria_2[i].second) * this->outer_dam.p_vals[i])
                                  + (1 - this->outer_dec) * std::abs(this->outer_solution->criteria_2[i].second) * this->outer_dam.q_vals[i]
                                  + this->inner_dec * inner_dam.s_vals[i].second;

        base_criteria_minmax[i] = this->net->minmax[i] * base_criteria[i].first;
    }

    // If we have the DCIP, update it specially
    if (config.dcip_index >= 0) {
        base_criteria[config.dcip_index].first +=
                outer_dec * pow(outer_solution->criteria_2[config.connectivity_index].first, 2.0);

        base_criteria[config.dcip_index].second += outer_dec * pow(outer_solution->criteria_2[config.connectivity_index].second, 2.0);

        base_criteria_minmax[config.dcip_index] = base_criteria[config.dcip_index].first;

        k_vals[config.dcip_index] = net->dcip_k;
    }

    for (int i = 0; i < net->num_criteria; i++) {
        base_criteria_all[i] = parentNode->node_data.r_all[i]
                               + this->outer_dec * (this->outer_dam.s_all[i] + this->outer_solution->criteria[i] * this->outer_dam.p_all[i])
                               + (1 - this->outer_dec) * this->outer_solution->criteria[i] * this->outer_dam.q_all[i]
                               + this->inner_dec * inner_dam.s_all[i];
    }
}

Transform::Transform(TreeNode* parentNode, TreeNode* inner_node, const Dam &inner_dam,
                           double inner_dec, int num_criteria, double k_epsilon,
                           const std::vector<std::vector<double> > &w, int batch_size, Network *net, DivideAndConquer *dnc,
                           bool verbose, Config config)
        : parentNode(parentNode), inner_node(inner_node), outer_solution(nullptr),
          inner_dam(inner_dam), outer_dec(-1), inner_dec(inner_dec), num_criteria(num_criteria),
          k_epsilon(k_epsilon), batch_size(batch_size), net(net),
          dnc(dnc), w(w), verbose(verbose), inferior(false), config(config) {

    decisionVec = {inner_dec, -1};

    for (int i = 0; i < num_criteria; i++) {
        built_criteria[i] = this->inner_dec * this->inner_dam.s_vals[i].first;

        k_vals[i] = this->parentNode->node_data.r_vals[i].first * k_epsilon + std::floor(
                (built_criteria[i] + std::numeric_limits<double>::epsilon()) / net->min_vals[i])
                        * net->min_vals[i] * k_epsilon / 2.0;

        base_criteria[i].first = this->parentNode->node_data.r_vals[i].first + this->inner_dec * inner_dam.s_vals[i].first;
        base_criteria[i].second = this->parentNode->node_data.r_vals[i].second + this->inner_dec * inner_dam.s_vals[i].second;

        base_criteria_minmax[i] = this->net->minmax[i] * base_criteria[i].first;
    }

    // If we have the DCIP, update it specially
    if (config.dcip_index >= 0) {
        base_criteria_minmax[config.dcip_index] = base_criteria[config.dcip_index].first;

        k_vals[config.dcip_index] = net->dcip_k;
    }

    for (int i = 0; i < net->num_criteria; i++) {
        base_criteria_all[i] = parentNode->node_data.r_all[i] + this->inner_dec * inner_dam.s_all[i];
    }
}

std::vector<Solution *> Transform::get_grouping(int j) {
    std::vector<Solution *> grouping;
    if (outer_dec == -1) {
        grouping = {inner_node->vec_frontier[j], nullptr};
    }
    else {
        grouping = {outer_solution, inner_node->vec_frontier[j]};
    }

    return grouping;
}

void Transform::generate_criteria(std::pair<double, double> *criteria, double *criteria_all, int j) {
    initialize_criteria(criteria, criteria_all);
    update_criteria(inner_node->vec_frontier[j], criteria, criteria_all);
    round_criteria(criteria); // Round the criteria values based on local ks
}

bool Transform::add_solution(Solution* adding, std::unordered_set<Solution *, SolutionHash, EqualSolutions> &solution_set,
                                int &num_generated, std::atomic<unsigned int> &batch_num) {
    auto find = solution_set.find(adding);
    if (find != solution_set.end()) { // Solution tie exists!
        int rand_num = rand() % 2; // Only want to do if epsilon != 0 i.e. we are rounding!!!!

        if (rand_num == 0) { // Remove the current and replace it with the new node
            Solution *old = *find;

            solution_set.erase(find);
            delete (old);
            solution_set.insert(adding);
        } else {
            delete (adding);
        }
        return false;
    } else {
        solution_set.insert(adding);
        num_generated++;

        return true;
    }
}

void Transform::initialize_criteria(std::pair<double, double> *criteria, double *criteria_all) {
    for (int i = 0; i < num_criteria; i++) {
        criteria[i] = base_criteria[i];
    }

    for (int i = 0; i < net->num_criteria; i++) {
        criteria_all[i] = base_criteria_all[i];
    }
}

void Transform::update_criteria(Solution* inner_solution, std::pair<double, double> *criteria,
                                double *criteria_all) {
    for (int i = 0; i < num_criteria; i++) {
        criteria[i].first += inner_dec * (inner_solution->criteria_2[i].first * inner_dam.p_vals[i]) +
                             (1 - inner_dec) * inner_solution->criteria_2[i].first * inner_dam.q_vals[i];

        criteria[i].second += inner_dec * (std::abs(inner_solution->criteria_2[i].second) * inner_dam.p_vals[i]) +
                              (1 - inner_dec) * std::abs(inner_solution->criteria_2[i].second) * inner_dam.q_vals[i];
    }

    for (int i = 0; i < net->num_criteria; i++) {
        criteria_all[i] += inner_dec * (inner_solution->criteria[i] * inner_dam.p_all[i]) +
                (1 - inner_dec) * inner_solution->criteria[i] * inner_dam.q_all[i];
    }

    // If we have the DCIP, update it specially
    if (config.dcip_index >= 0) {
        criteria[config.dcip_index].first +=
                inner_dec * pow(inner_solution->criteria_2[config.connectivity_index].first, 2.0);

        criteria[config.dcip_index].second += inner_dec * pow(inner_solution->criteria_2[config.connectivity_index].second, 2.0);

        criteria_all[net->num_criteria] = criteria[config.dcip_index].first;
    }
}

void Transform::round_criteria(std::pair<double, double> *criteria) {
    for (int i = 0; i < num_criteria; i++) {
        criteria[i].second = net->minmax[i] * round_with_k(criteria[i].second, k_vals[i]);
    }
}

double Transform::round_with_k(double value, double k) {
    if (k == 0) {
        return value;
    }

    return std::floor((value + std::numeric_limits<double>::epsilon()) / k) * k;
}