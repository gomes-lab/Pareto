//
// Created by Jonathan Gomes Selman on 8/6/17.
//

#include "Solution.h"

#include <utility>
#include <iostream>

//Solution::Solution(int node_id, const std::vector<std::pair<double, double>> &criteria,
//               const std::vector<double> &dam_decisions,
//               const std::vector<Solution *> &pareto_decisions) : node_id(node_id),
//                                                                criteria(criteria),
//                                                                dam_decisions(dam_decisions),
//                                                                pareto_decisions(pareto_decisions) {
//    solution_id = -1;
//    inferior = false;
//}

Solution::Solution(int node_id, int tree_node_id, const double *criteria_all, const std::pair<double, double> *criteria,
               std::vector<double> dam_decisions, std::vector<Solution *> pareto_decisions, int total_criteria,
               int num_criteria, bool use_linear_prefs, const std::vector<std::vector<double> > &w, std::array<double, MAX_CRITERIA> &minmax,
               std::array<double, MAX_CRITERIA> &normalizing_factors) :  node_id(node_id), tree_node_id(tree_node_id),
               dam_decisions(std::move(dam_decisions)), pareto_decisions(std::move(pareto_decisions)), num_criteria(num_criteria),
               minmax(minmax), normalizing_factors(normalizing_factors), w(w), total_criteria(total_criteria), use_linear_prefs(use_linear_prefs) {
    solution_id = -1;
    inferior = false;
    // Copy in criteria values
    for (int i = 0; i < num_criteria; i++) {
        this->criteria_2[i] = criteria[i];
    }
    for (int i = 0; i < total_criteria; i++) {
        this->criteria[i] = criteria_all[i];
    }

    update_linear_preference();

    // Compute individual hash values each criterion's rounded value
    // Hash approach courtesy of ---
    // http://stackoverflow.com/a/1646913/126995
    hash_val = 17;
    for (int i = 0; i < getNum_criteria(); i++) {
        hash_val = hash_val * 31 + std::hash<double>{}(this->criteria_2[i].second);
    }

    boost::dynamic_bitset<> rep(0, 0);
    representation = rep;
}

Solution::Solution(const Solution &src) {
    this->inferior = src.inferior;
    this->num_criteria = src.num_criteria;
    this->total_criteria = src.total_criteria;
    this->solution_id = src.solution_id;
    this->use_linear_prefs = src.use_linear_prefs;
    // Copy in criteria values
    for (int i = 0; i < num_criteria; i++) {
        this->criteria_2[i] = src.criteria_2[i];
    }
    for (int i = 0; i < total_criteria; i++) {
        this->criteria[i] = src.criteria[i];
    }
    this->node_id = src.node_id;
    this->tree_node_id = src.tree_node_id;
    this->dam_decisions = src.dam_decisions;
    this->pareto_decisions = src.pareto_decisions;
    this->minmax = src.minmax;
    this->normalizing_factors = src.normalizing_factors;
    this->w = src.w;
    this->representation = src.representation;

    update_linear_preference();
}

int Solution::is_dominant(Solution *compare) {
    bool dominate = true; // Flag to see if the current solution dominates the compared solution
    bool is_dominated = true; // Flag to see if the current solution is dominated by compare
    for (int i = 0; i < num_criteria; i++) {
        // Check to see if current solution dominates
        if (this->criteria_2[i].second > compare->criteria_2[i].second) { // We want to check some double equality issues!
            is_dominated = false;
        } else if (this->criteria_2[i].second < compare->criteria_2[i].second) {
            dominate = false;
        }
    }

    if (dominate && is_dominated) { // If both dominating and dominated then they are the same
        return 2;
    } else if(!dominate && !is_dominated) { // Neither dominates
        return 0;
    } else if (dominate) { // Current solution dominates
        return 1;
    } else { // Current solution is dominated
        return -1;
    }
}

int Solution::is_dominant_2(Solution *compare) {
    bool dominate = true; // Flag to see if the current solution dominates the compared solution
    bool is_dominated = true; // Flag to see if the current solution is dominated by compare
    for (int i = 0; i < num_criteria; i++) {
        // Check to see if current solution dominates
        if (this->criteria_2[i].second > compare->criteria_2[i].second) { // We want to check some double equality issues!
            is_dominated = false;
        } else if (this->criteria_2[i].second < compare->criteria_2[i].second) {
            dominate = false;
        }
    }

    if (dominate && is_dominated) { // If is both dominating and dominated then they are the same
        return 2;
    } else if(!dominate && !is_dominated) { // Neither dominates
        return 0;
    } else if (dominate) { // Current solution dominates
        return 1;
    } else { // Current solution is dominated
        return -1;
    }
}

int Solution::is_dominant_strict(Solution *compare) {
    bool dominant = true;
    bool equals = true;
    bool is_dominated = true;
    bool equal_in_one = false;
    // Check for strict domination
    for (int i = 0; i < num_criteria; i++) {
        if (this->criteria_2[i].second <= compare->criteria_2[i].second) {
            dominant = false;
        }

        if (this->criteria_2[i].second >= compare->criteria_2[i].second) {
            is_dominated = false;
        }
        // Check for equality
        if (this->criteria_2[i].second != compare->criteria_2[i].second) {
            equals = false;
        } else {
            equal_in_one = true;
        }
    }

    if (equals) {
        return 2;
    } else if (equal_in_one) {
        return 3;
    } else if (!dominant && !is_dominated) {
        return 0;
    } else if (dominant) {
        return 1;
    } else {
        return -1;
    }
}

int Solution::is_dominant_strong(Solution *compare) {
    bool equals = true;
    for (int i = 0; i < num_criteria; i++) {
        if (this->criteria_2[i].second != compare->criteria_2[i].second) {
            equals = false;
        }
    }
    if (equals) {
        return 2;
    }

    bool dominant = true;
    bool is_dominated = true;
    for (int i = 0; i < num_criteria; i++) {
        if (this->criteria_2[i].second == compare->criteria_2[i].second) {
            // Check the true values to break ties in random
            // Check to see if current solution dominates
            if (this->criteria_2[i].first > compare->criteria_2[i].first) {
                is_dominated = false;
            } else if (this->criteria_2[i].first < compare->criteria_2[i].first) {
                dominant = false;
            }
        } else if (this->criteria_2[i].second < compare->criteria_2[i].second) {
            dominant = false;
        } else if (this->criteria_2[i].second > compare->criteria_2[i].second) {
            is_dominated = false;
        }
    }

    if (!dominant && !is_dominated) {
        return 0;
    } else if (dominant) {
        return 1;
    } else {
        return -1;
    }
}

int Solution::is_dominant_true(Solution *compare) {
    bool dominate = true; // Flag to see if the current solution dominates the compared solution
    bool is_dominated = true; // Flag to see if the current solution is dominated by compare
    for (int i = 0; i < num_criteria; i++) {
        // Check to see if current solution dominates
        if (this->criteria_2[i].first > compare->criteria_2[i].first) {
            is_dominated = false;
        } else if (this->criteria_2[i].first < compare->criteria_2[i].first) {
            dominate = false;
        }
    }

    if (dominate && is_dominated) { // If is both dominating and dominated then they are the same
        return 2;
    } else if(!dominate && !is_dominated) { // Neither dominates
        return 0;
    } else if (dominate) {
        return 1;
    } else {
        return -1;
    }
}


bool Solution::operator==(const Solution &rhs) const {
    if (use_linear_prefs) {
        // Compare the 'rounded-value' for each linear weighted criterion
        for (int i = 0; i < w.size(); i++) {
            if (this->criteria_lp[i].second != rhs.criteria_lp[i].second) { // Watch out for doubles equals!!
                return false;
            }
        }
    }
    else {
        // Compare the 'rounded-value' for each criterion
        for (int i = 0; i < num_criteria; i++) {
            if (this->criteria_2[i].second != rhs.criteria_2[i].second) { // Watch out for doubles equals!!
                return false;
            }
        }
    }
    return true;
}

bool Solution::operator!=(const Solution &rhs) const {
    return !(rhs == *this);
}

int Solution::getNum_criteria() const {
    return num_criteria;
}

void Solution::update_linear_preference() {
    int idx_offset = 0;
    for (int idx = 0; idx < (int)w.size(); idx++) {
        double dot_product1 = 0, dot_product2 = 0;
        for (int i = idx_offset; i < idx_offset + (int)w[idx].size(); i++) {
            dot_product1 += (w[idx][i - idx_offset] * criteria_2[i].first) / normalizing_factors[i];
            dot_product2 += (w[idx][i - idx_offset] * criteria_2[i].second) / normalizing_factors[i];
        }
        idx_offset += (int)w[idx].size();
        this->criteria_lp[idx] = std::make_pair(dot_product1, dot_product2);
    }
}

Solution::~Solution() {
    pareto_decisions.clear();
}

void Solution::generate_representation() {
    int size = 0;
    for (auto p : pareto_decisions) {
        if (p != nullptr) {
            size += p->representation.size();
        }
    }

    for (auto d : dam_decisions) {
        if (d >= 0) {
            size++;
        }
    }

    boost::dynamic_bitset<> rep(size, 0);

    int j = 0;
    for (auto p : pareto_decisions) {
        if (p != nullptr) {
            for (int i = 0; i < p->representation.size(); i++) {
                rep[j] = p->representation[i];
                j++;
            }
        }
    }

    for (auto d : dam_decisions) {
        if (d >= 0) {
            rep[j] = static_cast<bool>(d);
            j++;
        }
    }

    representation = rep;
}




