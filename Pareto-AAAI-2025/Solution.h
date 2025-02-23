//
// Created by Jonathan Gomes Selman on 8/6/17.
//

#ifndef AMAZON_PROJECT_PARETO_SOLUTION_H
#define AMAZON_PROJECT_PARETO_SOLUTION_H


#include <utility>
#include <vector>
#include "string"
#include "Node.h"
#include <functional>
#include <cmath>
#include <boost/dynamic_bitset.hpp>

class Solution {
public:
    // Class variables
    int node_id;
    int tree_node_id;
    int solution_id;
    std::size_t hash_val;

    // Holds the values of the criteria for a given solution.
    // The order of the criteria is determined by the order given
    // by the Network.
    // Ex. if the algorithm is running on the criteria 'conn' and 'energy'
    //      criteria[0] --> connectivity
    //      criteria[1] --> energy
    // Note: std::pair is used to hold the true value and rounded value of a criteria.
    // The data is stored as <true_value, rounded_value>. If we choose not to round then
    // rounded_value is just equal to the true_value.
//    std::vector<std::pair<double, double> > criteria; // NOTE: to allow generic types we should use std::variant

    // Static array representation of the criteria values for a given solution.
    // A static array is initialized to the size of MAX_CRITERIA, even though
    // the number of filled indexes is equal to the number of criteria being
    // considered for a given experimental run. A static array is used to avoid
    // the overhead of dynamically allocated memory on the heap that is needed
    // when working with std::vector. The static array is allocated on the stack
    // and allows for much quicker allocation and data retrieval / modification.
    double criteria[MAX_CRITERIA];
    double criteria_normalized[MAX_CRITERIA];
    std::pair<double, double> criteria_2[MAX_CRITERIA];
    std::pair<double, double> criteria_lp[MAX_CRITERIA];
    std::array<double, MAX_CRITERIA> minmax{};
    std::array<double, MAX_CRITERIA> normalizing_factors{};
    std::vector<std::vector<double>> w;

    std::vector<double> dam_decisions; // Array of 0 and 1 showing if each dam was built or not

    boost::dynamic_bitset<> representation;

    bool inferior; // Used in the divide-and-conquer algorithm. Signals if the solution is in the inferior or superior category

    // Represents an array of pointers to the Pareto_Solutions that was used for each
    // child to generate the new Solution
    std::vector<Solution*> pareto_decisions;

    // End class variables

    /*
     * Constructor used for creating a solution utilizing the static array for storing criterion values
     */
    Solution(int node_id, int tree_node_id, const double criteria_all[], const std::pair<double, double> criteria[], std::vector<double> dam_decisions,
           std::vector<Solution *> pareto_decisions, int total_criteria, int num_criteria, bool use_linear_prefs, const std::vector<std::vector<double> > &w,
           std::array<double, MAX_CRITERIA> &minmax, std::array<double, MAX_CRITERIA> &normalizing_factors);

    /*
     * Copy constructor --- Copy all data except pointer to next
     */
    Solution(const Solution& src);

    ~Solution();

    /*
     * Used with the vector implementation for storing criteria
     * Compares to Pareto_Solutions to see if one dominates the other.
     * Specifically returns:
     *      1 --> if the current solutions dominates 'compare'
     *      2 --> if the solutions are equal
     *      0 --> if the solutions do not dominate eachother
     *      -1 --> if the current solution is dominated by 'compare'
     */
    int is_dominant(Solution* compare);

    /*
     * Used with the static array implementation for storing criteria
     * See documentation for 'is_dominant' to see implementation details
     */
    int is_dominant_2(Solution* compare);

    /*
     * See if the current solution is strictly dominated by compare.
     * Additional flag (3) returned if solutions are equal in one category.
     *
     */
    int is_dominant_strict(Solution* compare);

    /*
     * If ties exist between rounded criteria values, the true
     * solution values are compared.
     */
    int is_dominant_strong(Solution* compare);

    /*
     * Compares non-rounded values of solutions
     */
    int is_dominant_true(Solution* compare);

    /*
     * Compares solutions based on the criteria values. Specifically, compares
     * the 'rounded_value' for each criteria. Returns true is the compared solution
     * has equal values for each rounded criterion value.
     */
    bool operator==(const Solution &rhs) const;

    /*
     * Uses equal operator to evaluate not-equals
     */
    bool operator!=(const Solution &rhs) const;

    /*
     * Gets the number of criteria that is actually being considered
     */
    int getNum_criteria() const;

    void update_linear_preference();
    void generate_representation();

    int num_criteria;
    int total_criteria;
private:
    bool use_linear_prefs;
    static constexpr double EQUALITY_EPSILON = 0.0001;

};

/*
 * defines the hash protocol for a Solution object
 */
namespace std
{
    template <>
    struct hash<Solution* >
    {
        size_t operator()( const Solution* sol ) const
        {
            // Compute individual hash values each criterion's rounded value
            // Hash approach courtesy of ---
            // http://stackoverflow.com/a/1646913/126995
            size_t res = 17;
            for (int i = 0; i < sol->getNum_criteria(); i++) {
                res = res * 31 + std::hash<double>{}(sol->criteria_2[i].second);
            }
            return res;
        }
    };
}

/*
 * Used to compare solutions based on a given dimension.
 * Ties are broken by comparing solutions lexicographically.
 * Used to sort partial policies
 */
struct Compare_Dimensions {
    int dimension;
    int total_dimensions;
    bool use_linear_prefs;

    explicit Compare_Dimensions(int dimension, int total_dimensions, bool use_linear_prefs) {
        this->dimension = dimension;
        this->total_dimensions = total_dimensions;
        this->use_linear_prefs = use_linear_prefs;
    }

    bool operator () (Solution* sol1, Solution* sol2) const {
        std::pair<double, double>* crit1;
        std::pair<double, double>* crit2;
        if (use_linear_prefs) {
            crit1 = sol1->criteria_lp;
            crit2 = sol2->criteria_lp;
        }
        else {
            crit1 = sol1->criteria_2;
            crit2 = sol2->criteria_2;
        }

        // Because arrays are 0 based index we subtract from the dimension
        for (int i = dimension - 1; i >= 0; i--) {
            // If the solutions are equal in the given dimension, try next dimension
            if (crit1[i].second == crit2[i].second) {
                continue;
            }
            return crit1[i].second > crit2[i].second; // Sort in descending order
        }
        // We will never have solutions equal in all dimensions so now we compare across all dimensions
        // Resolve tie in all lower dimensions by comparing upper dimensions
        int tie = dimension;
        while (tie <= total_dimensions - 1) {
            if (crit1[tie].second == crit2[tie].second) {
                tie++;
            } else {
                return crit1[tie].second > crit2[tie].second;
            }
        }
        return false;
    }
};

struct EqualSolutions{
    bool operator()(Solution const* a, Solution const* b) const {
        return *a == *b;
    }
};

struct SolutionHash {
    std::size_t operator()(const Solution* solution) const {
        return solution->hash_val;
    }
};

#endif //AMAZON_PROJECT_PARETO_SOLUTION_H
