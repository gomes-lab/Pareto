//
// Created by marc on 7/22/21.
//

#include "DivideAndConquer.h"

DivideAndConquer::DivideAndConquer(int num_criteria, bool use_linear_preferences) :
    num_criteria(num_criteria), use_linear_preferences(use_linear_preferences) {
    for (int i = 0; i < num_criteria; i++) {
        sorts[i] = 0;
        divides[i] = 0;
        marries[i] = 0;
    }
}


std::vector<Solution *> DivideAndConquer::divide_and_conquer(std::vector<Solution *> &curr_set, int dimension, bool presorted,
                                                           unsigned long &num_comparisons) {
    std::vector<Solution* > r;

    // In case there are no valid solutions in the given set, just return an empty solution set
    if (curr_set.size() == 0) {
        return r;
    }

    if (!presorted) {
        num_comparisons += (unsigned long) curr_set.size() * (unsigned long) log2(curr_set.size());
    }

    if (dimension == 2) {
        return L2D(curr_set);
    }

    if (!presorted) {
        // Sort the vector by the first criteria and then call helper function
        std::sort(curr_set.begin(), curr_set.end(),
                  Compare_Dimensions(dimension, num_criteria, use_linear_preferences));
        sorts[dimension - 1]++;
    }

    //openmp doesn't work here since curr_set is reference
    r = divide_and_conquer_helper(curr_set, 0, curr_set.size() - 1, dimension);

    return r;
}

std::vector<Solution *> DivideAndConquer::divide_and_conquer_helper(std::vector<Solution *> &curr_set, unsigned long low, unsigned long high,
                                                          int dimension) {
    if (high - low == 0) {
        // We want to return the single value
        return std::vector<Solution* > (curr_set.begin() + low, curr_set.begin() + high + 1);
    }

    divides[dimension - 1]++;

    unsigned long midpoint = (high + low) / 2;

    std::vector<Solution* > superior = divide_and_conquer_helper(curr_set, low, midpoint, dimension);
    std::vector<Solution* > inferior = divide_and_conquer_helper(curr_set, midpoint + 1, high, dimension);

    std::vector<Solution* > sup_inf;
    sup_inf.reserve(superior.size() + inferior.size());

    // Copy those from superior
    for (Solution* sol: superior) {
        sol->inferior = false;
        sup_inf.push_back(sol);
    }
    // Copy from inferior
    for (Solution* sol: inferior) {
        sol->inferior = true;
        sup_inf.push_back(sol);
    }

    // Sort by d - 1 before calling marry
    std::sort(sup_inf.begin(), sup_inf.end(), Compare_Dimensions(dimension - 1, num_criteria, use_linear_preferences));
    sorts[dimension - 1]++;

    std::vector<Solution* > keep_inferior = marry(sup_inf, 0, sup_inf.size() - 1, dimension - 1);
    // Take the union of superior and the inferior points kept
    superior.insert(superior.end(), keep_inferior.begin(), keep_inferior.end());

    return superior;

}

std::vector<Solution *>
DivideAndConquer::marry(std::vector<Solution *> &curr_set, unsigned long low, unsigned long high, int dimension) {
    marries[dimension - 1]++;
    if (dimension == 2) { // We can use marry_2D to quickly find the optimal pairs
        return marry_2d(curr_set);
    }
    if (curr_set.empty()){
        return curr_set;
    }
    // Create divide plane and count the inferior and superior in each section
    unsigned long median = (low + high) / 2;
    // We need to keep the superior points from the first half for future use
    std::vector<Solution* > sup_first;
    // Lets count the number in each as a base case!
    unsigned long sup1 = 0, inf1 = 0, sup2 = 0, inf2 = 0;
    // First half
    for (unsigned long i = low; i <= median; i++) {

        if (curr_set[i]->inferior) {
            inf1++;
        } else {
            sup1++;
            sup_first.push_back(curr_set[i]);
        }

    }
    // Second half
    for (unsigned long i = median + 1; i <= high; i++) {
        if (curr_set[i]->inferior) {
            inf2++;
        }
        else {
            sup2++;
        }
    }

    // Check base cases!!
    if (inf1 + inf2 == high - low + 1) { // All inferior so we return all!
        return std::vector<Solution* > (curr_set.begin() + low, curr_set.begin() + high + 1);
    }
    else if (sup1 + sup2 == high - low + 1) { // All superior so return none
        return {};
    }
    else if (inf1 == 0 && sup2 == 0) { // All superior points remain superior and inf remain inferior -- drop dimension
        std::vector<Solution* > immediate_drop(curr_set.begin() + low, curr_set.begin() + high + 1);
        // Sort and drop dimension on curr set
        std::sort(immediate_drop.begin(), immediate_drop.end(), Compare_Dimensions(dimension - 1, num_criteria, use_linear_preferences));
        sorts[dimension - 1]++;
        return marry(immediate_drop, 0, immediate_drop.size() - 1, dimension - 1);
    }
    else if (sup1 == 0 && inf2 == 0){ // Inferior points are in first half and superior in second
        return std::vector<Solution* > (curr_set.begin() + low, curr_set.begin() + median + 1);
    }

    std::vector<Solution* > inf_first = marry(curr_set, low, median, dimension);
    std::vector<Solution* > inf_second = marry(curr_set, median + 1, high, dimension);

    // Now we must merge the sup_first with inf_second by looking at d - 1
    // Take the union of the two and sort
    sup_first.insert(sup_first.end(), inf_second.begin(), inf_second.end());
    // Sort by d - 1 before calling marry
    std::sort(sup_first.begin(), sup_first.end(), Compare_Dimensions(dimension - 1, num_criteria, use_linear_preferences));
    sorts[dimension - 1]++;
    std::vector<Solution* > keep_inf_second = marry(sup_first, 0, sup_first.size() - 1, dimension - 1);

    // Take union of inf_first and keep_int_second
    inf_first.insert(inf_first.end(), keep_inf_second.begin(), keep_inf_second.end());

    return inf_first;
}

std::vector<Solution *> DivideAndConquer::marry_2d(std::vector<Solution *> &curr_set) {
    // Solutions to be saved
    std::vector<Solution* > saved_solutions;

    // Used to track the max first criteria seen so far
    float max_first = std::numeric_limits<float>::lowest();
    for (Solution* sol: curr_set) {
        // compare max_first or the first criteria
        if (sol->criteria_2[0].second > max_first) {
            // Keep solution if is inferior
            if (sol->inferior) {
                saved_solutions.push_back(sol);
            } else {
                max_first = sol->criteria_2[0].second;
            }
        } else if (sol->inferior) { // Free the memory pointed to by sol because it is dominated
            delete(sol);
        }
    }

    return saved_solutions;
}

/**
 * Two dimensional sorting of Policies
 *
 * @param all_possible
 * @return
 */
std::vector<Solution *> DivideAndConquer::L2D(std::vector<Solution *> &all_possible) {
    // Sort by the 1st criteria (i.e. the criteria at index 1)
    std::sort(all_possible.begin(), all_possible.end(), Compare_Dimensions(2, num_criteria, use_linear_preferences));
    sorts[2]++;
    std::vector<Solution* > solutions;
    float second_criteria_max;

    // Put the first element into the solution set and update second_criteria_max;
    solutions.push_back(all_possible[0]);
    second_criteria_max = all_possible[0]->criteria_2[0].second; // Get the last criteria
    int deleted = 0;

    for (unsigned long curr_index = 1; curr_index < all_possible.size(); curr_index++) {
        // If the last criteria is larger than the max so far update and include solution
        if (all_possible[curr_index]->criteria_2[0].second > second_criteria_max) {
            {
                solutions.push_back(all_possible[curr_index]);
            }
            second_criteria_max = all_possible[curr_index]->criteria_2[0].second;
        } else { // Free the memory for the dominated solution
            delete(all_possible[curr_index]);

            deleted++;
        }
    }

    return solutions;
}