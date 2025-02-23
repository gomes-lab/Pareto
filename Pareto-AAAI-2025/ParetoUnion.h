//
// Created by Marc Grimson on 7/31/23
//

#ifndef AMAZON_PROJECT_DP_ALGORITHM_H
#define AMAZON_PROJECT_DP_ALGORITHM_H

#include <cstdlib>
#include <stack>
#include <cmath>
#include <thread>
#include <algorithm>
#include <mutex>
#include <iomanip>
#include <unordered_set>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include "Solution.h"
#include "Consts.h"
#include "DivideAndConquer.h"
#include <map>

class ParetoUnion {
public:
    ParetoUnion(std::vector<std::string> criteria, std::vector<std::string> paths, std::string name, unsigned long batch_size); // Constructor
    virtual ~ParetoUnion();

    void run();

private:
    std::fstream get_solution_file(const std::string path);
    std::fstream get_data_file(const std::string path);
    void process_path(std::string path, std::unordered_set<Solution*, SolutionHash, EqualSolutions> &solution_set, std::map<std::size_t, std::string> &sol_dams);
    bool add_solution(Solution* adding, std::unordered_set<Solution*, SolutionHash, EqualSolutions> &solution_set);

    std::vector<std::string> criteria;
    std::vector<std::string> criteria_all;
    std::vector<std::string> paths;
    std::string name;
    unsigned long batch_size;

    DivideAndConquer *dnc;
};

#endif //AMAZON_PROJECT_DP_ALGORITHM_H
