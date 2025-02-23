#include <iostream>
#include <vector>
#include <algorithm>
#include "ParetoUnion.h"
#include <ctime>
#include <cmath>
#include <cassert>
#include "Consts.h"

using namespace std;

vector<string> optimized_criteria;
vector<string> paths;
std::string name;
int connectivity_index = -1;
int dcip_index = -1;
unsigned long batch_size = 0;

bool parse_arguments(int argc, char **argv) {
    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "-criteria") {
            int num_criteria = stoi(std::string(argv[++i]), 0); // Number of criteria to read in
            optimized_criteria.clear();
            std::string crit;
            for (int num = 0; num < num_criteria; num++) {
                crit = argv[++i];
                if (crit == "dcip") {
                    // Haven't added connectivity yet, add it here
                    if (connectivity_index == -1) {
                        connectivity_index = optimized_criteria.size();

                        optimized_criteria.emplace_back("connectivity");
                    }

                    dcip_index = optimized_criteria.size();

                    optimized_criteria.emplace_back(crit);
                }
                else if (crit == "connectivity") {
                    // If already added via dcip, ignore
                    if (connectivity_index == -1) {
                        connectivity_index = optimized_criteria.size();

                        optimized_criteria.emplace_back(crit);
                    }
                }
                else {
                    optimized_criteria.emplace_back(crit);
                }
            }
        }
        else if (std::string(argv[i]) == "-path") {
            int num_paths = stoi(std::string(argv[++i]), 0);

            paths.clear();

            std::string path;
            for (int num = 0; num < num_paths; num++) {
                path = argv[++i];

                paths.emplace_back(path);
            }
        }
        else if (std::string(argv[i]) == "-name") {
            name = std::string(argv[++i]);
        }else if (std::string(argv[i]) == "-batch_size") {
            batch_size = std::stol(argv[++i], nullptr, 0);
        }
    }
    return argc > 1;
}

int main (int argc, char **argv) {
    parse_arguments(argc, argv);

    ParetoUnion *paretoUnion = new ParetoUnion(optimized_criteria, paths, name, batch_size);
    paretoUnion->run();

    exit(0);
}
