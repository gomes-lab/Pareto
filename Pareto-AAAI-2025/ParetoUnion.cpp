//
// Created by Marc Grimson on 7/31/23
//

#include "ParetoUnion.h"


ParetoUnion::ParetoUnion(std::vector<std::string> criteria, std::vector<std::string> paths, std::string name,
                         unsigned long batch_size)
    : criteria(criteria), paths(paths), name(name), batch_size(batch_size) {
    dnc = new DivideAndConquer(criteria.size(), false);
}

ParetoUnion::~ParetoUnion() {
    delete(dnc);
}

bool ParetoUnion::add_solution(Solution* adding, std::unordered_set<Solution*, SolutionHash, EqualSolutions> &solution_set) {
    bool added_new;
    auto find = solution_set.find(adding);
    if (find != solution_set.end()) { // Solution tie exists!
        int rand_num = rand() % 2; // Only want to do if epsilon != 0 i.e. we are rounding!!!!

        if (rand_num == 0) { // Remove the current and replace it with the new node
            Solution* old = *find;

            solution_set.erase(find);
            delete (old);
            solution_set.insert(adding);
        } else {
            delete(adding);
        }

        added_new = false;
    } else {
        solution_set.insert(adding);

        added_new = true;
    }

    return added_new;
}

std::fstream ParetoUnion::get_solution_file(const std::string path) {
    std::fstream file;
    std::fstream::openmode mode = std::fstream::out | std::fstream::in;
    file.open(path + "/solution.txt", mode);

    if (!file.is_open()) {
        throw FILE_UNOPENED;
    }
    return file;
}

std::fstream ParetoUnion::get_data_file(const std::string path) {
    std::fstream file;
    std::fstream::openmode mode = std::fstream::out | std::fstream::in;
    file.open(path + "/data.txt", mode);

    if (!file.is_open()) {
        throw FILE_UNOPENED;
    }
    return file;
}

void ParetoUnion::process_path(std::string path, std::unordered_set<Solution*, SolutionHash, EqualSolutions> &solution_set,
                               std::map<std::size_t, std::string> &sol_dams) {
    std::fstream file = get_solution_file(path);
    std::fstream dfile = get_data_file(path);
    std::string line;
    int i = 0;
    int num_solns = 0;
    std::vector<int> indices;
    int total_criteria = 0;
    std::array<double, MAX_CRITERIA> minmax;

    while (std::getline(dfile, line)) {
        if (line.size() >= 8 && line.substr(0, 9) == "criteria ") {
            std::string delimiter = " ";
            std::string crit_line = line.substr(9);

            size_t pos = 0;
            std::string token;
            int j = 0;
            while ((pos = crit_line.find(delimiter)) != std::string::npos) {
                token = crit_line.substr(0, pos);

                if (token.substr(0, 1) == "-") {
                    minmax[j] = -1.0;
                }
                else {
                    minmax[j] = 1.0;
                }

                crit_line.erase(0, pos + delimiter.length());
                j++;
            }
            if (crit_line.substr(0, 1) == "-") {
                minmax[j] = -1.0;
            }
            else {
                minmax[j] = 1.0;
            }
        }
    }

    while (std::getline(file, line)) {
        int num_criteria = 0;
        if (i < 13) {
            // Do nothing
        }
        // Criteria all line
        else if (i == 13) {
            // Skip "criteria all: "
            std::string crit_all = line.substr(14);
            criteria_all.clear();

            std::string delimiter = ", ";

            size_t pos = 0;
            std::string token;
            int j = 0;
            while ((pos = crit_all.find(delimiter)) != std::string::npos) {
                token = crit_all.substr(0, pos);
                for (auto c : criteria) {
                    if (token == c) {
                        indices.emplace_back(j);
                    }
                }
                crit_all.erase(0, pos + delimiter.length());
                criteria_all.emplace_back(token);
                j++;
                total_criteria++;
            }
            for (auto c : criteria) {
                if (crit_all == c) {
                    indices.emplace_back(j);
                }
            }
            criteria_all.emplace_back(crit_all);
            total_criteria++;
        }
        else {
            std::string solution = line;
            std::string delimiter = ", ";

            double criteria_all[MAX_CRITERIA];
            std::pair<double, double> crit[MAX_CRITERIA];

            size_t pos = 0;
            std::string token;
            int j = 0;
            while ((pos = solution.find(delimiter)) != std::string::npos) {
                token = solution.substr(0, pos);
                if (j < total_criteria) {
                    double val = std::stod(token);

                    for (auto k : indices) {
                        if (k == j) {
                            crit[num_criteria] = std::make_pair(val, minmax[k] * val);
                            num_criteria++;
                            break;
                        }
                    }

                    criteria_all[j] = std::stod(token);
                }

                solution.erase(0, pos + delimiter.length());
                j++;
            }

            num_solns++;

            std::vector<double> decision = {-1, -1};
            std::vector<Solution*> pareto_decisions;
            std::vector<std::vector<double> > w;
            std::array<double, MAX_CRITERIA> nf;
            Solution* soln = new Solution(-1, -1, criteria_all, crit, decision, pareto_decisions,
                                          total_criteria, num_criteria, false, w, minmax, nf);

            if (add_solution(soln, solution_set)) {
                sol_dams[soln->hash_val] = solution;
            }
        }

        i++;
    }

    std::cout << "Number of solutions read " << num_solns << std::endl;

    file.close();
}

void ParetoUnion::run() {
    // Create large list of solutions
    std::unordered_set<Solution*, SolutionHash, EqualSolutions> solution_set;
    std::map<std::size_t, std::string> sol_dams;

    unsigned long generated = 0;
    unsigned long last_size = 0;
    std::cout << "Beginning union of solutions" << std::endl;
    for (std::string path : paths) {
        std::cout << "Parsing solutions from " << path << std::endl;

        // Loop over each solution in the file and add to the set
        process_path(path, solution_set, sol_dams);

        generated += solution_set.size() - last_size;

        if (generated >= batch_size) {
            std::cout << "Checking Pareto dominance" <<  std::endl;
            // Run the DNC
            unsigned long num_comparisons = 0;
            std::vector<Solution* > batched_solutions(solution_set.begin(), solution_set.end());
            std::vector<Solution* > temp_answers = dnc->divide_and_conquer(batched_solutions, criteria.size(), false, num_comparisons);
            batched_solutions.clear();
            solution_set = std::unordered_set<Solution*, SolutionHash, EqualSolutions> (temp_answers.begin(), temp_answers.end());
            temp_answers.clear();

            std::cout << "New set has " << solution_set.size() << " solutions" << std::endl;

            generated = 0;
        }

        last_size = solution_set.size();
    }

    std::cout << "Checking Pareto dominance" <<  std::endl;
    // Run the DNC
    unsigned long num_comparisons = 0;
    std::vector<Solution* > batched_solutions(solution_set.begin(), solution_set.end());
    std::vector<Solution* > temp_answers = dnc->divide_and_conquer(batched_solutions, criteria.size(), false, num_comparisons);
    batched_solutions.clear();
    solution_set = std::unordered_set<Solution*, SolutionHash, EqualSolutions> (temp_answers.begin(), temp_answers.end());
    temp_answers.clear();

    std::cout << "There are " << solution_set.size() << " solutions" << std::endl;

    std::fstream file;
    std::fstream::openmode mode = std::fstream::out | std::fstream::in;
    mode = mode | std::fstream::trunc;
    file.unsetf(std::ios::floatfield);
    file.precision(16);
    file.open(name, mode);

    if (!file.is_open()) {
        throw FILE_UNOPENED;
    }

    std::cout << "Writing out file" << std::endl;
    file << "Union file" << std::endl;
    file << "Number of solutions: " << solution_set.size() << std::endl;
    file << "Criteria optimized: ";
    bool first = true;
    for (auto c : criteria) {
        if (first) {
            first = false;
        }
        else {
            file << ", ";
        }
        file << c;
    }
    file << std::endl;
    file << "Criteria all: ";
    first = true;
    for (auto c : criteria_all) {
        if (first) {
            first = false;
        }
        else {
            file << ", ";
        }
        file << c;
    }
    file << std::endl;
    for (auto soln : solution_set) {
        for (int i = 0; i < soln->total_criteria; i++) {
            file << soln->criteria[i] << ", ";
        }

        // Count number of dams
        int count = 1;
        std::string dam_list = sol_dams[soln->hash_val];
        for(int j = 0; j < dam_list.size(); j++)
        {
            if (dam_list[j] == ' ') {
                count++;
            }
        }

        file << count << ", " << dam_list << std::endl;
    }

    file.close();
}
