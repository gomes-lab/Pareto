#include <iostream>
#include <vector>
#include <algorithm>
#include "Network.h"
#include "ParetoDP.h"
#include <ctime>
#include <cmath>
#include <cassert>
#include "IOHandler.h"
#include "Config.h"

using namespace std;

Config config;

void writeMeta(IOHandler &io) {
    std::fstream meta_file;
    meta_file = io.get_file(true, META_FILE_NAME);

    meta_file << "n " << config.experiment_name << std::endl;
    meta_file << "e " << config.epsilon << std::endl;
    meta_file << "s " << config.seed << std::endl;
    meta_file << "f " << config.file_path << std::endl;
    meta_file << "r " << config.river_network << std::endl;
    meta_file << "b " << config.batch_size << std::endl;
    meta_file << "k " << config.num_threads << std::endl;
    meta_file << "c " << config.optimized_criteria.size();
    for (auto const &crit : config.optimized_criteria) {
        meta_file << " " << crit;
    }
    meta_file << std::endl;

    meta_file.close();
}

void writeRunInfo(IOHandler &io, RunInfo &run_info) {
    std::fstream run_file = io.get_file(true, RUN_FILE_NAME);

    run_file << run_info.num_policies_final << " " << run_info.num_policies_generated << " ";
    run_file << run_info.num_policies_pruned << " " << run_info.max_node_policies << " " << run_info.num_comparisons << std::endl;

    run_file << run_info.wall_time << " " << run_info.cpu_time << std::endl;

    run_file << run_info.transforms_considered << " " << run_info.transforms_pruned << " " << run_info.policies_avoided << std::endl;

    run_file.close();
}

void run() {
    IOHandler io(config);

    if (io.config.save_progress) {
        writeMeta(io);
    }

    std::fstream config_file = io.get_file(true, CONFIGURATION_FILE_NAME);
    config.save(config_file);
    config_file.close();

    std::cout << "Beginning run for " << config.experiment_name << std::endl;

    Network *net = new Network(config, io);

    // Sanity check parameters
    assert(config.optimized_criteria.size() <= net->num_criteria);

    std::cout << "Beggining DP" << std::endl;

    assert ((config.use_linear_preferences && (int)config.w.size() > 0) || (!config.use_linear_preferences && (int)config.w.size() == 0));
    ParetoDP dp(net, io, config, true);

    dp.run_brute_force();
    writeRunInfo(io, dp.run_info);
    delete(net);
}

int main (int argc, char **argv) {
    config.parse_arguments(argc, argv);
    config.set_name();

    // Set the seed!
    srand(config.seed);

    run();

    exit(0);
}
