//
// Created by Marc Grimson on 4/1/23.
//

#ifndef HYDRODAMDP_CONFIG_H
#define HYDRODAMDP_CONFIG_H

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <cassert>
#include <fstream>

#include "inipp.h"
#include "Consts.h"

using namespace std;

class Config
{
public:
    double epsilon;
    string file_path;
    string river_network;
    string out_dir;
    vector<string> optimized_criteria;
    vector<double> upper_bounds;
    vector<double> lower_bounds;
    int sort_type;
    unsigned int seed;
    int batch_size;
    size_t num_threads;
    string experiment_name;
    bool name_set;
    bool verbose;
    bool use_transforms;
    bool save_progress;
    bool save_results;
    bool use_linear_preferences;
    int compress_mode;
    double compress_gamma;
    int compress_distance_mode;
    int compress_linkage;
    int compression_depth;
    long compression_threshold;
    long compression_chunk_size;
    int compress_objectives;
    vector<vector<double> > w;

    int connectivity_index;
    int dcip_index;

    Config();

    bool parse_arguments(int argc, char **argv);
    void save(std::fstream &config_file);
    void set_name();
private:
    void print_help();
    void parse_ini(std::string ini_file);
};

#endif //HYDRODAMDP_CONFIG_H
