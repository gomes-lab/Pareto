//
// Created by Marc Grimson on 6/28/23.
//

#include "Config.h"

Config::Config() {
    // Set default parameters
    epsilon = 0;
    seed = (unsigned int) time(nullptr);
    river_network = "Amazon";
    batch_size = 100000;
    file_path = "dam_data/amazon/amazon_built.txt";
    out_dir = "results";
    optimized_criteria = {"energy","connectivity"};
    use_transforms = true;
    save_progress = false;
    save_results = true;
    verbose = true;
    connectivity_index = -1;
    dcip_index = -1;
    sort_type = SORT_TYPE_TREE;
    compress_mode = COMPRESS_MODE_NONE;
    compress_gamma = 0.1;
    compress_distance_mode = COMPRESS_DISTANCE_OBJECTIVE_SPACE;
    compress_linkage = COMPRESS_LINKAGE_AVERAGE;
    compress_objectives = COMPRESS_OBJECTIVES_OPTIMIZED;
    compression_depth = 0;
    compression_threshold = 10000;
    compression_chunk_size = 1000;

    use_linear_preferences = false;
    w = {};

    num_threads=6;
}

void Config::save(std::fstream &config_file) {
    config_file << "epsilon = " << epsilon << std::endl;
    config_file << "path = " << file_path << std::endl;
    config_file << "outdir = " << out_dir << std::endl;
    config_file << "seed = " << seed << std::endl;
    config_file << "batch = " << batch_size << std::endl;
    config_file << "threads = " << num_threads << std::endl;
    config_file << "river = " << river_network << std::endl;
    config_file << "name = " << experiment_name << std::endl;
    config_file << "verbose = " << verbose << std::endl;
    config_file << "transforms = " << use_transforms << std::endl;
    config_file << "save_progress = " << save_progress << std::endl;
    config_file << "save_results = " << save_results << std::endl;
    config_file << "sort = " << sort_type << std::endl;
    config_file << "compress_mode = " << compress_mode << std::endl;
    config_file << "compress_gamma = " << compress_gamma << std::endl;
    config_file << "compress_distance_mode = " << compress_distance_mode << std::endl;
    config_file << "compress_linkage = " << compress_linkage << std::endl;
    config_file << "compress_threshold = " << compression_threshold << std::endl;
    config_file << "compress_chunk_size = " << compression_chunk_size << std::endl;
    config_file << "compress_depth = " << compression_depth << std::endl;
    config_file << "compress_objectives = " << compress_objectives << std::endl;

    config_file << "criteria =";
    for (auto crit : optimized_criteria) {
        config_file << " " << crit;
    }
    config_file << std::endl;

    if (!upper_bounds.empty()) {
        config_file << "ub =";
        for (auto ub : upper_bounds) {
            config_file << " " << ub;
        }
        config_file << std::endl;
    }

    if (!lower_bounds.empty()) {
        config_file << "lb =";
        for (auto lb : lower_bounds) {
            config_file << " " << lb;
        }
        config_file << std::endl;
    }

    config_file << "lp = " << use_linear_preferences << std::endl;
    if (!w.empty()) {
        config_file << "w =";
        for (auto cur_w : w) {
            for (auto i : cur_w) {
                config_file << " " << i;
            }
            config_file << ",";
        }
        config_file << std::endl;
    }
}

void Config::parse_ini(std::string ini_file) {
    inipp::Ini<char> ini;
    std::ifstream is(ini_file);

    ini.parse(is);
    ini.strip_trailing_comments();
    ini.interpolate();

    auto sec = ini.sections[""];

    // Handle the simple cases first
    inipp::get_value(sec, "epsilon", epsilon);
    inipp::get_value(sec, "path", file_path);
    inipp::get_value(sec, "outdir", out_dir);
    inipp::get_value(sec, "seed", seed);
    inipp::get_value(sec, "batch", batch_size);
    inipp::get_value(sec, "threads", num_threads);
    inipp::get_value(sec, "river", river_network);
    inipp::get_value(sec, "name", experiment_name);
    inipp::get_value(sec, "verbose", verbose);
    inipp::get_value(sec, "transforms", use_transforms);
    inipp::get_value(sec, "save_progress", save_progress);
    inipp::get_value(sec, "save_results", save_results);
    inipp::get_value(sec, "sort", sort_type);
    inipp::get_value(sec, "compress_mode", compress_mode);
    inipp::get_value(sec, "compress_gamma", compress_gamma);
    inipp::get_value(sec, "compress_distance_mode", compress_distance_mode);
    inipp::get_value(sec, "compress_threshold", compression_threshold);
    inipp::get_value(sec, "compress_chunk_size", compression_chunk_size);
    inipp::get_value(sec, "compress_depth", compression_depth);
    inipp::get_value(sec, "compress_linkage", compress_linkage);
    inipp::get_value(sec, "compress_objectives", compress_objectives);

    // Handle more complex structures
    std::string criteria;
    if (inipp::get_value(sec, "criteria", criteria)) {
        optimized_criteria.clear();
        std::string crit;
        std::istringstream iss(criteria);

        while (std::getline(iss, crit, ' ')) {
            optimized_criteria.emplace_back(crit);
        }
    }

    std::string ubs;
    if (inipp::get_value(sec, "ub", ubs)) {
        upper_bounds.clear();

        std::string ub;
        std::istringstream iss(ubs);

        while (std::getline(iss, ub, ' ')) {
            upper_bounds.emplace_back(std::stod(ub));
        }
    }

    std::string lbs;
    if (inipp::get_value(sec, "lb", lbs)) {
        lower_bounds.clear();

        std::string lb;
        std::istringstream iss(lbs);

        while (std::getline(iss, lb, ' ')) {
            lower_bounds.emplace_back(std::stod(lb));
        }
    }

    inipp::get_value(sec, "lp", use_linear_preferences);
    std::string ws;
    if (inipp::get_value(sec, "w", ws)) {
        w.clear();

        std::string s;
        std::istringstream iss(ws);

        while (std::getline(iss, s, ',')) {
            std::vector<double> cur_w;

            std::string s2;
            std::istringstream iss2(s);

            while (std::getline(iss2, s2, ' ')) {
                cur_w.emplace_back(std::stod(s2));
            }

            w.emplace_back(cur_w);
        }
    }

    is.close();
}

bool Config::parse_arguments(int argc, char **argv) {
    // First check if there's an ini file provided and apply further defaults
    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "-ini") {
            std::string ini_file = argv[++i];
            parse_ini(ini_file);
            break;
        }
    }

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
        else if (std::string(argv[i]) == "-ub") {
            int num_ub = stoi(std::string(argv[++i]), 0);
            upper_bounds.clear();

            for (int num = 0; num < num_ub; num++) {
                upper_bounds.emplace_back(stod(std::string(argv[++i]), nullptr));
            }
        }
        else if (std::string(argv[i]) == "-lb") {
            int num_lb = stoi(std::string(argv[++i]), 0);
            lower_bounds.clear();

            for (int num = 0; num < num_lb; num++) {
                lower_bounds.emplace_back(stod(std::string(argv[++i]), nullptr));
            }
        }
        else if (std::string(argv[i]) == "-w") {
            // The total number of criteria we considered, note the first criterion should always be the energy
            int num_c = stoi(std::string(argv[++i]), 0);
            w.clear();
            for (int num = 0; num < num_c; ++num) {
                int num_w = stoi(std::string(argv[++i]), 0);
                std::vector<double> cur_w;
                for (int numw = 0; numw < num_w; ++numw) {
                    cur_w.push_back(stod(std::string(argv[++i]), 0));
                }
                w.push_back(cur_w);
            }
            assert ((int)w.size() == num_c);
            //TODO: for sanity check, can remove after develop
            for (const auto& x: w) {
                for (auto tt: x) {
                    std::cout << tt << " ";
                }
                std::cout << std::endl;
            }
        } else if (std::string(argv[i]) == "-basin") {
            river_network = std::string(argv[++i]);
        } else if (std::string(argv[i]) == "-path") {
            file_path = std::string(argv[++i]);
        } else if (std::string(argv[i]) == "-outdir") {
            out_dir = std::string(argv[++i]);
        } else if (std::string(argv[i]) == "-seed") {
            seed = static_cast<unsigned int>(stoul(argv[++i], nullptr, 0));
        } else if (std::string(argv[i]) == "-epsilon") {
            epsilon = stod(argv[++i], nullptr);
        } else if (std::string(argv[i]) == "-batch") {
            batch_size = stoi(argv[++i], nullptr, 0);
        } else if (std::string(argv[i]) == "-thread") {
            num_threads = stoi(argv[++i], nullptr);
        } else if (std::string(argv[i]) == "-name") {
            experiment_name = std::string(argv[++i]);
            name_set = true;
        } else if (std::string(argv[i]) == "-verbose") {
            verbose = true;
        } else if (std::string(argv[i]) == "-quiet") {
            verbose = false;
        } else if (std::string(argv[i]) == "-transforms") {
            use_transforms = true;
        } else if (std::string(argv[i]) == "-ntransforms") {
            use_transforms = false;
        } else if (std::string(argv[i]) == "-save") {
            save_progress = true;
        } else if (std::string(argv[i]) == "-nsave") {
            save_progress = false;
        } else if (std::string(argv[i]) == "-saveresults") {
            save_results = true;
        } else if (std::string(argv[i]) == "-nsaveresults") {
            save_results = false;
        } else if (std::string(argv[i]) == "-lp") {
            use_linear_preferences = true;
        } else if (std::string(argv[i]) == "-compressmode") {
            int mode = stoi(argv[++i], nullptr);
            switch(mode) {
                case COMPRESS_MODE_NONE:
                case COMPRESS_MODE_FINAL:
                case COMPRESS_MODE_DEPTH_STATIC:
                case COMPRESS_MODE_THRESHOLD:
                case COMPRESS_MODE_DEPTH_DYNAMIC:
                    compress_mode = mode;
                    break;
                default:
                    std::cerr << "Unknown compression mode " << mode << std::endl;
                    exit(1);
            }
        } else if (std::string(argv[i]) == "-compressobjectives") {
            int mode = stoi(argv[++i], nullptr);
            switch(mode) {
                case COMPRESS_OBJECTIVES_ALL:
                case COMPRESS_OBJECTIVES_OPTIMIZED:
                    compress_objectives = mode;
                    break;
                default:
                    std::cerr << "Unknown objective compression mode " << mode << std::endl;
                    exit(1);
            }
        } else if (std::string(argv[i]) == "-compressgamma") {
            compress_gamma = stod(argv[++i], nullptr);
        } else if (std::string(argv[i]) == "-compressdepth") {
            compression_depth = stoi(argv[++i], nullptr);
        } else if (std::string(argv[i]) == "-compressthresh") {
            compression_threshold = stol(argv[++i], nullptr);
        } else if (std::string(argv[i]) == "-compresschunksize") {
            compression_chunk_size = stol(argv[++i], nullptr);
        } else if (std::string(argv[i]) == "-compressdistancemode") {
            int mode = stoi(argv[++i], nullptr);
            switch(mode) {
                case COMPRESS_DISTANCE_DECISION_SPACE:
                case COMPRESS_DISTANCE_OBJECTIVE_SPACE:
                    compress_distance_mode = mode;
                    break;
                default:
                    std::cerr << "Unknown compression distance mode " << mode << std::endl;
                    exit(1);
            }
        } else if (std::string(argv[i]) == "-compresslinkage") {
            int mode = stoi(argv[++i], nullptr);
            switch(mode) {
                case COMPRESS_LINKAGE_SINGLE:
                case COMPRESS_LINKAGE_AVERAGE:
                case COMPRESS_LINKAGE_COMPLETE:
                    compress_linkage = mode;
                    break;
                default:
                    std::cerr << "Unknown compression linkage mode " << mode << std::endl;
                    exit(1);
            }
        } else if (std::string(argv[i]) == "-sort") {
            int method = stoi(argv[++i], nullptr);
            switch(method) {
                case SORT_TYPE_NONE:
                case SORT_TYPE_TREE:
                case SORT_TYPE_TREE_REVERSE:
                case SORT_TYPE_FRONTIER:
                case SORT_TYPE_FRONTIER_REVERSE:
                case SORT_TYPE_RANDOM:
                    sort_type = method;
                    break;
                default:
                    std::cerr << "Unknown sorting method " << method << std::endl;
                    exit(1);
            }
        }
        else if (std::string(argv[i]) == "-river") {
            river_network = argv[++i];
        } else if (std::string(argv[i]) == "-ini") {
            // Already handled, skip
            i++;
        } else {
            print_help();
        }
    }
    return argc > 1;
}

void Config::print_help() {
    std::cout << "General Parameters" << std::endl;
    std::cout << "-ini STR = ini file with default parameters" << std::endl;
    std::cout << "-criteria N STR STR = reads in N relevant criteria" << std::endl;
    std::cout << "-basin STR = the name of the river basin that we are working on" << std::endl;
    std::cout << "-path FILE = the file with the data for the river basin" << std::endl;
    std::cout << "-outdir DIR = the base directory to write output results to. Defaults to the local directory" << std::endl;
    std::cout << "-seed N = the random seed to use in the experiment (NOTE: if this is not given a random seed is chosen)" << std::endl;
    std::cout << "-epsilon N = the epsilon that will be used for rounding in the experiment" << std::endl;
    std::cout << "-batch N = The batch size used for dynamically processing partial policies in nlogn algorithm" << std::endl;
    std::cout << "-name NAME = name of the experiment run, can be used to restart previous runs that may have failed partway through" << std::endl;
    std::cout << "-transforms = use transforms for removing Pareto frontier search spaces" << std::endl;
    std::cout << "-ntransforms = don't use transforms for removing Pareto frontier search spaces" << std::endl;
    std::cout << "-save = save progress on each node" << std::endl;
    std::cout << "-nsave = don't save progress on each node - just save final results" << std::endl;
    std::cout << "-saveresults = save the final frontier, this should be true unless you are just looking for run metrics" << std::endl;
    std::cout << "-nsaveresults = don't save the final frontier, this should not be set unless you are just looking for run metrics" << std::endl;
    std::cout << "-sort METHOD = sorting method to use:" << std::endl;
    std::cout << "\t" << SORT_TYPE_NONE << ": no sort" << std::endl;
    std::cout << "\t" << SORT_TYPE_TREE << ": sort by tree size" << std::endl;
    std::cout << "\t" << SORT_TYPE_TREE_REVERSE << ": sort by tree size reverse" << std::endl;
    std::cout << "\t" << SORT_TYPE_FRONTIER << ": sort by frontier size" << std::endl;
    std::cout << "\t" << SORT_TYPE_FRONTIER_REVERSE << ": sort by frontier size reverse" << std::endl;
    std::cout << "\t" << SORT_TYPE_RANDOM << ": sort randomly" << std::endl;
    std::cout << "-lp = use linear preferences" << std::endl;
    std::cout << "-compressmode METHOD = mode for compression" << std::endl;
    std::cout << "\t" << COMPRESS_MODE_NONE << ": no compression" << std::endl;
    std::cout << "\t" << COMPRESS_MODE_FINAL << ": compress final solution set only" << std::endl;
    std::cout << "\t" << COMPRESS_MODE_DEPTH_STATIC << ": compress nodes at depth n or lower before binarized (use with -compressdepth <n>)" << std::endl;
    std::cout << "\t" << COMPRESS_MODE_THRESHOLD << ": compress nodes with > n solutions (use with -compressthresh <n>)" << std::endl;
    std::cout << "\t" << COMPRESS_MODE_DEPTH_DYNAMIC << ": compress nodes at depth n or lower after binarized (use with -compressdepth <n>)" << std::endl;
    std::cout << "-compressgamma N = compress using a gamma of N" << std::endl;
    std::cout << "-compressdepth N = compress starting at a certain depth of the tree" << std::endl;
    std::cout << "-compressthresh N = compress nodes with more than N solutions" << std::endl;
    std::cout << "-compresschunksize N = split solutions into chunks of size N for compression in parallel" << std::endl;
    std::cout << "-compressdistancemode METHOD = mode for distance calculation in clustering" << std::endl;
    std::cout << "\t" << COMPRESS_DISTANCE_OBJECTIVE_SPACE << ": compress based on decision space distances" << std::endl;
    std::cout << "\t" << COMPRESS_DISTANCE_DECISION_SPACE << ": compress based on objective space distances" << std::endl;
    std::cout << "-compresslinkage METHOD = mode for linkage compress in clustering" << std::endl;
    std::cout << "\t" << COMPRESS_LINKAGE_SINGLE << ": compress based on single link" << std::endl;
    std::cout << "\t" << COMPRESS_LINKAGE_AVERAGE << ": compress based on average link" << std::endl;
    std::cout << "\t" << COMPRESS_LINKAGE_COMPLETE << ": compress based on complete link" << std::endl;
    exit(-1);
}

void Config::set_name() {
    if (!name_set) {
        // Create file name
        std::ostringstream name;
        // Append tree structure
        name << "tree_";
        // Add the river network
        name << river_network;
        for (const std::string &str: optimized_criteria) {
            name << "_" << str;
        }
        // Add the epsilon
        name << "_e_" << epsilon;

        name << "_K_" << num_threads;

        // Add the batching value
        name << "_b_" << batch_size;

        // Get the time
        time_t rawtime;
        struct tm *timeinfo;
        char buffer[80];

        time(&rawtime);
        timeinfo = localtime(&rawtime);

        strftime(buffer, 80, "%c", timeinfo);
        std::string time_stamp = buffer;
        // Get rid of the spaces
        std::replace(time_stamp.begin(), time_stamp.end(), ' ', '_');
        std::replace(time_stamp.begin(), time_stamp.end(), ':', '_');
        name << "_" << time_stamp;

        experiment_name = name.str();
    }
}