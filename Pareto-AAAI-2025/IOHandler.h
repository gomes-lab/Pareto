//
// Created by marc on 1/18/21.
//

#ifndef AMAZON_PROJECT_IOHANDLER_H
#define AMAZON_PROJECT_IOHANDLER_H

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
//#include <filesystem>
#include <sys/stat.h>
#include <sys/types.h>
#include "Consts.h"
#include "Config.h"

inline static const std::string TREE_FILE_NAME = "tree.txt";
inline static const std::string META_FILE_NAME = "meta.txt";
inline static const std::string DATA_FILE_NAME = "data.txt";
inline static const std::string NODE_DIRECTORY_NAME = "node_data";
inline static const std::string NODE_FILE_NAME = "node.txt";
inline static const std::string DAM_FILE_NAME = "dam.txt";
inline static const std::string POLICY_FILE_NAME = "policies.txt";
inline static const std::string SOLUTION_FILE_NAME = "solution.txt";
inline static const std::string RUN_FILE_NAME = "run.txt";
inline static const std::string CONFIGURATION_FILE_NAME = "config.ini";

class IOHandler {
public:
    IOHandler(const Config &config);

    IOHandler(const IOHandler &handler);

    void check_directory();
    void make_node_dir(int node_id);
    void get_line(std::fstream &file, std::stringstream &reader, std::string &line);
    std::fstream get_file(bool create, std::string file_name);
    std::fstream get_file(bool create, std::string, int node_id);
    bool check_file(std::string file_name);
    bool check_file(std::string file_name, int node_id);

    Config config;

private:
//    std::filesystem::path get_node_path(int node_id);
    std::string get_node_path(int node_id);

//    std::filesystem::path file_path;
    std::string file_path;
};


#endif //AMAZON_PROJECT_IOHANDLER_H
