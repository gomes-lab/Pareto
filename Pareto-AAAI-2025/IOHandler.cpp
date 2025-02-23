//
// Created by marc on 1/18/21.
//

#include "IOHandler.h"

IOHandler::IOHandler(const Config &config)
        : config(config) {
    this->file_path = config.out_dir + "/" + config.experiment_name;

    check_directory();
}

void IOHandler::check_directory() {
    struct stat info;

    if (stat(this->file_path.c_str(), &info) == 0) {
        if (info.st_mode & S_IFDIR) {
            // At minimum we require the meta and data files and the node_data directory
            if (!(this->check_file(META_FILE_NAME) && this->check_file(DATA_FILE_NAME) && this->check_file(NODE_DIRECTORY_NAME))) {
                throw ERROR_MISSING_FILES;
            }
        }
        else {
            throw ERROR_NOT_DIRECTORY;
        }
    }
    else {
        mkdir(this->file_path.c_str(), 0755);
        mkdir((this->file_path + "/" + NODE_DIRECTORY_NAME).c_str(), 0755);

        // Copy the data file over
        std::ifstream src(this->config.file_path, std::ios::binary);
        std::ofstream dst(this->file_path + "/" + DATA_FILE_NAME, std::ios::binary);
        dst << src.rdbuf();
    }
}

std::fstream IOHandler::get_file(bool create, const std::string file_name, int node_id) {
    std::fstream file;
    std::fstream::openmode mode = std::fstream::out | std::fstream::in;
    if (create) {
        mode = mode | std::fstream::trunc;
        file.unsetf(std::ios::floatfield);
        file.precision(16);
    }
    file.open(this->get_node_path(node_id) + "/" + file_name, mode);

    if (!file.is_open()) {
        throw FILE_UNOPENED;
    }
    return file;
}

std::fstream IOHandler::get_file(bool create, const std::string file_name) {
    std::fstream file;
    std::fstream::openmode mode = std::fstream::out | std::fstream::in;
    if (create) {
        mode = mode | std::fstream::trunc;
        file.unsetf(std::ios::floatfield);
        file.precision(16);
    }
    file.open(this->file_path + "/" + file_name, mode);

    if (!file.is_open()) {
        throw FILE_UNOPENED;
    }
    return file;
}

bool IOHandler::check_file(const std::string file_name) {
    std::ifstream ifile;
    ifile.open(this->file_path + "/" + file_name);
    if (ifile) {
        ifile.close();
        return true;
    }
    else {
        return false;
    }
}

bool IOHandler::check_file(const std::string file_name, int node_id) {
    std::ifstream ifile;
    ifile.open(this->get_node_path(node_id) + "/" + file_name);
    if (ifile) {
        ifile.close();
        return true;
    }
    else {
        return false;
    }
}

void IOHandler::make_node_dir(int node_id) {
    mkdir(this->get_node_path(node_id).c_str(), 0755);
}

std::string IOHandler::get_node_path(int node_id) {
    std::string node_path = this->file_path + "/" + NODE_DIRECTORY_NAME + "/" + std::to_string(node_id);

    return node_path;
}


void IOHandler::get_line(std::fstream &file, std::stringstream &reader, std::string &line) {
    getline(file, line);
    reader.str(line);
    reader.clear();
}

IOHandler::IOHandler(const IOHandler &handler) : config(handler.config), file_path(handler.file_path) {
}
