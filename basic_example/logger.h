#pragma once
#include <fstream>
#include "llmodel.h"
class Logger{
    public:
    std::ofstream file;
    LL_model& model;
    Logger(std::string output_file, LL_model& model);
    ~Logger();
    void log_params();
};