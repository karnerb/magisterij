#include "llmodel.h"
#include "logger.h"

Logger::Logger(std::string filename, LL_model& model) : model(model){
    file.open(filename);
    file << "n: " << model.n << "\n";
    file << "beta" << " " << "cycle" << " " << "acceptance_rate" << " " << "Energy" << " "
        << "P2_00" << " " << "P2_01" << " " << "P2_02" << " "
        << "P2_10" << " " << "P2_11" << " " << "P2_12" << " "
        << "P2_20" << " " << "P2_21" << " " << "P2_22"
        << "\n"; 

}

Logger::~Logger(){
    file.close();
}

void Logger::log_params(){
    file << model.beta << " " << model.cycle << " " << model.acceptance_rate << " " << model.E << " "
    << model.P2[0][0] << " " << model.P2[0][1] << " " << model.P2[0][2] << " "
    << model.P2[1][0] << " " << model.P2[1][1] << " " << model.P2[1][2] << " "
    << model.P2[2][0] << " " << model.P2[2][1] << " " << model.P2[2][2]
    << "\n"; 
}