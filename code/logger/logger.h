#include <fstream>
#include <iostream>

template <typename T> class Logger{
    public:
    T& model;
    std::ofstream file;

    
    Logger(std::string filename, T& model) : model(model){

    file.open("output/" + filename);
    if (!file.is_open()) std::cout << "Error opening output file!\n";
    file << model.metadata();
    file << model.names();
    }

    ~Logger(){
        file.close();
    }

    void log_params(){
        file << model;
       /* file << model.beta << " " << model.cycle << " " << model.acceptance_rate << " " << model.E << " "
        << model.P2[0][0] << " " << model.P2[0][1] << " " << model.P2[0][2] << " "
        << model.P2[1][0] << " " << model.P2[1][1] << " " << model.P2[1][2] << " "
        << model.P2[2][0] << " " << model.P2[2][1] << " " << model.P2[2][2] << " "
        << model.polar_order
        << "\n"; */
    }
};
