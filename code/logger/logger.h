#include <fstream>
#include <iostream>

template <typename T> class Logger{
    public:
    T& model;
    std::ofstream file;
    int i_lattice=0;
    std::string filename;

    
    Logger(std::string name, T& model) : model(model), filename(name){
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

    void output_lattice_configuration(){
        std::ofstream lattice_configuration;
        std::string name = filename;
        lattice_configuration.open("output/lattice_config/" + name.erase(name.length()-4) + "_lattice_config_" + std::to_string(i_lattice) + ".txt");
        if (!lattice_configuration.is_open()) std::cout << "Error opening output file!\n";
        lattice_configuration << 1.0/model.beta << " " << model.n << " " 
                     << model.count_broken_molecules() << " " 
                     << model.broken_neighbors << " " 
                     << model.broken_cluster_count <<
                      "\n";
        lattice_configuration << "broken v1x v1y v1z v2x v2y v2z \n";
        for (int i=0; i<model.n*model.n*model.n; i++){
            lattice_configuration << model.broken[i] << " "
                                  << model.spins[i][0][0] << " " 
                                  << model.spins[i][0][1] << " "
                                  << model.spins[i][0][2] << " " 
                                  << model.spins[i][1][0] << " " 
                                  << model.spins[i][1][1] << " "
                                  << model.spins[i][1][2] << " " 
                                  << "\n";
        }

        lattice_configuration.close();
        i_lattice++;
    }


};
