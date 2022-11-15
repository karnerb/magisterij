#include <iostream>
#include "LC_Elastomer.h"

std::string LC_Elastomer::metadata(){
    std::string space = " ";
    std::string metadata = std::to_string(n) + space + std::to_string(count_broken_molecules()) + space + std::to_string(alpha) + "\n";
    return metadata;
}

std::string LC_Elastomer::names(){
        std::string space = " ";
        std::string names = 
        "beta" + space 
        + "cycle" + space 
        + "acceptance_rate" + space 
        + "EnergyLL" + space
        + "EnergyElastic" + space
        + "EnergyCoupling" + space
        + "P2_00" + space + "P2_01" + space + "P2_02" + space
        + "P2_10" + space + "P2_11" + space + "P2_12" + space
        + "P2_20" + space + "P2_21" + space + "P2_22" + space
        + "same_neighbors" + space
        + "broken_cluster_count" + space
        + "lambda" + space
        + "sigma" 
        + "\n";
        return names;
}

std::ostream& operator <<(std::ostream& os, const LC_Elastomer& model){
    os  << model.beta << " " 
        << model.cycle << " " 
        << model.acceptance_rate << " "
        << model.Energy_LL_interaction << " "
        << model.Energy_elastic << " "
        << model.Energy_coupling << " "
        << model.P2[0][0] << " " << model.P2[0][1] << " " << model.P2[0][2] << " "
        << model.P2[1][0] << " " << model.P2[1][1] << " " << model.P2[1][2] << " "
        << model.P2[2][0] << " " << model.P2[2][1] << " " << model.P2[2][2] << " "
        << model.broken_neighbors << " "
        << model.broken_cluster_count << " "
        << model.lambda << " " 
        << model.sigma  
        << "\n";
        return os;
};
