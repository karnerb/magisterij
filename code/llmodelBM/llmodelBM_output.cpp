#include <iostream>
#include "llmodelBM.h"

std::string LL_model_BM::metadata(){
    std::string space = " ";
    std::string metadata = std::to_string(n) + space + std::to_string(count_broken_molecules()) + "\n";
    return metadata;
}

std::string LL_model_BM::names(){
        std::string space = " ";
        std::string names = "beta" + space + "cycle" + space + "acceptance_rate" + space + "Energy" + space
        + "P2_00" + space + "P2_01" + space + "P2_02" + space
        + "P2_10" + space + "P2_11" + space + "P2_12" + space
        + "P2_20" + space + "P2_21" + space + "P2_22" + space
        + "polar_order" + space + "same_neighbors" + space
        + "swap_acceptance_rate" + space
        + "broken_cluster_count" 
        + "\n";
        return names;
}

std::ostream& operator <<(std::ostream& os, const LL_model_BM& model){
    os << model.beta << " " << model.cycle << " " << model.acceptance_rate << " " << model.E << " "
        << model.P2[0][0] << " " << model.P2[0][1] << " " << model.P2[0][2] << " "
        << model.P2[1][0] << " " << model.P2[1][1] << " " << model.P2[1][2] << " "
        << model.P2[2][0] << " " << model.P2[2][1] << " " << model.P2[2][2] << " "
        << model.polar_order << " " << model.broken_neighbors << " "
        << model.swap_acceptance_rate << " " 
        << model.broken_cluster_count
        << "\n";
        return os;
};
