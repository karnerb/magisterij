#include "LC_Elastomer.h"

bool LC_Elastomer::swap_move(){
    //select a random site and attempt to switch with a random neighbor site
    int I = random_I(generator);
    
    //choose a random neighbor
    
    int neighbor = random_neighbor(generator);
    //if (neighbor > 6 || neighbor < 0) std::cout << "BIG PROBLEM";
    int J = neighbors_list[I][neighbor];

    //store the current energy of the pair
    double Eold = E_neighbors(I) + E_neighbors(J);

    //switch the two sites
    swap_sites(I, J);

    //calculate the new energy of the pair
    double Enew = E_neighbors(I) + E_neighbors(J);

    double deltaE = (Enew-Eold);
    //accept if deltaE < 0 otherwise Metropolis criterion

    if (deltaE > 0){
        if (p(generator) > exp(-beta*deltaE)){
            //revert to previous config
            swap_sites(I, J);
            return false;
        }
    }
    return true;
}


void LC_Elastomer::swap_sites(int I, int J){
    
    double temp[2][3];
    bool temp2;

    temp2 = broken[I];
    broken[I]=broken[J];
    broken[J]=temp2;


    temp[0][0] = spins[I][0][0];
    temp[0][1] = spins[I][0][1];
    temp[0][2] = spins[I][0][2];

    temp[1][0] = spins[I][1][0];
    temp[1][1] = spins[I][1][1];
    temp[1][2] = spins[I][1][2];

    spins[I][0][0] = spins[J][0][0];
    spins[I][0][1] = spins[J][0][1];
    spins[I][0][2] = spins[J][0][2];

    spins[I][1][0] = spins[J][1][0];
    spins[I][1][1] = spins[J][1][1];
    spins[I][1][2] = spins[J][1][2];

    spins[J][0][0] = temp[0][0];
    spins[J][0][1] = temp[0][1];
    spins[J][0][2] = temp[0][2];

    spins[J][1][0] = temp[1][0];
    spins[J][1][1] = temp[1][1];
    spins[J][1][2] = temp[1][2];

}

void LC_Elastomer::swap_cycle(){
    int accepted = 0;
    for (int i=0; i<n*n*n; i++){
        if (swap_move()) accepted++;
    }
    swap_acceptance_rate = (double) accepted / (n*n*n);
    //std::cout << swap_acceptance_rate << "\r" << std::flush;
}
