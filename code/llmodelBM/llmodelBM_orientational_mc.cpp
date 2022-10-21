#include "llmodelBM.h"
#include <algorithm>

void LL_model_BM::shuffle_I(){
    //randomly shuffles elements of shuffled_I[n*n*n]
    /*    int temp, j;
    for (int i=n*n*n-1; i>0; i--){
        temp=shuffled_I[i];
        j = random_I(generator);
        shuffled_I[i]=shuffled_I[j];
        shuffled_I[j]=temp;
    }
    */
   std::random_shuffle(&shuffled_I[0], &shuffled_I[n*n*n-1]);
}



bool LL_model_BM::BarkerWatts_move(int I){    
    //save the previous config
    double temp[2][3];
    for (int i=0; i<3; i++){
        temp[0][i] = spins[I][0][i];
        temp[1][i] = spins[I][1][i];
    }
    
    //save the energy of current config
    double Eold = E_neighbors(I);

    //generate a random unit vector ((n0,n1,n2)) and angle of rotation (gama)
    double theta = acos(2*p(generator)-1);
    double phi = 2 * pi * p(generator);

    double n0 = sin(theta)*cos(phi);
    double n1 = sin(theta)*sin(phi);
    double n2 = cos(theta);
        
    double gama = rotation_angle*(p(generator)-0.5);
    
    //rotate the spin
    rotation3D(spins[I][0], n0, n1, n2, gama);
    rotation3D(spins[I][1], n0, n1, n2, gama);
        
    //calculate new energy 
    double Enew = E_neighbors(I);

    //Metropolis acceptance criterion
    //if Enew < Eold always accept 
    if (Enew > Eold){
        double acceptance_p = p(generator);
        //reject with probability according to Metropolis and keep the previous config
        if (acceptance_p > exp(-beta * (Enew-Eold))){
            spins[I][0][0] = temp[0][0];
            spins[I][0][1] = temp[0][1];
            spins[I][0][2] = temp[0][2];
            spins[I][1][0] = temp[1][0];
            spins[I][1][1] = temp[1][1];
            spins[I][1][2] = temp[1][2];
            return false;
            }
        } 
    return true;
    }

    void LL_model_BM::BarkerWatts_cycle(){
    //attempts a BarkerWatts move for each lattice site
    shuffle_I();
    int accepted = 0;
    for (int i=0; i<n*n*n; i++){
        if (BarkerWatts_move(shuffled_I[i])) accepted++;
    }
    acceptance_rate = (double) accepted / (n*n*n);
    adjust_rotation_angle();
    cycle++;
}

void LL_model_BM::thermalize_BarkerWatts(int N){
    //thermalizes the lattice using the BarkerWatts cycle
    for (int i = 0; i<N; i++){
        BarkerWatts_cycle();
    }
}

void LL_model_BM::adjust_rotation_angle(){
    //adjusts the rotation angle so that approximately half of proposed moves are accepted
    if (acceptance_rate>0.5){
        rotation_angle = rotation_angle + 0.1;
        if (rotation_angle > 2*pi) rotation_angle = 2*pi;
    }
    if (acceptance_rate<0.5){
        rotation_angle = rotation_angle - 0.1;
        if (rotation_angle<0) rotation_angle = 0.001;
    }
}
