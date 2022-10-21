#include "llmodelBM.h"
#include <iostream>

void LL_model_BM::break_molecules(int count){
    //breaks int count molecules
    shuffle_I();
    for (int i=0; i<count; i++){
        break_molecule(shuffled_I[i]);
    }
}

void LL_model_BM::initialize_lattice_random(){
    //sets the spins on the lattice in random directions
    for (int I=0; I<n*n*n; I++){
        double theta = acos(2*p(generator)-1);
        double phi = 2 * pi * p(generator);
        spins[I][0][0]=sin(theta)*cos(phi);
        spins[I][0][1]=sin(theta)*sin(phi);
        spins[I][0][2]=cos(theta);
        spins[I][1][0]=sin(theta)*cos(phi);
        spins[I][1][1]=sin(theta)*sin(phi);
        spins[I][1][2]=cos(theta);
    }
}

void LL_model_BM::initialize_lattice_parallel(){
    //sets the spins in the same random direction
    double theta = acos(2*p(generator)-1);
    double phi = 2 * pi * p(generator);
    double n1 = sin(theta)*cos(phi);
    double n2 = sin(theta)*sin(phi);
    double n3 = cos(theta);    
    for (int I=0; I<n*n*n; I++){
        spins[I][0][0]=n1;
        spins[I][0][1]=n2;
        spins[I][0][2]=n3;
        spins[I][1][0]=n1;
        spins[I][1][1]=n2;
        spins[I][1][2]=n3;

    }        
}




void LL_model_BM::break_molecule(int I){

    if (broken[I]) std::cout << "molecule already broken!";
    //choose a random perpendicular axis of rotation and rotate the two vectors by angle +-alpha/2

    //choose a random axis
    //first find a vector perpendicular to spins[I][0] using the cross product with x system axis
    double v[3];
    v[0] = 0;
    v[1] = spins[I][0][2] / sqrt(spins[I][0][2]*spins[I][0][2]+spins[I][0][1]*spins[I][0][1]);
    v[2] = -spins[I][0][1] / sqrt(spins[I][0][2]*spins[I][0][2]+spins[I][0][1]*spins[I][0][1]);

    //rotate the perpendicular vector by angle phi € [0,2pi) about spins[I][0] to get a random axis of rotation

    double n0 = spins[I][0][0];
    double n1 = spins[I][0][1];
    double n2 = spins[I][0][2];    
    double gama = 2*p(generator)*pi;

    rotation3D(v, n0, n1, n2, gama);

    
    //rotate the spins by +-alpha/2 about (v0, v1, v2)
    double alpha = pi / 3;
    n0 = v[0];
    n1 = v[1];
    n2 = v[2];        
    gama = alpha/2;

    rotation3D(spins[I][0], n0, n1, n2, gama);
    rotation3D(spins[I][1], n0, n1, n2, -gama);

    broken[I]=true;
   //std::cout << spins[I][1][0]*spins[I][0][0] + spins[I][1][1]*spins[I][0][1] + spins[I][1][2]*spins[I][0][2]; 

}

void LL_model_BM::align_molecule(int I){
    //aligns a broken molecule
    //sum the two vectors and normalize
    double temp[3];
    temp[0] = spins[I][0][0] + spins[I][1][0];
    temp[1] = spins[I][0][1] + spins[I][1][1];
    temp[2] = spins[I][0][2] + spins[I][1][2];

    double norm = sqrt(temp[0]*temp[0] + temp[1]*temp[1]+temp[2]*temp[2]);

    spins[I][0][0] = temp[0]/norm;
    spins[I][0][1] = temp[1]/norm;
    spins[I][0][2] = temp[2]/norm;

    spins[I][1][0] = temp[0]/norm;
    spins[I][1][1] = temp[1]/norm;
    spins[I][1][2] = temp[2]/norm;

    broken[I] = false;
}
