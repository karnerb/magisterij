#include "llmodelBM.h"
#include <iostream>
#include <cmath>

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


void LL_model_BM::make_one_cluster(int a){
    //creates a cube of volume a^3 of broken sites
    for (int I=0; I<n*n*n; I++) broken[I]=false;
    for (int k=0; k<a; k++){
        for (int j=0; j<a; j++){
            for (int i=0; i<a; i++){
                int I = i+j*n+k*n*n;
                broken[I]=true;
            }
        }
    }
}

void LL_model_BM::make_two_clusters(int a){
    //creates two cubes of volume a^3 of broken sites. a should be less or equal to size/2.
    for (int I=0; I<n*n*n; I++) broken[I]=false;
    for (int k=0; k<a; k++){
        for (int j=0; j<a; j++){
            for (int i=0; i<a; i++){
                int I = i+j*n+k*n*n;
                broken[I]=true;
            }
        }
    }
    
    for (int k=0; k<a; k++){
        for (int j=0; j<a; j++){
            for (int i=0; i<a; i++){
                int I = (n-1-i)+(n-1-j)*n+(n-1-k)*n*n;
                broken[I]=true;
            }
        }
    }
}

void LL_model_BM::checkerboard_pattern(){
    //creates a checkerboard pattern of broken sites
    for (int I=0; I<n*n*n; I++) broken[I] = false;
    for (int k=0; k<n; k++){
        for (int j=0; j<n; j++){
            for (int i=0; i<n; i=i+2){
                int I =  i+(j+k)%2 + n*j + n*n*k;
                broken[I]=true;
            }
        }
    }
}

void LL_model_BM::make_empty_cube_cluster(int a){
    //creates an "empty cube" with side a of broken sites
    for (int I=0; I<n*n*n; I++) broken[I] = false;

    for (int k=0; k<a; k++){
        for (int j=0; j<a; j++){
            for (int i=0; i<a; i++){
                int I = i+j*n+k*n*n;
                broken[I]=true;
            }
        }
    }

    for (int k=1; k<a-1; k++){
        for (int j=1; j<a-1; j++){
            for (int i=1; i<a-1; i++){
                int I = i+j*n+k*n*n;
                broken[I]=false;
            }
        }
    }
}

void LL_model_BM::make_L_shaped_cluster(int a, int b){
    //creates an L shaped cluster of broken sites with arm lengths a, b 
    for (int I=0; I<n*n*n; I++) broken[I] = false;
    for (int i=0; i<a; i++) broken[i] = true;
    for (int i=0; i<b; i++) broken[n*n*i] = true;
}