//LL model with additional bent molecules
#include <iostream>
#include <time.h>
#include <cmath>
#include "llmodelBM.h"


LL_model_BM::LL_model_BM(int n) : n(n){
    spins = new double** [n*n*n];
    P2 = new double* [3];
    broken = new bool [n*n*n];
    neighbors_list = new int*[n*n*n];
    broken_neighbors = 0;
    broken_cluster_count = 0;
    swap_acceptance_rate = 0;
    shuffled_I = new int[n*n*n];
    for (int i=0; i<3; i++) P2[i] = new double [3];
    for (int I=0; I<n*n*n; I++){
        
        shuffled_I[I]=I;
        broken[I] = false;
        spins[I] = new double* [2];
        

        //setup neighbors list
        neighbors_list[I]=new int [6];
        int k = I/(n*n);
        int j = (I%(n*n))/n;
        int i = (I%(n*n))%n;

        int up, down, north, south, west, east;

        up = n*n*((k+1)%n) + n*(j) + i;
        if (k>0) down = n*n*(k-1) + n*(j) + i;
        else down = n*n*(n-1) + n*(j) + i;
        north = n*n*k + n*((j+1)%n) + i;
        if (j>0) south = n*n*k + n*(j-1) + i;
        else south = n*n*k + n*(n-1) + i;
        east = n*n*k + n*j + (i + 1)%n;
        if (i>0) west = n*n*k + n*j + i - 1;
        else west = n*n*k + n*j + n - 1;
       
        neighbors_list[I][0]=up;
        neighbors_list[I][1]=down;
        neighbors_list[I][2]=north;
        neighbors_list[I][3]=south;
        neighbors_list[I][4]=east;
        neighbors_list[I][5]=west;
    
    }
    for (int I=0; I<n*n*n; I++){
        spins[I][0] = new double[3];
        spins[I][1] = new double[3];

    }
    
    generator.seed(time(NULL));
    random_I = std::uniform_int_distribution <int>(0,n*n*n-1);
    random_neighbor = std::uniform_int_distribution<int>(0, 5);
    p = std::uniform_real_distribution<double>(0,1);
}


int LL_model_BM::count_broken_molecules(){
    int count = 0;
    for (int i=0; i<n*n*n; i++) if (broken[i]) count++;
    return count;
}
void LL_model_BM::set_beta(double beta){
    this -> beta = beta;
}

double LL_model_BM::dot(int I, int i,  int J, int j){
    //dot product of spins[I][i] and spins[J][j]
    double sum = 0.0;
    for (int k=0; k<3; k++){
        sum+=spins[I][i][k]*spins[J][j][k];
    }
    return sum;
}

double LL_model_BM::E_IJ(int I, int J){
    //interaction energy of sites I, J
    double dot1 = dot(I,0,J,0);
    double dot2 = dot(I,0,J,1);
    double dot3 = dot(I,1,J,0);
    double dot4 = dot(I,1,J,1);
    return - ( 1.5 * 0.25 * (dot1*dot1+dot2*dot2+dot3*dot3+dot4*dot4) - 0.5);
}


double LL_model_BM::E_neighbors(int I){
    //sums the interaction energies with nearest neighbors
    double sum=0.0;
    //neighbors(I);
    for (int i=0; i<6; i++){
        sum+=E_IJ(I, neighbors_list[I][i]);
    }
    return sum;
}

/*void LL_model_BM::H(){
    //sums E_neighbors over all lattice sites
    double sum=0;
    for (int I=0;  I<n*n*n; I++){
        sum+=E_neighbors(I);
    }
    E = 0.5*sum/(n*n*n);
}*/

void LL_model_BM::H(){
    double sum=0;
    
    for (int k=0; k<n; k++){
        for (int j=0; j<n; j++){
            for (int i=0; i<n; i=i+2){
                int I =  i+(j+k)%2 + n*j + n*n*k; 
                sum+=E_neighbors(I);
            }
        }
    }
    E=sum/(n*n*n);
}

void LL_model_BM::calculate_P2(){
        //calculates the matrix elements of P2 matrix order parameter
       for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                P2[i][j] = 0.0;            
            }
        }

        for (int I = 0; I<n*n*n; I++){
            for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                P2[i][j] = P2[i][j] + 0.25*((spins[I][0][i]+spins[I][1][i])*(spins[I][0][j]+spins[I][1][j]));
                }
            }
        }
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                P2[i][j] = 1.5*P2[i][j] / (n*n*n);            
            }
        }
        P2[0][0] = P2[0][0] - 0.5;
        P2[1][1] = P2[1][1] - 0.5;
        P2[2][2] = P2[2][2] - 0.5;
}






