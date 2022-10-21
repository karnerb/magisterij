#include "llmodelBM.h"

void LL_model_BM::count_neighbors_of_same_kind(){
    int count=0;
    for (int i=0; i<n*n*n; i++){
        if (broken[i]){
            for (int j=0; j<6; j++){
                if (broken[neighbors_list[i][j]]){
                    count++;
            }
            }
        }
    }
    broken_neighbors=count/2;
}

void LL_model_BM::calculate_polar_order(){
    //calculates the polar order parameter of the broken molecules P = sum_i ()
    int broken_count=0;
    double vector[3];
    vector[0]=0;
    vector[1]=0;
    vector[2]=0;
    double temp[3];
    for (int I=0; I<n*n*n; I++){
        if (broken[I]){
            temp[0]=spins[I][0][0]+spins[I][1][0];
            temp[1]=spins[I][0][1]+spins[I][1][1];
            temp[2]=spins[I][0][2]+spins[I][1][2];
            vector[0]=vector[0]-0.5*temp[0]+spins[I][0][0];
            vector[1]=vector[1]-0.5*temp[1]+spins[I][0][1];
            vector[2]=vector[2]-0.5*temp[2]+spins[I][0][2];
            broken_count++;
        }
    if (broken_count>0) polar_order = sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]) / broken_count;
    else polar_order = 0.0;
    }
}
