#include "llmodelBM.h"
#include <iostream>

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

void LL_model_BM::cluster_count_and_size(){
    int cluster[n*n*n];
    int label=0;
    int properLabels[n*n*n];
    for (int i=0; i<n*n*n; i++) cluster[i]=-1;
    
    for (int I=0; I<n*n*n; I++){
        int minimal_neighboring_label = label;
        if (broken[I]){
            //visit neighbors of site I
            for (int j=0; j<6; j++){
                int J = neighbors_list[I][j];
                if (broken[J]){
                    if (cluster[J]!=-1){
                        //check if neighbor is of broken type and if it already is in a cluster
                        if (cluster[J]<minimal_neighboring_label) minimal_neighboring_label = cluster[J];
                    }
                }
            }
        //label the site I with the minimal_neighboring_label (if no neighbors in cluster, minmal_neighboring_label = label) 
        cluster[I]=minimal_neighboring_label;
        properLabels[minimal_neighboring_label]=minimal_neighboring_label;

        //visit neighbors of site I and set their properLabels to minimal_neighboring_label to connect clusters
        for (int j=0; j<6; j++){
            int J=neighbors_list[I][j];
            if (broken[J]){
                if (cluster[J]!=-1){
                    properLabels[cluster[J]]=minimal_neighboring_label;
                    }
                }
            }
            //if minimal_neighboring_label==label no neighbors are in cluster yet - start a new cluster 
            if (minimal_neighboring_label==label) label++;    
        }
    }
   //relabel cluster[] with properLabels
   for (int I=0; I<n*n*n; I++){
        if (broken[I]){
            int current_label = cluster[I];
            while (current_label!=properLabels[current_label]) current_label=properLabels[current_label];
            cluster[I] = current_label;
       }
   }
   //find the number of clusters
   int max_cluster_label=-1;
   for (int I=0; I<n*n*n; I++){
       if (cluster[I]>max_cluster_label) max_cluster_label = cluster[I];
   }
   //offset by one
   max_cluster_label++;
   //find the average cluster size
   broken_cluster_count = max_cluster_label;
   
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
