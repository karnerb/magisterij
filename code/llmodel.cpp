#include "llmodel.h"
#include <list>
#include <set>
#define dim 3
#define eps 1
#define n_neighbors 6
#define sc true
#define pi 3.141592653589793238462

LL_model::LL_model(int n) : n(n){
    spins = new double* [n*n*n];
    P2 = new double* [dim];
    neighbors_list = new int[n_neighbors];
    for (int i=0; i<dim; i++) P2[i] = new double [dim];
    for (int I=0; I<n*n*n; I++){
        spins[I] = new double [dim];
    }    
    generator.seed(123213);
    random_I = std::uniform_int_distribution <int>(0,n*n*n-1);
    p = std::uniform_real_distribution<double>(0, 1);
}

void LL_model::set_beta(double beta){
    this -> beta = beta;
}

double LL_model::dot(int I, int J){
    double sum = 0.0;
    for (int k=0; k<dim; k++){
        sum+=spins[I][k]*spins[J][k];
    }
    return sum;
}

double LL_model::E_IJ(int I, int J){
    return -eps * ( 1.5 * dot(I, J)*dot(I,J) - 0.5);
}

void LL_model::neighbors(int I){
    //define rules for neighbors here
    //updates neighbor_list with indeces
    //of neighbors of site I
    //simple cubic lattice with PBC here
    if (sc){
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
        neighbors_list[0]=up;
        neighbors_list[1]=down;
        neighbors_list[2]=north;
        neighbors_list[3]=south;
        neighbors_list[4]=east;
        neighbors_list[5]=west;
    }
}

double LL_model::E_neighbors(int I){
    double sum=0.0;
    neighbors(I);
    for (int i=0; i<n_neighbors; i++){
        sum+=E_IJ(I, neighbors_list[i]);
    }
    return sum;
}

void LL_model::H(){
    double sum=0;
    for (int I=0;  I<n*n*n; I++){
        sum+=E_neighbors(I);
    }
    E = 0.5*sum/(n*n*n);
}

void LL_model::calculate_P2(){

       for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                P2[i][j] = 0.0;            
            }
        }

        for (int I = 0; I<n*n*n; I++){
            for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                P2[i][j] = P2[i][j] + spins[I][i]*spins[I][j];
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

void LL_model::adjust_rotation_angle(){
    if (acceptance_rate>0.5){
        rotation_angle = rotation_angle + 0.1;
        if (rotation_angle > 2*pi) rotation_angle = 2*pi;
    }
    if (acceptance_rate<0.5){
        rotation_angle = rotation_angle - 0.1;
        if (rotation_angle<0) rotation_angle = 0.001;
    }
}

bool LL_model::BarkerWatts_move(){
    int I = random_I(generator); //random mesto
        
    //save the previous config
    double temp[dim];
    for (int i=0; i<dim; i++){
        temp[i] = spins[I][i];
    }
    
    //save the energy of current config
    double Eold = E_neighbors(I);

    //generate a random rotation and rotate vector at site I
    double theta = acos(2*p(generator)-1);
    double phi = 2 * pi * p(generator);

    double n0 = sin(theta)*cos(phi);
    double n1 = sin(theta)*sin(phi);
    double n2 = cos(theta);
        
    double gama = rotation_angle*(p(generator)-0.5);

    double R00 = (cos(gama)+n0*n0*(1-cos(gama)));
    double R01 = (n0*n1*(1-cos(gama)) - n2 * sin(gama));
    double R02 = (n0*n2*(1-cos(gama))+n1*sin(gama));

    double R10 = (n0*n1*(1-cos(gama))+n2*sin(gama));
    double R11 = (cos(gama) + n1*n1*(1-cos(gama)));
    double R12 = (n1*n2*(1-cos(gama))-n0*sin(gama)); 

    double R20 = (n0*n2*(1-cos(gama))-n1*sin(gama));
    double R21 = (n1*n2*(1-cos(gama))+ n0*sin(gama));
    double R22 = (cos(gama)+ n2*n2*(1-cos(gama)));

    //rotate the spin
    spins[I][0] = R00 * temp[0] +
                  R01 * temp[1] +
                  R02 * temp[2];

    spins[I][1] = R10 * temp[0] +
                  R11 * temp[1] +
                  R12 * temp[2];

    spins[I][2] = R20 * temp[0] +
                  R21 * temp[1] +   
                  R22 * temp[2];

    
        
    //calculate new energy 
    double Enew = E_neighbors(I);

    //Metropolis acceptance criterion
    //if Enew < Eold always accept 
    //else 
    if (Enew > Eold){
        double acceptance_p = p(generator);
        //reject with probability accoring to Metropolis
        //and keep the previous config
        if (acceptance_p > exp(-beta * (Enew-Eold))){
            spins[I][0] = temp[0];
            spins[I][1] = temp[1];
            spins[I][2] = temp[2];
            return false;
            }
        } 
    return true;
    }

void LL_model::initialize_lattice_random(){
    //sets the spins on the lattice in random directions
    for (int I=0; I<n*n*n; I++){
        double theta = acos(2*p(generator)-1);
        double phi = 2 * pi * p(generator);
        spins[I][0]=sin(theta)*cos(phi);
        spins[I][1]=sin(theta)*sin(phi);
        spins[I][2]=cos(theta);
    }
}

void LL_model::initialize_lattice_parallel(){
    //sets the spins in the same random direction
    double theta = acos(2*p(generator)-1);
    double phi = 2 * pi * p(generator);
    double n1 = sin(theta)*cos(phi);
    double n2 = sin(theta)*sin(phi);
    double n3 = cos(theta);    
    for (int I=0; I<n*n*n; I++){
        spins[I][0]=n1;
        spins[I][1]=n2;
        spins[I][2]=n3;
    }        
}

void LL_model::thermalize_BarkerWatts(int N){
    //thermalizes the lattice using the BarkerWatts mo
    for (int i = 0; i<N; i++){
        int acc=0;
        for (int j=0; j<n*n*n; j++){
            if(BarkerWatts_move()) acc++;
        }
        acceptance_rate = (double) acc / (n*n*n);
        adjust_rotation_angle();
    }
}

void LL_model::thermalize_cluster(int N){
    for (int i = 0; i<N; i++){
        cluster_move();
    }
}

void LL_model::cluster_move(){
    std::list<int> sites_to_visit_neighbors_of = {};
    std::set<int> cluster_members = {};
    int startI = random_I(generator);

    cluster_members.insert(startI);
    sites_to_visit_neighbors_of.push_back(startI);

    //reflect the initial particle
    double theta = acos(2*p(generator)-1);
    double phi = 2 * pi * p(generator);
    //random unit vector
    double n0 = sin(theta)*cos(phi);
    double n1 = sin(theta)*sin(phi);
    double n2 = cos(theta);
    double dot = n0*spins[startI][0] + n1*spins[startI][1] + n2*spins[startI][2];
    //reflection
    spins[startI][0] = -spins[startI][0] + 2 * dot * n0;
    spins[startI][1] = -spins[startI][1] + 2 * dot * n1;
    spins[startI][2] = -spins[startI][2] + 2 * dot * n2;


    while (!sites_to_visit_neighbors_of.empty()){
        
        int currentI = sites_to_visit_neighbors_of.back();
        sites_to_visit_neighbors_of.pop_back();
        //calculate neighbors of currentI
        neighbors(currentI);
        //loop over neighbors
        for (int l = 0; l<6; l++){
            //if neighbor not in cluster
            if (!cluster_members.count(neighbors_list[l])==1){

                double dot1 = spins[currentI][0] * spins[neighbors_list[l]][0] + 
                              spins[currentI][1] * spins[neighbors_list[l]][1] + 
                              spins[currentI][2] * spins[neighbors_list[l]][2]; 
                
                double dot3 = n0*spins[neighbors_list[l]][0] + n1*spins[neighbors_list[l]][1] + n2*spins[neighbors_list[l]][2];

                double t1 = -spins[neighbors_list[l]][0] + 2* dot3 * n0;
                double t2 = -spins[neighbors_list[l]][1] + 2* dot3 * n1;
                double t3 = -spins[neighbors_list[l]][2] + 2* dot3 * n2;

                double dot2 = spins[currentI][0] * t1 +
                              spins[currentI][1] * t2 +
                              spins[currentI][2] * t3;

                double exponent = 3.0/2.0 * beta * ( dot1*dot1 - dot2*dot2);
                if (exponent < 0){                
                //if visit - detailed balance
                    if (p(generator) < 1-exp(exponent)){
                    //add the site to cluster members and to list of sites_to_visit_neighbors_of
                    cluster_members.insert(neighbors_list[l]);
                    sites_to_visit_neighbors_of.push_back(neighbors_list[l]);
                    //flip
                    spins[neighbors_list[l]][0]=t1;
                    spins[neighbors_list[l]][1]=t2;
                    spins[neighbors_list[l]][2]=t3;
                    
                }
                }
            }

        }
        }
    }

