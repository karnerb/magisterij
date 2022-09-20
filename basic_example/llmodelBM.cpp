//LL model with additional bent molecules
#include <iostream>
#include "llmodelBM.h"
#include "rotation3d.cpp"
#define dim 3
#define eps 1
#define n_neighbors 6
#define sc true
#define pi 3.141592653589793238462

LL_model_BM::LL_model_BM(int n) : n(n){
    spins = new double** [n*n*n];
    P2 = new double* [dim];
    broken = new bool [n*n*n];
    neighbors_list = new int[n_neighbors];
    shuffled_I = new int[n*n*n];
    for (int i=0; i<dim; i++) P2[i] = new double [dim];
    for (int I=0; I<n*n*n; I++){
        shuffled_I[I]=I;
        broken[I] = false;
        spins[I] = new double* [2];
    }
    for (int I=0; I<n*n*n; I++){
        spins[I][0] = new double[dim];
        spins[I][1] = new double[dim];

    }
    //TODO:: add time dependent seeding
    generator.seed(123213);
    random_I = std::uniform_int_distribution <int>(0,n*n*n-1);
    p = std::uniform_real_distribution<double>(0,1);
}

void LL_model_BM::shuffle_I(){
    //randomly shuffles elements of shuffled_I[n*n*n]
    int temp, j;
    for (int i=n*n*n-1; i>0; i--){
        temp=shuffled_I[i];
        j = random_I(generator);
        shuffled_I[i]=shuffled_I[j];
        shuffled_I[j]=temp;
    }
}

void LL_model_BM::set_beta(double beta){
    this -> beta = beta;
}

double LL_model_BM::dot(int I, int i,  int J, int j){
    //dot product of spins[I][i] and spins[J][j]
    double sum = 0.0;
    for (int k=0; k<dim; k++){
        sum+=spins[I][i][k]*spins[J][j][k];
    }
    return sum;
}

double LL_model_BM::E_IJ(int I, int J){
    //interaction energy of sites I, J
    return - ( 1.5 * 0.25 * (dot(I,0,J,0)*dot(I,0,J,0) + dot(I,1,J,0)*dot(I,1,J,0) + dot(I,0,J,1)*dot(I,0,J,1) + dot(I,1,J,1)*dot(I,1,J,1)) - 0.5);
}

void LL_model_BM::neighbors(int I){
    //define rules for neighbors here
    //updates neighbor_list with indeces
    //of neighbors of site I
    //TODO: implement a neighbor list for all sites
    //simple cubic here
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

double LL_model_BM::E_neighbors(int I){
    //sums the interaction energies with nearest neighbors
    double sum=0.0;
    neighbors(I);
    for (int i=0; i<n_neighbors; i++){
        sum+=E_IJ(I, neighbors_list[i]);
    }
    return sum;
}

void LL_model_BM::H(){
    //sums E_neighbors over all lattice sites
    double sum=0;
    for (int I=0;  I<n*n*n; I++){
        sum+=E_neighbors(I);
    }
    E = 0.5*sum/(n*n*n);
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

bool LL_model_BM::BarkerWatts_move(int I){    
    //save the previous config
    double temp[2][dim];
    for (int i=0; i<dim; i++){
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

bool LL_model_BM::switch_move(){
    //select a random site and attempt to switch with a random neighbor site
    int I = random_I(generator);
    neighbors(I);
    //choose a random neighbor
    
    int neighbor = floor(p(generator)*6.0);
    if (neighbor > 6 || neighbor < 0) std::cout << "BIG PROBLEM";
    int J = neighbors_list[neighbor];

    //store the current energy of the pair
    double Eold = E_neighbors(I) + E_neighbors(J);

    //switch the two sites
    switch_sites(I, J);

    //calculate the new energy of the pair
    double Enew = E_neighbors(I) + E_neighbors(J);

    double deltaE = Enew-Eold;
    //accept if deltaE < 0 otherwise Metropolis criterion

    if (Enew > 0){
        if (p(generator) > exp(-beta*deltaE)){
            //revert to previous config
            switch_sites(I, J);
            return false;
        }
    }
    return true;
}

void LL_model_BM::switch_sites(int I, int J){
    
    double temp[2][3];

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

void LL_model_BM::thermalize_BarkerWatts(int N){
    //thermalizes the lattice using the BarkerWatts cycle
    for (int i = 0; i<N; i++){
        int acc=0;
        BarkerWatts_cycle();
        acceptance_rate = (double) acc / (n*n*n);
        adjust_rotation_angle();
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
    double alpha = pi / 6;
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





