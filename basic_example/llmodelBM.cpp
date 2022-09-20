//LL model with additional bent molecules
#include <iostream>
#include "llmodelBM.h"
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

    generator.seed(123213);
    random_I = std::uniform_int_distribution <int>(0,n*n*n-1);
    p = std::uniform_real_distribution<double>(0,1);
}

void LL_model_BM::shuffle_I(){
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
    double sum = 0.0;
    for (int k=0; k<dim; k++){
        sum+=spins[I][i][k]*spins[J][j][k];
    }
    return sum;
}

double LL_model_BM::E_IJ(int I, int J){
    return - ( 1.5 * 0.25 * (dot(I,0,J,0)*dot(I,0,J,0) + dot(I,1,J,0)*dot(I,1,J,0) + dot(I,0,J,1)*dot(I,0,J,1) + dot(I,1,J,1)*dot(I,1,J,1)) - 0.5);
}

void LL_model_BM::neighbors(int I){
    //define rules for neighbors here
    //updates neighbor_list with indeces
    //of neighbors of site I
    
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
    double sum=0.0;
    neighbors(I);
    for (int i=0; i<n_neighbors; i++){
        sum+=E_IJ(I, neighbors_list[i]);
    }
    return sum;
}

void LL_model_BM::H(){
    double sum=0;
    for (int I=0;  I<n*n*n; I++){
        sum+=E_neighbors(I);
    }
    E = 0.5*sum/(n*n*n);
}

void LL_model_BM::calculate_P2(){

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
    //int I = random_I(generator); //random mesto
        
    //save the previous config
    double temp[2][dim];
    for (int i=0; i<dim; i++){
        temp[0][i] = spins[I][0][i];
        temp[1][i] = spins[I][1][i];
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
    spins[I][0][0] = R00 * temp[0][0] +
                  R01 * temp[0][1] +
                  R02 * temp[0][2];

    spins[I][0][1] = R10 * temp[0][0] +
                  R11 * temp[0][1] +
                  R12 * temp[0][2];

    spins[I][0][2] = R20 * temp[0][0] +
                  R21 * temp[0][1] +   
                  R22 * temp[0][2];

    spins[I][1][0] = R00 * temp[1][0] +
                  R01 * temp[1][1] +
                  R02 * temp[1][2];

    spins[I][1][1] = R10 * temp[1][0] +
                  R11 * temp[1][1] +
                  R12 * temp[1][2];

    spins[I][1][2] = R20 * temp[1][0] +
                  R21 * temp[1][1] +   
                  R22 * temp[1][2];

    
        
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
    while (count > 0){
        int I = random_I(generator);
        if (!broken[I]){
            break_molecule(I);
            broken[I] = true;
            count = count - 1;
        }
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
    //thermalizes the lattice using the BarkerWatts mo
    for (int i = 0; i<N; i++){
        int acc=0;
        BarkerWatts_cycle();
        acceptance_rate = (double) acc / (n*n*n);
        adjust_rotation_angle();
    }
}

void LL_model_BM::break_molecule(int I){

    if (broken[I]) std::cout << "BIG PROBLEM!";
    //choose a random perpendicular axis of rotation and rotate the two vectors by angle +-alpha/2

    //choose a random axis
    //first find a vector perpendicular to spins[I][0] using the cross product with x system axis

    double v0 = 0;
    double v1 = spins[I][0][2] / sqrt(spins[I][0][2]*spins[I][0][2]+spins[I][0][1]*spins[I][0][1]);
    double v2 = -spins[I][0][1] / sqrt(spins[I][0][2]*spins[I][0][2]+spins[I][0][1]*spins[I][0][1]);

    //rotate the perpendicular vector by angle phi € [0,2pi) about spins[I][0] to get a random axis of rotation

    double n0 = spins[I][0][0];
    double n1 = spins[I][0][1];
    double n2 = spins[I][0][2];
        
    double gama = 2*p(generator)*pi;

    double R00 = (cos(gama)+n0*n0*(1-cos(gama)));
    double R01 = (n0*n1*(1-cos(gama)) - n2 * sin(gama));
    double R02 = (n0*n2*(1-cos(gama))+n1*sin(gama));

    double R10 = (n0*n1*(1-cos(gama))+n2*sin(gama));
    double R11 = (cos(gama) + n1*n1*(1-cos(gama)));
    double R12 = (n1*n2*(1-cos(gama))-n0*sin(gama)); 

    double R20 = (n0*n2*(1-cos(gama))-n1*sin(gama));
    double R21 = (n1*n2*(1-cos(gama))+ n0*sin(gama));
    double R22 = (cos(gama)+ n2*n2*(1-cos(gama)));

    double temp[3];
    temp[0]=v0;
    temp[1]=v1;
    temp[2]=v2;

    v0 = R00 * temp[0] +
         R01 * temp[1] +
         R02 * temp[2];

    v1 = R10 * temp[0] +
         R11 * temp[1] +
         R12 * temp[2];

    v2 = R20 * temp[0] +
         R21 * temp[1] +   
         R22 * temp[2];

    //rotate the spins by -+alpha/2 about (v0, v1, v2)

    double alpha = pi / 6;

    n0 = v0;
    n1 = v1;
    n2 = v2;
        
    gama = alpha/2;

    R00 = (cos(gama)+n0*n0*(1-cos(gama)));
    R01 = (n0*n1*(1-cos(gama)) - n2 * sin(gama));
    R02 = (n0*n2*(1-cos(gama))+n1*sin(gama));

    R10 = (n0*n1*(1-cos(gama))+n2*sin(gama));
    R11 = (cos(gama) + n1*n1*(1-cos(gama)));
    R12 = (n1*n2*(1-cos(gama))-n0*sin(gama)); 

    R20 = (n0*n2*(1-cos(gama))-n1*sin(gama));
    R21 = (n1*n2*(1-cos(gama))+ n0*sin(gama));
    R22 = (cos(gama)+ n2*n2*(1-cos(gama)));

    double temp1[3];
    temp1[0]=spins[I][0][0];
    temp1[1]=spins[I][0][1];
    temp1[2]=spins[I][0][2];

    //rotate the spin
    spins[I][0][0] = R00 * temp1[0] +
                     R01 * temp1[1] +
                     R02 * temp1[2];

    spins[I][0][1] = R10 * temp1[0] +
                     R11 * temp1[1] +
                     R12 * temp1[2];

    spins[I][0][2] = R20 * temp1[0] +
                     R21 * temp1[1] +   
                     R22 * temp1[2];


    gama = -alpha/2;

    R00 = (cos(gama)+n0*n0*(1-cos(gama)));
    R01 = (n0*n1*(1-cos(gama)) - n2 * sin(gama));
    R02 = (n0*n2*(1-cos(gama))+n1*sin(gama));

    R10 = (n0*n1*(1-cos(gama))+n2*sin(gama));
    R11 = (cos(gama) + n1*n1*(1-cos(gama)));
    R12 = (n1*n2*(1-cos(gama))-n0*sin(gama)); 

    R20 = (n0*n2*(1-cos(gama))-n1*sin(gama));
    R21 = (n1*n2*(1-cos(gama))+ n0*sin(gama));
    R22 = (cos(gama)+ n2*n2*(1-cos(gama)));

    double temp2[3];
    temp2[0]=spins[I][1][0];
    temp2[1]=spins[I][1][1];
    temp2[2]=spins[I][1][2];

    //rotate the spin
    spins[I][1][0] = R00 * temp2[0] +
                     R01 * temp2[1] +
                     R02 * temp2[2];

    spins[I][1][1] = R10 * temp2[0] +
                     R11 * temp2[1] +
                     R12 * temp2[2];

    spins[I][1][2] = R20 * temp2[0] +
                     R21 * temp2[1] +   
                     R22 * temp2[2];


   //std::cout << spins[I][1][0]*spins[I][0][0] + spins[I][1][1]*spins[I][0][1] + spins[I][1][2]*spins[I][0][2]; 

}

void LL_model_BM::align_molecule(int I){
    //aligns a broken molecule

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

    //std::cout << spins[I][1][0]*spins[I][0][0] + spins[I][1][1]*spins[I][0][1] + spins[I][1][2]*spins[I][0][2] << " "; 
}





