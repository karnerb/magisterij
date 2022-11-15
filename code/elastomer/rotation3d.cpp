#include <cmath>
#include "LC_Elastomer.h"

void LC_Elastomer::rotation3D(double* vector, double n0, double n1, double n2, double gama){
    //rotates double vector[3] by angle gama about a unit vector (n0, n1, n2)
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
    temp[0]=vector[0];
    temp[1]=vector[1];
    temp[2]=vector[2];

    vector[0] = R00 * temp[0] +
                R01 * temp[1] +
                R02 * temp[2];

    vector[1] = R10 * temp[0] +
                R11 * temp[1] +
                R12 * temp[2];

    vector[2] = R20 * temp[0] +
                R21 * temp[1] +   
                R22 * temp[2];

}
