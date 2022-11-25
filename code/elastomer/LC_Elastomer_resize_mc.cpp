#include "LC_Elastomer.h"
#include <algorithm>



bool LC_Elastomer::resize_move(){    
   double new_lambda = lambda + (p(generator)-0.5)*resize_step;
   while (new_lambda<0) new_lambda = lambda + (p(generator)-0.5)*resize_step;
   //calculate old energies
   H_elastic();
   H_coupling();
   double Eold = Energy_elastic + Energy_coupling;
   double old_lambda = lambda;
   lambda = new_lambda;

   //calculate new energies
   H_elastic();
   H_coupling();
   double Enew = Energy_elastic + Energy_coupling;
   double eps = 1.0;
   double K = Enew-Eold + n*n*n*eps*sigma*(old_lambda-new_lambda);
   if (K>0){
       if (p(generator) > exp(-beta * K)){
           lambda = old_lambda;
           return false;
           }
   }
   return true;
}


void LC_Elastomer::adjust_resize_step(){
    int accepted=0;
    int total = 100;
    for (int i=0; i<total; i++){
        if (resize_move()) accepted++;
    }
    if ((double)accepted/total > 0.5){
        resize_step = resize_step*1.05;
    }
    else{
        resize_step = resize_step*0.95;
    }
}