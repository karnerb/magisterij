#include "llmodel.h"
#include "logger.h"
#include <iostream>

/* void temperature_dependence(){
    
    LL_model model(10);
    std::cout << " asdjjas";
    Logger logger("output.txt", model);
    int Nthermal = 100;
    int N = 1000;
    
    double T0 = 0.1;
    double Tmax = 3;
    int nsteps = 100;
    std::cout << "thermalizing0";
    model.initialize_lattice_parallel();
    for (int i=0; i<nsteps; i++){
        double T = T0 + (Tmax-T0)*(double) i / nsteps;
        model.set_beta(1/T);
        std::cout << "thermalizing";
        model.thermalize_BarkerWatts(Nthermal*model.n*model.n*model.n);
        for (int j=0; j<N; j++){
            std::cout << "i: " << i << "\n";
            for (int k=0; k<(model.n*model.n*model.n); k++){
                model.BarkerWatts_move();
            }
            model.calculate_P2();
            model.H();
            logger.log_params();
        }
    }
} */

void temperature_dependence(int n, double Tstart, double Tend, std::string filename){
    LL_model model(n);
    Logger logger(filename, model);
    int Nthermal = 100;
    int N = 10000;
    model.initialize_lattice_parallel();
    
    for (int i=0; i<nsteps; i++){
        std::cout << "i: " << i << "\n" << std::flush;
        double T = T0 + (Tmax-T0)*(double) i / nsteps;
        model.set_beta(1/T);
        model.cycle=0;
        
        model.thermalize_BarkerWatts(Nthermal);
        for (int j=0; j<N; j++){
            for (int k=0; k<(model.n*model.n*model.n); k++){
                model.BarkerWatts_move();
            }
            model.cycle++;
            model.calculate_P2();
            model.H();
            //std::cout << model.E << "\n";
            logger.log_params();
        }
    }
}

int main(void){
    
    LL_model model(5);
    Logger logger("output.txt", model);
    int Nthermal = 1000;
    int N = 1000;
    
    double T0 = 0.1;
    double Tmax = 3;
    int nsteps = 100;
    
    model.initialize_lattice_parallel();
    for (int i=0; i<nsteps; i++){
        std::cout << "i: " << i << "\n" << std::flush;
        double T = T0 + (Tmax-T0)*(double) i / nsteps;
        model.set_beta(1/T);
        model.cycle=0;
        
        model.thermalize_BarkerWatts(Nthermal);
        for (int j=0; j<N; j++){
            for (int k=0; k<(model.n*model.n*model.n); k++){
                model.BarkerWatts_move();
            }
            model.cycle++;
            model.calculate_P2();
            model.H();
            //std::cout << model.E << "\n";
            logger.log_params();
        }
    }
//temperature_dependence();
    
}