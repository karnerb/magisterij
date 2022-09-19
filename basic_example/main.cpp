#include "llmodel.h"
#include "llmodelBM.h"
#include "logger.h"
#include <iostream>
#include <chrono>

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

void temperature_dependence_BW(int n, int nsteps, double Tstart, double Tend, std::string filename, std::string start_phase){
    LL_model model(n);
    Logger<LL_model> logger(filename, model);
    int Nthermal = 1000;
    int N = 1000;
    
    if (start_phase == "isotropic") model.initialize_lattice_random();
    else model.initialize_lattice_parallel();

    for (int i=0; i<nsteps; i++){
        std::cout << "i: " << i << "/" << nsteps << "\r" << std::flush;
        double T = Tstart + (Tend-Tstart)*(double) i / nsteps;
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

void temperature_dependence_cluster(int n, int nsteps, double Tstart, double Tend, std::string filename, std::string start_phase){
    LL_model model(n);
    Logger<LL_model> logger(filename, model);
    int Nthermal = 100;
    int N = 1000;
    
    if (start_phase == "isotropic") model.initialize_lattice_random();
    else model.initialize_lattice_parallel();
    
    model.thermalize_cluster(Nthermal);

    for (int i=0; i<nsteps; i++){
        std::cout << "i: " << i << "/" << nsteps << "\r" << std::flush;
        double T = Tstart + (Tend-Tstart)*(double) i / nsteps;
        model.set_beta(1/T);
        model.cycle=0;
        model.thermalize_cluster(Nthermal);
        for (int j=0; j<N; j++){
            for (int k=0; k<(model.n); k++){
                model.cluster_move();
            }
            model.cycle++;
            model.calculate_P2();
            model.H();
            //std::cout << model.E << "\n";
            logger.log_params();
        }
    }
}


void temperature_dependence_BM(int n, int nsteps, double Tstart, double Tend, std::string filename, std::string start_phase, double ratioBM, bool switchmove){
    LL_model_BM model(n);
    Logger<LL_model_BM> logger(filename, model);
    int Nthermal = 1000;
    int N = 5000;
    
    if (start_phase == "isotropic") model.initialize_lattice_random();
    else model.initialize_lattice_parallel();
    model.break_molecules((int)n*n*n*ratioBM);


    for (int i=0; i<nsteps; i++){
        int switch_moves=0;
        int successful_switch_moves=0;
        
        double T = Tstart + (Tend-Tstart)*(double) i / nsteps;
        model.set_beta(1/T);
        model.cycle=0;
        model.thermalize_BarkerWatts(Nthermal);
        for (int j=0; j<N; j++){
            for (int k=0; k<(model.n*model.n*model.n); k++){
                model.BarkerWatts_move();
            }
            if (switchmove){
            for (int k=0; k<model.n*model.n; k++){
                switch_moves++;
                if (model.switch_move()) successful_switch_moves++;
            }}
            model.cycle++;
            model.calculate_P2();
            model.H();
            //std::cout << model.E << "\n";
            logger.log_params();
        }
        std::cout << "i: " << i << "/" << nsteps << "switch moves: " << successful_switch_moves << "/" << switch_moves << "\r" << std::flush;
    }
    }

int main(void){
    
    //auto start = std::chrono::high_resolution_clock::now();
    //temperature_dependence_cluster(30, 50, 1.15, 1.1, "cooldown_zoom_cluster.txt", "isotropic");
    //auto end = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::minutes>(end-start);
    //std::cout << "time: " << duration.count() << "minutes\n";

    //temperature_dependence_BW(30, 100, 1.15, 1.1, "cooldown_zoom.txt", "isotropic");
    //temperature_dependence_cluster(20, 100, 2, 0.1, "cluster_cooldown.txt", "isotropic");
    //temperature_dependence_cluster(20, 100, 0.1, 2, "cluster_heatup.txt", "nematic");
    temperature_dependence_BM(10, 10, 1.0, 1.3, "00.txt", "nematic", 0.0, false);
    temperature_dependence_BM(10, 10, 1.0, 1.3, "01.txt", "nematic", 0.5, false);
    temperature_dependence_BM(10, 10, 1.0, 1.3, "05.txt", "nematic", 0.5, true);
    

    
}