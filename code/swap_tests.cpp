#include <thread>
#include "llmodel.h"
#include "llmodelBM.h"
#include "logger.h"
#define size 10
#define kilo 100

void ceplak_swap_heatup(double ratio, std::string filename, int sim_id){

    int temperatures_count = 70;

    double temperatures[temperatures_count] = {1.000000,  1.010000,  1.020000,  1.030000,  1.040000,  1.050000,
                                               1.051754,  1.053509,  1.055263,  1.057018,  1.058772,  1.060526,
                                               1.062281,  1.064035,  1.065789,  1.067544,  1.069298,  1.071053,
                                               1.072807,  1.074561,  1.076316,  1.078070,  1.079825,  1.081579,
                                               1.083333,  1.085088,  1.086842,  1.088596,  1.090351,  1.092105,
                                               1.093860,  1.095614,  1.097368,  1.099123,  1.100877,  1.102632,
                                               1.104386,  1.106140,  1.107895,  1.109649,  1.111404,  1.113158,
                                               1.114912,  1.116667,  1.118421,  1.120175,  1.121930,  1.123684,
                                               1.125439,  1.127193,  1.128947,  1.130702,  1.132456,  1.134211,
                                               1.135965,  1.137719,  1.139474,  1.141228,  1.142982,  1.144737,
                                               1.146491,  1.148246,  1.150000,  1.160000,  1.170000,  1.180000,
                                               1.190000,  1.200000,  1.250000,  1.300000};
    
    int cycles[temperatures_count] = {50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
                                      50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
                                      50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
                                      50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
                                      50, 50};
    
    for (int i=0; i<temperatures_count; i++) cycles[i]=cycles[i]*kilo;


    LL_model_BM model_heatup(size);
    Logger<LL_model_BM> logger_heatup(filename, model_heatup);
    model_heatup.initialize_lattice_parallel();
    model_heatup.break_molecules((int) (ratio*size*size*size));
    model_heatup.set_beta(1.0/temperatures[0]);
    model_heatup.thermalize_BarkerWatts(cycles[0]*5);


    for (int i=0; i<temperatures_count; i++){
        model_heatup.set_beta(1.0/temperatures[i]);
        model_heatup.thermalize_BarkerWatts(cycles[i]);
        model_heatup.cycle = 0;
        for (int j=0; j<cycles[i]; j++){
            model_heatup.BarkerWatts_cycle();
            model_heatup.calculate_P2();
            model_heatup.H();
            model_heatup.calculate_polar_order();
            logger_heatup.log_params();
        }
         if (sim_id==0) std::cout << "heatup: " << i << "/" << temperatures_count << "\r" << std::flush;
    }
    
}

void ceplak_swap_cooldown(double ratio, std::string filename, int sim_id){
    
    int temperatures_count = 70;

    double temperatures[temperatures_count] = {1.000000,  1.010000,  1.020000,  1.030000,  1.040000,  1.050000,
                                               1.051754,  1.053509,  1.055263,  1.057018,  1.058772,  1.060526,
                                               1.062281,  1.064035,  1.065789,  1.067544,  1.069298,  1.071053,
                                               1.072807,  1.074561,  1.076316,  1.078070,  1.079825,  1.081579,
                                               1.083333,  1.085088,  1.086842,  1.088596,  1.090351,  1.092105,
                                               1.093860,  1.095614,  1.097368,  1.099123,  1.100877,  1.102632,
                                               1.104386,  1.106140,  1.107895,  1.109649,  1.111404,  1.113158,
                                               1.114912,  1.116667,  1.118421,  1.120175,  1.121930,  1.123684,
                                               1.125439,  1.127193,  1.128947,  1.130702,  1.132456,  1.134211,
                                               1.135965,  1.137719,  1.139474,  1.141228,  1.142982,  1.144737,
                                               1.146491,  1.148246,  1.150000,  1.160000,  1.170000,  1.180000,
                                               1.190000,  1.200000,  1.250000,  1.300000};
    
    int cycles[temperatures_count] = {50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
                                      50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
                                      50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
                                      50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
                                      50, 50};

    for (int i=0; i<temperatures_count; i++) cycles[i]=cycles[i]*kilo;


    LL_model_BM model_cooldown(size);
    Logger<LL_model_BM> logger_cooldown(filename, model_cooldown);
    model_cooldown.initialize_lattice_random();
    model_cooldown.break_molecules((int) (ratio*size*size*size));
    model_cooldown.set_beta(1.0/temperatures[temperatures_count-1]);
    model_cooldown.thermalize_BarkerWatts(cycles[temperatures_count-1]*5);

    for (int i=temperatures_count-1; i>=0; i--){
        model_cooldown.set_beta(1.0/temperatures[i]);
        model_cooldown.thermalize_BarkerWatts(cycles[i]);
        model_cooldown.cycle=0;
        for (int j=0; j<cycles[i]; j++){
            model_cooldown.BarkerWatts_cycle();
            model_cooldown.calculate_P2();
            model_cooldown.H();
            model_cooldown.calculate_polar_order();
            logger_cooldown.log_params();
        }
        if (sim_id==0) std::cout << "cooldown: " << i << "/" << temperatures_count << "\r" << std::flush;
    }
    
}

void ceplak_with_swaps(){
    double ratios[12] = {0.0000, 0.0019, 0.0037,
                         0.0074, 0.0093, 0.0111,
                         0.0185, 0.0278, 0.0370,
                         0.0463, 0.0556, 0.0741};
    
    std::string filenames_heatup[12] = {"ceplak_heatup_r_0.00.txt", "ceplak_heatup_r_0.19.txt", "ceplak_heatup_r_0.37.txt",
                                        "ceplak_heatup_r_0.74.txt", "ceplak_heatup_r_0.93.txt", "ceplak_heatup_r_1.11.txt",
                                        "ceplak_heatup_r_1.85.txt", "ceplak_heatup_r_2.78.txt", "ceplak_heatup_r_3.70.txt",
                                        "ceplak_heatup_r_4.63.txt", "ceplak_heatup_r_5.56.txt", "ceplak_heatup_r_7.41.txt"};

    std::string filenames_cooldown[12] = {"ceplak_cooldown_r_0.00.txt", "ceplak_cooldown_r_0.19.txt", "ceplak_cooldown_r_0.37.txt",
                                          "ceplak_cooldown_r_0.74.txt", "ceplak_cooldown_r_0.93.txt", "ceplak_cooldown_r_1.11.txt",
                                          "ceplak_cooldown_r_1.85.txt", "ceplak_cooldown_r_2.78.txt", "ceplak_cooldown_r_3.70.txt",
                                          "ceplak_cooldown_r_4.63.txt", "ceplak_cooldown_r_5.56.txt", "ceplak_cooldown_r_7.41.txt"};

    #pragma omp parallel for
    for (int i=0; i<12; i++){
        ceplak_swap_heatup(ratios[i], filenames_heatup[i], i);
    }


    #pragma omp parallel for
    for (int i=0; i<12; i++){
        ceplak_swap_cooldown(ratios[i], filenames_cooldown[i], i);
    }
}


void swap_moves_thermalization(std::string filename){
    int temperatures_count = 5;
    LL_model_BM model(size);
    Logger<LL_model_BM> logger(filename, model);
    double ratio = 0.00;
    model.initialize_lattice_parallel();
    model.break_molecules((int) (ratio*size*size*size));
    double temperatures[6] = {1.0, 1.05, 1.0925, 1.1125, 1.20};
    model.set_beta(1.0/temperatures[0]);
    model.thermalize_BarkerWatts(1000);
    

    for (int i=0; i<temperatures_count; i++){
        model.set_beta(1.0/temperatures[i]);
        model.thermalize_BarkerWatts(1000);
        int swap_cycles = 1000;
        model.cycle = 0;
        for (int j=0; j<swap_cycles; j++){
        model.H();
        model.calculate_P2();
        logger.log_params();
        model.swap_cycle();
        model.cycle++;
        }
    }
}
    
