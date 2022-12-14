#include <thread>
#include "../llmodelBM/llmodelBM.h"
#include "../logger/logger.h"
#define size 30
#define kilo 1000

void fabri_zannoni_heatup(){
    
    int temperatures_count = 30;

    double temperatures[temperatures_count] = {1.0000, 1.0500, 1.1000, 1.1100, 1.1140, 1.1160,
                                               1.1180, 1.1200, 1.1210, 1.1220, 1.1230, 1.1235,
                                               1.1238, 1.1240, 1.1250, 1.1260, 1.1270, 1.1280,
                                               1.1290, 1.1300, 1.1320, 1.1340, 1.1360, 1.1380,
                                               1.1400, 1.1600, 1.1800, 1.2000, 1.2500, 1.3000};
    
    int cycles[temperatures_count] = {10, 10, 10, 10, 13, 14,
                                      10, 15, 20, 25, 50, 70,
                                      70, 70, 30, 25, 25, 25,
                                      25, 15, 15, 15, 15, 15,
                                      15, 15, 15, 15, 15, 15};
    
    for (int i=0; i<temperatures_count; i++) cycles[i]=cycles[i]*kilo;


    LL_model_BM model_heatup(size);
    Logger<LL_model_BM> logger_heatup("fabri_zannoni_test_heatup.txt", model_heatup);
    model_heatup.initialize_lattice_parallel();
    model_heatup.set_beta(1.0/temperatures[0]);
    model_heatup.thermalize_BarkerWatts(cycles[0]);


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
    }
    
}

void fabri_zannoni_cooldown(){
    
    int temperatures_count = 30;

    double temperatures[temperatures_count] = {1.0000, 1.0500, 1.1000, 1.1100, 1.1140, 1.1160,
                                               1.1180, 1.1200, 1.1210, 1.1220, 1.1230, 1.1235,
                                               1.1238, 1.1240, 1.1250, 1.1260, 1.1270, 1.1280,
                                               1.1290, 1.1300, 1.1320, 1.1340, 1.1360, 1.1380,
                                               1.1400, 1.1600, 1.1800, 1.2000, 1.2500, 1.3000};
    
    int cycles[temperatures_count] = {10, 10, 10, 10, 13, 14,
                                      10, 15, 20, 25, 50, 70,
                                      70, 70, 30, 25, 25, 25,
                                      25, 15, 15, 15, 15, 15,
                                      15, 15, 15, 15, 15, 15};
    
    for (int i=0; i<temperatures_count; i++) cycles[i]=cycles[i]*kilo;


    LL_model_BM model_cooldown(size);
    Logger<LL_model_BM> logger_cooldown("fabri_zannoni_test_cooldown.txt", model_cooldown);
    model_cooldown.initialize_lattice_random();
    model_cooldown.set_beta(1.0/temperatures[temperatures_count-1]);
    model_cooldown.thermalize_BarkerWatts(cycles[temperatures_count-1]);

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
    }
    
}

void run(){
    std::thread t1(fabri_zannoni_cooldown);
    std::thread t2(fabri_zannoni_heatup);
    t1.join();
    t2.join();
}

