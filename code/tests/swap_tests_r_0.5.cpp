#include <thread>
#include "../llmodelBM/llmodelBM.h"
#include "../logger/logger.h"
#define size 30
#define temperatures_count 30

void cooldown_swap(double ratio, double* temperatures, int* cycles, std::string filename, int sim_id){
    
    LL_model_BM model(size);
    model.initialize_lattice_parallel();
    model.break_molecules((int) size*size*size*ratio);
    Logger<LL_model_BM> logger(filename, model);
    model.set_beta(1.0/temperatures[temperatures_count-1]);
    model.thermalize_BarkerWatts(cycles[temperatures_count-1]*5);
    for (int i=0; i<cycles[temperatures_count-1]*5; i++) model.swap_cycle();
    
    for (int i=temperatures_count-1; i>=0; i--){
        model.set_beta(1.0/temperatures[i]);
        model.thermalize_BarkerWatts(cycles[i]);
        model.cycle=0;
        for (int j=0; j<cycles[i]; j++){
            model.BarkerWatts_cycle();
            model.swap_cycle();
            model.calculate_P2();
            model.H();
            model.calculate_polar_order();
            model.count_neighbors_of_same_kind();
            model.cluster_count_and_size();
            logger.log_params();

        }
        logger.output_lattice_configuration();
        if (sim_id==0) std::cout << "cooldown with swap: " << i << "/" << std::to_string(temperatures_count) << "\r" << std::flush;
    }
}

void cooldown(double ratio, double* temperatures, int* cycles, std::string filename, int sim_id){
    
    
    LL_model_BM model(size);
    model.initialize_lattice_parallel();
    model.break_molecules((int) size*size*size*ratio);
    Logger<LL_model_BM> logger(filename, model);
    model.set_beta(1.0/temperatures[temperatures_count-1]);
    model.thermalize_BarkerWatts(cycles[temperatures_count-1]*5);
    
    
    for (int i=temperatures_count-1; i>=0; i--){
        model.set_beta(1.0/temperatures[i]);
        model.thermalize_BarkerWatts(cycles[i]);
        model.cycle=0;
        for (int j=0; j<cycles[i]; j++){
            model.BarkerWatts_cycle();
            model.calculate_P2();
            model.H();
            model.calculate_polar_order();
            model.count_neighbors_of_same_kind();
            model.cluster_count_and_size();
            logger.log_params();
           
        }
        //logger.output_lattice_configuration();
        if (sim_id==0) std::cout << "cooldown with swap: " << i << "/" << std::to_string(temperatures_count) << "\r" << std::flush;
    }
}

void heatup_swap(double ratio, double* temperatures, int* cycles, std::string filename, int sim_id){
    
    LL_model_BM model(size);
    model.initialize_lattice_parallel();
    //model.break_molecules((int) size*size*size*ratio);
    model.break_molecules(ratio*size*size*size);
    Logger<LL_model_BM> logger(filename, model);
    model.set_beta(1.0/temperatures[0]);
    model.thermalize_BarkerWatts(cycles[0]*5);
    for (int i=0; i<cycles[0]*5; i++) model.swap_cycle();
    
    for (int i=0; i<temperatures_count; i++){
        model.set_beta(1.0/temperatures[i]);
        model.thermalize_BarkerWatts(cycles[i]);
        model.cycle=0;
        for (int j=0; j<cycles[i]; j++){
            model.BarkerWatts_cycle();
            model.swap_cycle();
            model.calculate_P2();
            model.H();
            model.calculate_polar_order();
            model.count_neighbors_of_same_kind();
            model.cluster_count_and_size();
            logger.log_params();
        }
        logger.output_lattice_configuration();
        if (sim_id==0) std::cout << "heatup with swap: " << i << "/" << std::to_string(temperatures_count) << "\r" << std::flush;
    }
}

void heatup(double ratio, double* temperatures, int* cycles, std::string filename, int sim_id){
    
    
    LL_model_BM model(size);
    model.initialize_lattice_parallel();
    model.break_molecules((int) size*size*size*ratio);
    Logger<LL_model_BM> logger(filename, model);
    model.set_beta(1.0/temperatures[0]);
    model.thermalize_BarkerWatts(cycles[0]*5);
    
    for (int i=0; i<temperatures_count; i++){
        model.set_beta(1.0/temperatures[i]);
        model.thermalize_BarkerWatts(cycles[i]);
        model.cycle=0;
        for (int j=0; j<cycles[i]; j++){
            model.BarkerWatts_cycle();
            model.calculate_P2();
            model.H();
            model.calculate_polar_order();
            model.count_neighbors_of_same_kind();
            model.cluster_count_and_size();
            logger.log_params();
        }
        //logger.output_lattice_configuration();
        if (sim_id==0) std::cout << "heatup: " << i << "/" << std::to_string(temperatures_count) << "\r" << std::flush;
    }
}


void cooldown_then_heatup_swap(double ratio, double* temperatures, int* cycles, std::string filename1, std::string filename2, int sim_id){
    
    LL_model_BM model(size);
    model.initialize_lattice_parallel();
    model.break_molecules((int) size*size*size*ratio);
    Logger<LL_model_BM> logger(filename1, model);
    Logger<LL_model_BM> logger2(filename2, model);
    model.set_beta(1.0/temperatures[temperatures_count-1]);
    model.thermalize_BarkerWatts(cycles[temperatures_count-1]*5);
    for (int i=0; i<cycles[temperatures_count-1]*5; i++) model.swap_cycle();
    
    for (int i=temperatures_count-1; i>=0; i--){
        model.set_beta(1.0/temperatures[i]);
        model.thermalize_BarkerWatts(cycles[i]);
        model.cycle=0;
        for (int j=0; j<cycles[i]; j++){
            model.BarkerWatts_cycle();
            model.swap_cycle();
            model.calculate_P2();
            model.H();
            model.calculate_polar_order();
            model.count_neighbors_of_same_kind();
            model.cluster_count_and_size();
            logger.log_params();

        }
        logger.output_lattice_configuration();
        if (sim_id==0) std::cout << "cooldown with swap: " << i << "/" << std::to_string(temperatures_count) << "\r" << std::flush;
    }
    for (int i=0; i<temperatures_count; i++){

        model.set_beta(1.0/temperatures[i]);
        model.thermalize_BarkerWatts(cycles[i]);
        model.cycle=0;
        for (int j=0; j<cycles[i]; j++){
            model.BarkerWatts_cycle();
            model.swap_cycle();
            model.calculate_P2();
            model.H();
            model.calculate_polar_order();
            model.count_neighbors_of_same_kind();
            model.cluster_count_and_size();
            logger2.log_params();

        }
        logger2.output_lattice_configuration();
        if (sim_id==0) std::cout << "heatup with swap: " << i << "/" << std::to_string(temperatures_count) << "\r" << std::flush;
    }
}


void cooldown_then_heatup(double ratio, double* temperatures, int* cycles, std::string filename1, std::string filename2, int sim_id){
    
    LL_model_BM model(size);
    model.initialize_lattice_parallel();
    model.break_molecules((int) size*size*size*ratio);
    Logger<LL_model_BM> logger(filename1, model);
    Logger<LL_model_BM> logger2(filename2, model);
    model.set_beta(1.0/temperatures[temperatures_count-1]);
    model.thermalize_BarkerWatts(cycles[temperatures_count-1]*5);
 //   for (int i=0; i<cycles[temperatures_count-1]*5; i++) model.swap_cycle();
    
    for (int i=temperatures_count-1; i>=0; i--){
        model.set_beta(1.0/temperatures[i]);
        model.thermalize_BarkerWatts(cycles[i]);
        model.cycle=0;
        for (int j=0; j<cycles[i]; j++){
            model.BarkerWatts_cycle();
   //         model.swap_cycle();
            model.calculate_P2();
            model.H();
            model.calculate_polar_order();
            model.count_neighbors_of_same_kind();
            model.cluster_count_and_size();
            logger.log_params();

        }
        logger.output_lattice_configuration();
        if (sim_id==0) std::cout << "cooldown with swap: " << i << "/" << std::to_string(temperatures_count) << "\r" << std::flush;
    }
    for (int i=0; i<temperatures_count; i++){

        model.set_beta(1.0/temperatures[i]);
        model.thermalize_BarkerWatts(cycles[i]);
        model.cycle=0;
        for (int j=0; j<cycles[i]; j++){
            model.BarkerWatts_cycle();
 //           model.swap_cycle();
            model.calculate_P2();
            model.H();
            model.calculate_polar_order();
            model.count_neighbors_of_same_kind();
            model.cluster_count_and_size();
            logger2.log_params();

        }
        logger2.output_lattice_configuration();
        if (sim_id==0) std::cout << "heatup with swap: " << i << "/" << std::to_string(temperatures_count) << "\r" << std::flush;
    }
}




void run(){
  
  /*
    std::thread t1(cooldown, 0.0, "c00.txt",0);
    std::thread t2(cooldown_swap, 0.0, "c00_swap.txt", 1);
    std::thread t3(heatup, 0.0, "h00.txt", 1);
    std::thread t4(heatup_swap, 0.0, "h00_swap.txt", 1);

    t1.join();
    t2.join();
    t3.join();
    t4.join();

    */
    double temperatures[temperatures_count] = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1 , 0.11,
       0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2 , 0.21, 0.22,
       0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3 };
    
    int cycles[temperatures_count] = {500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500,
       500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500,
       500, 500, 500, 500};


    

    std::thread t5(cooldown_then_heatup, 0.216, temperatures, cycles, "c02.txt", "h02.txt", 1);
    std::thread t6(cooldown_then_heatup_swap, 0.216, temperatures, cycles, "c02_swap.txt", "h02_swap.txt", 0);
    //std::thread t7(heatup, 0.216, temperatures, cycles, "h02.txt", 1);
    //std::thread t8(heatup_swap, 0.216, temperatures, cycles, "h02_swap.txt", 1);

    t5.join();
    t6.join();
    //t7.join();
    //t8.join();
/*

    std::thread t9(cooldown, 1.0, "c10.txt",0);
    std::thread t10(cooldown_swap, 1.0, "c10_swap.txt", 1);
    std::thread t11(heatup, 1.0, "h10.txt", 1);
    std::thread t12(heatup_swap, 1.0, "h10_swap.txt", 1);

    t9.join();
    t10.join();
    t11.join();
    t12.join();
*/
}


