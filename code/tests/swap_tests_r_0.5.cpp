#include <thread>
#include "../llmodelBM/llmodelBM.h"
#include "../logger/logger.h"
#define size 30
#define temperatures_count 100

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
   /* double temperatures[temperatures_count] = {0.1       , 0.11919192, 0.13838384, 0.15757576, 0.17676768,
       0.1959596 , 0.21515152, 0.23434343, 0.25353535, 0.27272727,
       0.29191919, 0.31111111, 0.33030303, 0.34949495, 0.36868687,
       0.38787879, 0.40707071, 0.42626263, 0.44545455, 0.46464646,
       0.48383838, 0.5030303 , 0.52222222, 0.54141414, 0.56060606,
       0.57979798, 0.5989899 , 0.61818182, 0.63737374, 0.65656566,
       0.67575758, 0.69494949, 0.71414141, 0.73333333, 0.75252525,
       0.77171717, 0.79090909, 0.81010101, 0.82929293, 0.84848485,
       0.86767677, 0.88686869, 0.90606061, 0.92525253, 0.94444444,
       0.96363636, 0.98282828, 1.0020202 , 1.02121212, 1.04040404,
       1.05959596, 1.07878788, 1.0979798 , 1.11717172, 1.13636364,
       1.15555556, 1.17474747, 1.19393939, 1.21313131, 1.23232323,
       1.25151515, 1.27070707, 1.28989899, 1.30909091, 1.32828283,
       1.34747475, 1.36666667, 1.38585859, 1.40505051, 1.42424242,
       1.44343434, 1.46262626, 1.48181818, 1.5010101 , 1.52020202,
       1.53939394, 1.55858586, 1.57777778, 1.5969697 , 1.61616162,
       1.63535354, 1.65454545, 1.67373737, 1.69292929, 1.71212121,
       1.73131313, 1.75050505, 1.76969697, 1.78888889, 1.80808081,
       1.82727273, 1.84646465, 1.86565657, 1.88484848, 1.9040404 ,
       1.92323232, 1.94242424, 1.96161616, 1.98080808, 2.        };
    */
    double temperatures[temperatures_count] = {0.01      , 0.01393939, 0.01787879, 0.02181818, 0.02575758,
       0.02969697, 0.03363636, 0.03757576, 0.04151515, 0.04545455,
       0.04939394, 0.05333333, 0.05727273, 0.06121212, 0.06515152,
       0.06909091, 0.0730303 , 0.0769697 , 0.08090909, 0.08484848,
       0.08878788, 0.09272727, 0.09666667, 0.10060606, 0.10454545,
       0.10848485, 0.11242424, 0.11636364, 0.12030303, 0.12424242,
       0.12818182, 0.13212121, 0.13606061, 0.14      , 0.14393939,
       0.14787879, 0.15181818, 0.15575758, 0.15969697, 0.16363636,
       0.16757576, 0.17151515, 0.17545455, 0.17939394, 0.18333333,
       0.18727273, 0.19121212, 0.19515152, 0.19909091, 0.2030303 ,
       0.2069697 , 0.21090909, 0.21484848, 0.21878788, 0.22272727,
       0.22666667, 0.23060606, 0.23454545, 0.23848485, 0.24242424,
       0.24636364, 0.25030303, 0.25424242, 0.25818182, 0.26212121,
       0.26606061, 0.27      , 0.27393939, 0.27787879, 0.28181818,
       0.28575758, 0.28969697, 0.29363636, 0.29757576, 0.30151515,
       0.30545455, 0.30939394, 0.31333333, 0.31727273, 0.32121212,
       0.32515152, 0.32909091, 0.3330303 , 0.3369697 , 0.34090909,
       0.34484848, 0.34878788, 0.35272727, 0.35666667, 0.36060606,
       0.36454545, 0.36848485, 0.37242424, 0.37636364, 0.38030303,
       0.38424242, 0.38818182, 0.39212121, 0.39606061, 0.4       };

    int cycles[temperatures_count] = {1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000};


    

    std::thread t1(cooldown_then_heatup, 0.1, temperatures, cycles, "c01.txt", "h01.txt", 1);
    std::thread t2(cooldown_then_heatup_swap, 0.1, temperatures, cycles, "c01_swap.txt", "h01_swap.txt", 0);
    std::thread t3(cooldown_then_heatup, 0.3, temperatures, cycles, "c03.txt", "h03.txt", 1);
    std::thread t4(cooldown_then_heatup_swap, 0.3, temperatures, cycles, "c03_swap.txt", "h03_swap.txt", 0);
    std::thread t5(cooldown_then_heatup, 0.0, temperatures, cycles, "c0.txt", "h0.txt", 1);
    std::thread t6(cooldown_then_heatup_swap, 0.0, temperatures, cycles, "c0_swap.txt", "h0_swap.txt", 0);
    std::thread t7(cooldown_then_heatup, 0.2, temperatures, cycles, "c02.txt", "h02.txt", 1);
    std::thread t8(cooldown_then_heatup_swap, 0.2, temperatures, cycles, "c02_swap.txt", "h02_swap.txt", 0);
    std::thread t9(cooldown_then_heatup, 0.5, temperatures, cycles, "c05.txt", "h05.txt", 1);
    std::thread t10(cooldown_then_heatup_swap, 0.5, temperatures, cycles, "c05_swap.txt", "h05_swap.txt", 0);
    std::thread t11(cooldown_then_heatup, 1.0, temperatures, cycles, "c1.txt", "h1.txt", 1);
    std::thread t12(cooldown_then_heatup_swap, 1.0, temperatures, cycles, "c1_swap.txt", "h1_swap.txt", 0);



    t5.join();
    t6.join();
    t7.join();
    t8.join();
    t9.join();
    t10.join();
    t11.join();
    t12.join();
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


