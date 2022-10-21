#include <thread>
#include "../llmodelBM/llmodelBM.h"
#include "../logger/logger.h"
#define size 10
#define kilo 10

void r05_swap_cooldown(double ratio, std::string filename, int sim_id){
    
    int temperatures_count = 52;

    double temperatures[temperatures_count] = {0.01      , 0.05882353, 0.10764706, 0.15647059, 0.20529412,
       0.25411765, 0.30294118, 0.35176471, 0.40058824, 0.44941176,
       0.49823529, 0.54705882, 0.59588235, 0.64470588, 0.69352941,
       0.74235294, 0.79117647, 0.84      , 0.88882353, 0.93764706,
       0.98647059, 1.03529412, 1.08411765, 1.13294118, 1.18176471,
       1.23058824, 1.27941176, 1.32823529, 1.37705882, 1.42588235,
       1.47470588, 1.52352941, 1.57235294, 1.62117647, 1.67      ,
       1.71882353, 1.76764706, 1.81647059, 1.86529412, 1.91411765,
       1.96294118, 2.01176471, 2.06058824, 2.10941176, 2.15823529,
       2.20705882, 2.25588235, 2.30470588, 2.35352941, 2.40235294,
       2.45117647, 2.5       };
    
    int cycles[temperatures_count] = {50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
       50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
       50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
       50};

    for (int i=0; i<temperatures_count; i++) cycles[i]=cycles[i]*kilo;


    LL_model_BM model_cooldown(size);
    model_cooldown.initialize_lattice_random();
    model_cooldown.break_molecules((int) (ratio*size*size*size));
    Logger<LL_model_BM> logger_cooldown(filename, model_cooldown);
    model_cooldown.set_beta(1.0/temperatures[temperatures_count-1]);
    model_cooldown.thermalize_BarkerWatts(cycles[temperatures_count-1]*5);

    for (int i=temperatures_count-1; i>=0; i--){
        model_cooldown.set_beta(1.0/temperatures[i]);
        model_cooldown.thermalize_BarkerWatts(cycles[i]);
        model_cooldown.cycle=0;
        for (int j=0; j<cycles[i]; j++){
            model_cooldown.BarkerWatts_cycle();
            model_cooldown.swap_cycle();
            model_cooldown.calculate_P2();
            model_cooldown.H();
            model_cooldown.calculate_polar_order();
            logger_cooldown.log_params();
            model_cooldown.cycle++;
        }
        if (sim_id==0) std::cout << "cooldown with swap: " << i << "/" << temperatures_count << "\r" << std::flush;
    }
    
}


void r05_cooldown(double ratio, std::string filename, int sim_id){
    
    int temperatures_count = 52;

    double temperatures[temperatures_count] = {0.01      , 0.05882353, 0.10764706, 0.15647059, 0.20529412,
       0.25411765, 0.30294118, 0.35176471, 0.40058824, 0.44941176,
       0.49823529, 0.54705882, 0.59588235, 0.64470588, 0.69352941,
       0.74235294, 0.79117647, 0.84      , 0.88882353, 0.93764706,
       0.98647059, 1.03529412, 1.08411765, 1.13294118, 1.18176471,
       1.23058824, 1.27941176, 1.32823529, 1.37705882, 1.42588235,
       1.47470588, 1.52352941, 1.57235294, 1.62117647, 1.67      ,
       1.71882353, 1.76764706, 1.81647059, 1.86529412, 1.91411765,
       1.96294118, 2.01176471, 2.06058824, 2.10941176, 2.15823529,
       2.20705882, 2.25588235, 2.30470588, 2.35352941, 2.40235294,
       2.45117647, 2.5        };
    
    int cycles[temperatures_count] = {50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
       50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
       50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
       50};

    for (int i=0; i<temperatures_count; i++) cycles[i]=cycles[i]*kilo;


    LL_model_BM model_cooldown(size);
    model_cooldown.initialize_lattice_random();
    model_cooldown.break_molecules((int) (ratio*size*size*size));
    Logger<LL_model_BM> logger_cooldown(filename, model_cooldown);
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
            model_cooldown.cycle++;
        }
        if (sim_id==0) std::cout << "cooldown without swap: " << i << "/" << temperatures_count << "\r" << std::flush;
    }
    
}


void r05_heatup(double ratio, std::string filename, int sim_id){
    
    int temperatures_count = 52;

    double temperatures[temperatures_count] = {0.01      , 0.05882353, 0.10764706, 0.15647059, 0.20529412,
       0.25411765, 0.30294118, 0.35176471, 0.40058824, 0.44941176,
       0.49823529, 0.54705882, 0.59588235, 0.64470588, 0.69352941,
       0.74235294, 0.79117647, 0.84      , 0.88882353, 0.93764706,
       0.98647059, 1.03529412, 1.08411765, 1.13294118, 1.18176471,
       1.23058824, 1.27941176, 1.32823529, 1.37705882, 1.42588235,
       1.47470588, 1.52352941, 1.57235294, 1.62117647, 1.67      ,
       1.71882353, 1.76764706, 1.81647059, 1.86529412, 1.91411765,
       1.96294118, 2.01176471, 2.06058824, 2.10941176, 2.15823529,
       2.20705882, 2.25588235, 2.30470588, 2.35352941, 2.40235294,
       2.45117647, 2.5         };
    
    int cycles[temperatures_count] = {50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
       50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
       50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
       50};

    for (int i=0; i<temperatures_count; i++) cycles[i]=cycles[i]*kilo;


    LL_model_BM model_cooldown(size);
    model_cooldown.initialize_lattice_parallel();
    model_cooldown.break_molecules((int) (ratio*size*size*size));
    Logger<LL_model_BM> logger_cooldown(filename, model_cooldown);
    model_cooldown.set_beta(1.0/temperatures[0]);
    model_cooldown.thermalize_BarkerWatts(cycles[0]*5);

    for (int i=0; i<temperatures_count; i++){
        model_cooldown.set_beta(1.0/temperatures[i]);
        model_cooldown.thermalize_BarkerWatts(cycles[i]);
        model_cooldown.cycle=0;
        for (int j=0; j<cycles[i]; j++){
            model_cooldown.BarkerWatts_cycle();
            model_cooldown.calculate_P2();
            model_cooldown.H();
            model_cooldown.calculate_polar_order();
            logger_cooldown.log_params();
            model_cooldown.cycle++;
        }
        if (sim_id==0) std::cout << "heatup without swap: " << i << "/" << temperatures_count << "\r" << std::flush;
    }
    
}


void r05_swap_heatup(double ratio, std::string filename, int sim_id){
    
    int temperatures_count = 52;

    double temperatures[temperatures_count] = {0.01      , 0.05882353, 0.10764706, 0.15647059, 0.20529412,
       0.25411765, 0.30294118, 0.35176471, 0.40058824, 0.44941176,
       0.49823529, 0.54705882, 0.59588235, 0.64470588, 0.69352941,
       0.74235294, 0.79117647, 0.84      , 0.88882353, 0.93764706,
       0.98647059, 1.03529412, 1.08411765, 1.13294118, 1.18176471,
       1.23058824, 1.27941176, 1.32823529, 1.37705882, 1.42588235,
       1.47470588, 1.52352941, 1.57235294, 1.62117647, 1.67      ,
       1.71882353, 1.76764706, 1.81647059, 1.86529412, 1.91411765,
       1.96294118, 2.01176471, 2.06058824, 2.10941176, 2.15823529,
       2.20705882, 2.25588235, 2.30470588, 2.35352941, 2.40235294,
       2.45117647, 2.5       };
    
    int cycles[temperatures_count] = {50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
       50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
       50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
       50};

    for (int i=0; i<temperatures_count; i++) cycles[i]=cycles[i]*kilo;


    LL_model_BM model_cooldown(size);
    model_cooldown.initialize_lattice_parallel();
    model_cooldown.break_molecules((int) (ratio*size*size*size));
    Logger<LL_model_BM> logger_cooldown(filename, model_cooldown);
    model_cooldown.set_beta(1.0/temperatures[0]);
    model_cooldown.thermalize_BarkerWatts(cycles[0]*5);

    for (int i=0; i<temperatures_count; i++){
        model_cooldown.set_beta(1.0/temperatures[i]);
        model_cooldown.thermalize_BarkerWatts(cycles[i]);
        model_cooldown.cycle=0;
        for (int j=0; j<cycles[i]; j++){
            model_cooldown.BarkerWatts_cycle();
            model_cooldown.swap_cycle();
            model_cooldown.calculate_P2();
            model_cooldown.H();
            model_cooldown.calculate_polar_order();
            logger_cooldown.log_params();
            model_cooldown.cycle++;
        }
        if (sim_id==0) std::cout << "heatup with swap: " << i << "/" << temperatures_count << "\r" << std::flush;
    }
    
}


void run_test_r05(){
    std::thread t1(r05_cooldown, 0.0, "test0c.txt",0);
    std::thread t2(r05_cooldown, 1.0, "test1c.txt", 1);
    std::thread t3(r05_heatup, 0.0, "test0h.txt", 1);
    std::thread t4(r05_heatup, 1.0, "test1h.txt", 1);
    t1.join();
    t2.join();
    t3.join();
    t4.join();
}


