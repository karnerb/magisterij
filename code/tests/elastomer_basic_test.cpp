#include <thread>
#include <chrono>
#include "../elastomer/LC_Elastomer.h"
#include "../logger/logger.h"
#include <cstdlib>
#include <unistd.h>
#define size 10

void stretch_and_compress(double temperature, std::string filename, int s){
    
    LC_Elastomer model(size);
    Logger<LC_Elastomer> logger(filename, model);

    model.initialize_lattice_random();
    model.lambda=1.0;
    model.sigma=0.0;
    model.set_beta(1.0/temperature);
    int steps = 30;
    int cycles = 10000;
    
    for (int i=0; i<cycles; i++){
        model.BarkerWatts_cycle();
    }
    for (int i=0; i<cycles; i++){
        model.BarkerWatts_cycle();
        if (i%50==0){
            model.adjust_resize_step();
        }
    }

    double sigma_max = 0.03;

    for (int i=0; i<steps; i++){        //swap MC
        bool swap_move();
        void swap_sites(int I, int J);
        void swap_cycle();

        for (int j=0; j<cycles; j++){
            model.BarkerWatts_cycle();
            if (j%10==0) model.adjust_resize_step();
        }
        if (s==0) std::cout << "resize step: " << model.resize_step << "\n";

        int accepted = 0;
        model.cycle=0;
        for (int j=0; j<cycles; j++){
            model.BarkerWatts_cycle();
            if (model.resize_move()) accepted++;
            if (j%5==0){
                model.H_neighbors();
                model.H_elastic();
                model.H_coupling();
                model.calculate_P2();
                logger.log_params();
            }
        }
        if (s==0) std::cout << "acc_rate: " << (double) accepted / cycles << "\n";
    }
}



void zero_stress_temperature_scan(std::string filename, double ratio, double angle, int s){
    int temperature_count = 50;
    double temperatures[temperature_count] = {1.6       , 1.56938776, 1.53877551, 1.50816327, 1.47755102,
       1.44693878, 1.41632653, 1.38571429, 1.35510204, 1.3244898 ,
       1.29387755, 1.26326531, 1.23265306, 1.20204082, 1.17142857,
       1.14081633, 1.11020408, 1.07959184, 1.04897959, 1.01836735,
       0.9877551 , 0.95714286, 0.92653061, 0.89591837, 0.86530612,
       0.83469388, 0.80408163, 0.77346939, 0.74285714, 0.7122449 ,
       0.68163265, 0.65102041, 0.62040816, 0.58979592, 0.55918367,
       0.52857143, 0.49795918, 0.46734694, 0.43673469, 0.40612245,
       0.3755102 , 0.34489796, 0.31428571, 0.28367347, 0.25306122,
       0.22244898, 0.19183673, 0.16122449, 0.13061224, 0.1};
    int cycles=10000;
    LC_Elastomer model(size);
    model.alpha = angle;
    Logger<LC_Elastomer> logger(filename, model);
    model.initialize_lattice_random();
    model.break_molecules((int) size*size*size*ratio);
    model.set_beta(1.0/temperatures[0]);
    model.sigma = 0.0;
    model.lambda=1.00;

    for (int i=0; i<cycles; i++){
        model.BarkerWatts_cycle();
        if (i%100==0){
            model.adjust_resize_step();
        }
    }
    for (int i=0; i<cycles; i++){
        model.BarkerWatts_cycle();
        model.resize_move();
    }

    
    for (int i=0; i<temperature_count; i++){        
        model.set_beta(1.0/temperatures[i]);
        for (int j=0; j<cycles; j++){
            model.BarkerWatts_cycle();
            if (j%10==0) model.adjust_resize_step();
        }

        if (s==0) std::cout << "i: " << i << " resize step: " << model.resize_step << "\n";

        int accepted = 0;
        model.cycle=0;
        for (int j=0; j<cycles; j++){
            model.BarkerWatts_cycle();
            if (model.resize_move()) accepted++;
            if (j%5==0){
                model.H_neighbors();
                model.H_elastic();
                model.H_coupling();
                model.calculate_P2();
                logger.log_params();
            }
        }
        if ((((double) accepted / cycles) < 0.25) || (((double) accepted / cycles) > 0.75)) std::cout << "need resize step adjustment!" << " " << s << "\n";
        if (s==0) std::cout << "acc_rate: " << (double) accepted / cycles << "\n";
    }
}


void speed(){
    LC_Elastomer model(size);
    model.initialize_lattice_random();
    model.set_beta(1.0/1.145);
    int cycles = 100;
    auto start = std::chrono::high_resolution_clock::now();
    for (int i=0; i<cycles; i++){
        model.BarkerWatts_cycle();
        model.resize_move();
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout << "100 MC cycle: " << duration.count() << "ms\n";
    std::cout << (double) 100 / duration.count() * 3600.0 << " cycles/h\n";

}

void run(){
    //speed();
    //speed();
    /*double temperatures[7] = {1.135 , 1.1375, 1.14, 1.1425, 1.145 , 1.1475, 1.15};
    std::string filenames[7] = {"elastomer1.135.txt", "elastomer1.1375.txt", "elastomer1.1400.txt", "elastomer1.1425.txt", "elastomer1.1450.txt", "elastomer1.1475.txt", "elastomer1.1500.txt"};
    
    std::thread t0(stretch, temperatures[0], filenames[0], 0);
    std::thread t1(stretch, temperatures[1], filenames[1], 1);
    std::thread t2(stretch, temperatures[2], filenames[2], 2);
    std::thread t3(stretch, temperatures[3], filenames[3], 3);
    std::thread t4(stretch, temperatures[4], filenames[4], 4);
    std::thread t5(stretch, temperatures[5], filenames[5], 5);
    std::thread t6(stretch, temperatures[6], filenames[6], 6);

    
    t0.join();
    t1.join();
    t2.join();
    t3.join();
    t4.join();
    t5.join();
    t6.join();*/

    std::string filenames[12] = {"temp_scan_r0.txt", "temp_scan_r001.txt", "temp_scan_r005.txt", "temp_scan_r01.txt",
                                "temp_scan_r05.txt", "temp_scan_r1.txt",
                                "temp_scan_r0_pravikot.txt", "temp_scan_r001_pravikot.txt", "temp_scan_r005_pravikot.txt", "temp_scan_r01_pravikot.txt",
                                "temp_scan_r05_pravikot.txt", "temp_scan_r1_pravikot.txt"};
    
    std::thread t0(zero_stress_temperature_scan,filenames[0], 0.0, pi/3.0, 0);
    sleep(1);
    std::thread t1(zero_stress_temperature_scan,filenames[1], 0.01, pi/3.0, 1);
    sleep(1);
    std::thread t2(zero_stress_temperature_scan,filenames[2], 0.05, pi/3.0, 2);
    sleep(1);
    std::thread t3(zero_stress_temperature_scan,filenames[3], 0.10, pi/3.0, 3);
    sleep(1);
    std::thread t4(zero_stress_temperature_scan,filenames[4], 0.5, pi/3.0, 4);
    sleep(1);
    std::thread t5(zero_stress_temperature_scan,filenames[5], 1.0, pi/3.0, 5);
    sleep(1);
    std::thread t6(zero_stress_temperature_scan,filenames[6], 0.0, pi/2.0, 6);
    sleep(1);
    std::thread t7(zero_stress_temperature_scan,filenames[7], 0.01, pi/2.0, 7);
    sleep(1);
    std::thread t8(zero_stress_temperature_scan,filenames[8], 0.05, pi/2.0, 8);
    sleep(1);
    std::thread t9(zero_stress_temperature_scan,filenames[9], 0.10, pi/2.0, 9);
    sleep(1);
    std::thread t10(zero_stress_temperature_scan,filenames[10], 0.5, pi/2.0, 10);
    sleep(1);
    std::thread t11(zero_stress_temperature_scan,filenames[11], 1.0, pi/2.0, 11);


    //
    t0.join();
    t1.join();
    t2.join();
    t3.join();
    t4.join();
    t5.join();
    t6.join();
    t7.join();
    t8.join();
    t9.join();
    t10.join();
    t11.join();
    //
}



    




