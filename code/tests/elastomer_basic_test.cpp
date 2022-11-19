#include <thread>
#include <chrono>
#include "../elastomer/LC_Elastomer.h"
#include "../logger/logger.h"
#define size 50

void stretch(double temperature, std::string filename, int s){
    LC_Elastomer model(size);
    Logger<LC_Elastomer> logger(filename, model);
    model.initialize_lattice_random();
    model.lambda=1.0;
    model.sigma=0.0;
    model.set_beta(1.0/temperature);
    int steps = 30;
    int cycles = 100;
    for (int i=0; i<cycles; i++){
        model.BarkerWatts_cycle();
    }
    for (int i=0; i<cycles; i++){
        model.BarkerWatts_cycle();
        model.resize_move();
    }
    double sigma_max = 0.03;
    for (int i=0; i<steps; i++){
        if (s==0) std::cout << "i: " << i << "\r" << std::flush;
        model.sigma=i*sigma_max/steps;
        for (int j=0; j<cycles; j++){
            model.BarkerWatts_cycle();
            model.resize_move();
        }
        for (int j=0; j<cycles; j++){
            model.BarkerWatts_cycle();
            model.resize_move();
            if (j%5==0){
                model.H_neighbors();
                model.H_elastic();
                model.H_coupling();
                model.calculate_P2();
                logger.log_params();
            }
        }
    }
}

void acceptance_rate(){
    LC_Elastomer model(size);
    double temperature = 1.1450;
    model.initialize_lattice_random();
    model.lambda=1.0;
    model.sigma=0.0;
    model.set_beta(1.0/temperature);
    int steps = 30;
    int cycles = 100;
    for (int i=0; i<cycles; i++){
        model.BarkerWatts_cycle();
    }
    int accepted = 0;
    for (int i=0; i<cycles; i++){
        model.BarkerWatts_cycle();
        if (model.resize_move()) accepted++;
    }
    std::cout << "accepted moves during thermalization: " << accepted << "\n";
    double sigma_max = 0.03;
    for (int i=0; i<steps; i++){
        model.sigma=i*sigma_max/steps;
        for (int j=0; j<cycles; j++){
            model.BarkerWatts_cycle();
            model.resize_move();
        }
        accepted=0;
        for (int j=0; j<cycles; j++){
            model.BarkerWatts_cycle();
            if (model.resize_move()) accepted++;
        }
    std::cout << "accepted moves step " << i << "acc: " << accepted << "\n";
    }
}
void speed(){
    LC_Elastomer model(size);
    double temperature = 1.1450;
    int cycles = 100;
    auto start = std::chrono::high_resolution_clock::now();
    for (int i=0; i<cycles; i++){
        model.BarkerWatts_cycle();
        model.resize_move();
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout << "100 MC cycle: " << duration.count() << "ms\n";

}

void run(){
    speed();
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

}



    




