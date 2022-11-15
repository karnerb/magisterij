#include "../elastomer/LC_Elastomer.h"
#include "../logger/logger.h"
#define size 50

void run(){
    double temperatures[7] = {1.135 , 1.1375, 1.14, 1.1425, 1.145 , 1.1475, 1.15};
    std::string filenames[7] = {"elastomer1.135.txt", "elastomer1.1375.txt", "elastomer1.1400.txt", "elastomer1.1425.txt", "elastomer1.1450.txt", "elastomer1.1475.txt", "elastomer1.1500.txt"};
    #pragma omp parallel for
    for (int s=0; s<7; s++){
        LC_Elastomer model(size);
        Logger<LC_Elastomer> logger(filenames[s], model);
        model.initialize_lattice_random();
        model.lambda=1.0;
        model.sigma=0.0;
        model.set_beta(1.0/temperatures[s]);
        int steps = 30;
        int cycles = 10;
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
                logger.log_params();
            }
        }
    }
}



    




