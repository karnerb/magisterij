#include <chrono>
#include <iostream>
#include "../llmodelBM/llmodelBM.h"

void performance_test(){
    //time individual functions
    LL_model_BM model(30);
    model.initialize_lattice_random();

    auto start = std::chrono::high_resolution_clock::now();
    for (int i=0; i<1000000; i++) model.dot(model.random_I(model.generator), 0, model.random_I(model.generator), 1);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout << "1000000 * dot product: " << duration.count() << "ms\n";

    start = std::chrono::high_resolution_clock::now();
    for (int i=0; i<1000000; i++) model.rotation3D(model.spins[model.random_I(model.generator)][0], 0.1,0.2,0.3,1.0);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout << "1000000 rotation3d time: " << duration.count() << "ms\n";

    start = std::chrono::high_resolution_clock::now();
    for (int i=0; i<1000; i++) model.H();
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout << "1000 H() time: " << duration.count() << "ms\n";

    start = std::chrono::high_resolution_clock::now();
    for (int i=0; i<10000; i++) model.calculate_P2();
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout << "10000 calculate_P2 time: " << duration.count() << "ms\n";


    start = std::chrono::high_resolution_clock::now();
    for (int i=0; i<100; i++) model.BarkerWatts_cycle();
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout << "100 BarkerWatts cycle time: " << duration.count() << "ms\n";

    start = std::chrono::high_resolution_clock::now();
    for (int i=0; i<1000; i++) model.shuffle_I();
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout << "1000 shuffleI cycle time: " << duration.count() << "ms\n";

    start = std::chrono::high_resolution_clock::now();
    for (int i=0; i<100; i++) model.swap_cycle();
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    std::cout << "100 swap_cycle time: " << duration.count() << "ms\n";




}