#include "llmodelBM.h"
#include "logger.h"
#include <iostream>
#include <chrono>
#include "swap_tests.cpp"

int main(void){
    
    
    auto start = std::chrono::high_resolution_clock::now();
    ceplak_with_swaps();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::minutes>(end-start);
    std::cout << "time: " << duration.count() << "minutes\n";

    return 0;

}