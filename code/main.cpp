#include <iostream>
#include <chrono>
#include "tests/performance.cpp"
//#include "tests/elastomer_basic_test.cpp"


int main(void){
    
    
    /*auto start = std::chrono::high_resolution_clock::now();
    run_test_r05();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::minutes>(end-start);
    std::cout << "time: " << duration.count() << "minutes\n";*/

    run();

    return 0;

}
