#include <random>
#include <iostream>

#include "low_discrepancy.hpp"

int main() {
    std::random_device rd;
    auto seed = rd();
    std::mt19937 gen(seed);
    std::uniform_real_distribution<> U(-1, 1);
    for (int n = 0; n < 10; ++n) {
        std::cout << U(gen) << std::endl;
    }
     std::cout << primes[100] << std::endl;
    
}
