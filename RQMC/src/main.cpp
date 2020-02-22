#include <random>
#include <iostream>

#include "low_discrepancy.hpp"

int main() {
    std::random_device rd;
    auto seed = rd();
    std::mt19937 gen(seed);
    std::uniform_real_distribution<> U(0, 1);
    for (int n = 0; n < 10; ++n) {
        std::cout << U(gen) << std::endl;
    }
     std::cout << primes[100] << std::endl;
    
     
     p_adic a(102,2);
     int res = a;
     std::cout<< res << std::endl; 
    
     int k = 100;
     halton p(2);
     for(int i=0; i<k; i++){
         std::vector<double> res2 =p();  
         std::cout<<res2[0]<<"\t"<<res2[1]<<std::endl;
    }
    
}
