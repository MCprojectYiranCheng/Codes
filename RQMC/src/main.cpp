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
    
     
    double x = 0.23;
    int dimension = 18;
    shifted_halton sh(x, dimension);
    std::cout<<"The shifted halton series of dimension "<< 8 <<" are as follows: "<< std::endl;
    for(int i=0; i<10; i++){
        std::vector<double> res =sh();  
        for (std::vector<double>::iterator j= res.begin(); j<res.end();j++){
            std::cout<<*j<<" ";
        }
        std::cout<<std::endl;
    }    
    
    
     
//      int k = 10;
//      halton p(2);
//      
//      for(int i=0; i<k; i++){
//          std::vector<double> res2 =p();  
//          std::cout<<res2[0]<<"\t"<<res2[1]<<std::endl;
//     }
//      
//      sobol sob(7);
//      std::cout<<"Sobol series are as follows: "<< std::endl;
//      for(int i=0; i<20; i++){
//        std::vector<double> res3 = sob();
//        for (std::vector<double>::iterator j= res3.begin(); j<res3.end();j++)
//          std::cout<<*j<<" ";
//        std::cout<<std::endl;
//        sob = sobol(sob);
//        
//      }
//      
     
    
}
