#include <random>
#include <iostream>
#include "low_discrepancy.hpp"
#include "Normal_distribution_LD.hpp"
#include <fstream>

int main() {
    std::random_device rd;
    auto seed = rd();
    std::mt19937 gen(seed);
    std::uniform_real_distribution<> U(0, 1);
//     for (int n = 0; n < 10; ++n) {
//         std::cout << U(gen) << std::endl;
//     }
     //std::cout << primes[100] << std::endl;
    
     
    double x = 0;//U(gen);
    int dimension = 2;
    shifted_halton sh(x, dimension+1);
   // std::cout<<"The shifted halton series of dimension "<< 8 <<" are as follows: "<< std::endl;
    
    std::vector<double> mean(dimension);
    
    
    int N = 10000;
    
    std::ofstream myfile;
    myfile.open ("../data/example.csv");
    myfile << "x1,x2\n";
    for(int i=0; i<N; i++){
         
        Gaussian_lds gaussian(dimension, sh());
        std::vector<double> res = gaussian();  
        myfile <<res[0]<<"," <<res[1]<<"\n";
    /*    
        for (std::vector<double>::iterator j= res.begin(); j<res.end();j++){
             std::cout<<*j<<" ";
         }*/
        //std::cout<<std::endl;
        std::transform (res.begin(), res.end(), mean.begin(), mean.begin(), std::plus<double>());
        
    }   
    
    myfile.close();
    
     for (std::vector<double>::iterator j= mean.begin(); j<mean.end();j++){
             std::cout<<*j/N<<" ";
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
