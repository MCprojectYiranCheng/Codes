#include <random>
#include <iostream>
#include "low_discrepancy.hpp"
#include "Normal_distribution_LD.hpp"
#include <fstream>
#include <string>

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
    int dimension = 100;
    shifted_halton sh(x, dimension+1);
   // std::cout<<"The shifted halton series of dimension "<< 8 <<" are as follows: "<< std::endl;
    
    std::vector<double> mean(dimension);
    
    
    int N = 10000;
    
    
    double* xg = new double[dimension];
    std::fill_n(xg,dimension,0);

    std::ofstream myfile;
    myfile.open ("../data/example.csv");
    //myfile << "x1,x2\n";
    
    for(int i=0; i<dimension-1; i++){
        myfile <<"x"<<i<<",";
    }
    myfile <<"x"<<dimension-1<< "\n";
    
    for(int i=0; i<N; i++){
        Gaussian_lds gaussian(dimension, sh());
        std::vector<double> res = gaussian();  
        for (int j=0; j< dimension; j++){
             xg[j]+= res[j];
             myfile << res[j];
             if(j< dimension-1)myfile <<",";
             else myfile <<"\n";
         }
        //std::cout<<std::endl;
        //std::transform (res.begin(), res.end(), mean.begin(), mean.begin(), std::plus<double>());
     }   
    
     myfile.close();
    
     std::cout<<"The mean of the generated guassian: " <<std::endl;
     for (int i=0; i< dimension; i++){
             std::cout<<xg[i]/N<<" ";
         }
     std::cout<<std::endl;
     
     delete[] xg;
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
