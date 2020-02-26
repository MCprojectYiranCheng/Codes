#include<vector>
#include<math.h>
#include <assert.h>
#include<numeric>
#include<iostream>
#include<random>
#include"asa241.hpp"
#include"importanceSamplingAndVectorTools.hpp"
#include"stratifiedMCClass.hpp"
#define I 10
using namespace std;
// type of random varaible
typedef double VAR;
// uniform distribution generator
default_random_engine generator;
uniform_real_distribution<double> distribution(0.0,1.0);
// transform uniform distribution to stratified gaussian
double transformToStratifiedNormal(double a,double b, double u){
    assert((b>a)&(a>=0)&(b<=1) );
    double v=a+u*(b-a);
    return r8_normal_01_cdf_inverse (v);
};
// generate stratified gaussian vector with respect to direction u  
VAR generateStratifiedGaussian(default_random_engine &generator,uniform_real_distribution<double> &distribution,double a,double b){
    
    //generate Z
    double Z=distribution(generator);
    Z=transformToStratifiedNormal(a,b,Z);
    return Z;
}

int testMCNormal(char method){
    // Initialize
    assert((method=='A')|(method=='B'));
    const double pi_ = 1.0/I;
    std::vector<double>p(I,pi_);
    std::vector<int> N={0,300,1300,11300,31300};
    assert(N[0]==0);
    // PAYOFF 
    const function<double(VAR)> payoff=[](VAR x){return x;};

    // Start MC
    cout<<"start MonteCarlo by Mi update method:"<<method<<":"<<endl;
    auto mcAS=AdaptiveStratifiedMC<VAR>(I,N,p,payoff,generateStratifiedGaussian,method);
    mcAS.montecarloStratified();
    
    // estimators:
    double sigmaStar=mcAS.calculateSigmaStar();
    double price=mcAS.calculateExpectation();


    cout<<"price:"<<price<<endl;
    cout<<"sigmaStar:"<<sigmaStar<<endl;
    cout<<"variance Star:"<<(sigmaStar*sigmaStar)<<endl;
    // show stratification!
    for(int i=0;i<I;++i){
        double q=(double)(mcAS.get_X()[i].size())/N[N.size()-1];
        cout<<"q"<<i<<" sample ratio:"<<q<<endl;
    }

    
    return 0;


}

int main(){
       
    testMCNormal('A');
    testMCNormal('B');
    return 0;


}






