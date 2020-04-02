#include<vector>
#include<math.h>
#include <assert.h>
#include<numeric>
#include<iostream>
#include<random>
#include<functional>
#include <fstream>

#include"asa241.hpp"
#include"importanceSamplingAndVectorTools.hpp"
#include"stratifiedMCClass.hpp"
#define I 100
using namespace std;
// type of random varaible
typedef vector<double> VAR;
typedef function<VAR(default_random_engine&,uniform_real_distribution<double>&,double,double)> StratifiedRandomGenerateAlgo;
// transform uniform distribution to stratified gaussian
double transformToStratifiedNormal(double a,double b, double u){
    assert((b>a)&(a>=0)&(b<=1) );
    double v=a+u*(b-a);
    return r8_normal_01_cdf_inverse (v);
};
// generate stratified gaussian vector with respect to direction u  
VAR generateStratifiedGaussianVector(default_random_engine& generator,uniform_real_distribution<double>& distribution,double a,double b,VAR u){
    int d=u.size();
    assert(d>0);
    double modU=sqrt(dot(u,u));
    // normalize
    u*=(1.0/modU);
    modU=dot(u,u);
    assert((modU-1.0<0.0000001)&(modU-1.0>-0.0000001));
    //generate Z
    double Z=distribution(generator);
    Z=transformToStratifiedNormal(a,b,Z);
    //generate Y
    vector<double> Y;
    for(int i=0;i<d;++i){
        double yi=transformToStratifiedNormal(0.0,1.0, distribution(generator));
        Y.push_back(yi);
    }
    vector<double> X=(u*Z);
    X+=Y;
    X-=(u*(dot(u,Y)));
    return X;
}


int testMC(const char& method,const bool CALL,const int& k,const int& d){
    // Initialize parameters 
    const double S0=50.,r=0.05,v=0.1,T=1.0;
        

    const double pi_ = 1.0/I;
    std::vector<double>p(I,pi_);
    vector< vector<VAR> > X(I);
    std::vector<int> N={0,100000,500000,1000000};
    assert(N[0]==0);
    // get optimal importance sampling shift
    VAR u=getOptimalDirection(S0,r,v,T,k,d,15,CALL);
    // cout << "Optimal Importance Sampling Shift u:" << endl; 
    // payoff function and Random generation Algo
    
    const function<double(VAR)>& payoff=bind(newPayoff,placeholders::_1,u,r,T,d,k,v,S0,CALL);
    const StratifiedRandomGenerateAlgo &generateRandomAlgo=bind(generateStratifiedGaussianVector,placeholders::_1,placeholders::_2,placeholders::_3,placeholders::_4,u);


    // for (const double &item:u){
    //     cout<<item<<' ';
    // }
    // cout<<endl;
    // MC
    auto asmc=AdaptiveStratifiedMC<VAR> (I,N,p,payoff,generateRandomAlgo,method);
    asmc.montecarloStratified();

    // estimators
    double sigmaStar=asmc.calculateSigmaStar();
    double price=asmc.calculateExpectation();
    
    if(CALL){
        cout<<"Call:";
    }
    else{
        cout<<"Put:";
    }
      std::ofstream myfile;
      if(CALL){
        myfile.open ("call_res.csv",std::ios::app);
      }
      else{
        myfile.open ("put_res.csv",std::ios::app);
      }
      myfile<<method<<","<<k<<","<<d<<","<<price<<","<<sigmaStar<<","<<(sigmaStar*sigmaStar)/N[N.size()-1]<<"\n";
      myfile.close();
    cout<<"__k:"<<k<<"__d:"<<d<<endl;
    cout<<"price estimator:"<<price<<endl;
    cout<<"sigmaStar:"<<sigmaStar<<endl;
    cout<<"variance of estimator:"<<(sigmaStar*sigmaStar)/N[N.size()-1]<<endl;

    
    return 0;


}
int main(){
std::vector<int> Ks{45,50,55};
std::vector<int> ds{16,64};
// Call option
for (auto &K : Ks){
    for (auto &d : ds){
        testMC('A',true,K,d);
        testMC('B',true,K,d);
    }

}
// Put option
for (auto &K : Ks){
    for (auto &d : ds){
        testMC('A',false,K,d);
        testMC('B',false,K,d);
    }

}
// testMC('A',false,55,16);
// testMC('B',false,55,16);
}





