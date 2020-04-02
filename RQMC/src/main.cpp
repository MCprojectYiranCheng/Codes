#include <random>
#include <iostream>
#include "low_discrepancy.hpp"
#include "Normal_distribution_LD.hpp"
#include <fstream>
#include <string>
#include "Randomized_MC.hpp"
#include <numeric>
#include <cassert>
#include <math.h>

template < class T >
std::ostream& operator << (std::ostream& os, const std::vector<T>& v)
{
    os << "[";
    for (typename std::vector<T>::const_iterator ii = v.begin(); ii != v.end(); ++ii)
    {
        os << " " << *ii;
    }
    os << "]";
    return os;
}

void to_csv(std::string dir, const RQMC_Gaussian rqmc){
  std::ofstream myfile;
  myfile.open (dir);
  myfile << rqmc <<"\n";
  myfile.close();
}

double f(std::vector<double> x){
  //std::cout << "input of f results is " << x << '\n';
  double res = std::accumulate(x.begin(), x.end(), 0.);
  //std::cout << "f results is " << res << '\n';
  return res;
}

double positive_part(double x){
  return (x>=0) ? x : 0;
}


double payoff_arithematic_Asian_option(std::vector<double> X, double T, int d, double s0, double K, double V, double r,bool is_call){
  double St,payoff=0.,dt;
  St= s0;
  dt = T/((double) d);
  assert(X.size()==d);
  for(int i=0; i<d;i++){
    St *= exp((r-0.5*V*V)*dt+V*sqrt(dt)*X[i]);
    payoff+= St;
  }

  if(is_call){
    payoff = exp(0.-r*T)*positive_part(payoff/((double)d)-K);
  }
  else{
    payoff = exp(0.-r*T)*positive_part(K-payoff/((double)d));
  }
  return payoff;
}


int main() {
  std::random_device rd;
  auto seed = rd();
  std::mt19937 gen(seed);

  std::vector<double> x;

  int N;
  std::cout << "Please input the N :" << '\n';
  std::cin >> N;
  // QMC_Gaussian rqmc_g(100,N,1,f,gen);
  // rqmc_g();
  //
  // std::cout << "RQMC results are: \n" << rqmc_g<<'\n';
  //
  // std::string out_dir = "../data/RQMC.csv";
  // to_csv(out_dir,rqmc_g);

  double T=1.0;
  //int d=16;
  double s0=50;
  //double K =45;
  double V=0.1;
  double r=0.05;
  bool is_call =true;

//  std::cout << "Please input the strik price K :" << '\n';
//  std::cin >> K;
  int d [] ={16,64};
  double K [] ={45,50,55};

  for(int i=0;i<2;i++){
    for(int j=0;j<3;j++){
        auto f2 = std::bind(payoff_arithematic_Asian_option,std::placeholders::_1,T,d[i],s0,K[j],V,r,is_call);
        RQMC_Gaussian pricing(100,N,d[i],f2,gen);
        pricing();
        std::cout << "Pricing results are: \n" << pricing <<'\n';
        std::string out_dir = "../data/Call_"+std::to_string(d[i])+"_"+std::to_string((int)K[j])+".csv";
        to_csv(out_dir,pricing);
    }
  }

  is_call = false;
    for(int i=0;i<2;i++){
      for(int j=0;j<3;j++){
          auto f2 = std::bind(payoff_arithematic_Asian_option,std::placeholders::_1,T,d[i],s0,K[j],V,r,is_call);
          RQMC_Gaussian pricing(100,N,d[i],f2,gen);
          pricing();
          std::cout << "Pricing results are: \n" << pricing <<'\n';
          std::string out_dir = "../data/put_"+std::to_string(d[i])+"_"+std::to_string((int)K[j])+".csv";
          to_csv(out_dir,pricing);
      }
    }

}
