#include "Randomized_MC.hpp"
#include <math.h>       /* sqrt */ /* pow */
//template<typename TGen>
RQMC_Gaussian::RQMC_Gaussian(int I, int N, int dimension,std::function<double(std::vector<double>)> f,std::mt19937 gen):
  I(I),N(N),dimension(dimension),f(f),U(0,1),gen(gen),estimator(),half_CI(){};
void RQMC_Gaussian::operator ()(){
    std::vector<Gaussian_SH> Gn;
    for(int i=0; i<I; i++){
      std::vector<double> x;
      for(int j=0;j<((dimension%2==0) ? dimension : dimension+1);j++){
        x.push_back(U(gen));
      }
      Gn.push_back(Gaussian_SH(dimension,x));
    }
    std::vector<double> S(I); // For saving the sum of results for each Gn[i]
    for(int n=0; n< N; n++){
       for(int i=0;i<I; i++){
         double f_res = f(Gn[i]());
         S[i] += f_res;
       }
      double mu = std::accumulate(S.begin(), S.end(), 0.)/((double)(I*(n+1)));
      double sigma2 = 0.0;
      double beta3 = 0.00;

      for( std::vector<double>::iterator it=S.begin();it!=S.end(); it++){
        double s_tmp = std::abs(((*it)/((double)(n+1))-mu));
        sigma2+= s_tmp*s_tmp;
        beta3 += s_tmp*s_tmp*s_tmp;
      }
      sigma2/=(double)I;
      beta3/= (double)I;
      double half_l =1.96*sqrt(sigma2/(double)(I))+0.3354*(beta3/pow(sigma2,1.5)+0.415)*sqrt(sigma2)/((double)I);
      estimator.push_back(mu);
      half_CI.push_back(half_l);
    }
  };

  std::ostream& operator << (std::ostream& os, const RQMC_Gaussian rqmc)
  {
      int N = (rqmc.estimator).size();
      os <<"estimator,half_CI" << std::endl;
      for (int i=0; i<N;i++)
      {
          os <<rqmc.estimator[i]<<","<<rqmc.half_CI[i]<< std::endl;
      }
      return os;
  }
