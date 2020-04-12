#ifndef NORMAL_DISTRIBUTION_LD_HPP
#define NORMAL_DISTRIBUTION_LD_HPP
#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <list>
#include <vector>
#include <numeric>

#include "low_discrepancy.hpp"


using d_vect = std::vector<double>;
// The class to generate gaussian with low discrepancy series by using Box_Muller method
class Gaussian_lds{
  public:
      int dimension;
      d_vect lds ; // Low discrepancy sequence
      Gaussian_lds (){};
      Gaussian_lds (int dim, d_vect lds): dimension(dim), lds(lds) {};
      //copy constructor
      Gaussian_lds& operator= (const Gaussian_lds& x) {
          this->dimension = x.dimension;
          this->lds = x.lds;
          return *this;
    }

      d_vect Box_Muller(double u1, double u2){
          d_vect BM;
          double R = std::sqrt(-2.0* std::log(u1));
          double theta = 2.0*M_PI*u2;
          BM.push_back(R* std::cos(theta));
          BM.push_back( R* std::sin(theta));
          return BM;
      }
      d_vect operator ()(){
          int n = dimension%2==0 ? dimension: dimension+1 ;
          d_vect res;
          for ( int i =0; i<n; i+=2){
              d_vect BM = Box_Muller(lds[i], lds[i+1]);
              res.push_back(BM[0]);
              res.push_back(BM[1]);
          }
          if(n > dimension)  res.pop_back();
          return res;
    }
};

//Gaussian generator by the shifted halton sequence
class Gaussian_SH{
    typedef std::vector<double> result_type;
    int dimension;
    //double shifted_x;
    shifted_halton sh; // shifted halton sequence
    Gaussian_lds value;

  public:
    //Gaussian_SH(){}; //Defaut constructor
    Gaussian_SH(int dimension, result_type x): dimension(dimension), sh(x, (dimension%2==0) ? dimension : dimension+1){};
    result_type operator ()(){
      result_type out_res;
      value = Gaussian_lds(dimension,sh());
      return value();
    }
};

#endif
