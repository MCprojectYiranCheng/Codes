#ifndef __STRATIFIEDMCCLASS__
#define __STRATIFIEDMCCLASS__ 
#include<vector>
#include<math.h>
#include <assert.h>
#include<numeric>
#include<iostream>
#include<random>
#include<functional>
using namespace std;


// Functions and templates for montecarlo with adaptive Stratification 
// Here typename VAR represents the type of random sample: could be double or vector of double

template <typename VAR> class AdaptiveStratifiedMC{
    private:
        // number of strata
        const int I_=0;
        // number of samples in each Iteration
        const vector<int> N;
        // probability for each strata
        const vector<double> p;
        // payoff function
        const function<double(VAR)> payoff;
        // function for stratified gaussian vector generation
        const function<VAR(default_random_engine&,uniform_real_distribution<double>&,double,double)> generateRandomAlgo;
        // Mi updating method='A' or 'B'
        const char method;
        // sigma for each strata, mi
        vector<double> sigma,mi;
        vector<int> Mi;
        // generated random samples for each strata
        vector< vector<VAR> > X;

        // Member Functions:
        //
        //


        // Declarations: adaptive calculation for mi by method a
        void getNewmi_double(int delta);
        // Declarations: upDate Mi= Number of samples to be generated in each strata for the next step
        void upDateMi();
        // Templates Declarations!
        // Calculate StandardDev for a given Strata with certain Payoff Function
        double standardDev(const vector<VAR>& x)const;
        //  adaptive calculation for mi by method b
        void getNewmi_double_b(int delta);

    public:
        //constructor
        AdaptiveStratifiedMC(const int &I,const vector<int> &N,const vector<double> &p,const function<double(VAR)> &payoff,const function<VAR(default_random_engine&,uniform_real_distribution<double>&,double,double)> &generateRandomAlgo,const char method);

        // Data member access
        const vector<double>& get_sigma()const{
            return this->sigma;
        }
        const vector<double>& get_mi()const{
            return this->mi;
        }
        const vector<int>& get_Mi()const{
            return this->Mi;
        }
        const vector<vector<VAR> >& get_X()const{
            return this->X;
        }

        // Declarations: compute the SigmaStar defined in the paper
        double calculateSigmaStar()const;
        // calculate the estimator for E(payoff(X)) based on all simulated data
        double calculateExpectation()const;
        // Run MC with Adaptive Stratification 
        void montecarloStratified();

};



#endif
