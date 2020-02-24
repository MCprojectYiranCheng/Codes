
#include<vector>
#include<math.h>
#include <assert.h>
#include<numeric>
#include<iostream>
#include<random>
#include<functional>
#include"gradientDescent.hpp"
#include "stratifiedMCTemplates.hpp"
// This file contains implementations for undefined functions in "stratifiedMCTemplates.hpp"!

// adaptive calculation for mi by method a
void getNewmi_double(const vector<double> &p, const vector<double> &sigma,int delta,vector<double> &m){
    double sum=0.0;
    const int I_=p.size();
    assert (I_==sigma.size());
    for (int i=0;i<I_;++i){
        sum+=p[i]*sigma[i];
        m[i]=p[i]*sigma[i];
    }
    for (int i=0;i<I_;++i){
        assert(m[i]>=0);
        m[i]/=sum;
        m[i]*=delta;
    }
    return;
    
}

//upDate Mi= Number of samples to be generated in each strata for the next step
void upDateMi(const vector<double>& mi,vector<int>& Mi){
    const int I_=mi.size();
    assert(I_==Mi.size());
    
    double Scurr=0.0;
    double Spast=0.0;
    for(int i=0;i<I_;++i){
        Scurr+=mi[i];
        Mi[i]=(int)(Scurr)-(int)(Spast)+1;
        Spast=Scurr;
    }
}

// compute the SigmStar defined in the paper
double calculateSigmaStar(const vector<double> & p, const vector<double> & sigma){
    return dot(p,sigma);
}

