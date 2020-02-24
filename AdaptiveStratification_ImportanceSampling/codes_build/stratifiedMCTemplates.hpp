#ifndef STRATIFIEDMCTEMPLATE
#define STRATIFIEDMCTEMPLATE

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


// Declarations: adaptive calculation for mi by method a
void getNewmi_double(const vector<double> &p, const vector<double> &sigma,int delta,vector<double> &m);

// Declarations: upDate Mi= Number of samples to be generated in each strata for the next step
void upDateMi(const vector<double>& mi,vector<int>& Mi);

// Declarations: compute the SigmStar defined in the paper
double calculateSigmaStar(const vector<double> & p, const vector<double> & sigma);






// Templates Definations!
// Calculate StandardDev for a given Strata with certain Payoff Function
template <typename VAR>
double standardDev(const vector<VAR>& x,const function<double(VAR) >& payoff){
    double m1=0.0;
    double m2=0.0;
    int N=x.size();
    assert(N>0);
    
    for (const VAR &item:x){
        double fx=payoff(item);
        //newPayoff(item,u,r,T,d,k,v,S0);
        m1+=fx;
        m2+=fx*fx;
    }
    m1/=N;
    m2/=N;
    return sqrt(m2-m1*m1);
    
}

//  adaptive calculation for mi by method b
template <typename VAR>
void getNewmi_double_b(const vector<double> &p, const vector<double> &sigma,int delta,vector<double> &mi,const vector<vector<VAR> > &X){

    const int I_=p.size();
    assert (I_==sigma.size());
    // define a_n
    vector<double> a_n;
    for (int i=0;i<I_;++i){
        a_n.push_back((X[i].size()+1)/(p[i]*sigma[i]));

    }

    //sorting and get odered index
    vector<int> newIndex;
    for(int i=0;i<I_;++i){
        newIndex.push_back(i);
    }
    sort( newIndex.begin(),newIndex.end(), [&](int i,int j){return a_n[i]<a_n[j];} );
    sort(a_n.begin(),a_n.end());

    // get Ik
    auto itrInf=find(a_n.begin(),a_n.end(),INFINITY);
    int Ik=distance(a_n.begin(), itrInf);
    assert((Ik>=0)&(Ik<=I_));

    // get b_n
    vector<double> b_n;
    vector<double> tempb;
    double s1=0.0,s2=0.0;

    for(int i=Ik-1;i>=0;--i){
        int i_origin=newIndex[i];
        tempb.push_back((delta+s1)/s2);
        s1+=X[i_origin].size();
        s2+=p[i_origin]*sigma[i_origin];
        
        
    }

    for(int i=Ik-1;i>=0;--i){
        b_n.push_back(tempb[i]);
    }
    // define b_n[-1]
    double b_n_=(delta+s1)/s2;
    // get I star
    int istar=Ik-1;
    while(istar>=0){
        if(a_n[istar]>=b_n[istar]){
            break;
        }
        --istar;    
    }

    assert(istar>=-1);

    // update mi
    // i>istar
    for(int i=istar+1;i<Ik;++i){
        int i_origin=newIndex[i];
        if(istar==-1){
            mi[i_origin]=p[i_origin]*sigma[i_origin]*b_n_-X[i_origin].size()-1;
        }
        else{
            mi[i_origin]=p[i_origin]*sigma[i_origin]*b_n[istar]-X[i_origin].size()-1;
        }

    }
    //i<=istar
    for(int i=0;i<=istar;++i){
        int i_origin=newIndex[i];
        mi[i_origin]=0;
    }
    // sigam[i]=0
    for(int i=Ik;i<I_;++i){
        int i_origin=newIndex[i];
        mi[i_origin]=0;
    }    
    
}

// calculate the estimator for E(payoff(X)) based on all simulated data
template <typename VAR>
double calculateExpectation(const vector<vector<VAR> >& X,const vector<double> & p,const function<double(VAR)>& payoff){
    double m1=0.0;
    const int I_=X.size();
    assert((I_>0)&(I_==p.size()));
    for(int i=0;i<I_;++i){
        int Ni=X[i].size();
        assert(Ni>0);
        double s=0.0;
        for(const VAR &xi:X[i]){
            s+=payoff(xi);
            //newPayoff(xi,u,r,T,d,k,v,S0);
        }
        s=s/Ni*p[i];
        m1+=s;
    }
    return m1;
    
}



// Run MC with Adaptive Stratification 
template <typename VAR>
void montecarloStratified(const vector<int> &N,const vector<double> p, vector<double> &sigma, vector<double> &mi, vector<int> &Mi, vector< vector<VAR> > &X,const function<double(VAR)>& payoff,const  function<VAR(default_random_engine&,uniform_real_distribution<double>&,double,double)> &generateRandomAlgo,const char &method){
    assert((method=='A')|(method=='B'));
    cout<<"current updating method is:"<<method<<endl;
    // uniform distribution generator
    default_random_engine generator;
    uniform_real_distribution<double> distribution(0.0,1.0);
    // start iterations
    const int I_=p.size();
    const double pi_=1.0/((double)I_);
    const int NumIteration= N.size();
    assert(NumIteration>0);
    for(int nItr=1;nItr<NumIteration;++nItr){

        if(method=='B'){
            getNewmi_double_b(p, sigma,N[nItr]-N[nItr-1]-I_,mi,X);
        }
        else{
            getNewmi_double(p, sigma,N[nItr]-N[nItr-1]-I_,mi);
        }
        upDateMi( mi,Mi);
        for(int i=0;i<I_;++i){
            for(int j=0;j<Mi[i];++j){
                double ai=pi_*i;
                double bi=pi_*(i+1);
                X[i].push_back(generateRandomAlgo(generator,distribution,ai,bi));   
            }
        }

        for(int i=0;i<I_;++i){
            sigma[i]=standardDev(X[i],payoff);
            assert(sigma[i]>=0);
        }
    }
}




#endif
