
#include<vector>
#include<math.h>
#include <assert.h>
#include<numeric>
#include<iostream>
#include<random>
#include<functional>
#include "importanceSamplingAndVectorTools.hpp"
#include "stratifiedMCClass.hpp"
//implementations of AdaptiveStratifiedMC Class 

// Constructor
template <typename VAR>  AdaptiveStratifiedMC<VAR>::AdaptiveStratifiedMC( const int &I,const vector<int> &N,const vector<double> &p,const function<double(VAR)> &payoff,const function<VAR(default_random_engine&,uniform_real_distribution<double>&,double,double)> &generateRandomAlgo,const char method):I_(I),N(N),p(p),payoff(payoff),generateRandomAlgo(generateRandomAlgo),method(method)
{
            assert(I>0);
            for(int i=0;i<this->I_;++i){
                this->sigma.push_back(1.0);
                this->mi.push_back(0.0);
                this->Mi.push_back(0.0);
                vector<VAR> xi;
                this->X.push_back(xi);

            }
            assert((this->method=='A')|(this->method=='B'));
            assert(this->I_>0);
            assert(this->N[0]==0);
            assert(this->X.size()==this->I_);
            assert(this->p.size()==this->I_);
            assert(this->sigma.size()==this->I_);
            assert(this->mi.size()==this->I_);
            assert(this->Mi.size()==this->I_);
        }

// adaptive calculation for mi by method a
template <typename VAR> void AdaptiveStratifiedMC<VAR>::getNewmi_double(int delta){
    
    const vector<double> &sigma=this->get_sigma();
    double sum=0.0;
    const int I_=p.size();
    assert (I_==sigma.size());
    for (int i=0;i<I_;++i){
        sum+=p[i]*sigma[i];
        mi[i]=p[i]*sigma[i];
    }
    for (int i=0;i<I_;++i){
        assert(mi[i]>=0);
        mi[i]/=sum;
        mi[i]*=delta;
    }
    return;
    
}


//upDate Mi= Number of samples to be generated in each strata for the next step
template <typename VAR> 
void AdaptiveStratifiedMC<VAR>::upDateMi(){

    assert(I_==Mi.size());
    const vector<double>& mi=this->get_mi();

    double Scurr=0.0;
    double Spast=0.0;
    for(int i=0;i<I_;++i){
        Scurr+=mi[i];
        Mi[i]=(int)(Scurr)-(int)(Spast)+1;
        Spast=Scurr;
    }
}
// Templates Definations!
// Calculate StandardDev for a given Strata with certain Payoff Function
//
template <typename VAR>
double AdaptiveStratifiedMC<VAR>::standardDev(const vector<VAR>& x)const{
    double m1=0.0;
    double m2=0.0;
    int N=x.size();
    assert(N>0);
    
    for (int i=0;i<N;++i){
        double fx=payoff(x[i]);
        m1=double(i)/(double(i)+1.0)*m1+fx/(double(i)+1.0);
        m2=double(i)/(double(i)+1.0)*m2+fx*fx/(double(i)+1.0);
       
    }

    return sqrt(m2-m1*m1);
    
}

//  adaptive calculation for mi by method b
template <typename VAR>
void AdaptiveStratifiedMC<VAR>::getNewmi_double_b(int delta){

    const vector<double> &sigma=this->get_sigma();
    const vector<vector<VAR> > &X=this->get_X();
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



// compute the SigmStar defined in the paper
template <typename VAR>
double AdaptiveStratifiedMC<VAR>::calculateSigmaStar()const{
    return dot(p,sigma);
}

// calculate the estimator for E(payoff(X)) based on all simulated data
template <typename VAR>
double AdaptiveStratifiedMC<VAR>::calculateExpectation()const{

    double m1=0.0;
    const int I_=X.size();
    assert((I_>0)&(I_==p.size()));
    for(int i=0;i<I_;++i){
        int Ni=X[i].size();
        assert(Ni>0);
        double s=0.0;
        for(int j=0;j<Ni;++j){
            s=double(j)/(1.0+double(j))*s+payoff(X[i][j])/(1.0+double(j));
            //newPayoff(xi,u,r,T,d,k,v,S0);
        }
        s=s*p[i];
        m1+=s;
    }
    return m1;
    
}

// Run MC with Adaptive Stratification 
template <typename VAR>
void AdaptiveStratifiedMC<VAR>::montecarloStratified(){

    cout<<"current updating method is:"<<method<<endl;
    assert(this->X.size()==this->I_);
    // uniform distribution generator
    random_device rd;
    default_random_engine generator(rd());
    uniform_real_distribution<double> distribution(0.0,1.0);
    // start iterations
    const double pi_=1.0/((double)I_);
    const int NumIteration= N.size();
    assert(NumIteration>0);
    cout<<"start iteration!"<<endl;

    for(int nItr=1;nItr<NumIteration;++nItr){

        if(method=='B'){
            getNewmi_double_b(N[nItr]-N[nItr-1]-I_);
        }
        else{
            getNewmi_double(N[nItr]-N[nItr-1]-I_);
        }
        upDateMi();
        cout<<"update finished, start generations for step:"<<nItr<<endl;

        for(int i=0;i<I_;++i){
            for(int j=0;j<Mi[i];++j){
                double ai=pi_*i;
                double bi=pi_*(i+1);
                X[i].push_back(generateRandomAlgo(generator,distribution,ai,bi));   
            }
        }

        for(int i=0;i<I_;++i){
            sigma[i]=standardDev(X[i]);
            assert(sigma[i]>=0);
        }
    }
}
template class AdaptiveStratifiedMC<double>;
template class AdaptiveStratifiedMC<vector<double> >;


