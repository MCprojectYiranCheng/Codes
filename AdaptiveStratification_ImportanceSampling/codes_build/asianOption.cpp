#include<vector>
#include<math.h>
#include <assert.h>
#include<numeric>
#include<iostream>
#include<random>
#include"asa241.hpp"
#include"gradientDescent.hpp"
#define I 100
using namespace std;
// type of random varaible
typedef vector<double> VAR;
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
VAR generateStratifiedGaussianVector(double a,double b,VAR u){
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

// Functions for montecarlo 
double standardDev(const vector<VAR>& x,const VAR& u,double r,double T,int d,double k,double v, double S0){
    double m1=0.0;
    double m2=0.0;
    int N=x.size();
    assert(N>0);
    
    for (const VAR &item:x){
        double fx=newPayoff(item,u,r,T,d,k,v,S0);
        m1+=fx;
        m2+=fx*fx;
    }
    m1/=N;
    m2/=N;
    return sqrt(m2-m1*m1);
    
}
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
double calculateExpectation(const vector<vector<VAR> >& X,const vector<double> & p,const VAR& u,double r,double T,int d,double k, double v, double S0 ){
    double m1=0.0;
    const int I_=X.size();
    assert((I_>0)&(I_==p.size()));
    for(int i=0;i<I_;++i){
        int Ni=X[i].size();
        assert(Ni>0);
        double s=0.0;
        for(const VAR &xi:X[i]){
            s+=newPayoff(xi,u,r,T,d,k,v,S0);
        }
        s=s/Ni*p[i];
        m1+=s;
    }
    return m1;
    
}



double calculateSigmaStar(const vector<double> & p, const vector<double> & sigma){

    return dot(p,sigma);
    
}
void montecarlo(const vector<int> &N,const vector<double> p,const VAR &u, double S0,double r, double v, double T, double k, int d, vector<double> &sigma, vector<double> &mi, vector<int> &Mi, vector< vector<VAR> > &X){
    const int I_=p.size();
    const double pi_=1.0/((double)I_);
    const int NumIteration= N.size();
    assert(NumIteration>0);
    for(int nItr=1;nItr<NumIteration;++nItr){
        getNewmi_double(p, sigma,N[nItr]-N[nItr-1]-I_,mi);
        upDateMi( mi,Mi);
        for(int i=0;i<I_;++i){
            for(int j=0;j<Mi[i];++j){
                double ai=pi_*i;
                double bi=pi_*(i+1);
                X[i].push_back(generateStratifiedGaussianVector(ai,bi,u));   
            }
        }

        for(int i=0;i<I_;++i){
            sigma[i]=standardDev(X[i],u,r,T,d,k,v,S0);
            assert(sigma[i]>=0);
        }
    }
}


int main(){
    double S0=50.,r=0.05,v=0.1,T=1.0,k=45.;
    int d=16;
    const double pi_ = 1.0/I;
    std::vector<double>p(I,pi_);
    vector< vector<VAR> > X(I);
    std::vector<int> N={0,100000,500000,1000000};
    assert(N[0]==0);
    vector<double> sigma(I,1.0);
    vector<double> mi(I,0.0);
    vector<int> Mi(I,0);
    VAR u=getOptimalDirection(S0,r,v,T,k,d,10);
    cout << "u:" << endl; 
    for (const double &item:u){
        cout<<item<<' ';
    }
    cout<<endl;


    montecarlo(N,p,u,S0,r,v,T,k,d,sigma,mi,Mi,X);
    double sigmaStar=calculateSigmaStar(p,sigma);
    double price=calculateExpectation(X,p,u,r,T,d,k,v,S0);


    cout<<"price:"<<price<<endl;
    cout<<"sigmaStar:"<<sigmaStar<<endl;
    cout<<"variance Star:"<<(sigmaStar*sigmaStar)<<endl;

    
    return 0;


}





