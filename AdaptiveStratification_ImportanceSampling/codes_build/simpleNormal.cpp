#include<vector>
#include<math.h>
#include <assert.h>
#include<numeric>
#include<iostream>
#include<random>
#include"asa241.hpp"
#include"gradientDescent.hpp"
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
VAR generateStratifiedGaussian(double a,double b){
    
    //generate Z
    double Z=distribution(generator);
    Z=transformToStratifiedNormal(a,b,Z);
    return Z;
}

// Functions for montecarlo 
double standardDev(const vector<VAR>& x){
    double m1=0.0;
    double m2=0.0;
    int N=x.size();
    assert(N>0);
    
    for (const VAR &item:x){
        double fx=item;
        m1+=fx;
        m2+=fx*fx;
    }
    m1/=N;
    m2/=N;
    return sqrt(m2-m1*m1);
    
}
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
//  adaptive calculation for mi by method b
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
double calculateExpectation(const vector<vector<VAR> >& X,const vector<double> & p){
    double m1=0.0;
    const int I_=X.size();
    assert((I_>0)&(I_==p.size()));
    for(int i=0;i<I_;++i){
        int Ni=X[i].size();
        assert(Ni>0);
        double s=0.0;
        for(const VAR &xi:X[i]){
            s+=xi;
        }
        s=s/Ni*p[i];
        m1+=s;
    }
    return m1;
    
}



double calculateSigmaStar(const vector<double> & p, const vector<double> & sigma){

    return dot(p,sigma);
    
}
void montecarlo(const vector<int> &N,const vector<double> p,vector<double> &sigma, vector<double> &mi, vector<int> &Mi, vector< vector<VAR> > &X,char method){
    assert((method=='A')|(method=='B'));

    const int I_=p.size();
    const double pi_=1.0/((double)I_);
    const int NumIteration= N.size();
    assert(NumIteration>0);
    for(int nItr=1;nItr<NumIteration;++nItr){
        if(method=='A'){
             getNewmi_double(p, sigma,N[nItr]-N[nItr-1]-I_,mi);
        }
        else{
            getNewmi_double_b(p,sigma,N[nItr]-N[nItr-1]-I_,mi,X);            
        }

        upDateMi( mi,Mi);
        for(int i=0;i<I_;++i){
            for(int j=0;j<Mi[i];++j){
                double ai=pi_*i;
                double bi=pi_*(i+1);
                X[i].push_back(generateStratifiedGaussian(ai,bi));   
            }
        }

        for(int i=0;i<I_;++i){
            sigma[i]=standardDev(X[i]);
            assert(sigma[i]>=0);
        }
    }
}


int testMCNormal(char method){
    assert((method=='A')|(method=='B'));
    const double pi_ = 1.0/I;
    std::vector<double>p(I,pi_);
    vector< vector<VAR> > X(I);
    std::vector<int> N={0,300,1300,11300,31300};
    assert(N[0]==0);
    vector<double> sigma(I,1.0);
    vector<double> mi(I,0.0);
    vector<int> Mi(I,0);
   
    cout<<"start MonteCarlo by Mi update method:"<<method<<":"<<endl;
    montecarlo(N,p,sigma,mi,Mi,X,method);
    double sigmaStar=calculateSigmaStar(p,sigma);
    double price=calculateExpectation(X,p);


    cout<<"price:"<<price<<endl;
    cout<<"sigmaStar:"<<sigmaStar<<endl;
    cout<<"variance Star:"<<(sigmaStar*sigmaStar)<<endl;
    // show stratification!
    for(int i=0;i<I;++i){
        double q=(double)(X[i].size())/N[N.size()-1];
        cout<<"q"<<i<<" sample ratio:"<<q<<endl;
    }

    
    return 0;


}

int main(){
        
    
    testMCNormal('A');
    testMCNormal('B');
    return 0;


}






