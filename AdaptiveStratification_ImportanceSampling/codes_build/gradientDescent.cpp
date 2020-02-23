#include <iostream>
#include<vector>
#include<math.h>
#include <assert.h>
#include<numeric>
#include "gradientDescent.hpp"
using namespace std;
vector<double> operator+(const vector<double>& a, const vector<double>& b){
    int N=a.size();
    assert(a.size()==b.size());
    vector<double> res;
    for(int i=0;i<N;++i){
        res.push_back(a[i]+b[i]);
    }
    return res;

}
vector<double> operator*(const vector<double>& a, const double& b){
    vector<double> res;
    int N=a.size();
    assert(N>0);
    for(int i=0;i<N;++i){
        res.push_back(a[i]*b);
    }
    return res;

}



vector<double>& operator+=(vector<double>& a, const vector<double>& b)
{
    assert((a.size()>0)&(a.size()==b.size()));
    transform(a.begin(), a.end(), b.begin(),
               a.begin(), plus<double>());
    return a;
};
vector<double>& operator-=(vector<double>& a, const vector<double>& b)
{
    assert((a.size()>0)&(a.size()==b.size()));
    transform(a.begin(), a.end(), b.begin(),
               a.begin(), minus<double>());
    return a;
};



vector<double>& operator+=(vector<double>& a, const double& k)
{   int N=a.size();
    assert(N>0);
    for(int i=0;i<N;++i){
        a[i]+=k;
    };
    return a;
};

vector<double>& operator*=(vector<double>& a, const double& k)
{   int N=a.size();
    assert((N>0));
    for(int i=0;i<N;++i){
        a[i]*=k;
    };
    return a;
};

vector<double> cumsum(const vector<double>& a)
{   vector<double> sum;
    sum.reserve(a.size());
    partial_sum(a.begin(), a.end(), back_inserter(sum), plus<double>());
    return sum;
};
double sum(const vector<double>& a)
{   double sum=0;
    int N=a.size();
    assert(N>0);
    for(int i=0;i<N;++i){
        sum+=a[i];
    }
    return sum;
};
double dot(const vector<double>& a,const vector<double>& b)
{   double sum=0;
    int N=a.size();
    assert((N>0)&(b.size()==N));
    for(int i=0;i<N;++i){
        sum+=a[i]*b[i];
    }
    return sum;
};





vector<double> calculateS(double S0,double r,double v,double T,int d,const vector<double> &X){
    vector<double> S;
    double delta=T/d;
    double a=r-0.5*v*v;
    double b=v*sqrt(delta);
    int N=X.size();
    assert(N>0); 
    double cumsumX=0.0;
    for(int m=0;m<N;++m){
        cumsumX+=X[m];
        S.push_back(S0*exp(a*(m+1)*delta+b*cumsumX));
    }
    return S;
};

vector<double> initialize(double S0,double k, double T, double v,int d, double r){

    double delta=T/d;
    double a=v*sqrt(delta);
    double b=(r-v*v*0.5)*delta;
    double c=log(1.5*k/S0);
    double C_i_1=0.0;
    vector<double>X;
    for(int i=0;i<d;++i){
        double temp=C_i_1;
        double Ci=(c-b*(i+1))/a;
        C_i_1=Ci;
        X.push_back(Ci-temp);
            }
    return X;
}

double g(double r,double T,int d,double k,const vector<double>& S){
    double a=sum(S);
    a=(a/d-k)*exp(-1.0*r*T);
    return a;
}

vector<double> gradient(double T,double r,double v,int d,const vector<double>& S){
    double a=v*sqrt(T/d)/d*exp(-1.0*r*T);
    vector<double> gradtemp;
    vector<double> grad;
    double s=0.0;
    int N=S.size();
    for(int i=N-1;i>=0;--i){
        s+=S[i];
        gradtemp.push_back(s*a);
    }
    
    for(int i=N-1;i>=0;--i){
        grad.push_back(gradtemp[i]);
    }
    return grad;
}

void updateOneStep(const double& gx,const vector<double>&grad,vector<double>& X){
    int N=X.size();
    assert(grad.size()==N);
    double graddotX=dot(grad,X);
    double squaredGrad=dot(grad,grad);
    double up=-1.0*(gx-graddotX)+sqrt((gx-graddotX)*(gx-graddotX)+4*squaredGrad);
    double down=2*squaredGrad;
    double beta=up/down;
    for(int i=0;i<N;++i){
        X[i]=beta*grad[i];
    }
    return;
}

vector<double> getOptimalDirection(double S0,double r,double v, double T, double k,int d, int numItr){
    vector<double> X=initialize(S0,k,T,v,d,r);
    for(int i=0;i<numItr;++i){
        auto S=calculateS(S0,r,v,T,d,X);
        double gx=g(r,T,d,k,S);
        auto grad=gradient(T,r,v,d,S);
        updateOneStep(gx,grad,X);
    }
    return X;
}
double newPayoff(const vector<double>& X, const vector<double>& u,double r, double T,int d,double k,double v,double S0){
    vector<double> Xshifted=X+u;
    vector<double> S =calculateS(S0,r,v,T,d,Xshifted);
    double gXshifted=g(r,T,d,k,S);
    if(gXshifted<=0){
        return 0;
    }
    double temp= dot(u,X)+0.5*dot(u,u);
    return gXshifted*exp(-1.0*temp);


}

int test() {
  
    double S0=50.,r=0.05,v=0.1,T=1.0,k=45.;
    int d=16;
    vector<double> X0=initialize(S0,k,T,v,d,r);
    

    int numItr=5;
    vector<double> u=getOptimalDirection(S0,r,v,T,k,d,numItr);
    cout << "u:" << endl; 
    for (const double &item:u){
        cout<<item<<' ';
    }
    cout<<endl;
    
    cout<<"newPayOff:"<<newPayoff(X0,u,r,T,d,k,v,S0);


    return 0;
}
