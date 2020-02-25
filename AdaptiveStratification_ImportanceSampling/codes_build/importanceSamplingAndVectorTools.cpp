#include <iostream>
#include<vector>
#include<math.h>
#include <assert.h>
#include<numeric>
#include "importanceSamplingAndVectorTools.hpp" 
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




// From standard gaussian X to asset price S
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


// Call/put option payoff function when in the money
double g(double S0,double r,double v,double T,int d,double k,const vector<double>& X,bool CALL){
    const vector<double>S=calculateS(S0,r,v,T,d,X);
    assert(S.size()==d);
    double a=sum(S);
    if(CALL){
        a=(a/d-k)*exp(-1.0*r*T);
    }
    else{
        a=(k-a/d)*exp(-1.0*r*T);
    }
    return a;
}


// Initialize a X0 making Call/Put Option  in the money
vector<double> initialize(double S0,double k, double T, double v,int d, double r, bool CALL){
    // initial value of Sm = ratio * K 
    double ratio=0.0;
    if(CALL){
        ratio=1.5;
    }
    else{
        ratio=0.5;
    }
    double delta=T/d;
    double a=v*sqrt(delta);
    double b=(r-v*v*0.5)*delta;
    double c=log(ratio*k/S0);
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

vector<double> gradient(double T,double r,double v,int d,const vector<double>& S, bool CALL){
    double a=v*sqrt(T/d)/d*exp(-1.0*r*T);
    if(!CALL){
        a*=-1.0;
    }
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

vector<double> getOptimalDirection(double S0,double r,double v, double T, double k,int d, int numItr, bool CALL){

    vector<double> X=initialize(S0,k,T,v,d,r,CALL);
    for(int i=0;i<numItr;++i){
        auto S=calculateS(S0,r,v,T,d,X);
        double gx=g(S0,r,v,T,d,k,X,CALL);
        auto grad=gradient(T,r,v,d,S,CALL);
        updateOneStep(gx,grad,X);
    }
    return X;
}
double newPayoff(const vector<double>& X, const vector<double>& u,double r, double T,int d,double k,double v,double S0,bool CALL){
    vector<double> Xshifted=X+u;
    double gXshifted=g(S0,r,v,T,d,k,Xshifted,CALL);
    if(gXshifted<=0){
        return 0;
    }
    double temp= dot(u,X)+0.5*dot(u,u);
    return gXshifted*exp(-1.0*temp);


}

int test() {
  
    double S0=50.,r=0.05,v=0.1,T=1.0,k=45.;
    int d=16;
    bool CALL=false;
    
    vector<double> X0=initialize(S0,k,T,v,d,r,CALL);

    cout << "X0:" <<endl;
    for (const double &item:X0){
        cout<<item<<' ';
    }
    cout<<endl;

    int numItr=10;
    vector<double> u=getOptimalDirection(S0,r,v,T,k,d,numItr,CALL);
    cout << "u:" << endl; 
    for (const double &item:u){
        cout<<item<<' ';
    }
    cout<<endl;
    cout<<"newPayOff:"<<newPayoff(X0,u,r,T,d,k,v,S0,CALL)<<endl;
    return 0;
}

