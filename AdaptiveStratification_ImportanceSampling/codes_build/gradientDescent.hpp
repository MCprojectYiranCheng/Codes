#ifndef GRADIENTDESCENT
#define GRADIENTDESCENT

#include<vector>
#include<math.h>
#include <assert.h>
#include<numeric>
using namespace std;
// Basic Vector Operations
vector<double> operator+(const vector<double>& a, const vector<double>& b);
vector<double> operator*(const vector<double>& a, const double& b);

vector<double>& operator+=(vector<double>& a, const vector<double>& b);
vector<double>& operator-=(vector<double>& a, const vector<double>& b);
vector<double>& operator+=(vector<double>& a, const double& k);
vector<double>& operator*=(vector<double>& a, const double& k);
vector<double> cumsum(const vector<double>& a);
double sum(const vector<double>& a);
double dot(const vector<double>& a,const vector<double>& b);
// Calculate Asset Price S(X) on m timestamps based on X
vector<double> calculateS(double S0,double r,double v,double T,int d,const vector<double> &X);

// Initialize X0 such that g(X0)>0
vector<double> initialize(double S0,double k, double T, double v,int d, double r);

// calculate payoff function g(S) based on asset price S(X) on m timestamps
double g(double r,double T,int d,double k,const vector<double>& S);

// calculate gradient of g(S(X)) with respect to X
vector<double> gradient(double T,double r,double v,int d,const vector<double>& S);

// iteration methods for solving maxization problem of max(G(X)-0.5*X.dot(X), X in D) proposed in 
// ASYMPTOTICALLY OPTIMAL IMPORTANCE SAMPLING AND STRATIFICATION FOR PRICING PATH-DEPENDENT OPTIONS
void updateOneStep(const double& gx,const vector<double>&grad,vector<double>& X);

// solve the maximization problem by numItr of iterations (convergence could be achieved within 3 steps, so we suggest
// to take numItr=5)
vector<double> getOptimalDirection(double S0,double r,double v, double T, double k,int d, int numItr);

//new Payoff function f(x,u) in the new measure
double newPayoff(const vector<double>& X, const vector<double>& u,double r, double T,int d,double k,double v,double S0);




#endif
