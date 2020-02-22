#include "low_discrepancy.hpp"

// constructeur a partir d'un entier
p_adic::p_adic(int n, unsigned p) : p(p) {
    int puiss = 1;
    while (n > 0) {
        ak.push_back(n % p);
        pk.push_back(puiss);
        puiss *= p;
        n -= ak.back();
        n /= p;
    }
    pk.push_back(puiss);
};

// constructeur a partir d'un reel de [0,1] !
p_adic::p_adic(double x, unsigned p) : p(p) {
    double intpart;
    double puiss = 1;
    x *= p;
    while (puiss*p < INT_MAX && modf(x, &intpart) != 0) {
        ak.push_back((int) intpart);
        pk.push_back(puiss);
        puiss *= p;
        x -= ak.back();
        x *= p;
    }
    ak.push_back((int) intpart);
    pk.push_back(puiss);
    pk.push_back(puiss*p);
};

// fonction d'increment 
void p_adic::increment() {
    coeff::iterator i = ak.begin();
    while ((i != ak.end()) && ((*i)+1 == p)) { (*i) = 0; i++; }
    if (i == ak.end()) {
        ak.push_back(1);
        pk.push_back(pk.back()*p);
    }
    else (*i) += 1;
};

