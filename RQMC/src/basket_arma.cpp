#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <numeric>
#include "compose.hpp"
#include "mc.hpp"

#include <armadillo>

template <int d>
using mat = arma::mat::fixed<d, d>;

template <int d>
using vec = arma::vec::fixed<d>;

using vecMC = std::vector<double>; 
void affiche(std::string const & s, vecMC const & v) {
    //std::cout << setiosflags(std::ios_base::fixed)
    std::cout << s << ": \t" << v[0] << "\t" << v[1] << "\t" << v[2] << std::endl;
};

template <int d>
struct normal_dist_md {
    normal_dist_md() = default;
    normal_dist_md(mat<d> const & L) : rac_covar(L), G(0,1) {};
    template <typename TGen>
    vec<d> operator()(TGen & gen) {
        vec<d> x, y;
        for (auto & xk : x) { xk = G(gen); }
        y = rac_covar * x;
        return y;
    };
    private:
        mat<d> rac_covar;
        std::normal_distribution<double> G;
};

template <int d>
struct black_scholes_md {
    black_scholes_md() = default;
    black_scholes_md(vec<d> const & x0, double r,
                     vec<d> const & s, double T)
        : x0(x0), mu(r - 0.5 * s % s), sigma(s), T(T) {};
    vec<d> operator()(vec<d> const & g) const {
        vec<d> result = x0 % exp(mu * T + sqrt(T) * sigma % g);
        return result;
    };
    private:
        vec<d> x0, mu, sigma;
        double T;
};

template <int d>
struct basket {
    basket() = default;
    basket(vec<d> const & alpha, double K)
        : alpha(alpha), K(K) {};
    double operator()(vec<d> const & x) const {
        double y = dot(x, alpha);
        return (y > K) ? (y - K) : 0;
    };
    private:
        vec<d> alpha;
        double K;
};

template <int d>
struct basket_cv {
    basket_cv() = default;
    basket_cv(vec<d> const & alpha_x, vec<d> const & x0, double K, double I0, double price)
        : alpha_x(alpha_x), x0(x0), K(K), I0(I0), price(price) {};
    double operator()(vec<d> const & x) const {
        double y = I0 * exp(dot(alpha_x, log(x / x0)) / I0);
        return (y > K) ? ((y - K) - price) : - price;
    };
    private:
        vec<d> alpha_x, x0;
        double K, I0, price;
};

template <int d>
struct tneg {
    vec<d> operator()(vec<d> x) const { return -x; };
};

double cdf_normal(double x) {
    return 0.5*(1 + std::erf(x / std::sqrt(2)));
};

using namespace std;
int main() {
    random_device rd;
    auto seed = rd();
    mt19937_64 gen(seed);

    mat<40> correl;
    correl.fill(0.5);
    for (int i = 0; i < 40; ++i) correl(i, i) = 1;
    normal_dist_md<40> G(chol(correl, "lower"));
    
    vec<40> x0, sigma; 
    x0.fill(100);
    sigma.fill(0.3);
    double r = 0.1;
    double T = 1;
    black_scholes_md<40> bs(x0, r, sigma, T);
    
    vec<40> alpha;
    alpha.fill(1./40.);
    double K = 100;
    basket<40> payoff(alpha, K);
    
    vec<40> alpha_x = alpha % x0;
    vec<40> J = alpha_x % sigma;
    double I0 = sum(alpha_x);
    double mean_Z_T = (r - 0.5 * dot(alpha_x, sigma % sigma) / I0) * T;
    double var_Z_T = T / (I0 * I0) * dot(J, correl * J);
    double d = (mean_Z_T - log(K / I0)) / sqrt(var_Z_T);
    double price_cv = I0 * (exp(mean_Z_T + 0.5 * var_Z_T)
            * cdf_normal(d + sqrt(var_Z_T)) - (K / I0) * cdf_normal(d));
    basket_cv<40> payoff_cv(alpha_x, x0, K, I0, price_cv);
    
    cout << "Option Basket, n = 1e5 fixé" << endl;
    {
        // Monte Carlo avec n = 1e5 fixé (bof)
        auto Y = compose(payoff, bs, G);
        auto res1 = monte_carlo(Y, gen, 1e5);
        cout << "Plain MC:\t" << res1 << endl;

        auto phi_anti = antithetic(compose(payoff, bs), tneg<40>());
        auto Y_anti = compose(phi_anti, G);
        auto res2 = monte_carlo(Y_anti, gen, 1e5);
        cout << "Antithetic:\t" << res2 << endl;

        auto phi_cv = compose(control_variate(payoff, payoff_cv), bs);
        auto Y_cv = compose(phi_cv, G);
        auto res3 = monte_carlo(Y_cv, gen, 1e5);
        cout << "Control Var.:\t" << res3 << endl;
        
        auto phi_anticv = compose(control_variate(payoff, payoff_cv), bs);
        auto Y_anticv = compose(antithetic(phi_anticv, tneg<40>()), G);
        auto res4 = monte_carlo(Y_anticv, gen, 1e5);
        cout << "Ant + CV:\t" << res4 << endl;
    }        
    cout << endl << "Option Basket, précision epsilon = 0.01 fixée" << endl;
    {
        double epsilon = 0.01;
        // Monte Carlo avec précision epsilon fixée
        auto Y = compose(payoff, bs, G);
        auto res1 = monte_carlo(Y, gen, 1e3, epsilon);
        cout << "Plain MC:\t" << res1 << endl;

        auto phi_anti = antithetic(compose(payoff, bs), tneg<40>());
        auto Y_anti = compose(phi_anti, G);
        auto res2 = monte_carlo(Y_anti, gen, 1e3, epsilon);
        cout << "Antithetic:\t" << res2 << endl;

        auto phi_cv = compose(control_variate(payoff, payoff_cv), bs);
        auto Y_cv = compose(phi_cv, G);
        auto res3 = monte_carlo(Y_cv, gen, 1e3, epsilon);
        cout << "Control Var.:\t" << res3 << endl;
        
        auto phi_anticv = compose(control_variate(payoff, payoff_cv), bs);
        auto Y_anticv = compose(antithetic(phi_anticv, tneg<40>()), G);
        auto res4 = monte_carlo(Y_anticv, gen, 1e3, epsilon);
        cout << "Ant + CV:\t" << res4 << endl;
    }
    return 0;
};
