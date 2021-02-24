#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <numeric>
#include "monte-carlo.hpp"
#include "var_alea.hpp"
#include "fcts.hpp"

#include <tvmet/Matrix.h>
#include <tvmet/Vector.h>
using namespace tvmet;

void affiche(std::string s, vect v) {
	std::cout << setiosflags(std::ios_base::fixed)
	<< s << ": \t" << v[0] << "\t" << v[1] << "\t" << v[2] << std::endl;
};

template <int d>
struct gauss_md : public var_alea< Vector<double, d> >{
	gauss_md(Matrix<double, d, d> L) : rac_covar(L), G(0,1) {};
	Vector<double, d> operator()() {
		Vector<double, d> x, y;
		std::generate(x.begin(), x.end(), G);
		y = rac_covar * x;
		return y;
	};
	private:
		gaussian G;
		Matrix<double, d, d> rac_covar;
};

template <int d>
struct black_scholes_md 
	: public std::unary_function< Vector<double,d>, Vector<double,d> > {
	black_scholes_md(const Vector<double, d> &x0, double r, 
				     const Vector<double, d> &s, double T)
		: x0(x0), mu(r-0.5*s*s), sigma(s), T(T) {};
	Vector<double, d> operator()(Vector<double, d> g) const {
		Vector<double, d> result(x0*exp(mu*T + sigma*sqrt(T)*g));
		return result;
	};
	private:
		Vector<double, d> x0, mu, sigma;
		double T;
};

template <int d>
struct basket : public std::unary_function< Vector<double, d>, double> {
	basket(const Vector<double, d> &alpha, double K)
		: alpha(alpha), K(K) {};
	double operator()(const Vector<double, d> &x) const {
		double y = dot(x, alpha);
		return (y > K) ? (y - K) : 0;
	};
	private:
		Vector<double, d> alpha;
		double K;
};

template <int d>
struct basket_vc : public std::unary_function< Vector<double, d>, double> {
	basket_vc(const Vector<double, d> &alpha, double K, double esp)
		: alpha(alpha), K(K), esp(esp) {};
	double operator()(const Vector<double, d> &x) const {
		double y = exp(dot(alpha, log(x)));
		return (y > K) ? (y - K) - esp : - esp;
	};
	private:
		Vector<double, d> alpha;
		double K, esp;
};

template <int d>
struct transfo : public std::unary_function<Vector<double, d>, Vector<double, d> >{
	Vector<double, d> operator()(Vector<double, d> x) const {
		Vector<double, d> y(-x);
		return y;
	};
};

using namespace std;
int main() {
	init_alea();

	Matrix<double, 40, 40> L(0);
	gaussian G1(0, 1);
	for (int i = 0; i < 40; i++) {
		double norm = 0;
		for (int j = 0; j <= i; j++) {
			L(i, j) = G1(); 
			norm += G1.current() * G1.current(); 
		}
		for (int j = 0; j <= i; j++) {
			L(i, j) /= sqrt(norm); 
		}
	}
	Matrix<double, 40, 40> correl(L*trans(L));
	gauss_md<40> G(L);
	
	Vector<double, 40> x0(100);
	Vector<double, 40> sigma(0.2);
	Matrix<double, 40, 40> mat_sig(0);
	for (int i = 0; i < 40; i++) mat_sig(i,i) = sigma(i);
	double r = 0.040;
	double T = 1;
	black_scholes_md<40> bs(x0, r, sigma, T);
	
	Vector<double, 40> alpha(1./40.);
	double K = 100;
	basket<40> payoff(alpha, K);
	
	Matrix<double, 40, 40> covar(mat_sig*correl*mat_sig);
	double new_x0 = exp(dot(alpha, log(x0))); 
	double new_sigma = sqrt(dot(alpha, covar*alpha)*T);
	double new_r = r-0.5*dot(alpha, sigma*sigma)+0.5*new_sigma*new_sigma; 
	double prix_vc = call_black_scholes(new_x0, K, new_r, new_sigma, T);
	basket_vc<40> payoff_vc(alpha, K, exp(new_r*T)*prix_vc);

	vect res1 = monte_carlo(1e5, payoff, bs, G);
	affiche("Sans RV", res1);
	vect res2 = monte_carlo(1e5, var_control(payoff, payoff_vc), bs, G);
	affiche("VarCont", res2);
	vect res3 = monte_carlo(1e5/2, antithet(compo_f(payoff, bs), transfo<40>()), G);
	affiche("Antith", res3);
	vect res4 = monte_carlo(1e5/2, antithet(compo_f(
		var_control(payoff, payoff_vc), bs), transfo<40>()), G);
	affiche("Ant+VC", res4);

	return 0;
};
