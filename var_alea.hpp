#ifndef VAR_ALEA_HPP
#define VAR_ALEA_HPP
#include <cmath>
#include "mt19937.h"

void init_alea(unsigned seed = static_cast<unsigned>(std::time(0))) {
	init_genrand(seed);
};

template <typename T>
struct var_alea {
    typedef T result_type;
    var_alea() : value(0) {};
    var_alea(T value) : value(value) {};
    virtual ~var_alea() {}; 
    virtual T operator()() = 0;
    T current() const { return value; };
    protected: 
    	T value;
};

struct uniform : public var_alea<double>
{
    uniform(double left = 0, double right = 1)
        : left(left), size(right-left), genrand(genrand_real3) {};
    double operator()() {
        return value = left + size * genrand();
    };
    private: 
        double left, size;
        double (*genrand)(void);
};

struct expo : public var_alea<double>
{
	expo(double lambda) : inv_lambda(1./lambda), U(0,1) {};
	double operator()() {
		return value = - inv_lambda * log(U());
	};
	private:
		double inv_lambda;
		uniform U;
};

struct gaussian : public var_alea<double>
{
    gaussian(double mean = 0, double std = 1) 
        : mean(mean), std(std), flag(true), unif(-1,1) {};
    double operator()() {
        flag = !flag;
        if (!flag) {
            do {
                U = unif(); V = unif();
                R2 = U*U + V*V;
            } while (R2 > 1);
            rac = sqrt(-2 * log(R2) / R2);
            return value = mean + std * U * rac;
        } else
            return value = mean + std * V * rac;
    };
    private: 
        double mean, std, U, V, R2, rac;
        uniform unif;
        bool flag;
};

struct chi_deux : public var_alea<double>
{
	chi_deux(int n) : n(n), G(0,1) {};
	double operator()() {
		value = 0;
		for (int j = 0; j < n; j++) value += G()*G.current(); 
		return value;
	};
	private:
		int n;
		gaussian G;
};

struct inverse_gaussian : public var_alea<double>
{
	inverse_gaussian(double lambda, double mu) 
		: lambda(lambda), mu(mu), Y(1), U(0,1) {};
	double operator()() {
		double Z = mu + 0.5*mu*mu/lambda*Y();
		double rac = sqrt(Z*Z - mu*mu);
		return value = (U() < mu/(mu+Z+rac)) ? Z+rac : Z-rac; 
	};
	private: 
		double lambda, mu;
		chi_deux Y;
		uniform U;
};

struct normal_inverse_gaussian : public var_alea<double>
{
	normal_inverse_gaussian(double alpha, double beta, 
							double mu, double delta) 
		: alpha(alpha), beta(beta), mu(mu), delta(delta), G(0,1), 
		  Y(delta/sqrt(alpha*alpha-beta*beta), delta*delta)  {};
	double operator()() {
		double y_ = Y();
		return value = mu + beta*y_ * sqrt(y_) * G();
	};
	private: 
		double alpha, beta, mu, delta;
		gaussian G;
		inverse_gaussian Y;
};

#endif
