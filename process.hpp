// process.hpp 
// Simulation (exacte) de processus aleatoires

#include <iostream>
#include <list>
#include <utility>
#include <var_alea.hpp>

// Classe abstraite processus 
template <typename T>
struct processus
{
  typedef std::pair<double, T> state;
  typedef std::list<state> result_type;
  typedef typename result_type::iterator iter;
  typedef typename result_type::const_iterator cst_iter;
  processus(int size = 0) : value(size) {};
  virtual result_type operator()() = 0;
  result_type current() const { return value; };
  template <typename S>
  friend std::ostream& operator<<(std::ostream &o,
  const processus<S> &p);
protected:
  result_type value;
};

// Surcharge de l'operateur << 
template <typename T>
std::ostream& operator <<(std::ostream &o, const processus<T> &p){
  typename processus<T>::cst_iter i;
  for(i = p.value.begin(); i != p.value.end(); ++i)
  o << (*i).first << "\t" << (*i).second << std::endl;
  return o;
}

// Brownian
struct brownian : public processus<double>{
  brownian(int n, double T=1):processus<double>(pow(2,n)+1), n(n), T(T), h(T/pow(2., n)), G(0,sqrt(h)) {};
  result_type operator()();
  result_type affine();
  friend struct black_scholes;
protected:
  int n;
  double h, T;
  gaussian G;
};


// Construction directe (exacte aux points discretisation)
result_type brownian::operator()(){
  value.clear();
  state val_k(0,0);
  value.push_back(val_k);
  do {
    val_k.first += h;
    val_k.second += G();
    value.push_back(val_k);
  } while (val_k.first < T);
  return value;
};
// Exercice: ecriture en 5 ligne sans la variable val_k

// Definition de la methode affine (Construction de Levy)
result_type brownian::affine() {
  n++; h *= 0.5;
  G = gaussian(0, sqrt(0.5*h));
  iter precedent = value.begin(), current = ++value.begin();
  while (current != value.end()) {
    value.insert(current,
    state(0.5*((*current).first+(*precedent).first),
    0.5*((*current).second+(*precedent).second)+G()));
    precedent = current;
    current++;
  }
  return value;
};


struct brownian_bridge : public brownian{
  brownian_bridge(int n, double T = 1, double B_T = 0): brownian(n, T), B_T(B_T) {};
  result_type operator()() {
    value.clear();
    value.push_back(state(0,0));
    value.push_back(state(T, B_T));
    int n_tmp = n; n = 0; h = T;
    for (int j = 0; j < n_tmp; j++) affine();
    n = n_tmp;
    return value;
  };
private:
  double B_T;
};
