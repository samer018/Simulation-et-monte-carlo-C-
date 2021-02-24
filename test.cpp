// Test des classes de simulation var_alea.hpp
// Compilation: g++ -c mt19937.c
//              g++ mt19937.o test.cpp -o prog

#include <iostream>
#include "var_alea.hpp"

int main() {
	init_alea();

    std::cout << "Simulation d'une uniforme sur [0,1]" << std::endl;
	uniform U(0,1);
	std::cout << U() << std::endl;
	
	std::cout << "Simulation d'une exponentielle de parametre 1" << std::endl;
	expo E(1);
	std::cout << E() << std::endl;
	
	std::cout << "Simulation d'une gaussienne standard" << std::endl;
	gaussian G(0,1);
	std::cout << G() << std::endl;
	
	std::cout << "Simulation d'une chi deux a une dimension" << std::endl;
	chi_deux X(1);
	std::cout << X() << std::endl;
	
	std::cout << "Simulation d'une gausienne inverse (0.5,1)" << std::endl;
	inverse_gaussian Y(0.5,1);
	std::cout << Y() << std::endl;
	
	std::cout << "Simulation d'une gausienne inverse (0.8,0.1,0,2)" << std::endl;
	normal_inverse_gaussian Z(0.8,0.1,0,2);
	std::cout << Z() << std::endl;

	return 0;
};
