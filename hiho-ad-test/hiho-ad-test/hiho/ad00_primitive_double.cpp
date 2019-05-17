#include "ad00_primitive_double.hpp"
#include "ad_american.hpp"

#include <iostream>
#include <iomanip>
#include <limits>
#include "timer.hpp"

void hiho::ad00_primitive_double(double s, double sigma, double k, double r, double t, int simulation)
{
	auto func = [&]() { return american(s, sigma, k, r, t, simulation); };
	auto time = hiho::measureTime(func);
	auto value = func();

	auto diff = value - hiho::american(s, sigma, k, r, t, simulation);
	std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
	std::cout.setf(std::ios::left);
	std::cout << __func__ << " : " << value
		<< ", time : " << time << " msec" 
		<< ", delta : " << (hiho::american(s + 0.0001, sigma, k, r, t, simulation) - hiho::american(s - 0.0001, sigma, k, r, t, simulation)) / 0.0001 / 2
		<< ", vega : " << (hiho::american(s, sigma + 0.0001, k, r, t, simulation) - hiho::american(s, sigma - 0.0001, k, r, t, simulation)) / 0.0001 / 2
		<< ", theta : " << (hiho::american(s, sigma, k, r, t + 0.0001, simulation) - hiho::american(s, sigma, k, r, t - 0.0001, simulation)) / 0.0001 / 2
		<< std::endl;
}
