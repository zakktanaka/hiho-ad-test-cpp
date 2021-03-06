#include "ad00_primitive_double.hpp"
#include "ad_american.hpp"

#include "timer.hpp"
#include "iohelper.hpp"

void hiho::ad00_primitive_double(double s, double sigma, double k, double r, double t, int simulation)
{
	auto func = [&]() { return american(s, sigma, k, r, t, simulation); };
	auto time = hiho::measureTime(func);
	auto value = func();

	auto delta = (hiho::american(s + 0.0001, sigma, k, r, t, simulation) - hiho::american(s - 0.0001, sigma, k, r, t, simulation)) / 0.0001 / 2;
	auto vega = (hiho::american(s, sigma + 0.0001, k, r, t, simulation) - hiho::american(s, sigma - 0.0001, k, r, t, simulation)) / 0.0001 / 2;
	auto theta = (hiho::american(s, sigma, k, r, t + 0.0001, simulation) - hiho::american(s, sigma, k, r, t - 0.0001, simulation)) / 0.0001 / 2;
	auto gamma = ((hiho::american(s + 0.0002, sigma, k, r, t, simulation) - hiho::american(s, sigma, k, r, t, simulation)) / 0.0001 / 2)
					- ((hiho::american(s, sigma, k, r, t, simulation) - hiho::american(s - 0.0002, sigma, k, r, t, simulation)) / 0.0001 / 2);

	HIHO_IO_MAX_LEN_DOUBLE_LSHOW;
	HIHO_IO_LEFT_COUT
		<< HIHO_IO_FUNC_WIDTH << __func__ << " ( " << simulation << " )";
	HIHO_IO_RIGHT_COUT
		<< ", " << HIHO_IO_VALUE(value)
		<< ", " << HIHO_IO_TIME(time)
		<< ", " << HIHO_IO_VALUE(delta)
		<< ", " << HIHO_IO_VALUE(vega)
		<< ", " << HIHO_IO_VALUE(theta)
		<< ", " << HIHO_IO_VALUE(gamma)
		<< std::endl;
}
