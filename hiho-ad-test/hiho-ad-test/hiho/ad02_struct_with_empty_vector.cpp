#include "ad02_struct_with_empty_vector.hpp"

#include <cmath>
#include <vector>
#include <utility>

#include <iostream>
#include <iomanip>
#include <limits>

#include "ad_american.hpp"
#include "timer.hpp"

namespace {

	namespace math {
		using Term       = std::pair<double, double>;
		using Polynomial = std::vector<Term>;


		struct Number {
			double v;
			Polynomial polynomial;

			Number(double vv) : v{ vv } {}

			Number operator-() const { return Number{ -v }; }
		};

		Number operator+(const Number& l, const Number& r) { return Number{ l.v + r.v }; }
		Number operator-(const Number& l, const Number& r) { return Number{ l.v - r.v }; }
		Number operator*(const Number& l, const Number& r) { return Number{ l.v * r.v }; }
		Number operator/(const Number& l, const Number& r) { return Number{ l.v / r.v }; }
		Number operator+(const Number& l, double r) { return Number{ l.v + r }; }
		Number operator-(const Number& l, double r) { return Number{ l.v - r }; }
		Number operator*(const Number& l, double r) { return Number{ l.v * r }; }
		Number operator/(const Number& l, double r) { return Number{ l.v / r }; }
		Number operator+(double l, const Number& r) { return Number{ l + r.v }; }
		Number operator-(double l, const Number& r) { return Number{ l - r.v }; }
		Number operator*(double l, const Number& r) { return Number{ l * r.v }; }
		Number operator/(double l, const Number& r) { return Number{ l / r.v }; }
		bool operator>(const Number& l, const Number& r) { return l.v > r.v; }
		Number exp(const Number& l) { return Number{ std::exp(l.v) }; }
		Number sqrt(const Number& l) { return Number{ std::sqrt(l.v) }; }
		Number pow(const Number& l, double r) { return Number{ std::pow(l.v, r) }; }

		using std::exp;
		using std::sqrt;
		using std::pow;
	}

	using Real = math::Number;

	inline Real putAmericanOption(const Real& s, const Real& sigma, const Real& k, const Real& r, const Real& t, int simulation) {

		auto dt = t / simulation;
		auto up = math::exp(sigma * math::sqrt(dt));

		auto p0 = (up - math::exp(-r * dt)) / (up * up - 1);
		auto p1 = math::exp(-r * dt) - p0;

		std::vector<Real> p;
		for (int i = 0; i != simulation; ++i) {
			auto pp = k - s * math::pow(up, 2.0 * i - simulation);
			pp = pp > 0.0 ? pp : 0.0;
			p.push_back(pp);
		}

		for (int j = simulation - 1; j != 0; --j) {
			for (int i = 0; i != j; ++i) {
				p[i] = p0 * p[i + 1] + p1 * p[i];    // binomial value
				auto exercise = k - s * math::pow(up, 2.0 * i - j);  // exercise value
				p[i] = p[i] > exercise ? p[i] : exercise;
			}
		}

		return p[0];
	}
}

void hiho::ad02_struct_with_empty_vector(double s, double sigma, double k, double r, double t, int simulation)
{
	auto func = [&]() { return putAmericanOption(s, sigma, k, r, t, simulation); };
	auto time = hiho::measureTime(func);
	auto value = func();

	auto diff = value.v - hiho::american(s, sigma, k, r, t, simulation);
	std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
	std::cout.setf(std::ios::left);
	std::cout << std::setw(30) << __func__ << " ( " << simulation << " )"
		<< ", diff : " << diff
		<< ", time : " << time << " msec"
		<< std::endl;
}
