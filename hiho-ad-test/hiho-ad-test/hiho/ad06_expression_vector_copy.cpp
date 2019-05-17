#include "ad06_expression_vector_copy.hpp"

#include <cmath>
#include <vector>
#include <utility>

#include <iostream>
#include <iomanip>
#include <limits>

#include "ad_american.hpp"
#include "timer.hpp"

namespace {

	size_t indexer = 0;
	using ValueType = double;

	namespace math {

		struct Expression {
			using Term       = std::pair<ValueType, Expression>;
			using Polynomial = std::vector<Term>;

			size_t     index;
			Polynomial polynomial;

			Expression() : index(indexer++), polynomial{} {}
			Expression(ValueType cof, const Expression& expr) :
				index(indexer++), polynomial{ {cof, expr} } {}
			Expression(ValueType lcof, const Expression& lhs,
				ValueType rcof, const Expression& rhs) :
				index(indexer++), polynomial{ {lcof, lhs},{rcof, rhs} } {}

			ValueType d(const Expression& expr) const {
				if (index == expr.index) {
					return 1;
				}

				ValueType dx = 0;
				for (auto& term : polynomial) {
					dx += term.first * term.second.d(expr);
				}
				return dx;
			}
		};



		struct Number {
			ValueType  v;
			Expression expression;

			Number() : v{ 0 }, expression{} {}
			Number(ValueType vv) : v{ vv }, expression{}  {}
			Number(ValueType vv, const Expression& expr) : v{ vv }, expression{ expr }  {}
			Number(ValueType vv, Expression&& expr) : v{ vv }, expression{ expr }  {}

			ValueType d(const Number& x) const {
				return expression.d(x.expression);
			}

			Number operator-() const { return Number{ -v }; }
		};

		Number operator+(const Number& l, const Number& r) { return Number{ l.v + r.v, Expression{1, l.expression, 1, r.expression} }; }
		Number operator-(const Number& l, const Number& r) { return Number{ l.v - r.v, Expression{1, l.expression, -1, r.expression} }; }
		Number operator*(const Number& l, const Number& r) { return Number{ l.v * r.v, Expression{r.v, l.expression, l.v, r.expression, } }; }
		Number operator/(const Number& l, const Number& r) {
			auto ll = l.v;
			auto rr = r.v;
			return Number{ l.v / r.v, Expression{1.0 / rr, l.expression, -ll / (rr * rr), r.expression, } };
		}
		Number operator+(const Number& l, ValueType r) { return Number{ l.v + r, Expression{1, l.expression, } }; }
		Number operator-(const Number& l, ValueType r) { return Number{ l.v - r, Expression{1, l.expression, } }; }
		Number operator*(const Number& l, ValueType r) { return Number{ l.v * r, Expression{r, l.expression, } }; }
		Number operator/(const Number& l, ValueType r) { return Number{ l.v / r, Expression{1.0 / r, l.expression, } }; }
		Number operator+(ValueType l, const Number& r) { return Number{ l + r.v, Expression{1, r.expression, } }; }
		Number operator-(ValueType l, const Number& r) { return Number{ l - r.v, Expression{-1, r.expression, } }; }
		Number operator*(ValueType l, const Number& r) { return Number{ l * r.v, Expression{l, r.expression, } }; }
		Number operator/(ValueType l, const Number& r) { return Number{ l / r.v, Expression{-l / (r.v * r.v), r.expression, } }; }
		bool operator>(const Number& l, const Number& r) { return l.v > r.v; }
		Number exp(const Number& l) {
			auto ll = std::exp(l.v);
			return Number{ ll, Expression{ll, l.expression, } };
		}
		Number sqrt(const Number& l) {
			auto ll = std::sqrt(l.v);
			return Number{ ll, Expression{0.5 / ll, l.expression, } };
		}
		Number pow(const Number& l, ValueType r) {
			auto ll = std::pow(l.v, r);
			return Number{ ll, Expression{r * ll / l.v, l.expression, } };
		}

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


void hiho::ad06_expression_vector_copy(double s, double sigma, double k, double r, double t, int simulation)
{
	Real rs{ s };
	Real rsigma{ sigma };
	Real rr{ r };
	Real rt{ t };

	auto timer = hiho::newTimer(
		[&]() { return putAmericanOption(rs, rsigma, k, rr, rt, simulation); }
	);

	auto diff = timer.value.v - hiho::american(s, sigma, k, r, t, simulation);
	std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
	std::cout.setf(std::ios::left);
	std::cout << std::setw(30)
		<< __func__ << " ( " << simulation << " )"
		<< ", diff : " << diff
		<< ", time : " << timer.duration() << " msec "
		<< ", delta : " << timer.value.d(rs)
		<< ", vega : " << timer.value.d(rsigma)
		<< ", theta : " << timer.value.d(rt)
		<< std::endl;
}
