#include "ad07_expr_vec_tape_vec.hpp"

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

	struct CountUp {
		size_t& counter;
		CountUp(size_t& c) : counter(c) {}
		~CountUp() { ++counter; }
	};

	namespace math {

		struct Expression {
			using Term       = std::pair<ValueType, size_t>;
			using Polynomial = std::vector<Term>;

			static inline size_t counter = 0;
			static inline std::vector<Expression> expressions{};

			static size_t newExpression() { CountUp cp(counter); expressions.emplace_back(counter); return counter; }
			static size_t newExpression(ValueType  cof, size_t  exp) { CountUp cp(counter); expressions.emplace_back(counter, cof, exp); return counter; }
			static size_t newExpression(ValueType lcof, size_t lexp, ValueType rcof, size_t rexp) {
				CountUp cp(counter); 
				expressions.emplace_back(counter, lcof, lexp, rcof, rexp);
				return counter; 
			}
			static const Expression& getExpression(size_t index) {
				return expressions[index];
			}

			size_t index;
			Polynomial polynomial;

			Expression(size_t indx) : index(indx), polynomial{} {}
			Expression(size_t indx, 
				ValueType cof, size_t expr
			) : index(indx), polynomial{ {cof, expr} } {}
			Expression(size_t indx, 
				ValueType lcof, size_t lhs,
				ValueType rcof, size_t rhs
			) : index(indx), polynomial{ {lcof, lhs},{rcof, rhs} } {}

			ValueType d(const Expression& expr) const {
				if (index == expr.index) {
					return 1;
				}

				ValueType dx = 0;
				for (auto& term : polynomial) {
					dx += term.first * getExpression(term.second).d(expr);
				}
				return dx;
			}
		};

		struct Number {
			ValueType  v;
			size_t expression;

			Number() : v{ 0 }, expression{} {}
			Number(ValueType vv) : v{ vv }, expression{Expression::newExpression()}  {}
			Number(ValueType vv, size_t expr) : v{ vv }, expression{ expr }  {}

			ValueType d(const Number& x) const {
				return Expression::getExpression(expression).d(Expression::getExpression(x.expression));
			}

			Number operator-() const { return Number{ -v, Expression::newExpression(-1, expression) }; }
		};

		Number operator+(const Number& l, const Number& r) { return Number{ l.v + r.v, Expression::newExpression(1, l.expression, 1, r.expression) }; }
		Number operator-(const Number& l, const Number& r) { return Number{ l.v - r.v, Expression::newExpression(1, l.expression, -1, r.expression) }; }
		Number operator*(const Number& l, const Number& r) { return Number{ l.v * r.v, Expression::newExpression(r.v, l.expression, l.v, r.expression) }; }
		Number operator/(const Number& l, const Number& r) {
			auto ll = l.v;
			auto rr = r.v;
			return Number{ l.v / r.v, Expression::newExpression(1.0 / rr, l.expression, -ll / (rr * rr), r.expression) };
		}
		Number operator+(const Number& l, ValueType r) { return Number{ l.v + r, Expression::newExpression(1, l.expression) }; }
		Number operator-(const Number& l, ValueType r) { return Number{ l.v - r, Expression::newExpression(1, l.expression) }; }
		Number operator*(const Number& l, ValueType r) { return Number{ l.v * r, Expression::newExpression(r, l.expression) }; }
		Number operator/(const Number& l, ValueType r) { return Number{ l.v / r, Expression::newExpression(1.0 / r, l.expression) }; }
		Number operator+(ValueType l, const Number& r) { return Number{ l + r.v, Expression::newExpression(1, r.expression) }; }
		Number operator-(ValueType l, const Number& r) { return Number{ l - r.v, Expression::newExpression(-1, r.expression) }; }
		Number operator*(ValueType l, const Number& r) { return Number{ l * r.v, Expression::newExpression(l, r.expression) }; }
		Number operator/(ValueType l, const Number& r) { return Number{ l / r.v, Expression::newExpression(-l / (r.v * r.v), r.expression) }; }
		bool operator>(const Number& l, const Number& r) { return l.v > r.v; }
		Number exp(const Number& l) {
			auto ll = std::exp(l.v);
			return Number{ ll, Expression::newExpression(ll, l.expression) };
		}
		Number sqrt(const Number& l) {
			auto ll = std::sqrt(l.v);
			return Number{ ll, Expression::newExpression(0.5 / ll, l.expression) };
		}
		Number pow(const Number& l, ValueType r) {
			auto ll = std::pow(l.v, r);
			return Number{ ll, Expression::newExpression(r * ll / l.v, l.expression) };
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

void hiho::ad07_expr_vec_tape_vec(double s, double sigma, double k, double r, double t, int simulation)
{
	Real rs{ s };
	Real rsigma{ sigma };
	Real rr{ r };
	Real rt{ t };

	auto func = [&]() { return putAmericanOption(rs, rsigma, k, rr, rt, simulation); };
	auto timer = hiho::newTimer(func);
	auto time = timer.duration();
	auto& value = timer.value;

	auto diff = value.v - hiho::american(s, sigma, k, r, t, simulation);
	std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
	std::cout.setf(std::ios::left);
	std::cout << std::setw(30) << __func__ << " ( " << simulation << " )";
	std::cout.setf(std::ios::right);
	std::cout
		<< ", diff : " << diff
		<< ", time : " << std::setw(6) << time << " msec"
		<< ", greeks calculation is too late"
		//<< ", delta : " << value.d(rs)
		//<< ", vega : " << value.d(rsigma)
		//<< ", theta : " << value.d(rt)
		<< std::endl;

	math::Expression::counter = 0;
	math::Expression::expressions = {};
}
