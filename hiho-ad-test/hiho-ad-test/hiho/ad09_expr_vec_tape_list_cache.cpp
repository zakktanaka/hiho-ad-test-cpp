#include "ad09_expr_vec_tape_list_cache.hpp"

#include <cmath>
#include <list>
#include <vector>
#include <utility>
#include <unordered_map>

#include <iostream>
#include <iomanip>
#include <limits>

#include "ad_american.hpp"
#include "timer.hpp"

namespace {

	using ValueType = double;

	namespace math {

		struct Expression;

		using Cache = std::unordered_map<const Expression*, ValueType>;

		struct Expression {
			using Term       = std::pair<ValueType, Expression*>;
			using Polynomial = std::vector<Term>;

			static inline std::list<Expression> expressions{};

			static Expression* newExpression() { expressions.emplace_back(); return &*std::rbegin(expressions); }
			static Expression* newExpression(ValueType  cof, Expression*  exp) { 
				expressions.emplace_back(cof, exp); 
				return &*std::rbegin(expressions); 
			}
			static Expression* newExpression(ValueType lcof, Expression* lexp, ValueType rcof, Expression* rexp) {
				expressions.emplace_back(lcof, lexp, rcof, rexp);
				return &*std::rbegin(expressions); ;
			}

			Polynomial polynomial;

			Expression() : polynomial{} {}
			Expression(ValueType cof, Expression* expr) : polynomial{ {cof, expr} } {}
			Expression(
				ValueType lcof, Expression* lhs,
				ValueType rcof, Expression* rhs
			) : polynomial{ {lcof, lhs},{rcof, rhs} } {}

			ValueType d(const Expression* expr, Cache& cache) const {
				if (expr == this) {
					return 1;
				}

				auto it = cache.find(this);
				if (it != std::end(cache)) {
					return it->second;
				}

				ValueType dx = 0;
				for (auto& term : polynomial) {
					dx += term.first * term.second->d(expr, cache);
				}
				cache.emplace(this, dx);
				return dx;
			}
		};

		struct Number {
			ValueType  v;
			Expression* expression;

			Number(ValueType vv) : v{ vv }, expression{ Expression::newExpression() }  {}
			Number(ValueType vv, Expression* expr) : v{ vv }, expression{ expr }  {}

			ValueType d(const Number& x) const {
				Cache cache;
				return expression->d(x.expression, cache);
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

void hiho::ad09_expr_vec_tape_list_cache(double s, double sigma, double k, double r, double t, int simulation)
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
		<< ", delta : " << value.d(rs)
		<< ", vega : " << value.d(rsigma)
		<< ", theta : " << value.d(rt)
		<< std::endl;

	math::Expression::expressions = {};
}
