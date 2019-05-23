#include "ad12_expr_vec_tape_vec_lazy.hpp"

#include <cmath>
#include <utility>

#include <iostream>
#include <iomanip>
#include <limits>

#include "ad_american.hpp"
#include "timer.hpp"

namespace {

	size_t indexer = 0;
	using ValueType = double;
	using Cache = std::unordered_map<size_t, ValueType>;

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
			static Expression& getExpression(size_t index) {
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

			ValueType d(const Expression& expr, Cache& cache) const {
				if (index == expr.index) {
					return 1;
				}

				auto it = cache.find(index);
				if (it != std::end(cache)) {
					return it->second;
				}

				ValueType dx = 0;
				for (auto& term : polynomial) {
					dx += term.first * getExpression(term.second).d(expr, cache);
				}
				cache[index] = dx;
				return dx;
			}

			void addTerm(const Term& term) {
				for (auto& tm : polynomial) {
					if (tm.second == term.second) {
						tm.first += term.first;
						return;
					}
				}
				polynomial.emplace_back(term);
			}
		};

		struct INumber {
			virtual ValueType  v() const = 0;
			virtual ~INumber() {}
			virtual void update(Expression&, ValueType) const = 0;
			
		};


		struct Number : public INumber {
			ValueType  v_;
			size_t expr_;

			Number(ValueType vv) : v_{ vv }, expr_{ Expression::newExpression() }  {}
			Number(ValueType vv, size_t expr) : v_{ vv }, expr_{ expr }  {}
			Number(const Number& other) : v_{ other.v_ }, expr_{ other.expr_ } {};
			Number(Number&& other) : v_{ other.v_ }, expr_{ other.expr_ } { };
			~Number() {}

			ValueType v() const override { return v_; }
			void update(Expression& updated, ValueType coef) const override {
				for (auto& term : expression().polynomial) {
					updated.addTerm({coef * term.first, term.second});
				}
			}

			ValueType d(const Number& x) const {
				Cache cache;
				return expression().d(x.expression(), cache);
			}

			Expression& expression() const { return Expression::getExpression(expr_); }

			Number operator-() const { return Number{ -v_, Expression::newExpression(-1, expr_) }; }
			Number& operator=(const Number& other) {
				this->v_ = other.v_;
				this->expr_ = other.expr_;
				return *this;
			}
		};

		Number operator+(const Number& l, const Number& r) { return Number{ l.v() + r.v(), Expression::newExpression(1, l.expr_, 1, r.expr_) }; }
		Number operator-(const Number& l, const Number& r) { return Number{ l.v() - r.v(), Expression::newExpression(1, l.expr_, -1, r.expr_) }; }
		Number operator*(const Number& l, const Number& r) { return Number{ l.v() * r.v(), Expression::newExpression(r.v_, l.expr_, l.v_, r.expr_) }; }
		Number operator/(const Number& l, const Number& r) {
			auto ll = l.v();
			auto rr = r.v();
			return Number{ l.v() / r.v(), Expression::newExpression(1.0 / rr, l.expr_, -ll / (rr * rr), r.expr_) };
		}
		Number operator+(const Number& l, ValueType r) { return Number{ l.v() + r, Expression::newExpression(1, l.expr_) }; }
		Number operator-(const Number& l, ValueType r) { return Number{ l.v() - r, Expression::newExpression(1, l.expr_) }; }
		Number operator*(const Number& l, ValueType r) { return Number{ l.v() * r, Expression::newExpression(r, l.expr_) }; }
		Number operator/(const Number& l, ValueType r) { return Number{ l.v() / r, Expression::newExpression(1.0 / r, l.expr_) }; }
		Number operator+(ValueType l, const Number& r) { return Number{ l + r.v(), Expression::newExpression(1, r.expr_) }; }
		Number operator-(ValueType l, const Number& r) { return Number{ l - r.v(), Expression::newExpression(-1, r.expr_) }; }
		Number operator*(ValueType l, const Number& r) { return Number{ l * r.v(), Expression::newExpression(l, r.expr_) }; }
		Number operator/(ValueType l, const Number& r) { return Number{ l / r.v(), Expression::newExpression(-l / (r.v() * r.v()), r.expr_) }; }
		bool operator>(const Number& l, const Number& r) { return l.v() > r.v(); }
		Number exp(const Number& l) {
			auto ll = std::exp(l.v());
			return Number{ ll, Expression::newExpression(ll, l.expr_) };
		}
		Number sqrt(const Number& l) {
			auto ll = std::sqrt(l.v_);
			return Number{ ll, Expression::newExpression(0.5 / ll, l.expr_) };
		}
		Number pow(const Number& l, ValueType r) {
			auto ll = std::pow(l.v_, r);
			return Number{ ll, Expression::newExpression(r * ll / l.v(), l.expr_) };
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

void hiho::ad12_expr_vec_tape_vec_lazy(double s, double sigma, double k, double r, double t, int simulation)
{
	{
		Real rs{ s };
		Real rsigma{ sigma };
		Real rr{ r };
		Real rt{ t };

		auto func = [&]() { return putAmericanOption(rs, rsigma, k, rr, rt, simulation); };
		auto timer = hiho::newTimer(func);
		auto time = timer.duration();
		auto& value = timer.value;

		auto diff = value.v() - hiho::american(s, sigma, k, r, t, simulation);
		std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
		std::cout.setf(std::ios::left);
		std::cout << std::setw(30)
			<< __func__ << " ( " << simulation << " )"
			<< ", diff : " << diff
			<< ", time : " << time << " msec "
			<< ", delta : " << value.d(rs)
			<< ", vega : " << value.d(rsigma)
			<< ", theta : " << value.d(rt)
			<< std::endl;

	}
	math::Expression::counter = 0;
	math::Expression::expressions = {};
}
