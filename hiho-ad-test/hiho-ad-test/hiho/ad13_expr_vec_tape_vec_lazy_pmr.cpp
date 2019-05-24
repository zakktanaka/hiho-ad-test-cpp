#include "ad13_expr_vec_tape_vec_lazy_pmr.hpp"

#include <cmath>
#include <utility>

#include <iostream>
#include <iomanip>
#include <limits>

#include <memory_resource>

#include "ad_american.hpp"
#include "timer.hpp"

namespace {

	namespace pmr = std::pmr;

	size_t indexer = 0;
	using ValueType = double;
	using Cache = pmr::unordered_map<size_t, ValueType>;

	struct CountUp {
		size_t& counter;
		CountUp(size_t& c) : counter(c) {}
		~CountUp() { ++counter; }
	};

	namespace math {

		struct Expression {
			using Term       = std::pair<ValueType, size_t>;
			using Polynomial = pmr::vector<Term>;

			static inline size_t counter = 0;
			static inline pmr::vector<Expression> expressions{};

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

		struct Unary : INumber {
			ValueType v_;
			ValueType coef_;
			const INumber& num_;

			Unary(ValueType vv, ValueType coef, const INumber& num) :
				v_(vv), coef_(coef), num_(num) {}

			ValueType v() const override { return v_; }
			void update(Expression& updated, ValueType coef) const override {
				num_.update(updated, coef * coef_);
			}
		};

		struct Binary : INumber {
			ValueType v_;
			ValueType lcoef_;
			ValueType rcoef_;
			const INumber& lnum_;
			const INumber& rnum_;

			Binary(
				ValueType vv,
				ValueType lcoef, const INumber& lnum,
				ValueType rcoef, const INumber& rnum) :
				v_(vv), lcoef_(lcoef), rcoef_(rcoef), lnum_(lnum), rnum_(rnum) {}

			ValueType v() const override { return v_; }
			void update(Expression& updated, ValueType coef) const override {
				lnum_.update(updated, coef * lcoef_);
				rnum_.update(updated, coef * rcoef_);
			}
		};

		struct Number : public INumber {
			ValueType  v_;
			size_t expr_;

			Number(ValueType vv) : v_{ vv }, expr_{ Expression::newExpression() }  {}
			Number(ValueType vv, size_t expr) : v_{ vv }, expr_{ expr }  {}
			Number(const INumber& other) : v_{ other.v() }, expr_{ Expression::newExpression() } { other.update(expression(), 1); };
			~Number() {}

			ValueType v() const override { return v_; }
			void update(Expression& updated, ValueType coef) const override {
				updated.addTerm({ coef, expr_ });
			}

			ValueType d(const Number& x) const {
				Cache cache;
				return expression().d(x.expression(), cache);
			}

			Expression& expression() const { return Expression::getExpression(expr_); }

			Number operator-() const { return Number{ -v_, Expression::newExpression(-1, expr_) }; }
			Number& operator=(const INumber& other) {
				Number newone{ other.v() };
				other.update(newone.expression(), 1);
				this->v_ = newone.v_;
				this->expr_ = newone.expr_;
				return *this;
			}
		};

		Binary operator+(const INumber& l, const INumber& r) { return Binary{ l.v() + r.v(), 1, l, 1, r }; }
		Binary operator-(const INumber& l, const INumber& r) { return Binary{ l.v() - r.v(), 1, l, -1, r }; }
		Binary operator*(const INumber& l, const INumber& r) { return Binary{ l.v() * r.v(), r.v(), l, l.v(), r }; }
		Binary operator/(const INumber& l, const INumber& r) {
			auto ll = l.v();
			auto rr = r.v();
			return Binary{ l.v() / r.v(), 1.0 / rr, l, -ll / (rr * rr), r };
		}
		Unary operator+(const INumber& l, ValueType r) { return Unary{ l.v() + r, 1, l }; }
		Unary operator-(const INumber& l, ValueType r) { return Unary{ l.v() - r, 1, l }; }
		Unary operator*(const INumber& l, ValueType r) { return Unary{ l.v() * r, r, l }; }
		Unary operator/(const INumber& l, ValueType r) { return Unary{ l.v() / r, 1.0 / r, l }; }
		Unary operator+(ValueType l, const INumber& r) { return Unary{ l + r.v(), 1, r }; }
		Unary operator-(ValueType l, const INumber& r) { return Unary{ l - r.v(), -1, r }; }
		Unary operator*(ValueType l, const INumber& r) { return Unary{ l * r.v(), l, r }; }
		Unary operator/(ValueType l, const INumber& r) { return Unary{ l / r.v(), -l / (r.v() * r.v()), r }; }
		bool operator>(const INumber& l, const INumber& r) { return l.v() > r.v(); }
		bool operator>(const INumber& l, ValueType r) { return l.v() > r; }
		Unary exp(const INumber& l) {
			auto ll = std::exp(l.v());
			return Unary{ ll, ll, l };
		}
		Unary sqrt(const INumber& l) {
			auto ll = std::sqrt(l.v());
			return Unary{ ll, 0.5 / ll, l };
		}
		Unary pow(const INumber& l, ValueType r) {
			auto ll = std::pow(l.v(), r);
			return Unary{ ll, r * ll / l.v(), l };
		}

		using std::exp;
		using std::sqrt;
		using std::pow;
	}

	using Real = math::Number;

	inline Real putAmericanOption(const Real& s, const Real& sigma, const Real& k, const Real& r, const Real& t, int simulation) {

		Real dt = t / simulation;
		Real up = math::exp(sigma * math::sqrt(dt));

		Real p0 = (up - math::exp(-r * dt)) / (up * up - 1);
		Real p1 = math::exp(-r * dt) - p0;

		std::vector<Real> p;
		for (int i = 0; i != simulation; ++i) {
			Real pp = k - s * math::pow(up, 2.0 * i - simulation);
			pp = pp > 0.0 ? pp : 0.0;
			p.push_back(pp);
		}

		for (int j = simulation - 1; j != 0; --j) {
			for (int i = 0; i != j; ++i) {
				p[i] = p0 * p[i + 1] + p1 * p[i];    // binomial value
				Real exercise = k - s * math::pow(up, 2.0 * i - j);  // exercise value
				p[i] = p[i] > exercise ? p[i] : exercise;
			}
		}

		return p[0];
	}
}

void hiho::ad13_expr_vec_tape_vec_lazy_pmr(double s, double sigma, double k, double r, double t, int simulation)
{
	auto dfpm = pmr::get_default_resource();
	auto pm = pmr::monotonic_buffer_resource();
	pmr::set_default_resource(&pm);
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
	pmr::set_default_resource(dfpm);
}