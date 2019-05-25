#include "ad10_expr_vec_tape_vec_shrink.hpp"

#include <cmath>
#include <vector>
#include <utility>
#include <unordered_map>

#include "ad_american.hpp"
#include "timer.hpp"
#include "iohelper.hpp"

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
				const static size_t dmmy = std::numeric_limits<size_t>::max();
				static Expression dummy{ dmmy };
				if (index == dmmy) { return dummy; }
				return expressions[index];
			}
			static void shrinkExpressions() {
				for (auto& expr : expressions) {
					expr.shrink();
				}
			}

			size_t index;
			size_t ref;
			Polynomial polynomial;

			Expression(size_t indx) : index(indx), ref(0), polynomial{} {}
			Expression(size_t indx,
				ValueType cof, size_t expr
			) : index(indx), ref(0), polynomial{ {cof, expr} } {}
			Expression(size_t indx,
				ValueType lcof, size_t lhs,
				ValueType rcof, size_t rhs
			) : index(indx), ref(0), polynomial{ {lcof, lhs},{rcof, rhs} } {}

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

			static void appendTerm(Polynomial& polynomial, const Term& term) {
				for (auto& t : polynomial) {
					if (t.second == term.second) {
						t.first += term.first;
						return;
					}
				}
				polynomial.emplace_back(term);
			}

			void shrink() {
				Polynomial newpoly;

				for (auto& t0 : polynomial) {
					auto& t0exp = getExpression(t0.second);
					if (t0exp.ref != 0) {
						appendTerm(newpoly, t0);
					} else {
						for (auto& t1 : t0exp.polynomial) {
							if (getExpression(t1.second).ref != 0) {
								appendTerm(newpoly, {t0.first * t1.first, t1.second});
							}
						}
					}
				}

				polynomial = std::move(newpoly);
			}

			void reference() { ++ref; }
			void dreference() { --ref; }
		};

		struct Number {
			ValueType  v;
			size_t expr_;

			Number() : v{ 0 }, expr_{ std::numeric_limits<size_t>::max() }  {}
			Number(ValueType vv) : v{ vv }, expr_{ Expression::newExpression() }  {expression().reference(); }
			Number(ValueType vv, size_t expr) : v{ vv }, expr_{ expr }  {expression().reference(); }
			Number(const Number& other) : v { other.v }, expr_{other.expr_} { expression().reference(); };
			Number(Number&& other) : v{ other.v }, expr_{ other.expr_ } { expression().reference(); };
			~Number() { expression().dreference(); }

			ValueType d(const Number& x) const {
				Cache cache;
				return expression().d(x.expression(), cache);
			}

			Expression& expression() const { return Expression::getExpression(expr_); }

			Number operator-() const { return Number{ -v, Expression::newExpression(-1, expr_) }; }
			Number& operator=(const Number& other) {
				this->expression().dreference();
				this->v = other.v;
				this->expr_ = other.expr_;
				this->expression().reference();
				return *this;
			}
		};

		Number operator+(const Number& l, const Number& r) { return Number{ l.v + r.v, Expression::newExpression(1, l.expr_, 1, r.expr_) }; }
		Number operator-(const Number& l, const Number& r) { return Number{ l.v - r.v, Expression::newExpression(1, l.expr_, -1, r.expr_) }; }
		Number operator*(const Number& l, const Number& r) { return Number{ l.v * r.v, Expression::newExpression(r.v, l.expr_, l.v, r.expr_) }; }
		Number operator/(const Number& l, const Number& r) {
			auto ll = l.v;
			auto rr = r.v;
			return Number{ l.v / r.v, Expression::newExpression(1.0 / rr, l.expr_, -ll / (rr * rr), r.expr_) };
		}
		Number operator+(const Number& l, ValueType r) { return Number{ l.v + r, Expression::newExpression(1, l.expr_) }; }
		Number operator-(const Number& l, ValueType r) { return Number{ l.v - r, Expression::newExpression(1, l.expr_) }; }
		Number operator*(const Number& l, ValueType r) { return Number{ l.v * r, Expression::newExpression(r, l.expr_) }; }
		Number operator/(const Number& l, ValueType r) { return Number{ l.v / r, Expression::newExpression(1.0 / r, l.expr_) }; }
		Number operator+(ValueType l, const Number& r) { return Number{ l + r.v, Expression::newExpression(1, r.expr_) }; }
		Number operator-(ValueType l, const Number& r) { return Number{ l - r.v, Expression::newExpression(-1, r.expr_) }; }
		Number operator*(ValueType l, const Number& r) { return Number{ l * r.v, Expression::newExpression(l, r.expr_) }; }
		Number operator/(ValueType l, const Number& r) { return Number{ l / r.v, Expression::newExpression(-l / (r.v * r.v), r.expr_) }; }
		bool operator>(const Number& l, const Number& r) { return l.v > r.v; }
		Number exp(const Number& l) {
			auto ll = std::exp(l.v);
			return Number{ ll, Expression::newExpression(ll, l.expr_) };
		}
		Number sqrt(const Number& l) {
			auto ll = std::sqrt(l.v);
			return Number{ ll, Expression::newExpression(0.5 / ll, l.expr_) };
		}
		Number pow(const Number& l, ValueType r) {
			auto ll = std::pow(l.v, r);
			return Number{ ll, Expression::newExpression(r * ll / l.v, l.expr_) };
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

		math::Expression::shrinkExpressions();

		for (int j = simulation - 1; j != 0; --j) {
			for (int i = 0; i != j; ++i) {
				p[i] = p0 * p[i + 1] + p1 * p[i];    // binomial value
				auto exercise = k - s * math::pow(up, 2.0 * i - j);  // exercise value
				p[i] = p[i] > exercise ? p[i] : exercise;
			}
		}

		math::Expression::shrinkExpressions();

		return p[0];
	}
}

void hiho::ad10_expr_vec_tape_vec_shrink(double s, double sigma, double k, double r, double t, int simulation)
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
		math::Expression::shrinkExpressions();

		auto diff = value.v - hiho::american(s, sigma, k, r, t, simulation);

		auto delta = hiho::newTimer([&]() {return value.d(rs); });
		auto vega = hiho::newTimer([&]() {return value.d(rsigma); });
		auto theta = hiho::newTimer([&]() {return value.d(rt); });

		HIHO_IO_MAX_LEN_DOUBLE_LSHOW;
		HIHO_IO_LEFT_COUT
			<< HIHO_IO_FUNC_WIDTH << __func__ << " ( " << simulation << " )";
		HIHO_IO_RIGHT_COUT
			<< ", " << HIHO_IO_VALUE(diff)
			<< ", " << HIHO_IO_TIME(time)
			<< ", " << HIHO_IO_VALUE_TIME(delta)
			<< ", " << HIHO_IO_VALUE_TIME(vega)
			<< ", " << HIHO_IO_VALUE_TIME(theta)
			<< std::endl;
	}
	math::Expression::counter = 0;
	math::Expression::expressions = {};
}
