#include "ad21_exprump_tpv_lazy_pmr_shrnk.hpp"

#include <cmath>
#include <utility>

#include <memory_resource>

#include "ad_american.hpp"
#include "timer.hpp"
#include "iohelper.hpp"

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
			using Polynomial = pmr::unordered_map<size_t, ValueType>;

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

			size_t ref;
			size_t index;
			Polynomial polynomial;

			Expression(size_t indx) : ref(0), index(indx), polynomial{} {}
			Expression(size_t indx,
				ValueType cof, size_t expr
			) : ref(0), index(indx), polynomial{} {
				polynomial[expr] = cof;
			}
			Expression(size_t indx,
				ValueType lcof, size_t lhs,
				ValueType rcof, size_t rhs
			) : ref(0), index(indx), polynomial{} {
				polynomial[lhs] = lcof;
				polynomial[rhs] = rcof;
			}

			void reference() { ++ref; }
			void dreference() { --ref; }

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
					dx += term.second * getExpression(term.first).d(expr, cache);
				}
				cache[index] = dx;
				return dx;
			}

			void addTerm(size_t expr, ValueType val) {
				auto it = polynomial.find(expr);
				if (it == std::end(polynomial)) {
					polynomial[expr] = val;
				} else {
					it->second += val;
				}
			}

			static void appendTerm(Polynomial& polynomial, size_t expr, ValueType val) {
				auto it = polynomial.find(expr);
				if (it == std::end(polynomial)) {
					polynomial[expr] = val;
				} else {
					it->second += val;
				}
			}

			void shrink() {
				Polynomial newpoly;

				for (auto& t0 : polynomial) {
					auto& t0exp = getExpression(t0.first);
					if (t0exp.ref != 0) {
						appendTerm(newpoly, t0.first, t0.second);
					}
					else {
						for (auto& t1 : t0exp.polynomial) {
							if (getExpression(t1.first).ref != 0) {
								appendTerm(newpoly, t1.first, t0.second * t1.second);
							}
						}
					}
				}

				polynomial = std::move(newpoly);
			}

			static void shrinkExpressions() {
				for (auto& expr : expressions) {
					expr.shrink();
				}
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

			Number(ValueType vv) : v_{ vv }, expr_{ Expression::newExpression() }  { expression().reference(); }
			Number(ValueType vv, size_t expr) : v_{ vv }, expr_{ expr }  {expression().reference(); }
			Number(const Number& other) : v_{ other.v_ }, expr_{ other.expr_ } {expression().reference(); }
			Number(const INumber& other) : v_{ other.v() }, expr_{ Expression::newExpression() } {
				expression().reference();
				other.update(expression(), 1);
			};
			~Number() { expression().dreference(); }

			ValueType v() const override { return v_; }
			void update(Expression& updated, ValueType coef) const override {
				updated.addTerm(expr_, coef);
			}

			ValueType d(const Number& x) const {
				Cache cache;
				return expression().d(x.expression(), cache);
			}

			Expression& expression() const { return Expression::getExpression(expr_); }

			Number operator-() const { return Number{ -v_, Expression::newExpression(-1, expr_) }; }
			Number& operator=(const Number& other) {
				this->expression().dreference();
				this->v_ = other.v_;
				this->expr_ = other.expr_;
				this->expression().reference();
				return *this;
			}
			Number& operator=(const INumber& other) {
				this->expression().dreference();
				Number newone{ other.v() };
				other.update(newone.expression(), 1);
				this->v_ = newone.v_;
				this->expr_ = newone.expr_;
				this->expression().reference();
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
			Real ppp = pp > 0.0 ? pp : 0.0;
			p.emplace_back(ppp);
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

void hiho::ad21_exprump_tpv_lazy_pmr_shrnk(double s, double sigma, double k, double r, double t, int simulation)
{
	auto dfpm = pmr::get_default_resource();
	auto pm = pmr::unsynchronized_pool_resource();
	pmr::set_default_resource(&pm);
	{
		auto func = [&]() {
			Real rs{ s };
			Real rsigma{ sigma };
			Real rr{ r };
			Real rt{ t };
			auto value = putAmericanOption(rs, rsigma, k, rr, rt, simulation);
			math::Expression::shrinkExpressions();
			return pmr::vector<Real>{ value, rs, rsigma, rr, rt };
		};
		auto postprocess = []() {
			//auto c = 0; for (auto& e : math::Expression::expressions) { c += e.ref != 0 ? 1 : 0; }
			//std::cout << math::Expression::counter << ", " << c << std::endl;

			math::Expression::counter = 0;
			math::Expression::expressions = {};
		};

		{
			auto time = hiho::measureTime<3>(func, postprocess);
			auto vv = func();
			auto& value = vv[0];
			auto& rs = vv[1];
			auto& rsigma = vv[2];
			auto& rr = vv[3];
			auto& rt = vv[4];

			auto diff = value.v() - hiho::american(s, sigma, k, r, t, simulation);

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
		postprocess();
	}
	pmr::set_default_resource(dfpm);
}
