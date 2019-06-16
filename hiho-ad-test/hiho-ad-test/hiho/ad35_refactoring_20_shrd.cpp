#include "ad35_refactoring_20_shrd.hpp"

#include <cmath>
#include <utility>
#include <memory>

#include <memory_resource>

#include "ad_american.hpp"
#include "timer.hpp"
#include "iohelper.hpp"

namespace {

	namespace pmr = std::pmr;

	using ValueType = double;

	namespace math {

		struct Expression {
			using Cache      = std::unordered_map<const void*, ValueType>;
			using ExpPtr     = std::shared_ptr<math::Expression>;
			using Term       = std::pair<ExpPtr, ValueType>;
			using Polynomial = pmr::vector<Term>;

			bool marked;
			Polynomial polynomial;

			Expression() : marked{ false }, polynomial{} {}
			Expression(
				ValueType cof, ExpPtr expr
			) : marked{ false }, polynomial{ {expr, cof} } {}
			Expression(
				ValueType lcof, ExpPtr lhs,
				ValueType rcof, ExpPtr rhs
			) : marked{ false }, polynomial{ {lhs, lcof},{rhs, rcof} } {}

			void mark() { marked = true; }

			ValueType d(ExpPtr& expr, Cache& cache) const {
				if (this == expr.get()) {
					return 1;
				}

				auto it = cache.find(this);
				if (it != std::end(cache)) {
					return it->second;
				}

				ValueType dx = 0;
				for (auto& term : polynomial) {
					dx += term.second * term.first->d(expr, cache);
				}
				cache[this] = dx;
				return dx;
			}

			void addTerm(const Term& term) {
				for (auto& tm : polynomial) {
					if (tm.first == term.first) {
						tm.second += term.second;
						return;
					}
				}
				polynomial.emplace_back(term);
			}

			void addExpresison(ValueType coef, const ExpPtr& other) {
				if (other->marked) {
					addTerm({ other, coef});
				}
				else {
					for (auto& otm : other->polynomial) {
						auto c = coef * otm.second;
						auto& e = otm.first;
						addExpresison(c, e);
					}
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
			using Cache  = typename Expression::Cache;
			using ExpPtr = typename Expression::ExpPtr;
			ValueType  v_;
			ExpPtr expr_;

			Number(ValueType vv) : v_{ vv }, expr_{ std::make_shared<Expression>() }  {}
			Number(ValueType vv, ExpPtr expr) : v_{ vv }, expr_{ expr }  {}
			Number(const INumber& other) : v_{ other.v() }, expr_{ std::make_shared<Expression>() } { other.update(expression(), 1); };
			~Number() {}

			void mark() { expr_->mark(); }

			ValueType v() const override { return v_; }
			void update(Expression& updated, ValueType coef) const override {
				updated.addExpresison(coef, expr_);
			}

			ValueType d(Number& x) const {
				Cache cache;
				return expression().d(x.expr_, cache);
			}

			Expression& expression() const { return *expr_; }

			Number operator-() const { return Number{ -v_, std::make_shared<Expression>(-1, expr_) }; }
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
		
		using std::exp;
		using std::sqrt;
		using std::pow;

		Unary exp(const INumber& l) {
			auto ll = exp(l.v());
			return Unary{ ll, ll, l };
		}
		Unary sqrt(const INumber& l) {
			auto ll = sqrt(l.v());
			return Unary{ ll, 0.5 / ll, l };
		}
		Unary pow(const INumber& l, ValueType r) {
			auto ll = pow(l.v(), r);
			return Unary{ ll, r * ll / l.v(), l };
		}

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

void hiho::ad35_refactoring_20_shrd(double s, double sigma, double k, double r, double t, int simulation)
{
	auto dpm = pmr::get_default_resource();
	auto pm = pmr::unsynchronized_pool_resource();
	pmr::set_default_resource(&pm);

	{
		Real rs{ s }; rs.mark();
		Real rsigma{ sigma }; rsigma.mark();
		Real rr{ r }; rr.mark();
		Real rt{ t }; rt.mark();

		auto func = [&]() {
			return putAmericanOption(rs, rsigma, k, rr, rt, simulation);
		};
		auto time = hiho::measureTime<3>(func);
		auto value = func();

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

	pmr::set_default_resource(dpm);
}
