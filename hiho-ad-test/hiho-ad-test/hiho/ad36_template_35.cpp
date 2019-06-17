#include "ad36_template_35.hpp"

#include <cmath>
#include <utility>
#include <memory>

#include <memory_resource>

#include "ad_american.hpp"
#include "timer.hpp"
#include "iohelper.hpp"

namespace {

	namespace pmr = std::pmr;

	namespace math {

		template<typename EV>
		class Expression {
		public:
			using ExprValueType = EV;
			using ThisType      = Expression<ExprValueType>;
			using Cache         = std::unordered_map<const void*, ExprValueType>;
			using ExpPtr        = std::shared_ptr<ThisType>;
			using Term          = std::pair<ExpPtr, ExprValueType>;
			using Polynomial    = pmr::vector<Term>;

		private:
			bool marked;
			Polynomial polynomial;

			void addTerm(ExprValueType coef, const ExpPtr& other) {
				for (auto& tm : polynomial) {
					if (tm.first == other) {
						tm.second += coef;
						return;
					}
				}
				polynomial.emplace_back(other, coef);
			}

		public:
			Expression() : marked{ false }, polynomial{} {}

			void mark() { marked = true; }

			ExprValueType d(ExpPtr& expr, Cache& cache) const {
				if (this == expr.get()) {
					return 1;
				}

				auto it = cache.find(this);
				if (it != std::end(cache)) {
					return it->second;
				}

				ExprValueType dx = 0;
				for (auto& term : polynomial) {
					dx += term.second * term.first->d(expr, cache);
				}
				cache[this] = dx;
				return dx;
			}

			void addExpresison(ExprValueType coef, const ExpPtr& other) {
				if (other->marked) {
					addTerm(coef, other);
				} else {
					for (auto& otm : other->polynomial) {
						auto c = coef * otm.second;
						auto& e = otm.first;
						addExpresison(c, e);
					}
				}
			}
		};

		template<typename V, typename E>
		struct INumber {
			using ValueType     = V;
			using Expr          = E;
			using ExprValueType = typename E::ExprValueType;
			using ThisType      = INumber<ValueType, Expr>;

			virtual ~INumber() {}
			virtual ValueType  v() const = 0;
			virtual void update(Expr&, const ExprValueType&) const = 0;

		};

		template<typename V, typename E>
		struct Unary : INumber<V, E> {
			using ValueType     = V;
			using Expr          = E;
			using ExprValueType = typename E::ExprValueType;
			using BaseType      = INumber<ValueType, Expr>;
			using ThisType      = Unary<ValueType, Expr>;

			ValueType       v_;
			ExprValueType   coef_;
			const BaseType& num_;

			Unary(ValueType vv, const ExprValueType& coef, const BaseType& num) :
				v_(vv), coef_(coef), num_(num) {}

			ValueType v() const override { return v_; }

			void update(Expr& updated, const ExprValueType& coef) const override {
				num_.update(updated, coef * coef_);
			}
		};

		template<typename V, typename E>
		struct Binary : INumber<V, E> {
			using ValueType     = V;
			using Expr          = E;
			using ExprValueType = typename E::ExprValueType;
			using BaseType      = INumber<ValueType, Expr>;
			using ThisType      = Binary<ValueType, Expr>;

			ValueType      v_;
			ExprValueType  lcoef_;
			ExprValueType  rcoef_;
			const BaseType& lnum_;
			const BaseType& rnum_;

			Binary(
				ValueType vv,
				const ExprValueType& lcoef, const BaseType& lnum,
				const ExprValueType& rcoef, const BaseType& rnum) :
				v_(vv), lcoef_(lcoef), rcoef_(rcoef), lnum_(lnum), rnum_(rnum) {}

			ValueType v() const override { return v_; }

			void update(Expr& updated, const ExprValueType& coef) const override {
				lnum_.update(updated, coef * lcoef_);
				rnum_.update(updated, coef * rcoef_);
			}
		};

		template<typename V, typename E>
		class Number : public INumber<V, E> {
		public :
			using ValueType     = V;
			using Expr          = E;
			using ExprValueType = typename E::ExprValueType;
			using BaseType      = INumber<ValueType, Expr>;
			using ThisType      = Number<ValueType, Expr>;

		private:
			using Cache  = typename Expr::Cache;
			using ExpPtr = typename Expr::ExpPtr;

			ValueType  v_;
			ExpPtr     expr_;

			Number(ValueType vv, ExpPtr expr) : v_{ vv }, expr_{ expr }  {}

		public:
			Number(ValueType vv) :
				Number{ vv , std::make_shared<Expr>() } {}
			Number(const BaseType& other) :
				Number{ other.v(), std::make_shared<Expr>() } {
				other.update(*expr_, 1);
			};
			Number() : Number(ValueType(1)) {}
			~Number() {}

			void mark() { expr_->mark(); }

			ValueType v() const override { return v_; }

			void update(Expr& updated, const ExprValueType& coef) const override {
				updated.addExpresison(coef, expr_);
			}

			ExprValueType d(ThisType& x) const {
				Cache cache;
				return expr_->d(x.expr_, cache);
			}

			ThisType operator-() const {
				auto expr = std::make_shared<Expr>();
				expr->addExpresison(-1, expr_);
				return Number{ -v_, expr };
			}
			ThisType& operator=(const BaseType& other) {
				ThisType newone{ other.v() };
				other.update(*(newone.expr_), 1);
				this->v_ = newone.v_;
				this->expr_ = newone.expr_;
				return *this;
			}
			ThisType& operator+=(const ThisType& other) {
				ThisType newone = *this + other;
				this->v_ = newone.v_;
				this->expr_ = newone.expr_;
				return *this;
			}
		};

		template<typename V, typename E> Binary<V, E> operator+(const INumber<V, E>& l, const INumber<V, E>& r) { return Binary<V, E>{ l.v() + r.v(), 1, l, 1, r }; }
		template<typename V, typename E> Binary<V, E> operator-(const INumber<V, E>& l, const INumber<V, E>& r) { return Binary<V, E>{ l.v() - r.v(), 1, l, -1, r }; }
		template<typename V, typename E> Binary<V, E> operator*(const INumber<V, E>& l, const INumber<V, E>& r) { return Binary<V, E>{ l.v() * r.v(), r.v(), l, l.v(), r }; }
		template<typename V, typename E> Binary<V, E> operator/(const INumber<V, E>& l, const INumber<V, E>& r) {
			auto ll = l.v();
			auto rr = r.v();
			return Binary<V, E>{ l.v() / r.v(), 1.0 / rr, l, -ll / (rr * rr), r };
		}
		template<typename V, typename E> Unary<V, E> operator+(const INumber<V, E>& l, typename INumber<V, E>::ValueType r) { return Unary<V, E>{ l.v() + r, 1, l }; }
		template<typename V, typename E> Unary<V, E> operator-(const INumber<V, E>& l, typename INumber<V, E>::ValueType r) { return Unary<V, E>{ l.v() - r, 1, l }; }
		template<typename V, typename E> Unary<V, E> operator*(const INumber<V, E>& l, typename INumber<V, E>::ValueType r) { return Unary<V, E>{ l.v() * r, r, l }; }
		template<typename V, typename E> Unary<V, E> operator/(const INumber<V, E>& l, typename INumber<V, E>::ValueType r) { return Unary<V, E>{ l.v() / r, 1.0 / r, l }; }
		template<typename V, typename E> Unary<V, E> operator+(typename INumber<V, E>::ValueType l, const INumber<V, E>& r) { return Unary<V, E>{ l + r.v(), 1, r }; }
		template<typename V, typename E> Unary<V, E> operator-(typename INumber<V, E>::ValueType l, const INumber<V, E>& r) { return Unary<V, E>{ l - r.v(), -1, r }; }
		template<typename V, typename E> Unary<V, E> operator*(typename INumber<V, E>::ValueType l, const INumber<V, E>& r) { return Unary<V, E>{ l * r.v(), l, r }; }
		template<typename V, typename E> Unary<V, E> operator/(typename INumber<V, E>::ValueType l, const INumber<V, E>& r) { return Unary<V, E>{ l / r.v(), -l / (r.v() * r.v()), r }; }
		template<typename V, typename E> bool operator>(const INumber<V, E>& l, const INumber<V, E>& r) { return l.v() > r.v(); }
		template<typename V, typename E> bool operator>(const INumber<V, E>& l, typename INumber<V, E>::ValueType r) { return l.v() > r; }

		using std::exp;
		using std::sqrt;
		using std::pow;

		template<typename V, typename E>
		Unary<V, E> exp(const INumber<V, E>& l) {
			auto ll = exp(l.v());
			return Unary{ ll, ll, l };
		}
		template<typename V, typename E>
		Unary<V, E> sqrt(const INumber<V, E>& l) {
			auto ll = sqrt(l.v());
			return Unary{ ll, 0.5 / ll, l };
		}
		template<typename V, typename E>
		Unary<V, E> pow(const INumber<V, E>& l, typename INumber<V, E>::ValueType r) {
			auto ll = pow(l.v(), r);
			return Unary{ ll, r * ll / l.v(), l };
		}

	}

	using Number  = math::Number<double, math::Expression<double>>;

	using Real = Number;

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
				p[i] = p0 * p[i + int(1)] + p1 * p[i];    // binomial value
				Real exercise = k - s * math::pow(up, 2.0 * i - j);  // exercise value
				p[i] = p[i] > exercise ? p[i] : exercise;
			}
		}

		return p[0];
	}
}

void hiho::ad36_template_35(double s, double sigma, double k, double r, double t, int simulation)
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
