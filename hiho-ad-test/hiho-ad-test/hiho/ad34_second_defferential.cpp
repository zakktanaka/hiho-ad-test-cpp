#include "ad34_second_defferential.hpp"

#include <cmath>
#include <utility>

#include <memory_resource>
#include <array>

#include "ad_american.hpp"
#include "timer.hpp"
#include "iohelper.hpp"
#include <cassert>

namespace {

	namespace pmr = std::pmr;

	namespace math {

		template<typename T>
		class Expression {
		private:
			using ValueType  = typename T;
			using ExprPtr    = Expression *;
			using Term       = std::pair<ExprPtr, ValueType>;
			using Polynomial = pmr::vector<Term>;


			bool marked;
			size_t ref;
			Polynomial polynomial;

		public:
			using Cache = pmr::unordered_map<const void*, ValueType>;

			Expression() : marked{ false }, ref(0), polynomial{} {}
			Expression(const Expression&) = default;
			Expression(Expression&&) = default;

			Expression& operator=(const Expression&) = default;
			Expression& operator=(Expression&&) = default;

			void mark() { marked = true; }
			void reference() { ++ref; }
			void dereference() { --ref; }
			bool referred() const { return ref != 0; }

			ValueType d(const ExprPtr expr, Cache& cache) const {
				if (this == expr) {
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

			void append(const ExprPtr expr, ValueType coef) {
				if (expr->marked) {
					for (auto& tm : polynomial) {
						if (tm.first == expr) {
							tm.second += coef;
							return;
						}
					}
					polynomial.emplace_back(expr, coef);
				}
				else {
					for (auto& tt : expr->polynomial) {
						append(tt.first, coef * tt.second);
					}
				}
			}

			void clear() {
				marked = false;
				ref = 0;
				polynomial = {};
			}
		};

		template<typename T>
		class Store {
		private:
			static constexpr size_t bulk = 100000;
			using ValueType = typename T;
			using Datum     = Expression<ValueType>;

			struct Element {
				Element* next;
				Datum* datum;
			};

			using Data     = pmr::vector<Datum>;
			using Elements = pmr::vector<Element>;

			struct DataElement {
				Elements elements;
				Data     data;
			};

			Element* remains;
			Element* unassigned;
			size_t   capacity;
			size_t   used;
			Datum    defaultDatum;
			pmr::vector<DataElement> reserved;

			void resize() {
				reserved.emplace_back(DataElement{ Elements{bulk, Element{nullptr, nullptr}}, Data{bulk, defaultDatum} });
				auto& elems = *std::rbegin(reserved);

				auto current = &remains;
				for (size_t i = 0; i < bulk; ++i) {
					auto& elem = elems.elements[i];
					*current = &elem;
					current = &(elem.next);
					elem.datum = &(elems.data[i]);
				}

				capacity += bulk;
			}

		public:
			Store() : remains{ nullptr }, unassigned{ nullptr }, capacity{ 0 }, used{ 0 }, defaultDatum{}, reserved{} {};

			Datum* getDatum() {
				if (remains == nullptr) { resize(); }

				auto current = remains;
				remains = current->next;

				auto d = current->datum;
				current->datum = nullptr;

				current->next = unassigned;
				unassigned = current;

				++used;
				return d;
			}

			void returnBack(Datum* datum) {
				auto current = unassigned;
				unassigned = current->next;

				*datum = defaultDatum;

				current->datum = datum;
				current->next = remains;

				remains = current;
				--used;
			}
		};

		template<typename T>
		Store<T>* repository;

		template<typename T>
		struct INumber {
			using ValueType  = typename T;
			using Expr     = Expression<ValueType>;
			using ExprPtr  = Expr *;

			virtual ~INumber() {}
			virtual ValueType v() const = 0;
			virtual void update(ExprPtr, ValueType) const = 0;

		};

		template<typename T>
		struct Unary : INumber<T> {
			using ValueType  = typename T;
			using Expr     = Expression<ValueType>;
			using ExprPtr  = Expr *;
			using BaseType = INumber<ValueType>;
			using ThisType = Unary<ValueType>;

			ValueType v_;
			ValueType coef_;
			const BaseType& num_;

			Unary(ValueType vv, ValueType coef, const BaseType& num) :
				v_(vv), coef_(coef), num_(num) {}

			ValueType v() const override { return v_; }
			void update(ExprPtr updated, ValueType coef) const override {
				num_.update(updated, coef * coef_);
			}
		};

		template<typename T>
		struct Binary : INumber<T> {
			using ValueType  = typename T;
			using Expr     = Expression<ValueType>;
			using ExprPtr  = Expr *;
			using BaseType = INumber<ValueType>;
			using ThisType = Unary<ValueType>;

			ValueType v_;
			ValueType lcoef_;
			ValueType rcoef_;
			const BaseType& lnum_;
			const BaseType& rnum_;

			Binary(
				ValueType vv,
				ValueType lcoef, const BaseType& lnum,
				ValueType rcoef, const BaseType& rnum) :
				v_(vv), lcoef_(lcoef), rcoef_(rcoef), lnum_(lnum), rnum_(rnum) {}

			ValueType v() const override { return v_; }
			void update(ExprPtr updated, ValueType coef) const override {
				lnum_.update(updated, coef * lcoef_);
				rnum_.update(updated, coef * rcoef_);
			}
		};

		template<typename T>
		class Number : public INumber<T> {
		private:
			using ValueType  = typename T;
			using Expr     = Expression<ValueType>;
			using ExprPtr  = Expr *;
			using BaseType = INumber<ValueType>;
			using ThisType = Number<ValueType>;
			using Cache    = typename Expr::Cache;

			ValueType  v_;
			ExprPtr expr_;

			void reference() { expr_->reference(); }
			void dereference() {
				expr_->dereference();
				if (!expr_->referred()) {
					repository<ValueType>->returnBack(expr_);
				}
			}

		public:
			Number(ValueType vv, ExprPtr expr) : v_{ vv }, expr_{ expr } { reference(); }
			Number(ValueType vv) : Number{ vv ,       repository<ValueType>->getDatum() } { }
			Number(const ThisType& other) : Number{ other.v_ , other.expr_ } { }
			Number(const BaseType& other) : Number{ other.v(), repository<ValueType>->getDatum() } { other.update(expr_, 1); }
			~Number() { dereference(); }

			void mark() { expr_->mark(); }

			ValueType v() const override { return v_; }
			void update(ExprPtr updated, ValueType coef) const override {
				updated->append(expr_, coef);
			}

			ValueType d(const ThisType& x) const {
				Cache cache;
				return expr_->d(x.expr_, cache);
			}

			ThisType operator-() const {
				auto expr = repository<ValueType>->getDatum();
				expr->append(expr_, -1);

				return ThisType{ -v_, expr };
			}

			ThisType& operator=(const ThisType& other) {
				if (this == &other) {
					return *this;
				}

				ThisType newone{ other.v() };
				other.update(newone.expr_, 1);

				dereference();
				v_ = newone.v_;
				expr_ = newone.expr_;
				reference();

				return *this;
			}

			ThisType& operator=(const BaseType& other) {
				if (this == &other) {
					return *this;
				}

				ThisType newone{ other.v() };
				other.update(newone.expr_, 1);

				dereference();
				v_ = newone.v_;
				expr_ = newone.expr_;
				reference();

				return *this;
			}

			ThisType& operator+=(const ThisType& other) {
				ThisType newone = *this + other;
				
				dereference();
				v_ = newone.v_;
				expr_ = newone.expr_;
				reference();

				return *this;
			}

			ThisType& operator+=(ValueType other) {
				ThisType newone = *this + other;

				dereference();
				v_ = newone.v_;
				expr_ = newone.expr_;
				reference();

				return *this;
			}

		};

		template<typename T> Binary<T> operator+(const INumber<T>& l, const INumber<T>& r) { return Binary<T>{ l.v() + r.v(), 1, l, 1, r }; }
		template<typename T> Binary<T> operator-(const INumber<T>& l, const INumber<T>& r) { return Binary<T>{ l.v() - r.v(), 1, l, -1, r }; }
		template<typename T> Binary<T> operator*(const INumber<T>& l, const INumber<T>& r) { return Binary<T>{ l.v()* r.v(), r.v(), l, l.v(), r }; }
		template<typename T> Binary<T> operator/(const INumber<T>& l, const INumber<T>& r) {
			auto ll = l.v();
			auto rr = r.v();
			return Binary<T>{ l.v() / r.v(), 1.0 / rr, l, -ll / (rr * rr), r };
		}
		template<typename T> Unary<T> operator+(const INumber<T>& l, typename INumber<T>::ValueType r) { return Unary<T>{ l.v() + r, 1, l }; }
		template<typename T> Unary<T> operator-(const INumber<T>& l, typename INumber<T>::ValueType r) { return Unary<T>{ l.v() - r, 1, l }; }
		template<typename T> Unary<T> operator*(const INumber<T>& l, typename INumber<T>::ValueType r) { return Unary<T>{ l.v()* r, r, l }; }
		template<typename T> Unary<T> operator/(const INumber<T>& l, typename INumber<T>::ValueType r) { return Unary<T>{ l.v() / r, 1.0 / r, l }; }
		template<typename T> Unary<T> operator+(typename INumber<T>::ValueType l, const INumber<T>& r) { return Unary<T>{ l + r.v(), 1, r }; }
		template<typename T> Unary<T> operator-(typename INumber<T>::ValueType l, const INumber<T>& r) { return Unary<T>{ l - r.v(), -1, r }; }
		template<typename T> Unary<T> operator*(typename INumber<T>::ValueType l, const INumber<T>& r) { return Unary<T>{ l* r.v(), l, r }; }
		template<typename T> Unary<T> operator/(typename INumber<T>::ValueType l, const INumber<T>& r) { return Unary<T>{ l / r.v(), -l / (r.v() * r.v()), r }; }
		template<typename T> bool operator>(const INumber<T>& l, const INumber<T>& r) { return l.v() > r.v(); }
		template<typename T> bool operator>(const INumber<T>& l, typename INumber<T>::ValueType r) { return l.v() > r; }

		using std::exp;
		using std::sqrt;
		using std::pow;

		template<typename T> Unary<T> exp(const INumber<T>& l) {
			using ValueType = typename INumber<T>::ValueType;

			ValueType ll = exp(l.v());
			return Unary<T>{ ll, ll, l };
		}
		template<typename T> Unary<T> sqrt(const INumber<T>& l) {
			using ValueType = typename INumber<T>::ValueType;

			ValueType ll = sqrt(l.v());
			return Unary<T>{ ll, 0.5 / ll, l };
		}
		template<typename T> Unary<T> pow(const INumber<T>& l, typename INumber<T>::ValueType r) {
			using ValueType = typename INumber<T>::ValueType;

			ValueType ll = 1; // math::pow(l.v(), r);
			return Unary<T>{ ll, r* ll / l.v(), l };
		}
	}

	using Real = math::Number < math::Number<double>>;

	inline Real putAmericanOption(Real s, Real sigma, Real k, Real r, Real t, int simulation) {

		Real dt = t / simulation;
		Real up = math::exp(sigma * math::sqrt(dt));

		Real p0 = (up - math::exp(-r * dt)) / (up * up - 1);
		Real p1 = math::exp(-r * dt) - p0;

		std::vector<Real> p;
		for (int i = 0; i != simulation; ++i) {
			Real pp = k - s * math::pow(up, 2.0 * i - simulation);
			Real ppp = pp > 0.0 ? pp : math::Number<double>{ 0.0 };
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

void hiho::ad34_second_defferential(double s, double sigma, double k, double r, double t, int simulation)
{
	//auto dfpm = pmr::get_default_resource();
	//auto pm = pmr::unsynchronized_pool_resource();
	//pmr::set_default_resource(&pm);
	//{
	//	math::Store<double> rep0;
	//	math::repository<double> = &rep0;
	//	math::Store<math::Number<double>> rep1;
	//	math::repository<math::Number<double>> = &rep1;

	//	auto func = [&]() {
	//		Real rs{ s }; rs.mark();
	//		Real rsigma{ sigma }; rsigma.mark();
	//		Real rk{ k };
	//		Real rr{ r }; rr.mark();
	//		Real rt{ t }; rt.mark();
	//		auto value = putAmericanOption(rs, rsigma, rk, rr, rt, simulation);
	//		return pmr::vector<Real>{ value, rs, rsigma, rr, rt };
	//	};

	//	{
	//		auto time = hiho::measureTime<3>(func);
	//		auto vv = func();
	//		auto& value = vv[0];
	//		auto& rs = vv[1];
	//		auto& rsigma = vv[2];
	//		auto& rr = vv[3];
	//		auto& rt = vv[4];

	//		auto diff = value.v().v() - hiho::american(s, sigma, k, r, t, simulation);

	//		auto delta = hiho::newTimer([&]() {return value.d(rs).v(); });
	//		auto vega = hiho::newTimer([&]() {return value.d(rsigma).v(); });
	//		auto theta = hiho::newTimer([&]() {return value.d(rt).v(); });

	//		HIHO_IO_MAX_LEN_DOUBLE_LSHOW;
	//		HIHO_IO_LEFT_COUT
	//			<< HIHO_IO_FUNC_WIDTH << __func__ << " ( " << simulation << " )";
	//		HIHO_IO_RIGHT_COUT
	//			<< ", " << HIHO_IO_VALUE(diff)
	//			<< ", " << HIHO_IO_TIME(time)
	//			<< ", " << HIHO_IO_VALUE_TIME(delta)
	//			<< ", " << HIHO_IO_VALUE_TIME(vega)
	//			<< ", " << HIHO_IO_VALUE_TIME(theta)
	//			<< std::endl;
	//	}
	//}
	//pmr::set_default_resource(dfpm);
}
