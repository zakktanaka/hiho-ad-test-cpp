#include "ad25_exprump_tp_pclass_pzypmrmrk.hpp"

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
	using Cache = pmr::unordered_map<const void*, ValueType>;

	namespace math {

		struct Expression {
			using ExprPtr = Expression*;
			using Term       = std::pair<ValueType, ExprPtr>;
			using Polynomial = pmr::vector<Term>;

			bool marked;
			size_t ref;
			Polynomial polynomial;

			Expression() : marked{false}, ref(0), polynomial{} {}

			void mark() { marked = true; }
			void reference() { ++ref; }
			void dreference() { --ref; }

			ValueType d(const Expression& expr, Cache& cache) const {
				if (this == &expr) {
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
				cache[this] = dx;
				return dx;
			}

			ExprPtr addTerm(const Term& term) {
				if (term.second->marked) {
					for (auto& tm : polynomial) {
						if (tm.second == term.second) {
							tm.first += term.first;
							return this;
						}
					}
					polynomial.emplace_back(term);
					return this;
				} else {
					for (auto& tt : term.second->polynomial) {
						addTerm({term.first * tt.first, tt.second});
					}
					return this;
				}

			}

		};

		struct DataRepository {
			const static size_t bulk = 100000;
			using Datum = Expression;

			struct Element {
				Element* next;
				Datum* datum;
			};

			using Data = pmr::vector<Datum>;
			using Elements = pmr::vector<Element>;

			struct DataElement {
				Elements elements;
				Data data;
			};

			Element* remains;
			Element* unassigned;
			size_t capacity;
			size_t used;
			Datum defaultDatum;
			pmr::vector<DataElement> reserved;

			DataRepository() : remains { nullptr }, unassigned{ nullptr }, capacity{ 0 }, used{ 0 }, defaultDatum{}, reserved{} {};

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

			Datum* datum() {
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

			void dispose(Datum* datum) {
				auto current = unassigned;
				unassigned = current->next;

				*datum = defaultDatum;

				current->datum = datum;
				current->next = remains;

				remains = current;
				--used;
			}
		};

		DataRepository* repository;

		using ExprPtr = Expression*;

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
			ExprPtr expr_;

			Number(ValueType vv) : v_{ vv }, expr_{ repository->datum() }  { expression().reference(); }
			Number(ValueType vv, ExprPtr expr) : v_{ vv }, expr_{ expr }  {expression().reference(); }
			Number(const Number& other) : v_{ other.v_ }, expr_{ other.expr_ } {expression().reference(); }
			Number(const INumber& other) : v_{ other.v() }, expr_{ repository->datum() } {
				expression().reference();
				other.update(expression(), 1);
			};
			~Number() { 
				expression().dreference(); 
				if (expression().ref == 0) { repository->dispose(expr_); }
			}

			void mark() { expression().mark(); }

			ValueType v() const override { return v_; }
			void update(Expression& updated, ValueType coef) const override {
				updated.addTerm({ coef, expr_ });
			}

			ValueType d(const Number& x) const {
				Cache cache;
				return expression().d(x.expression(), cache);
			}

			Expression& expression() const { return *expr_; }

			Number operator-() const { return Number{ -v_, repository->datum()->addTerm({-1, expr_}) }; }
			Number& operator=(const Number& other) {
				if (this == &other) {
					return *this;
				}

				Number newone{ other.v() };
				other.update(newone.expression(), 1);
				this->v_ = newone.v_;
				this->expression().dreference();
				if (this->expression().ref == 0) { repository->dispose(expr_); }
				this->expr_ = newone.expr_;
				this->expression().reference();
				return *this;
			}
			Number& operator=(const INumber& other) {
				if (this == &other) {
					return *this;
				}

				Number newone{ other.v() };
				other.update(newone.expression(), 1);
				this->v_ = newone.v_;
				this->expression().dreference();
				if (this->expression().ref == 0) { repository->dispose(expr_); }
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

void hiho::ad25_exprump_tp_pclass_pzypmrmrk(double s, double sigma, double k, double r, double t, int simulation)
{
	auto dfpm = pmr::get_default_resource();
	auto pm = pmr::unsynchronized_pool_resource();
	pmr::set_default_resource(&pm);
	{
		math::DataRepository rep;
		math::repository = &rep;

		auto func = [&]() {
			Real rs{ s }; rs.mark();
			Real rsigma{ sigma }; rsigma.mark();
			Real rr{ r }; rr.mark();
			Real rt{ t }; rt.mark();
			auto value = putAmericanOption(rs, rsigma, k, rr, rt, simulation);
			return pmr::vector<Real>{ value, rs, rsigma, rr, rt };
		};

		{
			auto time = hiho::measureTime<3>(func);
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

			std::cout << math::repository->capacity << std::endl;
			std::cout << math::repository->used << std::endl;

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
	}
	pmr::set_default_resource(dfpm);
}
