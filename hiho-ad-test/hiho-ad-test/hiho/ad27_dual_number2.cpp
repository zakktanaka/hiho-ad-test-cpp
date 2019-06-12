#include "ad27_dual_number2.hpp"

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
	using Differentials = pmr::unordered_map<size_t, ValueType>;

	namespace math {

		struct INumber {
			virtual ValueType  v() const = 0;
			virtual ~INumber() {}
			virtual void update(Differentials&, ValueType) const = 0;

		};

		struct Unary : INumber {
			ValueType v_;
			ValueType coef_;
			const INumber& num_;

			Unary(ValueType vv, ValueType coef, const INumber& num) :
				v_(vv), coef_(coef), num_(num) {}

			ValueType v() const override { return v_; }
			void update(Differentials& updated, ValueType coef) const override {
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
			void update(Differentials& updated, ValueType coef) const override {
				lnum_.update(updated, coef * lcoef_);
				rnum_.update(updated, coef * rcoef_);
			}
		};

		struct Number : public INumber {
			size_t id_;
			ValueType  v_;
			Differentials differentials_;

			Number() : id_{ indexer++ }, v_{ 0 }, differentials_{}  {}
			Number(ValueType vv) : id_{ indexer++ }, v_{ vv }, differentials_{}  { differentials_[id_] = 1; }
			Number(const Number&) = default;
			Number(Number&&) noexcept = default;
			Number(const INumber& other) : id_{ indexer++ }, v_{ other.v() }, differentials_{} {
				other.update(differentials_, 1);
			};

			ValueType v() const override { return v_; }
			void update(Differentials& updated, ValueType coef) const override {
				for (auto& p : differentials_) {
					auto it = updated.find(p.first);
					if (it != std::end(updated)) {
						it->second += coef * p.second;
					}
					else {
						updated[p.first] = coef * p.second;
					}
				}
			}

			Number& update(ValueType coef) {
				for (auto& p : differentials_) {
					p.second *= coef;
				}
				return *this;
			}

			ValueType d(const Number& x) const {
				auto it = differentials_.find(x.id_);
				if (it != std::end(differentials_)) {
					return it->second;
				}
				else {
					return 0;
				}
			}

			Number operator-() const { return Number{ *this }.update(-1); }
			Number& operator=(Number&&) = default;
			Number& operator=(const Number& other) {
				if (this == &other) {
					return *this;
				}
				Number newone;
				newone.v_ = other.v();
				other.update(newone.differentials_, 1);
				*this = std::move(newone);
				return *this;
			}
			Number& operator=(const INumber& other) {
				if (this == &other) {
					return *this;
				}
				Number newone;
				newone.v_ = other.v();
				other.update(newone.differentials_, 1);
				*this = std::move(newone);
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

void hiho::ad27_dual_number2(double s, double sigma, double k, double r, double t, int simulation)
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
			return pmr::vector<Real>{ value, rs, rsigma, rr, rt };
		};

		{
			auto time = hiho::measureTime<1>(func);
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
	}
	pmr::set_default_resource(dfpm);
}
