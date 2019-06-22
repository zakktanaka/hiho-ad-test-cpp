#include <Test.hpp>

#include <cmath>
#include <memory>
#include <vector>
#include <unordered_map>

namespace { namespace sandbox {

	using ValueType = double;

	class Expression {
	public:
		using Cache      = std::unordered_map<const void*, ValueType>;
		using ExpPtr     = std::shared_ptr<Expression>;
		using Term       = std::pair<ExpPtr, ValueType>;
		using Polynomial = std::vector<Term>;

	private:
		bool marked;
		Polynomial polynomial;

		void addTerm(ValueType coef, const ExpPtr& other) {
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

		void addExpresison(ValueType coef, const ExpPtr& other) {
			if (other->marked) {
				addTerm(coef, other);
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
		virtual ~INumber() {}
		virtual ValueType  v() const = 0;
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

	class Number : public INumber {
	private:
		using Cache  = typename Expression::Cache;
		using ExpPtr = typename Expression::ExpPtr;

		ValueType  v_;
		ExpPtr     expr_;

		Number(ValueType vv, ExpPtr expr) : v_{ vv }, expr_{ expr }  {}

	public:
		Number(ValueType vv) :
			Number{ vv , std::make_shared<Expression>() } {}
		Number(const INumber& other) :
			Number{ other.v(), std::make_shared<Expression>() } {
			other.update(*expr_, 1);
		};
		~Number() {}

		void mark() { expr_->mark(); }

		ValueType v() const override { return v_; }

		void update(Expression& updated, ValueType coef) const override {
			updated.addExpresison(coef, expr_);
		}

		ValueType d(Number& x) const {
			Cache cache;
			return expr_->d(x.expr_, cache);
		}

		Number& operator=(const INumber& other) {
			Number newone{ other.v() };
			other.update(*(newone.expr_), 1);
			this->v_ = newone.v_;
			this->expr_ = newone.expr_;
			return *this;
		}
	};

	Binary operator*(const INumber& l, const INumber& r) { return Binary{ l.v() * r.v(), r.v(), l, l.v(), r }; }

	using std::exp;

	Unary exp(const INumber& l) {
		auto ll = exp(l.v());
		return Unary{ ll, ll, l };
	}

} }

namespace {
	using namespace sandbox;

	using Number = sandbox::Number;
}

TEST(Sandbox, TestName) {
	using Real = Number;

	constexpr double err  = 1e-12;
	constexpr int    loop = 1000000;
	constexpr double xx   = 2.0 / loop;

	Real x{ xx }; x.mark();

	auto func = [&]() {
		Real ans = exp(Real{ -1 });
		for (int i = 0; i < loop; ++i) {
			ans = ans * exp(x);
		}

		return ans;
	};
	auto actual = hiho::test::newTimer(func);

	auto expectfunc = [&]() {
		auto ans = exp(-1);
		for (int i = 0; i < loop; ++i) {
			ans = ans * exp(xx);
		}
		return ans;
	};
	auto expected = hiho::test::newTimer(expectfunc);

	EXPECT_LE  (actual.duration(),     expected.duration() * 400);
	EXPECT_NEAR(expected.value,        actual.value.v(),  err);
	EXPECT_NEAR(expected.value * loop, actual.value.d(x), err * loop);
}