#include <Test.hpp>

#include <cmath>
#include <memory>
#include <vector>
#include <unordered_map>

namespace { namespace sandbox {

	template<typename T>
	class Expression {
	public:
		using ExprValueType = T;
		using ThisType = Expression<ExprValueType>;
		using Cache    = std::unordered_map<const void*, ExprValueType>;
		using ExpPtr   = std::shared_ptr<ThisType>;

	private:
		using Term       = std::pair<ExpPtr, ExprValueType>;
		using Polynomial = std::vector<Term>;

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
			}
			else {
				for (auto& otm : other->polynomial) {
					ExprValueType c = coef * otm.second;
					auto& e = otm.first;
					addExpresison(c, e);
				}
			}
		}
	};

	template<typename T>
	struct INumber {
		using ValueType = T;
		using Expr      = Expression<ValueType>;
		using ThisType  = INumber<ValueType>;

		virtual ~INumber() {}
		virtual ValueType  v() const = 0;
		virtual void update(Expr&, ValueType) const = 0;
	};

	template<typename T>
	struct Unary : INumber<T> {
		using BaseType  = INumber<T>;
		using ValueType = typename BaseType::ValueType;
		using Expr      = typename BaseType::Expr;
		using ThisType  = Unary<ValueType>;

		ValueType v_;
		ValueType coef_;
		const BaseType& num_;

		Unary(ValueType vv, ValueType coef, const BaseType& num) :
			v_(vv), coef_(coef), num_(num) {}

		ValueType v() const override { return v_; }

		void update(Expr& updated, ValueType coef) const override {
			num_.update(updated, coef * coef_);
		}
	};

	template<typename T>
	struct Binary : INumber<T> {
		using BaseType  = INumber<T>;
		using ValueType = typename BaseType::ValueType;
		using Expr      = typename BaseType::Expr;
		using ThisType  = Binary<ValueType>;

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

		void update(Expr& updated, ValueType coef) const override {
			lnum_.update(updated, coef * lcoef_);
			rnum_.update(updated, coef * rcoef_);
		}
	};

	template<typename T>
	class Number : public INumber<T> {
		using BaseType  = INumber<T>;
		using ValueType = typename BaseType::ValueType;
		using Expr      = typename BaseType::Expr;
		using ThisType  = Number<ValueType>;

	private:
		using Cache  = typename Expr::Cache;
		using ExpPtr = typename Expr::ExpPtr;

		ValueType  v_;
		ExpPtr     expr_;

		Number(ValueType vv, ExpPtr expr) : v_{ vv }, expr_{ expr }  {}

	public:
		Number(ValueType vv) :
			ThisType{ vv , std::make_shared<Expr>() } {}
		Number(const BaseType& other) :
			ThisType{ other.v(), std::make_shared<Expr>() } {
			other.update(*expr_, 1);
		};
		~Number() {}

		void mark() { expr_->mark(); }

		ValueType v() const override { return v_; }

		void update(Expr& updated, ValueType coef) const override {
			updated.addExpresison(coef, expr_);
		}

		ValueType d(ThisType& x) const {
			Cache cache;
			return expr_->d(x.expr_, cache);
		}

		ThisType& operator=(const BaseType& other) {
			ThisType newone{ other.v() };
			static ValueType one{ 1 };
			other.update(*(newone.expr_), one);
			this->v_ = newone.v_;
			this->expr_ = newone.expr_;
			return *this;
		}
	};

	template<typename T> Binary<T> operator*(const INumber<T>& l, const INumber<T>& r) { return Binary<T>{ l.v() * r.v(), r.v(), l, l.v(), r }; }

	using std::exp;

	template<typename T>
	Unary<T> exp(const INumber<T>& l) {
		auto ll = exp(l.v());
		return Unary<T>{ ll, ll, l };
	}

} }

TEST(Sandbox, TestName) {
	using namespace sandbox;
	using Real = Number<double>;

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