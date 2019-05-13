#pragma once

#include <chrono>
#include <functional>

namespace hiho {

	template<typename T>
	struct Timer {
		using Real      = T;
		using Clock     = std::chrono::system_clock;
		using TimePoint = std::chrono::time_point<Clock>;

		TimePoint end;
		TimePoint start;
		Real      value;

		Timer(std::function<Real()> f) : end(), start(Clock::now()), value(f()) {
			end = Clock::now();
		}

		long long duration() {
			return std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		}
	};

	template<>
	struct Timer<void> {
		using Clock     = std::chrono::system_clock;
		using TimePoint = std::chrono::time_point<Clock>;

		TimePoint start;
		TimePoint end;

		Timer(std::function<void()> f) {
			start = Clock::now();
			f();
			end = Clock::now();
		}

		long long duration() {
			return std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		}
	};

	template<typename FUNC>
	inline auto newTimer(FUNC f) { return Timer<decltype(f())>(f); }

}