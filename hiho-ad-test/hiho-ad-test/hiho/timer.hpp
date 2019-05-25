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

	namespace timer {
		template<int TIMES, typename FUNC>
		inline long long measureTime(FUNC f, std::function<void()> preprocess, std::function<void()> postprocess) {
			constexpr int loop = TIMES;

			long long time = 0;
			long long min = std::numeric_limits<long long>::max();
			long long max = 0;
			for (int i = 0; i != loop + 2; ++i) {
				preprocess();
				auto tm = hiho::newTimer(f).duration();
				postprocess();
				if (tm < min) { min = tm; }
				if (tm > max) { max = tm; }
				time += tm;
			}
			time -= min;
			time -= max;
			time /= loop;

			return time;
		}
	}

	template<typename FUNC>
	inline long long measureTime(FUNC f) {
		return timer::measureTime<10>(f, []() {}, []() {});
	}

	template<typename FUNC>
	inline long long measureTime(FUNC f, std::function<void()> postprocess) {
		return timer::measureTime<10>(f, []() {}, postprocess);
	}

	template<int TIMES, typename FUNC>
	inline long long measureTime(FUNC f) {
		return timer::measureTime<TIMES>(f, []() {}, []() {});
	}

	template<int TIMES, typename FUNC>
	inline long long measureTime(FUNC f, std::function<void()> postprocess = []() {}) {
		return timer::measureTime<TIMES>(f, []() {}, postprocess);
	}

}