#include <Test.hpp>

int main(int argc, char** argv) {

	auto logger = spdlog::stdout_color_mt("LOG");
	spdlog::set_default_logger(logger);

	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}