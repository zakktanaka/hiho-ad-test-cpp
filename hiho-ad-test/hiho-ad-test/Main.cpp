#include"hiho/ad.hpp"

using namespace hiho;

int main() {
	ad00_primitive_double(100, 0.2, 100, 0.005, 3, 1000);
	ad01_double_struct   (100, 0.2, 100, 0.005, 3, 1000);
}