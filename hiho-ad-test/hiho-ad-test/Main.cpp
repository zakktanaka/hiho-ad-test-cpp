#include"hiho/ad.hpp"

using namespace hiho;

int main() {
	ad00_primitive_double           (100, 0.2, 100, 0.005, 3, 1000);
	ad01_double_struct              (100, 0.2, 100, 0.005, 3, 1000);
	ad02_struct_with_empty_vector   (100, 0.2, 100, 0.005, 3, 1000);
	ad03_struct_with_two_elem_vector(100, 0.2, 100, 0.005, 3, 1000);
	ad04_struct_with_two_elem_array (100, 0.2, 100, 0.005, 3, 1000);
	ad05_expression_vector          (100, 0.2, 100, 0.005, 3, 1000);
	ad06_expression_vector_copy     (100, 0.2, 100, 0.005, 3,   10);
	ad07_expr_vec_tape_vec          (100, 0.2, 100, 0.005, 3, 1000);
	ad08_expr_vec_tape_vec_cache    (100, 0.2, 100, 0.005, 3, 1000);
}