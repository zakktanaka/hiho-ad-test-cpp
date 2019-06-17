#include"hiho/ad.hpp"
#include <string>


#define FUNC(fnum, func) if(allexec || funcnum == fnum) { func; }

using namespace hiho;

int main(int argc, char* argv[]) {
	auto allexec = argc == 1 ? true : false;
	auto funcnum = argc == 1 ? -1 : std::stoi(argv[1]);

	FUNC(  0, ad00_primitive_double           (100, 0.2, 100, 0.005, 3, 1000));
	FUNC(  1, ad01_double_struct              (100, 0.2, 100, 0.005, 3, 1000));
	FUNC(  2, ad02_struct_with_empty_vector   (100, 0.2, 100, 0.005, 3, 1000));
	FUNC(  3, ad03_struct_with_two_elem_vector(100, 0.2, 100, 0.005, 3, 1000));
	FUNC(  4, ad04_struct_with_two_elem_array (100, 0.2, 100, 0.005, 3, 1000));
	FUNC(  5, ad05_expression_vector          (100, 0.2, 100, 0.005, 3, 1000));
	FUNC(  6, ad06_expression_vector_copy     (100, 0.2, 100, 0.005, 3,   10));
	FUNC(  7, ad07_expr_vec_tape_vec          (100, 0.2, 100, 0.005, 3, 1000));
	FUNC(  8, ad08_expr_vec_tape_vec_cache    (100, 0.2, 100, 0.005, 3, 1000));
	FUNC(  9, ad09_expr_vec_tape_list_cache   (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 10, ad10_expr_vec_tape_vec_shrink   (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 11, ad11_expr_vec_tape_vec_sh_pmr   (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 12, ad12_expr_vec_tape_vec_lazy     (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 13, ad13_expr_vec_tape_vec_lazy_pmr (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 14, ad14_expr_vec_shared            (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 15, ad15_expr_vec_pmr_shared        (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 16, ad16_exprv_tapv_lazy_pmr_shrnk  (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 17, ad17_expr_vec_shared_shrink     (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 18, ad18_expr_vec_shared_pmr_shrink (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 19, ad19_expr_vec_shared_mark       (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 20, ad20_expr_vec_shrd_mrk_pmr      (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 21, ad21_exprump_tpv_lazy_pmr_shrnk (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 22, ad22_expr_umap_shrd_mrk_pmr     (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 23, ad23_exprump_tpv_lazy_pmr_mrk   (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 24, ad24_exprv_tapv_lazy_pmr_mrk    (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 25, ad25_exprump_tp_pclass_pzypmrmrk(100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 26, ad26_dual_number                (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 27, ad27_dual_number2               (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 28, ad28_struct_with_empty_umap     (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 29, ad29_struct_with_empty_map      (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 30, ad30_dual_number3               (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 31, ad31_exprvec_classic_pzypmrmrk  (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 32, ad32_refactoring_31             (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 33, ad33_template                   (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 35, ad35_refactoring_20_shrd        (100, 0.2, 100, 0.005, 3, 1000));
	FUNC( 36, ad36_template_35                (100, 0.2, 100, 0.005, 3, 1000));

	FUNC(300, ad00_primitive_double           (100, 0.2, 100, 0.005, 3, 3000));
	FUNC(308, ad08_expr_vec_tape_vec_cache    (100, 0.2, 100, 0.005, 3, 3000));
	FUNC(309, ad09_expr_vec_tape_list_cache   (100, 0.2, 100, 0.005, 3, 3000));
	FUNC(310, ad10_expr_vec_tape_vec_shrink   (100, 0.2, 100, 0.005, 3, 3000));
	FUNC(311, ad11_expr_vec_tape_vec_sh_pmr   (100, 0.2, 100, 0.005, 3, 3000));
	FUNC(312, ad12_expr_vec_tape_vec_lazy     (100, 0.2, 100, 0.005, 3, 3000));
	FUNC(313, ad13_expr_vec_tape_vec_lazy_pmr (100, 0.2, 100, 0.005, 3, 3000));
	FUNC(314, ad14_expr_vec_shared            (100, 0.2, 100, 0.005, 3, 3000));
	FUNC(315, ad15_expr_vec_pmr_shared        (100, 0.2, 100, 0.005, 3, 3000));
	FUNC(316, ad16_exprv_tapv_lazy_pmr_shrnk  (100, 0.2, 100, 0.005, 3, 3000));
	FUNC(317, ad17_expr_vec_shared_shrink     (100, 0.2, 100, 0.005, 3, 3000));
	FUNC(318, ad18_expr_vec_shared_pmr_shrink (100, 0.2, 100, 0.005, 3, 3000));
	FUNC(319, ad19_expr_vec_shared_mark       (100, 0.2, 100, 0.005, 3, 3000));
	FUNC(320, ad20_expr_vec_shrd_mrk_pmr      (100, 0.2, 100, 0.005, 3, 3000));
	//FUNC(321, ad21_exprump_tpv_lazy_pmr_shrnk (100, 0.2, 100, 0.005, 3, 3000)); // too late
	FUNC(322, ad22_expr_umap_shrd_mrk_pmr     (100, 0.2, 100, 0.005, 3, 3000));
	//FUNC(323, ad23_exprump_tpv_lazy_pmr_mrk   (100, 0.2, 100, 0.005, 3, 3000)); // too late
	FUNC(324, ad24_exprv_tapv_lazy_pmr_mrk    (100, 0.2, 100, 0.005, 3, 3000));
	//FUNC(325, ad25_exprump_tp_pclass_pzypmrmrk(100, 0.2, 100, 0.005, 3, 3000)); 
	FUNC(326, ad26_dual_number                (100, 0.2, 100, 0.005, 3, 3000));
	//FUNC(327, ad27_dual_number2               (100, 0.2, 100, 0.005, 3, 3000)); // too late 
	FUNC(331, ad31_exprvec_classic_pzypmrmrk  (100, 0.2, 100, 0.005, 3, 3000));
	FUNC(332, ad32_refactoring_31             (100, 0.2, 100, 0.005, 3, 3000));
	FUNC(333, ad33_template                   (100, 0.2, 100, 0.005, 3, 3000));
	FUNC(335, ad35_refactoring_20_shrd        (100, 0.2, 100, 0.005, 3, 3000));
}