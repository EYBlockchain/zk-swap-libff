#ifndef TEST_CURVE_PP_HPP_
#define TEST_CURVE_PP_HPP_
#include <libff/algebra/curves/test_curve/test_curve_g1.hpp>
#include <libff/algebra/curves/test_curve/test_curve_g2.hpp>
#include <libff/algebra/curves/test_curve/test_curve_init.hpp>
#include <libff/algebra/curves/test_curve/test_curve_pairing.hpp>
#include <libff/algebra/curves/public_params.hpp>

namespace libff {

class test_curve_pp {
public:
    typedef test_curve_Fr Fp_type;
    typedef test_curve_G1 G1_type;
    typedef test_curve_G2 G2_type;
    typedef test_curve_G1_precomp G1_precomp_type;
    typedef test_curve_G2_precomp G2_precomp_type;
    typedef test_curve_Fq Fq_type;
    typedef test_curve_Fq2 Fqe_type;
    typedef test_curve_Fq12 Fqk_type;
    typedef test_curve_GT GT_type;

    static const bool has_affine_pairing = false;

    static void init_public_params();
    static test_curve_GT final_exponentiation(const test_curve_Fq12 &elt);
    static test_curve_G1_precomp precompute_G1(const test_curve_G1 &P);
    static test_curve_G2_precomp precompute_G2(const test_curve_G2 &Q);
    static test_curve_Fq12 miller_loop(const test_curve_G1_precomp &prec_P,
                                      const test_curve_G2_precomp &prec_Q);
    static test_curve_Fq12 double_miller_loop(const test_curve_G1_precomp &prec_P1,
                                             const test_curve_G2_precomp &prec_Q1,
                                             const test_curve_G1_precomp &prec_P2,
                                             const test_curve_G2_precomp &prec_Q2);
    static test_curve_Fq12 pairing(const test_curve_G1 &P,
                                  const test_curve_G2 &Q);
    static test_curve_Fq12 reduced_pairing(const test_curve_G1 &P,
                                          const test_curve_G2 &Q);
};

} // libff

#endif // TEST_CURVE_PP_HPP_
