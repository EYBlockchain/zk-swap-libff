#include <libff/algebra/curves/test_curve/test_curve_pp.hpp>

namespace libff {

void test_curve_pp::init_public_params()
{
    init_test_curve_params();
}

test_curve_GT test_curve_pp::final_exponentiation(const test_curve_Fq12 &elt)
{
    return test_curve_final_exponentiation(elt);
}

test_curve_G1_precomp test_curve_pp::precompute_G1(const test_curve_G1 &P)
{
    return test_curve_precompute_G1(P);
}

test_curve_G2_precomp test_curve_pp::precompute_G2(const test_curve_G2 &Q)
{
    return test_curve_precompute_G2(Q);
}

test_curve_Fq12 test_curve_pp::miller_loop(const test_curve_G1_precomp &prec_P,
                                         const test_curve_G2_precomp &prec_Q)
{
    return test_curve_miller_loop(prec_P, prec_Q);
}

test_curve_Fq12 test_curve_pp::double_miller_loop(const test_curve_G1_precomp &prec_P1,
                                                const test_curve_G2_precomp &prec_Q1,
                                                const test_curve_G1_precomp &prec_P2,
                                                const test_curve_G2_precomp &prec_Q2)
{
    return test_curve_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

test_curve_Fq12 test_curve_pp::pairing(const test_curve_G1 &P,
                                     const test_curve_G2 &Q)
{
    return test_curve_pairing(P, Q);
}

test_curve_Fq12 test_curve_pp::reduced_pairing(const test_curve_G1 &P,
                                             const test_curve_G2 &Q)
{
    return test_curve_reduced_pairing(P, Q);
}

} // libff
