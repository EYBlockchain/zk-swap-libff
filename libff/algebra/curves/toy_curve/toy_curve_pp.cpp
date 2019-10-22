#include <libff/algebra/curves/toy_curve/toy_curve_pp.hpp>

namespace libff {

void toy_curve_pp::init_public_params()
{
    init_toy_curve_params();
}

toy_curve_GT toy_curve_pp::final_exponentiation(const toy_curve_Fq12 &elt)
{
    return toy_curve_final_exponentiation(elt);
}

toy_curve_G1_precomp toy_curve_pp::precompute_G1(const toy_curve_G1 &P)
{
    return toy_curve_precompute_G1(P);
}

toy_curve_G2_precomp toy_curve_pp::precompute_G2(const toy_curve_G2 &Q)
{
    return toy_curve_precompute_G2(Q);
}

toy_curve_Fq12 toy_curve_pp::miller_loop(const toy_curve_G1_precomp &prec_P,
                                         const toy_curve_G2_precomp &prec_Q)
{
    return toy_curve_miller_loop(prec_P, prec_Q);
}

toy_curve_Fq12 toy_curve_pp::double_miller_loop(const toy_curve_G1_precomp &prec_P1,
                                                const toy_curve_G2_precomp &prec_Q1,
                                                const toy_curve_G1_precomp &prec_P2,
                                                const toy_curve_G2_precomp &prec_Q2)
{
    return toy_curve_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

toy_curve_Fq12 toy_curve_pp::pairing(const toy_curve_G1 &P,
                                     const toy_curve_G2 &Q)
{
    return toy_curve_pairing(P, Q);
}

toy_curve_Fq12 toy_curve_pp::reduced_pairing(const toy_curve_G1 &P,
                                             const toy_curve_G2 &Q)
{
    return toy_curve_reduced_pairing(P, Q);
}

} // libff
