#include <libff/algebra/curves/sw6_bis/sw6_bis_pp.hpp>

namespace libff {

void sw6_bis_pp::init_public_params()
{
    init_sw6_bis_params();
}

sw6_bis_GT sw6_bis_pp::final_exponentiation(const sw6_bis_Fq6 &elt)
{
    return sw6_bis_final_exponentiation(elt);
}

sw6_bis_G1_precomp sw6_bis_pp::precompute_G1(const sw6_bis_G1 &P)
{
    return sw6_bis_precompute_G1(P);
}

sw6_bis_G2_precomp sw6_bis_pp::precompute_G2(const sw6_bis_G2 &Q)
{
    return sw6_bis_precompute_G2(Q);
}


sw6_bis_Fq6 sw6_bis_pp::miller_loop(const sw6_bis_G1_precomp &prec_P,
                              const sw6_bis_G2_precomp &prec_Q)
{
    return sw6_bis_miller_loop(prec_P, prec_Q);
}

sw6_bis_affine_ate_G1_precomputation sw6_bis_pp::affine_ate_precompute_G1(const sw6_bis_G1 &P)
{
    return sw6_bis_affine_ate_precompute_G1(P);
}

sw6_bis_affine_ate_G2_precomputation sw6_bis_pp::affine_ate_precompute_G2(const sw6_bis_G2 &Q)
{
    return sw6_bis_affine_ate_precompute_G2(Q);
}

sw6_bis_Fq6 sw6_bis_pp::affine_ate_miller_loop(const sw6_bis_affine_ate_G1_precomputation &prec_P,
                                         const sw6_bis_affine_ate_G2_precomputation &prec_Q)
{
    return sw6_bis_affine_ate_miller_loop(prec_P, prec_Q);
}

sw6_bis_Fq6 sw6_bis_pp::double_miller_loop(const sw6_bis_G1_precomp &prec_P1,
                                     const sw6_bis_G2_precomp &prec_Q1,
                                     const sw6_bis_G1_precomp &prec_P2,
                                     const sw6_bis_G2_precomp &prec_Q2)
{
    return sw6_bis_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

sw6_bis_Fq6 sw6_bis_pp::affine_ate_e_over_e_miller_loop(const sw6_bis_affine_ate_G1_precomputation &prec_P1,
                                                  const sw6_bis_affine_ate_G2_precomputation &prec_Q1,
                                                  const sw6_bis_affine_ate_G1_precomputation &prec_P2,
                                                  const sw6_bis_affine_ate_G2_precomputation &prec_Q2)
{
    return sw6_bis_affine_ate_miller_loop(prec_P1, prec_Q1) * sw6_bis_affine_ate_miller_loop(prec_P2, prec_Q2).unitary_inverse();
}

sw6_bis_Fq6 sw6_bis_pp::affine_ate_e_times_e_over_e_miller_loop(const sw6_bis_affine_ate_G1_precomputation &prec_P1,
                                                          const sw6_bis_affine_ate_G2_precomputation &prec_Q1,
                                                          const sw6_bis_affine_ate_G1_precomputation &prec_P2,
                                                          const sw6_bis_affine_ate_G2_precomputation &prec_Q2,
                                                          const sw6_bis_affine_ate_G1_precomputation &prec_P3,
                                                          const sw6_bis_affine_ate_G2_precomputation &prec_Q3)
{
    return ((sw6_bis_affine_ate_miller_loop(prec_P1, prec_Q1) * sw6_bis_affine_ate_miller_loop(prec_P2, prec_Q2)) *
            sw6_bis_affine_ate_miller_loop(prec_P3, prec_Q3).unitary_inverse());
}

sw6_bis_Fq6 sw6_bis_pp::pairing(const sw6_bis_G1 &P,
                          const sw6_bis_G2 &Q)
{
    return sw6_bis_pairing(P, Q);
}

sw6_bis_Fq6 sw6_bis_pp::reduced_pairing(const sw6_bis_G1 &P,
                                  const sw6_bis_G2 &Q)
{
    return sw6_bis_reduced_pairing(P, Q);
}

sw6_bis_Fq6 sw6_bis_pp::affine_reduced_pairing(const sw6_bis_G1 &P,
                                         const sw6_bis_G2 &Q)
{
    return sw6_bis_affine_reduced_pairing(P, Q);
}

} // libff
