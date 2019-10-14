/** @file
 *****************************************************************************

 Implementation of interfaces for public parameters of SW6.

 See sw6_pp.hpp .

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/sw6/sw6_pp.hpp>

namespace libff {

void sw6_pp::init_public_params()
{
    init_sw6_params();
}

sw6_GT sw6_pp::final_exponentiation(const sw6_Fq6 &elt)
{
    return sw6_final_exponentiation(elt);
}

sw6_G1_precomp sw6_pp::precompute_G1(const sw6_G1 &P)
{
    return sw6_precompute_G1(P);
}

sw6_G2_precomp sw6_pp::precompute_G2(const sw6_G2 &Q)
{
    return sw6_precompute_G2(Q);
}


sw6_Fq6 sw6_pp::miller_loop(const sw6_G1_precomp &prec_P,
                              const sw6_G2_precomp &prec_Q)
{
    return sw6_miller_loop(prec_P, prec_Q);
}

sw6_affine_ate_G1_precomputation sw6_pp::affine_ate_precompute_G1(const sw6_G1 &P)
{
    return sw6_affine_ate_precompute_G1(P);
}

sw6_affine_ate_G2_precomputation sw6_pp::affine_ate_precompute_G2(const sw6_G2 &Q)
{
    return sw6_affine_ate_precompute_G2(Q);
}

sw6_Fq6 sw6_pp::affine_ate_miller_loop(const sw6_affine_ate_G1_precomputation &prec_P,
                                         const sw6_affine_ate_G2_precomputation &prec_Q)
{
    return sw6_affine_ate_miller_loop(prec_P, prec_Q);
}

sw6_Fq6 sw6_pp::double_miller_loop(const sw6_G1_precomp &prec_P1,
                                     const sw6_G2_precomp &prec_Q1,
                                     const sw6_G1_precomp &prec_P2,
                                     const sw6_G2_precomp &prec_Q2)
{
    return sw6_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

sw6_Fq6 sw6_pp::affine_ate_e_over_e_miller_loop(const sw6_affine_ate_G1_precomputation &prec_P1,
                                                  const sw6_affine_ate_G2_precomputation &prec_Q1,
                                                  const sw6_affine_ate_G1_precomputation &prec_P2,
                                                  const sw6_affine_ate_G2_precomputation &prec_Q2)
{
    return sw6_affine_ate_miller_loop(prec_P1, prec_Q1) * sw6_affine_ate_miller_loop(prec_P2, prec_Q2).unitary_inverse();
}

sw6_Fq6 sw6_pp::affine_ate_e_times_e_over_e_miller_loop(const sw6_affine_ate_G1_precomputation &prec_P1,
                                                          const sw6_affine_ate_G2_precomputation &prec_Q1,
                                                          const sw6_affine_ate_G1_precomputation &prec_P2,
                                                          const sw6_affine_ate_G2_precomputation &prec_Q2,
                                                          const sw6_affine_ate_G1_precomputation &prec_P3,
                                                          const sw6_affine_ate_G2_precomputation &prec_Q3)
{
    return ((sw6_affine_ate_miller_loop(prec_P1, prec_Q1) * sw6_affine_ate_miller_loop(prec_P2, prec_Q2)) *
            sw6_affine_ate_miller_loop(prec_P3, prec_Q3).unitary_inverse());
}

sw6_Fq6 sw6_pp::pairing(const sw6_G1 &P,
                          const sw6_G2 &Q)
{
    return sw6_pairing(P, Q);
}

sw6_Fq6 sw6_pp::reduced_pairing(const sw6_G1 &P,
                                  const sw6_G2 &Q)
{
    return sw6_reduced_pairing(P, Q);
}

sw6_Fq6 sw6_pp::affine_reduced_pairing(const sw6_G1 &P,
                                         const sw6_G2 &Q)
{
    return sw6_affine_reduced_pairing(P, Q);
}

} // libff
