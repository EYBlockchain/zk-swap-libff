/** @file
 *****************************************************************************

 Implementation of interfaces for public parameters of MNT6.

 See pendulum_pp.hpp .

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/pendulum/pendulum_pp.hpp>

namespace libff {

void pendulum_pp::init_public_params()
{
    init_pendulum_params();
}

pendulum_GT pendulum_pp::final_exponentiation(const pendulum_Fq6 &elt)
{
    return pendulum_final_exponentiation(elt);
}

pendulum_G1_precomp pendulum_pp::precompute_G1(const pendulum_G1 &P)
{
    return pendulum_precompute_G1(P);
}

pendulum_G2_precomp pendulum_pp::precompute_G2(const pendulum_G2 &Q)
{
    return pendulum_precompute_G2(Q);
}


pendulum_Fq6 pendulum_pp::miller_loop(const pendulum_G1_precomp &prec_P,
                              const pendulum_G2_precomp &prec_Q)
{
    return pendulum_miller_loop(prec_P, prec_Q);
}

pendulum_affine_ate_G1_precomputation pendulum_pp::affine_ate_precompute_G1(const pendulum_G1 &P)
{
    return pendulum_affine_ate_precompute_G1(P);
}

pendulum_affine_ate_G2_precomputation pendulum_pp::affine_ate_precompute_G2(const pendulum_G2 &Q)
{
    return pendulum_affine_ate_precompute_G2(Q);
}

pendulum_Fq6 pendulum_pp::affine_ate_miller_loop(const pendulum_affine_ate_G1_precomputation &prec_P,
                                         const pendulum_affine_ate_G2_precomputation &prec_Q)
{
    return pendulum_affine_ate_miller_loop(prec_P, prec_Q);
}

pendulum_Fq6 pendulum_pp::double_miller_loop(const pendulum_G1_precomp &prec_P1,
                                     const pendulum_G2_precomp &prec_Q1,
                                     const pendulum_G1_precomp &prec_P2,
                                     const pendulum_G2_precomp &prec_Q2)
{
    return pendulum_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

pendulum_Fq6 pendulum_pp::affine_ate_e_over_e_miller_loop(const pendulum_affine_ate_G1_precomputation &prec_P1,
                                                  const pendulum_affine_ate_G2_precomputation &prec_Q1,
                                                  const pendulum_affine_ate_G1_precomputation &prec_P2,
                                                  const pendulum_affine_ate_G2_precomputation &prec_Q2)
{
    return pendulum_affine_ate_miller_loop(prec_P1, prec_Q1) * pendulum_affine_ate_miller_loop(prec_P2, prec_Q2).unitary_inverse();
}

pendulum_Fq6 pendulum_pp::affine_ate_e_times_e_over_e_miller_loop(const pendulum_affine_ate_G1_precomputation &prec_P1,
                                                          const pendulum_affine_ate_G2_precomputation &prec_Q1,
                                                          const pendulum_affine_ate_G1_precomputation &prec_P2,
                                                          const pendulum_affine_ate_G2_precomputation &prec_Q2,
                                                          const pendulum_affine_ate_G1_precomputation &prec_P3,
                                                          const pendulum_affine_ate_G2_precomputation &prec_Q3)
{
    return ((pendulum_affine_ate_miller_loop(prec_P1, prec_Q1) * pendulum_affine_ate_miller_loop(prec_P2, prec_Q2)) *
            pendulum_affine_ate_miller_loop(prec_P3, prec_Q3).unitary_inverse());
}

pendulum_Fq6 pendulum_pp::pairing(const pendulum_G1 &P,
                          const pendulum_G2 &Q)
{
    return pendulum_pairing(P, Q);
}

pendulum_Fq6 pendulum_pp::reduced_pairing(const pendulum_G1 &P,
                                  const pendulum_G2 &Q)
{
    return pendulum_reduced_pairing(P, Q);
}

pendulum_Fq6 pendulum_pp::affine_reduced_pairing(const pendulum_G1 &P,
                                         const pendulum_G2 &Q)
{
    return pendulum_affine_reduced_pairing(P, Q);
}

} // libff
