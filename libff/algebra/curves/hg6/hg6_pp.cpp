#include <libff/algebra/curves/hg6/hg6_pp.hpp>

namespace libff {

void hg6_pp::init_public_params()
{
    init_hg6_params();
}

hg6_GT hg6_pp::final_exponentiation(const hg6_Fq6 &elt)
{
    return hg6_final_exponentiation(elt);
}

hg6_G1_precomp hg6_pp::precompute_G1(const hg6_G1 &P)
{
    return hg6_precompute_G1(P);
}

hg6_G2_precomp hg6_pp::precompute_G2(const hg6_G2& Q, const bigint<hg6_Fq::num_limbs> &loop_count)
{
    return hg6_precompute_G2(Q, loop_count);
}

hg6_Fq6 hg6_pp::miller_loop(const hg6_G1_precomp &prec_P,
                              const hg6_G2_precomp &prec_Q_1,
                              const hg6_G2_precomp &prec_Q_2)
{
    return hg6_miller_loop(prec_P, prec_Q_1, prec_Q_2);
}

/*
hg6_affine_ate_G1_precomputation hg6_pp::affine_ate_precompute_G1(const hg6_G1 &P)
{
    return hg6_affine_ate_precompute_G1(P);
}

hg6_affine_ate_G2_precomputation hg6_pp::affine_ate_precompute_G2(const hg6_G2 &Q, const bigint<hg6_Fq::num_limbs> &loop_count)
{
    return hg6_affine_ate_precompute_G2(Q, loop_count);
}

hg6_Fq6 hg6_pp::affine_ate_miller_loop(const hg6_affine_ate_G1_precomputation &prec_P,
                                         const hg6_affine_ate_G2_precomputation &prec_Q_1,
                                         const hg6_affine_ate_G2_precomputation &prec_Q_2)
{
    return hg6_affine_ate_miller_loop(prec_P, prec_Q_1, prec_Q_2);
}

hg6_Fq6 hg6_pp::double_miller_loop(const hg6_G1_precomp &prec_P1,
                                     const hg6_G2_precomp &prec_Q1,
                                     const hg6_G1_precomp &prec_P2,
                                     const hg6_G2_precomp &prec_Q2)
{
    return hg6_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

hg6_Fq6 hg6_pp::affine_ate_e_over_e_miller_loop(const hg6_affine_ate_G1_precomputation &prec_P1,
                                                  const hg6_affine_ate_G2_precomputation &prec_Q1,
                                                  const hg6_affine_ate_G1_precomputation &prec_P2,
                                                  const hg6_affine_ate_G2_precomputation &prec_Q2)
{
    return hg6_affine_ate_miller_loop(prec_P1, prec_Q1) * hg6_affine_ate_miller_loop(prec_P2, prec_Q2).unitary_inverse();
}

hg6_Fq6 hg6_pp::affine_ate_e_times_e_over_e_miller_loop(const hg6_affine_ate_G1_precomputation &prec_P1,
                                                          const hg6_affine_ate_G2_precomputation &prec_Q1,
                                                          const hg6_affine_ate_G1_precomputation &prec_P2,
                                                          const hg6_affine_ate_G2_precomputation &prec_Q2,
                                                          const hg6_affine_ate_G1_precomputation &prec_P3,
                                                          const hg6_affine_ate_G2_precomputation &prec_Q3)
{
    return ((hg6_affine_ate_miller_loop(prec_P1, prec_Q1) * hg6_affine_ate_miller_loop(prec_P2, prec_Q2)) *
            hg6_affine_ate_miller_loop(prec_P3, prec_Q3).unitary_inverse());
}
*/

hg6_Fq6 hg6_pp::pairing(const hg6_G1 &P,
                          const hg6_G2 &Q)
{
    return hg6_pairing(P, Q);
}

hg6_Fq6 hg6_pp::reduced_pairing(const hg6_G1 &P,
                                  const hg6_G2 &Q)
{
    return hg6_reduced_pairing(P, Q);
}

/*
hg6_Fq6 hg6_pp::affine_reduced_pairing(const hg6_G1 &P,
                                         const hg6_G2 &Q)
{
    return hg6_affine_reduced_pairing(P, Q);
}
*/

} // libff
