/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/bw12_446/bw12_446_pp.hpp>

namespace libff {

void bw12_446_pp::init_public_params()
{
    init_bw12_446_params();
}

bw12_446_GT bw12_446_pp::final_exponentiation(const bw12_446_Fq12 &elt)
{
    return bw12_446_final_exponentiation(elt);
}

bw12_446_G1_precomp bw12_446_pp::precompute_G1(const bw12_446_G1 &P)
{
    return bw12_446_precompute_G1(P);
}

bw12_446_G2_precomp bw12_446_pp::precompute_G2(const bw12_446_G2 &Q)
{
    return bw12_446_precompute_G2(Q);
}

bw12_446_Fq12 bw12_446_pp::miller_loop(const bw12_446_G1_precomp &prec_P,
                                         const bw12_446_G2_precomp &prec_Q)
{
    return bw12_446_miller_loop(prec_P, prec_Q);
}

bw12_446_Fq12 bw12_446_pp::double_miller_loop(const bw12_446_G1_precomp &prec_P1,
                                                const bw12_446_G2_precomp &prec_Q1,
                                                const bw12_446_G1_precomp &prec_P2,
                                                const bw12_446_G2_precomp &prec_Q2)
{
    return bw12_446_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

bw12_446_Fq12 bw12_446_pp::pairing(const bw12_446_G1 &P,
                                     const bw12_446_G2 &Q)
{
    return bw12_446_pairing(P, Q);
}

bw12_446_Fq12 bw12_446_pp::reduced_pairing(const bw12_446_G1 &P,
                                             const bw12_446_G2 &Q)
{
    return bw12_446_reduced_pairing(P, Q);
}

} // libff
