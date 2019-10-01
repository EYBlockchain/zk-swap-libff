/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/bls12_377/bls12_377_pp.hpp>

namespace libff {

void bls12_377_pp::init_public_params()
{
    init_bls12_377_params();
}

bls12_377_GT bls12_377_pp::final_exponentiation(const bls12_377_Fq12 &elt)
{
    return bls12_377_final_exponentiation(elt);
}

bls12_377_G1_precomp bls12_377_pp::precompute_G1(const bls12_377_G1 &P)
{
    return bls12_377_precompute_G1(P);
}

bls12_377_G2_precomp bls12_377_pp::precompute_G2(const bls12_377_G2 &Q)
{
    return bls12_377_precompute_G2(Q);
}

bls12_377_Fq12 bls12_377_pp::miller_loop(const bls12_377_G1_precomp &prec_P,
                                         const bls12_377_G2_precomp &prec_Q)
{
    return bls12_377_miller_loop(prec_P, prec_Q);
}

bls12_377_Fq12 bls12_377_pp::double_miller_loop(const bls12_377_G1_precomp &prec_P1,
                                                const bls12_377_G2_precomp &prec_Q1,
                                                const bls12_377_G1_precomp &prec_P2,
                                                const bls12_377_G2_precomp &prec_Q2)
{
    return bls12_377_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

bls12_377_Fq12 bls12_377_pp::pairing(const bls12_377_G1 &P,
                                     const bls12_377_G2 &Q)
{
    return bls12_377_pairing(P, Q);
}

bls12_377_Fq12 bls12_377_pp::reduced_pairing(const bls12_377_G1 &P,
                                             const bls12_377_G2 &Q)
{
    return bls12_377_reduced_pairing(P, Q);
}

} // libff
