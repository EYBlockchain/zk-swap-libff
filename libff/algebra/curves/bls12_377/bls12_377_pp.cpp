/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/bls12_377/bls12_377_pp.hpp>

#include <libff/algebra/curves/bls12.tcc>


namespace libff {


void bls12_377_pp::init_public_params()
{
    init_bls12_377_params();
}


bls12_377_pp::G1_precomp_type bls12_377_pp::precompute_G1(const bls12_377_G1 &P)
{
    return P;
}


bls12_377_pp::G2_precomp_type bls12_377_pp::precompute_G2(const bls12_377_G2 &Q)
{
    return bls12::G2Prepared<bls12_377_pp>(Q);
}


bls12_377_Fq12 bls12_377_pp::miller_loop(const bls12_377_pp::G1_precomp_type &prec_P,
                                         const bls12_377_pp::G2_precomp_type &prec_Q)
{
    return bls12::miller_loop<bls12_377_pp>({
        bls12::PreparedPair<bls12_377_pp>(prec_P, prec_Q)
    });
}


bls12_377_Fq12 bls12_377_pp::double_miller_loop(const G1_precomp_type &prec_P1,
                                                const G2_precomp_type &prec_Q1,
                                                const G1_precomp_type &prec_P2,
                                                const G2_precomp_type &prec_Q2)
{
    return bls12::miller_loop<bls12_377_pp>({
        bls12::PreparedPair<bls12_377_pp>(prec_P1, prec_Q1),
        bls12::PreparedPair<bls12_377_pp>(prec_P2, prec_Q2)
    });
}


bls12_377_Fq12 bls12_377_pp::final_exponentiation(const bls12_377_Fq12 &elt)
{
    return bls12::final_exponentiation<bls12_377_pp>(elt);
}


bls12_377_Fq12 bls12_377_pp::pairing(const bls12_377_G1 &P,
                                     const bls12_377_G2 &Q)
{
    return bls12::miller_loop<bls12_377_pp>(P, Q);
}


bls12_377_Fq12 bls12_377_pp::reduced_pairing(const bls12_377_G1 &P,
                                             const bls12_377_G2 &Q)
{
    return final_exponentiation(pairing(P, Q));
}


} // libff
