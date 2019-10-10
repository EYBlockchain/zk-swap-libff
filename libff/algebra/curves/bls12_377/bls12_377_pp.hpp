/** @file
*****************************************************************************
* @author     This file is part of libff, developed by SCIPR Lab
*             and contributors (see AUTHORS).
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/

#ifndef BLS12_377_PP_HPP_
#define BLS12_377_PP_HPP_
#include <libff/algebra/curves/bls12_377/bls12_377_g1.hpp>
#include <libff/algebra/curves/bls12_377/bls12_377_g2.hpp>
#include <libff/algebra/curves/bls12_377/bls12_377_init.hpp>
#include <libff/algebra/curves/bls12.hpp>
#include <libff/algebra/curves/public_params.hpp>

namespace libff {

class bls12_377_pp {
public:
    typedef bls12_377_Fr Fp_type;
    typedef bls12_377_G1 G1_type;
    typedef bls12_377_G2 G2_type;
    typedef G1_type G1_precomp_type;
    typedef bls12::G2Prepared<bls12_377_pp> G2_precomp_type;
    typedef bls12_377_Fq Fq_type;
    typedef bls12_377_Fq2 Fqe_type;
    typedef bls12_377_Fq12 Fqk_type;
    typedef bls12_377_GT GT_type;

    static const bool has_affine_pairing = false;

    /* 
     * https://eprint.iacr.org/2017/1174.pdf
     * ate_loop_count=t where t in the value chosen in q(t) and r(t) parameterization of BLS curve.
     * for BLS12_377 t=3·2^46·(7·13·499)+1 s.t. t=1 (mod 3·2^46) to have a high 2-adicity for q and r. 
    */
    static constexpr bls12::TwistType TWIST_TYPE = bls12::TwistType::D;
    static constexpr uint64_t X = 0x8508c00000000001;
    static constexpr auto X_HIGHEST_BIT = bls12::FindMSB<X>::MSB;
    static constexpr auto X_NUM_ONES = bls12::CountOnes<X>::n;
    static constexpr bool X_IS_NEG = false;

    static Fqe_type TWIST_COEFF_B;

    static void init_public_params();

    static G1_precomp_type precompute_G1(const G1_type &P);
    static G2_precomp_type precompute_G2(const G2_type &Q);
    static Fqk_type miller_loop(const G1_precomp_type &prec_P,
                                const G2_precomp_type &prec_Q);
    static Fqk_type double_miller_loop(const G1_precomp_type &prec_P1,
                                       const G2_precomp_type &prec_Q1,
                                       const G1_precomp_type &prec_P2,
                                       const G2_precomp_type &prec_Q2);
    static Fqk_type final_exponentiation(const Fqk_type &elt);
    static Fqk_type pairing(const G1_type &P, const G2_type &Q);
    static Fqk_type reduced_pairing(const G1_type &P, const G2_type &Q);
};

} // libff

#endif // BLS12_377_PP_HPP_
