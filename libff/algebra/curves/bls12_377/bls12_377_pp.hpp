/** @file
*****************************************************************************
* @author     This file is part of libff, developed by SCIPR Lab
*             and contributors (see AUTHORS).
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/

#ifndef BLS12_377_PP_HPP_
#define BLS12_377_PP_HPP_

#include <libff/algebra/curves/bls12.hpp>
#include <libff/algebra/curves/sw.hpp>
#include <libff/algebra/fields/fp12_2over3over2.hpp>

namespace libff {


const mp_size_t bls12_377_r_bitcount = 253;
const mp_size_t bls12_377_q_bitcount = 377;

const mp_size_t bls12_377_r_limbs = (bls12_377_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t bls12_377_q_limbs = (bls12_377_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<bls12_377_r_limbs> bls12_377_modulus_r;
extern bigint<bls12_377_q_limbs> bls12_377_modulus_q;

typedef Fp_model<bls12_377_r_limbs, bls12_377_modulus_r> bls12_377_Fr;
typedef Fp_model<bls12_377_q_limbs, bls12_377_modulus_q> bls12_377_Fq;
typedef Fp2_model<bls12_377_q_limbs, bls12_377_modulus_q> bls12_377_Fq2;
typedef Fp6_3over2_model<bls12_377_q_limbs, bls12_377_modulus_q> bls12_377_Fq6;
typedef Fp12_2over3over2_model<bls12_377_q_limbs, bls12_377_modulus_q> bls12_377_Fq12;
typedef bls12_377_Fq12 bls12_377_GT;


class bls12_377_G1;
class bls12_377_G2;


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


using bls12_377_G1_base = sw::SWJacobianPoint<bls12_377_pp, bls12_377_pp::Fq_type, bls12_377_pp::Fp_type, bls12_377_G1>;

class bls12_377_G1 : public bls12_377_G1_base
{
    using bls12_377_G1_base::SWJacobianPoint;

public:
    // Pass-thru Copy constructor from base type
    bls12_377_G1( const bls12_377_G1_base &base )
    : bls12_377_G1_base(base) {}

    int sign_bit() const {
        return Y.as_bigint().data[0] & 1;
    }

    // Define storage in our own class, otherwise it gets messy
    static base_field coeff_b;
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static bls12_377_G1 _zero;
    static bls12_377_G1 _one;
};


using bls12_377_G2_base = sw::SWJacobianPoint<bls12_377_pp, bls12_377_pp::Fqe_type, bls12_377_pp::Fp_type, bls12_377_G2>;

class bls12_377_G2 : public bls12_377_G2_base
{
    using bls12_377_G2_base::SWJacobianPoint;

public:
    // Pass-thru Copy constructor from base type
    bls12_377_G2( const bls12_377_G2_base &base )
    : bls12_377_G2_base(base) {}

    int sign_bit() const {
        return Y.c0.as_bigint().data[0] & 1;
    }

    // Define storage in our own class, otherwise it gets messy
    static base_field coeff_b;
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static bls12_377_G2 _zero;
    static bls12_377_G2 _one;
};


} // libff


#endif // BLS12_377_PP_HPP_
