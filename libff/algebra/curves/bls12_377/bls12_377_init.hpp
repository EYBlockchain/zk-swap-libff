/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BLS12_377_INIT_HPP_
#define BLS12_377_INIT_HPP_
#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/fields/fp.hpp>
#include <libff/algebra/fields/fp12_2over3over2.hpp>
#include <libff/algebra/fields/fp2.hpp>
#include <libff/algebra/fields/fp6_3over2.hpp>

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

// parameters for Barreto--Naehrig curve E/Fq : y^2 = x^3 + b
extern bls12_377_Fq bls12_377_coeff_b;
// parameters for twisted Barreto--Naehrig curve E'/Fq2 : y^2 = x^3 + b/xi
extern bls12_377_Fq2 bls12_377_twist;
extern bls12_377_Fq2 bls12_377_twist_coeff_b;
extern bls12_377_Fq bls12_377_twist_mul_by_b_c0;
extern bls12_377_Fq bls12_377_twist_mul_by_b_c1;
extern bls12_377_Fq2 bls12_377_twist_mul_by_q_X;
extern bls12_377_Fq2 bls12_377_twist_mul_by_q_Y;

// parameters for pairing
extern bigint<bls12_377_q_limbs> bls12_377_ate_loop_count;
extern bool bls12_377_ate_is_loop_count_neg;
extern bigint<12*bls12_377_q_limbs> bls12_377_final_exponent;
extern bigint<bls12_377_q_limbs> bls12_377_final_exponent_z;
extern bool bls12_377_final_exponent_is_z_neg;

void init_bls12_377_params();

class bls12_377_G1;
class bls12_377_G2;

} // libff
#endif // BLS12_377_INIT_HPP_
