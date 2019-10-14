/** @file
 *****************************************************************************

 Declaration of interfaces for initializing PENDULUM.

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef PENDULUM_INIT_HPP_
#define PENDULUM_INIT_HPP_

#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/fields/fp.hpp>
#include <libff/algebra/fields/fp3.hpp>
#include <libff/algebra/fields/fp6_2over3.hpp>

namespace libff {

const mp_size_t pendulum_r_bitcount = 298;
const mp_size_t pendulum_q_bitcount = 613;

const mp_size_t pendulum_r_limbs = (pendulum_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t pendulum_q_limbs = (pendulum_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<pendulum_r_limbs> pendulum_modulus_r;
extern bigint<pendulum_q_limbs> pendulum_modulus_q;

typedef Fp_model<pendulum_r_limbs, pendulum_modulus_r> pendulum_Fr;
typedef Fp_model<pendulum_q_limbs, pendulum_modulus_q> pendulum_Fq;
typedef Fp3_model<pendulum_q_limbs, pendulum_modulus_q> pendulum_Fq3;
typedef Fp6_2over3_model<pendulum_q_limbs, pendulum_modulus_q> pendulum_Fq6;
typedef pendulum_Fq6 pendulum_GT;

// parameters for twisted short Weierstrass curve E'/Fq3 : y^2 = x^3 + (a * twist^2) * x + (b * twist^3)
extern pendulum_Fq3 pendulum_twist;
extern pendulum_Fq3 pendulum_twist_coeff_a;
extern pendulum_Fq3 pendulum_twist_coeff_b;
extern pendulum_Fq pendulum_twist_mul_by_a_c0;
extern pendulum_Fq pendulum_twist_mul_by_a_c1;
extern pendulum_Fq pendulum_twist_mul_by_a_c2;
extern pendulum_Fq pendulum_twist_mul_by_b_c0;
extern pendulum_Fq pendulum_twist_mul_by_b_c1;
extern pendulum_Fq pendulum_twist_mul_by_b_c2;
extern pendulum_Fq pendulum_twist_mul_by_q_X;
extern pendulum_Fq pendulum_twist_mul_by_q_Y;

// parameters for pairing
extern bigint<pendulum_q_limbs> pendulum_ate_loop_count;
extern bool pendulum_ate_is_loop_count_neg;
extern bigint<6*pendulum_q_limbs> pendulum_final_exponent;
extern bigint<pendulum_q_limbs> pendulum_final_exponent_last_chunk_abs_of_w0;
extern bool pendulum_final_exponent_last_chunk_is_w0_neg;
extern bigint<pendulum_q_limbs> pendulum_final_exponent_last_chunk_w1;

void init_pendulum_params();

class pendulum_G1;
class pendulum_G2;

} // libff

#endif // PENDULUM_INIT_HPP_
