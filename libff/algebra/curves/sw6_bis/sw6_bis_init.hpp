#ifndef SW6_BIS_INIT_HPP_
#define SW6_BIS_INIT_HPP_

#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/fields/fp.hpp>
#include <libff/algebra/fields/fp3.hpp>
#include <libff/algebra/fields/fp6_2over3.hpp>

namespace libff {

const mp_size_t sw6_bis_r_bitcount = 377;
const mp_size_t sw6_bis_q_bitcount = 761;

const mp_size_t sw6_bis_r_limbs = (sw6_bis_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t sw6_bis_q_limbs = (sw6_bis_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<sw6_bis_r_limbs> sw6_bis_modulus_r;
extern bigint<sw6_bis_q_limbs> sw6_bis_modulus_q;

typedef Fp_model<sw6_bis_r_limbs, sw6_bis_modulus_r> sw6_bis_Fr;
typedef Fp_model<sw6_bis_q_limbs, sw6_bis_modulus_q> sw6_bis_Fq;
typedef Fp3_model<sw6_bis_q_limbs, sw6_bis_modulus_q> sw6_bis_Fq3;
typedef Fp6_2over3_model<sw6_bis_q_limbs, sw6_bis_modulus_q> sw6_bis_Fq6;
typedef sw6_bis_Fq6 sw6_bis_GT;

// parameters for twisted short Weierstrass curve E'/Fq3 : y^2 = x^3 + (a * twist^2) * x + (b * twist^3)
extern sw6_bis_Fq3 sw6_bis_twist;
extern sw6_bis_Fq3 sw6_bis_twist_coeff_a;
extern sw6_bis_Fq3 sw6_bis_twist_coeff_b;
extern sw6_bis_Fq sw6_bis_twist_mul_by_a_c0;
extern sw6_bis_Fq sw6_bis_twist_mul_by_a_c1;
extern sw6_bis_Fq sw6_bis_twist_mul_by_a_c2;
extern sw6_bis_Fq sw6_bis_twist_mul_by_b_c0;
extern sw6_bis_Fq sw6_bis_twist_mul_by_b_c1;
extern sw6_bis_Fq sw6_bis_twist_mul_by_b_c2;
extern sw6_bis_Fq sw6_bis_twist_mul_by_q_X;
extern sw6_bis_Fq sw6_bis_twist_mul_by_q_Y;

// parameters for pairing
extern bigint<sw6_bis_q_limbs> sw6_bis_ate_loop_count;
extern bool sw6_bis_ate_is_loop_count_neg;
extern bigint<6*sw6_bis_q_limbs> sw6_bis_final_exponent;
extern bigint<sw6_bis_q_limbs> sw6_bis_final_exponent_z;
extern bool sw6_bis_final_exponent_is_z_neg;

// old
extern bigint<sw6_bis_q_limbs> sw6_bis_final_exponent_last_chunk_abs_of_w0;
extern bool sw6_bis_final_exponent_last_chunk_is_w0_neg;
extern bigint<sw6_bis_q_limbs> sw6_bis_final_exponent_last_chunk_w1;

void init_sw6_bis_params();

class sw6_bis_G1;
class sw6_bis_G2;

} // libff

#endif // SW6_BIS_INIT_HPP_
