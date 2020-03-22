#ifndef BW12_446_INIT_HPP_
#define BW12_446_INIT_HPP_
#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/fields/fp.hpp>
#include <libff/algebra/fields/fp12_2over3over2.hpp>
#include <libff/algebra/fields/fp2.hpp>
#include <libff/algebra/fields/fp6_3over2.hpp>

namespace libff {

const mp_size_t bw12_446_r_bitcount = 296;
const mp_size_t bw12_446_q_bitcount = 446;

const mp_size_t bw12_446_r_limbs = (bw12_446_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t bw12_446_q_limbs = (bw12_446_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<bw12_446_r_limbs> bw12_446_modulus_r;
extern bigint<bw12_446_q_limbs> bw12_446_modulus_q;

typedef Fp_model<bw12_446_r_limbs, bw12_446_modulus_r> bw12_446_Fr;
typedef Fp_model<bw12_446_q_limbs, bw12_446_modulus_q> bw12_446_Fq;
typedef Fp2_model<bw12_446_q_limbs, bw12_446_modulus_q> bw12_446_Fq2;
typedef Fp6_3over2_model<bw12_446_q_limbs, bw12_446_modulus_q> bw12_446_Fq6;
typedef Fp12_2over3over2_model<bw12_446_q_limbs, bw12_446_modulus_q> bw12_446_Fq12;
typedef bw12_446_Fq12 bw12_446_GT;

// parameters for curve E/Fq : y^2 = x^3 + b
extern bw12_446_Fq bw12_446_coeff_b;
// parameters for twisted curve E'/Fq2 : y^2 = x^3 + b/xi
extern bw12_446_Fq2 bw12_446_twist;
extern bw12_446_Fq2 bw12_446_twist_coeff_b;
extern bw12_446_Fq bw12_446_twist_mul_by_b_c0;
extern bw12_446_Fq bw12_446_twist_mul_by_b_c1;
extern bw12_446_Fq2 bw12_446_twist_mul_by_q_X;
extern bw12_446_Fq2 bw12_446_twist_mul_by_q_Y;

// parameters for pairing
extern bigint<bw12_446_q_limbs> bw12_446_ate_loop_count;
extern bool bw12_446_ate_is_loop_count_neg;
extern bigint<12*bw12_446_q_limbs> bw12_446_final_exponent;
extern bigint<bw12_446_q_limbs> bw12_446_final_exponent_z;
extern bool bw12_446_final_exponent_is_z_neg;

void init_bw12_446_params();

class bw12_446_G1;
class bw12_446_G2;

} // libff
#endif // BW12_446_INIT_HPP_
