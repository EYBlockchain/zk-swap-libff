#ifndef TOY_CURVE_INIT_HPP_
#define TOY_CURVE_INIT_HPP_
#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/fields/fp.hpp>
#include <libff/algebra/fields/fp12_2over3over2.hpp>
#include <libff/algebra/fields/fp2.hpp>
#include <libff/algebra/fields/fp6_3over2.hpp>

namespace libff {

const mp_size_t toy_curve_r_bitcount = 81;
const mp_size_t toy_curve_q_bitcount = 81;

const mp_size_t toy_curve_r_limbs = (toy_curve_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t toy_curve_q_limbs = (toy_curve_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<toy_curve_r_limbs> toy_curve_modulus_r;
extern bigint<toy_curve_q_limbs> toy_curve_modulus_q;

typedef Fp_model<toy_curve_r_limbs, toy_curve_modulus_r> toy_curve_Fr;
typedef Fp_model<toy_curve_q_limbs, toy_curve_modulus_q> toy_curve_Fq;
typedef Fp2_model<toy_curve_q_limbs, toy_curve_modulus_q> toy_curve_Fq2;
typedef Fp6_3over2_model<toy_curve_q_limbs, toy_curve_modulus_q> toy_curve_Fq6;
typedef Fp12_2over3over2_model<toy_curve_q_limbs, toy_curve_modulus_q> toy_curve_Fq12;
typedef toy_curve_Fq12 toy_curve_GT;

// parameters for Barreto--Naehrig curve E/Fq : y^2 = x^3 + b
extern toy_curve_Fq toy_curve_coeff_b;
// parameters for twisted Barreto--Naehrig curve E'/Fq2 : y^2 = x^3 + b/xi
extern toy_curve_Fq2 toy_curve_twist;
extern toy_curve_Fq2 toy_curve_twist_coeff_b;
extern toy_curve_Fq toy_curve_twist_mul_by_b_c0;
extern toy_curve_Fq toy_curve_twist_mul_by_b_c1;
extern toy_curve_Fq2 toy_curve_twist_mul_by_q_X;
extern toy_curve_Fq2 toy_curve_twist_mul_by_q_Y;

// parameters for pairing
extern bigint<toy_curve_q_limbs> toy_curve_ate_loop_count;
extern bool toy_curve_ate_is_loop_count_neg;
extern bigint<12*toy_curve_q_limbs> toy_curve_final_exponent;
extern bigint<toy_curve_q_limbs> toy_curve_final_exponent_z;
extern bool toy_curve_final_exponent_is_z_neg;

void init_toy_curve_params();

class toy_curve_G1;
class toy_curve_G2;

} // libff
#endif // TOY_CURVE_INIT_HPP_
