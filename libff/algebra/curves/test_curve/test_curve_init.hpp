#ifndef TEST_CURVE_INIT_HPP_
#define TEST_CURVE_INIT_HPP_
#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/fields/fp.hpp>
#include <libff/algebra/fields/fp12_2over3over2.hpp>
#include <libff/algebra/fields/fp2.hpp>
#include <libff/algebra/fields/fp6_3over2.hpp>

namespace libff {

const mp_size_t test_curve_r_bitcount = 253;
const mp_size_t test_curve_q_bitcount = 377;

const mp_size_t test_curve_r_limbs = (test_curve_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t test_curve_q_limbs = (test_curve_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<test_curve_r_limbs> test_curve_modulus_r;
extern bigint<test_curve_q_limbs> test_curve_modulus_q;

typedef Fp_model<test_curve_r_limbs, test_curve_modulus_r> test_curve_Fr;
typedef Fp_model<test_curve_q_limbs, test_curve_modulus_q> test_curve_Fq;
typedef Fp2_model<test_curve_q_limbs, test_curve_modulus_q> test_curve_Fq2;
typedef Fp6_3over2_model<test_curve_q_limbs, test_curve_modulus_q> test_curve_Fq6;
typedef Fp12_2over3over2_model<test_curve_q_limbs, test_curve_modulus_q> test_curve_Fq12;
typedef test_curve_Fq12 test_curve_GT;

// parameters for Barreto--Naehrig curve E/Fq : y^2 = x^3 + b
extern test_curve_Fq test_curve_coeff_b;
// parameters for twisted Barreto--Naehrig curve E'/Fq2 : y^2 = x^3 + b/xi
extern test_curve_Fq2 test_curve_twist;
extern test_curve_Fq2 test_curve_twist_coeff_b;
extern test_curve_Fq test_curve_twist_mul_by_b_c0;
extern test_curve_Fq test_curve_twist_mul_by_b_c1;
extern test_curve_Fq2 test_curve_twist_mul_by_q_X;
extern test_curve_Fq2 test_curve_twist_mul_by_q_Y;

// parameters for pairing
extern bigint<test_curve_q_limbs> test_curve_ate_loop_count;
extern bool test_curve_ate_is_loop_count_neg;
extern bigint<12*test_curve_q_limbs> test_curve_final_exponent;
extern bigint<test_curve_q_limbs> test_curve_final_exponent_z;
extern bool test_curve_final_exponent_is_z_neg;

void init_test_curve_params();

class test_curve_G1;
class test_curve_G2;

} // libff
#endif // TEST_CURVE_INIT_HPP_
