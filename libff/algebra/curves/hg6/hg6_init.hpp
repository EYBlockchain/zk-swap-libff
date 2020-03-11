#ifndef SW6_BIS_INIT_HPP_
#define SW6_BIS_INIT_HPP_

#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/fields/fp.hpp>
#include <libff/algebra/fields/fp3.hpp>
#include <libff/algebra/fields/fp6_2over3.hpp>

namespace libff {

const mp_size_t hg6_r_bitcount = 377;
const mp_size_t hg6_q_bitcount = 761;

const mp_size_t hg6_r_limbs = (hg6_r_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;
const mp_size_t hg6_q_limbs = (hg6_q_bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

extern bigint<hg6_r_limbs> hg6_modulus_r;
extern bigint<hg6_q_limbs> hg6_modulus_q;

typedef Fp_model<hg6_r_limbs, hg6_modulus_r> hg6_Fr;
typedef Fp_model<hg6_q_limbs, hg6_modulus_q> hg6_Fq;
typedef Fp3_model<hg6_q_limbs, hg6_modulus_q> hg6_Fq3;
typedef Fp6_2over3_model<hg6_q_limbs, hg6_modulus_q> hg6_Fq6;
typedef hg6_Fq6 hg6_GT;

// parameters for D-twisted short Weierstrass curve E'/Fq : y^2 = x^3 + (b / xi)
extern hg6_Fq hg6_twist;
extern hg6_Fq hg6_twist_coeff_a;
extern hg6_Fq hg6_twist_coeff_b;

// parameters for pairing
extern bigint<hg6_q_limbs> hg6_ate_loop_count1;
extern bigint<hg6_q_limbs> hg6_ate_loop_count2;
extern bool hg6_ate_is_loop_count_neg;
extern bigint<hg6_q_limbs> hg6_final_exponent_z;
extern bool hg6_final_exponent_is_z_neg;

void init_hg6_params();

class hg6_G1;
class hg6_G2;

} // libff

#endif // SW6_BIS_INIT_HPP_
