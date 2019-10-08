/** @file
*****************************************************************************
* @author     This file is part of libff, developed by SCIPR Lab
*             and contributors (see AUTHORS).
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/

#ifndef BLS12_381_PP_HPP_
#define BLS12_381_PP_HPP_
#include <libff/algebra/curves/bls12_381/bls12_381_g1.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_g2.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_init.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_pairing.hpp>
#include <libff/algebra/curves/public_params.hpp>

template <uint64_t x>
struct CountOnes {
    static constexpr unsigned int n = (x & 1) + CountOnes<(x>>1)>::n;
};
template <>
struct CountOnes<0> {
    static constexpr unsigned int n = 0;
};

template<uint64_t N, uint64_t B=std::numeric_limits<uint64_t>::digits-1>
struct FindMSB {
    static constexpr unsigned int MSB = (N&(1ul<<B)) ? B: FindMSB<N,B-1>::MSB;
};
template<uint64_t N>
struct FindMSB<N,0> {
    static constexpr unsigned int MSB = (N==0) ? -1 : 0;
};


namespace libff {

class bls12_381_pp {
public:
    typedef bls12_381_Fr Fp_type;
    typedef bls12_381_G1 G1_type;
    typedef bls12_381_G2 G2_type;
    //typedef bls12_381_G1_precomp G1_precomp_type;
    //typedef bls12_381_G2_precomp G2_precomp_type;
    typedef bls12_381_Fq Fq_type;
    typedef bls12_381_Fq2 Fqe_type;
    typedef bls12_381_Fq12 Fqk_type;
    typedef bls12_381_GT GT_type;

    static const bool has_affine_pairing = false;

    /*
    static constexpr unsigned int bls12_381_x_num_set_bits = 6;
    static constexpr unsigned int bls12_381_x_highest_set_bit = 63;
    */

    static constexpr BlsTwistType TWIST_TYPE = BlsTwistType::M;
    static constexpr uint64_t X = 0xd201000000010000;
    static constexpr auto X_HIGHEST_BIT = FindMSB<X>::MSB;
    static constexpr auto X_NUM_BITS = CountOnes<X>::n;
    static constexpr bool X_IS_NEG = true;

    static Fqe_type TWIST_COEFF_B;

    static void init_public_params();
    /*
    static bls12_381_GT final_exponentiation(const bls12_381_Fq12 &elt);
    static bls12_381_G1_precomp precompute_G1(const bls12_381_G1 &P);
    static bls12_381_G2_precomp precompute_G2(const bls12_381_G2 &Q);
    static bls12_381_Fq12 miller_loop(const bls12_381_G1_precomp &prec_P,
                                      const bls12_381_G2_precomp &prec_Q);
    static bls12_381_Fq12 double_miller_loop(const bls12_381_G1_precomp &prec_P1,
                                             const bls12_381_G2_precomp &prec_Q1,
                                             const bls12_381_G1_precomp &prec_P2,
                                             const bls12_381_G2_precomp &prec_Q2);
    */
    static bls12_381_Fq12 pairing(const bls12_381_G1 &P,
                                  const bls12_381_G2 &Q);
    static bls12_381_Fq12 reduced_pairing(const bls12_381_G1 &P,
                                          const bls12_381_G2 &Q);
};

} // libff

#endif // BLS12_381_PP_HPP_
