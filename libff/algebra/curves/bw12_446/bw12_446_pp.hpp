/** @file
*****************************************************************************
* @author     This file is part of libff, developed by SCIPR Lab
*             and contributors (see AUTHORS).
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/

#ifndef BW12_446_PP_HPP_
#define BW12_446_PP_HPP_
#include <libff/algebra/curves/bw12_446/bw12_446_g1.hpp>
#include <libff/algebra/curves/bw12_446/bw12_446_g2.hpp>
#include <libff/algebra/curves/bw12_446/bw12_446_init.hpp>
#include <libff/algebra/curves/bw12_446/bw12_446_pairing.hpp>
#include <libff/algebra/curves/public_params.hpp>

namespace libff {

class bw12_446_pp {
public:
    typedef bw12_446_Fr Fp_type;
    typedef bw12_446_G1 G1_type;
    typedef bw12_446_G2 G2_type;
    typedef bw12_446_G1_precomp G1_precomp_type;
    typedef bw12_446_G2_precomp G2_precomp_type;
    typedef bw12_446_Fq Fq_type;
    typedef bw12_446_Fq2 Fqe_type;
    typedef bw12_446_Fq12 Fqk_type;
    typedef bw12_446_GT GT_type;

    static const bool has_affine_pairing = false;

    static void init_public_params();
    static bw12_446_GT final_exponentiation(const bw12_446_Fq12 &elt);
    static bw12_446_G1_precomp precompute_G1(const bw12_446_G1 &P);
    static bw12_446_G2_precomp precompute_G2(const bw12_446_G2 &Q);
    static bw12_446_Fq12 miller_loop(const bw12_446_G1_precomp &prec_P,
                                      const bw12_446_G2_precomp &prec_Q);
    static bw12_446_Fq12 double_miller_loop(const bw12_446_G1_precomp &prec_P1,
                                             const bw12_446_G2_precomp &prec_Q1,
                                             const bw12_446_G1_precomp &prec_P2,
                                             const bw12_446_G2_precomp &prec_Q2);
    static bw12_446_Fq12 pairing(const bw12_446_G1 &P,
                                  const bw12_446_G2 &Q);
    static bw12_446_Fq12 reduced_pairing(const bw12_446_G1 &P,
                                          const bw12_446_G2 &Q);
};

} // libff

#endif // BW12_446_PP_HPP_
