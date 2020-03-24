#ifndef SW6_BIS_PP_HPP_
#define SW6_BIS_PP_HPP_

#include <libff/algebra/curves/bw6_761/bw6_761_g1.hpp>
#include <libff/algebra/curves/bw6_761/bw6_761_g2.hpp>
#include <libff/algebra/curves/bw6_761/bw6_761_init.hpp>
#include <libff/algebra/curves/bw6_761/bw6_761_pairing.hpp>
#include <libff/algebra/curves/public_params.hpp>

namespace libff {

class bw6_761_pp {
public:
    typedef bw6_761_Fr Fp_type;
    typedef bw6_761_G1 G1_type;
    typedef bw6_761_G2 G2_type;
    /*
    typedef bw6_761_affine_ate_G1_precomputation affine_ate_G1_precomp_type;
    typedef bw6_761_affine_ate_G2_precomputation affine_ate_G2_precomp_type;
    */
    typedef bw6_761_G1_precomp G1_precomp_type;
    typedef bw6_761_G2_precomp G2_precomp_type;
    typedef bw6_761_Fq Fq_type;
    typedef bw6_761_Fq3 Fqe_type;
    typedef bw6_761_Fq6 Fqk_type;
    typedef bw6_761_GT GT_type;

    // static const bool has_affine_pairing = true;

    static void init_public_params();
    static bw6_761_GT final_exponentiation(const bw6_761_Fq6 &elt);
    static bw6_761_G1_precomp precompute_G1(const bw6_761_G1 &P);
    static bw6_761_G2_precomp precompute_G2(const bw6_761_G2& Q, const bigint<bw6_761_Fq::num_limbs> &loop_count);
    static bw6_761_Fq6 miller_loop(const bw6_761_G1_precomp &prec_P,
                                const bw6_761_G2_precomp &prec_Q_1,
                                const bw6_761_G2_precomp &prec_Q_2);
    /*
    static bw6_761_affine_ate_G1_precomputation affine_ate_precompute_G1(const bw6_761_G1 &P);
    static bw6_761_affine_ate_G2_precomputation affine_ate_precompute_G2(const bw6_761_G2 &Q, const bigint<bw6_761_Fq::num_limbs> &loop_count);
    static bw6_761_Fq6 affine_ate_miller_loop(const bw6_761_affine_ate_G1_precomputation &prec_P,
                                           const bw6_761_affine_ate_G2_precomputation &prec_Q_1,
                                           const bw6_761_affine_ate_G2_precomputation &prec_Q_2);
    static bw6_761_Fq6 affine_ate_e_over_e_miller_loop(const bw6_761_affine_ate_G1_precomputation &prec_P1,
                                                    const bw6_761_affine_ate_G2_precomputation &prec_Q1,
                                                    const bw6_761_affine_ate_G1_precomputation &prec_P2,
                                                    const bw6_761_affine_ate_G2_precomputation &prec_Q2);
    static bw6_761_Fq6 affine_ate_e_times_e_over_e_miller_loop(const bw6_761_affine_ate_G1_precomputation &prec_P1,
                                                            const bw6_761_affine_ate_G2_precomputation &prec_Q1,
                                                            const bw6_761_affine_ate_G1_precomputation &prec_P2,
                                                            const bw6_761_affine_ate_G2_precomputation &prec_Q2,
                                                            const bw6_761_affine_ate_G1_precomputation &prec_P3,
                                                            const bw6_761_affine_ate_G2_precomputation &prec_Q3);
    static bw6_761_Fq6 double_miller_loop(const bw6_761_G1_precomp &prec_P1,
                                       const bw6_761_G2_precomp &prec_Q1,
                                       const bw6_761_G1_precomp &prec_P2,
                                       const bw6_761_G2_precomp &prec_Q2);
    */

    /* the following are used in test files */
    static bw6_761_Fq6 pairing(const bw6_761_G1 &P,
                            const bw6_761_G2 &Q);
    static bw6_761_Fq6 reduced_pairing(const bw6_761_G1 &P,
                                    const bw6_761_G2 &Q);
    /*
    static bw6_761_Fq6 affine_reduced_pairing(const bw6_761_G1 &P,
                                           const bw6_761_G2 &Q);
    */
};

} // libff

#endif // SW6_BIS_PP_HPP_
