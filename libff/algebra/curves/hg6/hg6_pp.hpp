#ifndef SW6_BIS_PP_HPP_
#define SW6_BIS_PP_HPP_

#include <libff/algebra/curves/hg6/hg6_g1.hpp>
#include <libff/algebra/curves/hg6/hg6_g2.hpp>
#include <libff/algebra/curves/hg6/hg6_init.hpp>
#include <libff/algebra/curves/hg6/hg6_pairing.hpp>
#include <libff/algebra/curves/public_params.hpp>

namespace libff {

class hg6_pp {
public:
    typedef hg6_Fr Fp_type;
    typedef hg6_G1 G1_type;
    typedef hg6_G2 G2_type;
    /*
    typedef hg6_affine_ate_G1_precomputation affine_ate_G1_precomp_type;
    typedef hg6_affine_ate_G2_precomputation affine_ate_G2_precomp_type;
    */
    typedef hg6_G1_precomp G1_precomp_type;
    typedef hg6_G2_precomp G2_precomp_type;
    typedef hg6_Fq Fq_type;
    typedef hg6_Fq3 Fqe_type;
    typedef hg6_Fq6 Fqk_type;
    typedef hg6_GT GT_type;

    // static const bool has_affine_pairing = true;

    static void init_public_params();
    static hg6_GT final_exponentiation(const hg6_Fq6 &elt);
    static hg6_G1_precomp precompute_G1(const hg6_G1 &P);
    static hg6_G2_precomp precompute_G2(const hg6_G2& Q, const bigint<hg6_Fq::num_limbs> &loop_count);
    static hg6_Fq6 miller_loop(const hg6_G1_precomp &prec_P,
                                const hg6_G2_precomp &prec_Q_1,
                                const hg6_G2_precomp &prec_Q_2);
    /*
    static hg6_affine_ate_G1_precomputation affine_ate_precompute_G1(const hg6_G1 &P);
    static hg6_affine_ate_G2_precomputation affine_ate_precompute_G2(const hg6_G2 &Q, const bigint<hg6_Fq::num_limbs> &loop_count);
    static hg6_Fq6 affine_ate_miller_loop(const hg6_affine_ate_G1_precomputation &prec_P,
                                           const hg6_affine_ate_G2_precomputation &prec_Q_1,
                                           const hg6_affine_ate_G2_precomputation &prec_Q_2);
    static hg6_Fq6 affine_ate_e_over_e_miller_loop(const hg6_affine_ate_G1_precomputation &prec_P1,
                                                    const hg6_affine_ate_G2_precomputation &prec_Q1,
                                                    const hg6_affine_ate_G1_precomputation &prec_P2,
                                                    const hg6_affine_ate_G2_precomputation &prec_Q2);
    static hg6_Fq6 affine_ate_e_times_e_over_e_miller_loop(const hg6_affine_ate_G1_precomputation &prec_P1,
                                                            const hg6_affine_ate_G2_precomputation &prec_Q1,
                                                            const hg6_affine_ate_G1_precomputation &prec_P2,
                                                            const hg6_affine_ate_G2_precomputation &prec_Q2,
                                                            const hg6_affine_ate_G1_precomputation &prec_P3,
                                                            const hg6_affine_ate_G2_precomputation &prec_Q3);
    static hg6_Fq6 double_miller_loop(const hg6_G1_precomp &prec_P1,
                                       const hg6_G2_precomp &prec_Q1,
                                       const hg6_G1_precomp &prec_P2,
                                       const hg6_G2_precomp &prec_Q2);
    */

    /* the following are used in test files */
    static hg6_Fq6 pairing(const hg6_G1 &P,
                            const hg6_G2 &Q);
    static hg6_Fq6 reduced_pairing(const hg6_G1 &P,
                                    const hg6_G2 &Q);
    /*
    static hg6_Fq6 affine_reduced_pairing(const hg6_G1 &P,
                                           const hg6_G2 &Q);
    */
};

} // libff

#endif // SW6_BIS_PP_HPP_
