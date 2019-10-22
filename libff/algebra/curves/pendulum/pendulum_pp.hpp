#ifndef PENDULUM_PP_HPP_
#define PENDULUM_PP_HPP_

#include <libff/algebra/curves/pendulum/pendulum_g1.hpp>
#include <libff/algebra/curves/pendulum/pendulum_g2.hpp>
#include <libff/algebra/curves/pendulum/pendulum_init.hpp>
#include <libff/algebra/curves/pendulum/pendulum_pairing.hpp>
#include <libff/algebra/curves/public_params.hpp>

namespace libff {

class pendulum_pp {
public:
    typedef pendulum_Fr Fp_type;
    typedef pendulum_G1 G1_type;
    typedef pendulum_G2 G2_type;
    typedef pendulum_affine_ate_G1_precomputation affine_ate_G1_precomp_type;
    typedef pendulum_affine_ate_G2_precomputation affine_ate_G2_precomp_type;
    typedef pendulum_G1_precomp G1_precomp_type;
    typedef pendulum_G2_precomp G2_precomp_type;
    typedef pendulum_Fq Fq_type;
    typedef pendulum_Fq3 Fqe_type;
    typedef pendulum_Fq6 Fqk_type;
    typedef pendulum_GT GT_type;

    static const bool has_affine_pairing = true;

    static void init_public_params();
    static pendulum_GT final_exponentiation(const pendulum_Fq6 &elt);
    static pendulum_G1_precomp precompute_G1(const pendulum_G1 &P);
    static pendulum_G2_precomp precompute_G2(const pendulum_G2 &Q);
    static pendulum_Fq6 miller_loop(const pendulum_G1_precomp &prec_P,
                                const pendulum_G2_precomp &prec_Q);
    static pendulum_affine_ate_G1_precomputation affine_ate_precompute_G1(const pendulum_G1 &P);
    static pendulum_affine_ate_G2_precomputation affine_ate_precompute_G2(const pendulum_G2 &Q);
    static pendulum_Fq6 affine_ate_miller_loop(const pendulum_affine_ate_G1_precomputation &prec_P,
                                           const pendulum_affine_ate_G2_precomputation &prec_Q);
    static pendulum_Fq6 affine_ate_e_over_e_miller_loop(const pendulum_affine_ate_G1_precomputation &prec_P1,
                                                    const pendulum_affine_ate_G2_precomputation &prec_Q1,
                                                    const pendulum_affine_ate_G1_precomputation &prec_P2,
                                                    const pendulum_affine_ate_G2_precomputation &prec_Q2);
    static pendulum_Fq6 affine_ate_e_times_e_over_e_miller_loop(const pendulum_affine_ate_G1_precomputation &prec_P1,
                                                            const pendulum_affine_ate_G2_precomputation &prec_Q1,
                                                            const pendulum_affine_ate_G1_precomputation &prec_P2,
                                                            const pendulum_affine_ate_G2_precomputation &prec_Q2,
                                                            const pendulum_affine_ate_G1_precomputation &prec_P3,
                                                            const pendulum_affine_ate_G2_precomputation &prec_Q3);
    static pendulum_Fq6 double_miller_loop(const pendulum_G1_precomp &prec_P1,
                                       const pendulum_G2_precomp &prec_Q1,
                                       const pendulum_G1_precomp &prec_P2,
                                       const pendulum_G2_precomp &prec_Q2);

    /* the following are used in test files */
    static pendulum_Fq6 pairing(const pendulum_G1 &P,
                            const pendulum_G2 &Q);
    static pendulum_Fq6 reduced_pairing(const pendulum_G1 &P,
                                    const pendulum_G2 &Q);
    static pendulum_Fq6 affine_reduced_pairing(const pendulum_G1 &P,
                                           const pendulum_G2 &Q);
};

} // libff

#endif // PENDULUM_PP_HPP_
