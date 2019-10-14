#ifndef SW6_BIS_PP_HPP_
#define SW6_BIS_PP_HPP_

#include <libff/algebra/curves/sw6_bis/sw6_bis_g1.hpp>
#include <libff/algebra/curves/sw6_bis/sw6_bis_g2.hpp>
#include <libff/algebra/curves/sw6_bis/sw6_bis_init.hpp>
#include <libff/algebra/curves/sw6_bis/sw6_bis_pairing.hpp>
#include <libff/algebra/curves/public_params.hpp>

namespace libff {

class sw6_bis_pp {
public:
    typedef sw6_bis_Fr Fp_type;
    typedef sw6_bis_G1 G1_type;
    typedef sw6_bis_G2 G2_type;
    typedef sw6_bis_affine_ate_G1_precomputation affine_ate_G1_precomp_type;
    typedef sw6_bis_affine_ate_G2_precomputation affine_ate_G2_precomp_type;
    typedef sw6_bis_G1_precomp G1_precomp_type;
    typedef sw6_bis_G2_precomp G2_precomp_type;
    typedef sw6_bis_Fq Fq_type;
    typedef sw6_bis_Fq3 Fqe_type;
    typedef sw6_bis_Fq6 Fqk_type;
    typedef sw6_bis_GT GT_type;

    static const bool has_affine_pairing = true;

    static void init_public_params();
    static sw6_bis_GT final_exponentiation(const sw6_bis_Fq6 &elt);
    static sw6_bis_G1_precomp precompute_G1(const sw6_bis_G1 &P);
    static sw6_bis_G2_precomp precompute_G2(const sw6_bis_G2 &Q);
    static sw6_bis_Fq6 miller_loop(const sw6_bis_G1_precomp &prec_P,
                                const sw6_bis_G2_precomp &prec_Q);
    static sw6_bis_affine_ate_G1_precomputation affine_ate_precompute_G1(const sw6_bis_G1 &P);
    static sw6_bis_affine_ate_G2_precomputation affine_ate_precompute_G2(const sw6_bis_G2 &Q);
    static sw6_bis_Fq6 affine_ate_miller_loop(const sw6_bis_affine_ate_G1_precomputation &prec_P,
                                           const sw6_bis_affine_ate_G2_precomputation &prec_Q);
    static sw6_bis_Fq6 affine_ate_e_over_e_miller_loop(const sw6_bis_affine_ate_G1_precomputation &prec_P1,
                                                    const sw6_bis_affine_ate_G2_precomputation &prec_Q1,
                                                    const sw6_bis_affine_ate_G1_precomputation &prec_P2,
                                                    const sw6_bis_affine_ate_G2_precomputation &prec_Q2);
    static sw6_bis_Fq6 affine_ate_e_times_e_over_e_miller_loop(const sw6_bis_affine_ate_G1_precomputation &prec_P1,
                                                            const sw6_bis_affine_ate_G2_precomputation &prec_Q1,
                                                            const sw6_bis_affine_ate_G1_precomputation &prec_P2,
                                                            const sw6_bis_affine_ate_G2_precomputation &prec_Q2,
                                                            const sw6_bis_affine_ate_G1_precomputation &prec_P3,
                                                            const sw6_bis_affine_ate_G2_precomputation &prec_Q3);
    static sw6_bis_Fq6 double_miller_loop(const sw6_bis_G1_precomp &prec_P1,
                                       const sw6_bis_G2_precomp &prec_Q1,
                                       const sw6_bis_G1_precomp &prec_P2,
                                       const sw6_bis_G2_precomp &prec_Q2);

    /* the following are used in test files */
    static sw6_bis_Fq6 pairing(const sw6_bis_G1 &P,
                            const sw6_bis_G2 &Q);
    static sw6_bis_Fq6 reduced_pairing(const sw6_bis_G1 &P,
                                    const sw6_bis_G2 &Q);
    static sw6_bis_Fq6 affine_reduced_pairing(const sw6_bis_G1 &P,
                                           const sw6_bis_G2 &Q);
};

} // libff

#endif // SW6_BIS_PP_HPP_
