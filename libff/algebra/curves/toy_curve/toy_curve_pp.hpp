#ifndef TOY_CURVE_PP_HPP_
#define TOY_CURVE_PP_HPP_
#include <libff/algebra/curves/toy_curve/toy_curve_g1.hpp>
#include <libff/algebra/curves/toy_curve/toy_curve_g2.hpp>
#include <libff/algebra/curves/toy_curve/toy_curve_init.hpp>
#include <libff/algebra/curves/toy_curve/toy_curve_pairing.hpp>
#include <libff/algebra/curves/public_params.hpp>

namespace libff {

class toy_curve_pp {
public:
    typedef toy_curve_Fr Fp_type;
    typedef toy_curve_G1 G1_type;
    typedef toy_curve_G2 G2_type;
    typedef toy_curve_G1_precomp G1_precomp_type;
    typedef toy_curve_G2_precomp G2_precomp_type;
    typedef toy_curve_Fq Fq_type;
    typedef toy_curve_Fq2 Fqe_type;
    typedef toy_curve_Fq12 Fqk_type;
    typedef toy_curve_GT GT_type;

    static const bool has_affine_pairing = false;

    static void init_public_params();
    static toy_curve_GT final_exponentiation(const toy_curve_Fq12 &elt);
    static toy_curve_G1_precomp precompute_G1(const toy_curve_G1 &P);
    static toy_curve_G2_precomp precompute_G2(const toy_curve_G2 &Q);
    static toy_curve_Fq12 miller_loop(const toy_curve_G1_precomp &prec_P,
                                      const toy_curve_G2_precomp &prec_Q);
    static toy_curve_Fq12 double_miller_loop(const toy_curve_G1_precomp &prec_P1,
                                             const toy_curve_G2_precomp &prec_Q1,
                                             const toy_curve_G1_precomp &prec_P2,
                                             const toy_curve_G2_precomp &prec_Q2);
    static toy_curve_Fq12 pairing(const toy_curve_G1 &P,
                                  const toy_curve_G2 &Q);
    static toy_curve_Fq12 reduced_pairing(const toy_curve_G1 &P,
                                          const toy_curve_G2 &Q);
};

} // libff

#endif // TOY_CURVE_PP_HPP_
