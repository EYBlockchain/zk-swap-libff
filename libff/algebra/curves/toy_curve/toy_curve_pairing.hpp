#ifndef TOY_CURVE_PAIRING_HPP_
#define TOY_CURVE_PAIRING_HPP_
#include <vector>

#include <libff/algebra/curves/toy_curve/toy_curve_init.hpp>

namespace libff {

/* final exponentiation */

toy_curve_GT toy_curve_final_exponentiation(const toy_curve_Fq12 &elt);

/* ate pairing */

struct toy_curve_ate_G1_precomp {
    toy_curve_Fq PX;
    toy_curve_Fq PY;

    bool operator==(const toy_curve_ate_G1_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const toy_curve_ate_G1_precomp &prec_P);
    friend std::istream& operator>>(std::istream &in, toy_curve_ate_G1_precomp &prec_P);
};

struct toy_curve_ate_ell_coeffs {
    toy_curve_Fq2 ell_0;
    toy_curve_Fq2 ell_VW;
    toy_curve_Fq2 ell_VV;

    bool operator==(const toy_curve_ate_ell_coeffs &other) const;
    friend std::ostream& operator<<(std::ostream &out, const toy_curve_ate_ell_coeffs &dc);
    friend std::istream& operator>>(std::istream &in, toy_curve_ate_ell_coeffs &dc);
};

struct toy_curve_ate_G2_precomp {
    toy_curve_Fq2 QX;
    toy_curve_Fq2 QY;
    std::vector<toy_curve_ate_ell_coeffs> coeffs;

    bool operator==(const toy_curve_ate_G2_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const toy_curve_ate_G2_precomp &prec_Q);
    friend std::istream& operator>>(std::istream &in, toy_curve_ate_G2_precomp &prec_Q);
};

toy_curve_ate_G1_precomp toy_curve_ate_precompute_G1(const toy_curve_G1& P);
toy_curve_ate_G2_precomp toy_curve_ate_precompute_G2(const toy_curve_G2& Q);

toy_curve_Fq12 toy_curve_ate_miller_loop(const toy_curve_ate_G1_precomp &prec_P,
                              const toy_curve_ate_G2_precomp &prec_Q);
toy_curve_Fq12 toy_curve_ate_double_miller_loop(const toy_curve_ate_G1_precomp &prec_P1,
                                     const toy_curve_ate_G2_precomp &prec_Q1,
                                     const toy_curve_ate_G1_precomp &prec_P2,
                                     const toy_curve_ate_G2_precomp &prec_Q2);

toy_curve_Fq12 toy_curve_ate_pairing(const toy_curve_G1& P,
                          const toy_curve_G2 &Q);
toy_curve_GT toy_curve_ate_reduced_pairing(const toy_curve_G1 &P,
                                 const toy_curve_G2 &Q);

/* choice of pairing */

typedef toy_curve_ate_G1_precomp toy_curve_G1_precomp;
typedef toy_curve_ate_G2_precomp toy_curve_G2_precomp;

toy_curve_G1_precomp toy_curve_precompute_G1(const toy_curve_G1& P);

toy_curve_G2_precomp toy_curve_precompute_G2(const toy_curve_G2& Q);

toy_curve_Fq12 toy_curve_miller_loop(const toy_curve_G1_precomp &prec_P,
                          const toy_curve_G2_precomp &prec_Q);

toy_curve_Fq12 toy_curve_double_miller_loop(const toy_curve_G1_precomp &prec_P1,
                                 const toy_curve_G2_precomp &prec_Q1,
                                 const toy_curve_G1_precomp &prec_P2,
                                 const toy_curve_G2_precomp &prec_Q2);

toy_curve_Fq12 toy_curve_pairing(const toy_curve_G1& P,
                      const toy_curve_G2 &Q);

toy_curve_GT toy_curve_reduced_pairing(const toy_curve_G1 &P,
                             const toy_curve_G2 &Q);

toy_curve_GT toy_curve_affine_reduced_pairing(const toy_curve_G1 &P,
                                    const toy_curve_G2 &Q);

} // libff
#endif // TOY_CURVE_PAIRING_HPP_
