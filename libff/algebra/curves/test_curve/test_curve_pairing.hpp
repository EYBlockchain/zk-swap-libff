#ifndef TEST_CURVE_PAIRING_HPP_
#define TEST_CURVE_PAIRING_HPP_
#include <vector>

#include <libff/algebra/curves/test_curve/test_curve_init.hpp>

namespace libff {

/* final exponentiation */

test_curve_GT test_curve_final_exponentiation(const test_curve_Fq12 &elt);

/* ate pairing */

struct test_curve_ate_G1_precomp {
    test_curve_Fq PX;
    test_curve_Fq PY;

    bool operator==(const test_curve_ate_G1_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const test_curve_ate_G1_precomp &prec_P);
    friend std::istream& operator>>(std::istream &in, test_curve_ate_G1_precomp &prec_P);
};

struct test_curve_ate_ell_coeffs {
    test_curve_Fq2 ell_0;
    test_curve_Fq2 ell_VW;
    test_curve_Fq2 ell_VV;

    bool operator==(const test_curve_ate_ell_coeffs &other) const;
    friend std::ostream& operator<<(std::ostream &out, const test_curve_ate_ell_coeffs &dc);
    friend std::istream& operator>>(std::istream &in, test_curve_ate_ell_coeffs &dc);
};

struct test_curve_ate_G2_precomp {
    test_curve_Fq2 QX;
    test_curve_Fq2 QY;
    std::vector<test_curve_ate_ell_coeffs> coeffs;

    bool operator==(const test_curve_ate_G2_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const test_curve_ate_G2_precomp &prec_Q);
    friend std::istream& operator>>(std::istream &in, test_curve_ate_G2_precomp &prec_Q);
};

test_curve_ate_G1_precomp test_curve_ate_precompute_G1(const test_curve_G1& P);
test_curve_ate_G2_precomp test_curve_ate_precompute_G2(const test_curve_G2& Q);

test_curve_Fq12 test_curve_ate_miller_loop(const test_curve_ate_G1_precomp &prec_P,
                              const test_curve_ate_G2_precomp &prec_Q);
test_curve_Fq12 test_curve_ate_double_miller_loop(const test_curve_ate_G1_precomp &prec_P1,
                                     const test_curve_ate_G2_precomp &prec_Q1,
                                     const test_curve_ate_G1_precomp &prec_P2,
                                     const test_curve_ate_G2_precomp &prec_Q2);

test_curve_Fq12 test_curve_ate_pairing(const test_curve_G1& P,
                          const test_curve_G2 &Q);
test_curve_GT test_curve_ate_reduced_pairing(const test_curve_G1 &P,
                                 const test_curve_G2 &Q);

/* choice of pairing */

typedef test_curve_ate_G1_precomp test_curve_G1_precomp;
typedef test_curve_ate_G2_precomp test_curve_G2_precomp;

test_curve_G1_precomp test_curve_precompute_G1(const test_curve_G1& P);

test_curve_G2_precomp test_curve_precompute_G2(const test_curve_G2& Q);

test_curve_Fq12 test_curve_miller_loop(const test_curve_G1_precomp &prec_P,
                          const test_curve_G2_precomp &prec_Q);

test_curve_Fq12 test_curve_double_miller_loop(const test_curve_G1_precomp &prec_P1,
                                 const test_curve_G2_precomp &prec_Q1,
                                 const test_curve_G1_precomp &prec_P2,
                                 const test_curve_G2_precomp &prec_Q2);

test_curve_Fq12 test_curve_pairing(const test_curve_G1& P,
                      const test_curve_G2 &Q);

test_curve_GT test_curve_reduced_pairing(const test_curve_G1 &P,
                             const test_curve_G2 &Q);

test_curve_GT test_curve_affine_reduced_pairing(const test_curve_G1 &P,
                                    const test_curve_G2 &Q);

} // libff
#endif // TEST_CURVE_PAIRING_HPP_
