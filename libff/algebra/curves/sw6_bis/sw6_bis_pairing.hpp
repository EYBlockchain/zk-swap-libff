#ifndef SW6_BIS_PAIRING_HPP_
#define SW6_BIS_PAIRING_HPP_

#include <vector>

#include <libff/algebra/curves/sw6_bis/sw6_bis_init.hpp>

namespace libff {

/* final exponentiation */

sw6_bis_Fq6 sw6_bis_final_exponentiation_last_chunk(const sw6_bis_Fq6 &elt);
sw6_bis_Fq6 sw6_bis_final_exponentiation_first_chunk(const sw6_bis_Fq6 &elt,
                                               const sw6_bis_Fq6 &elt_inv);
sw6_bis_GT sw6_bis_final_exponentiation(const sw6_bis_Fq6 &elt);

/* affine ate miller loop */

struct sw6_bis_affine_ate_G1_precomputation {
    sw6_bis_Fq PX;
    sw6_bis_Fq PY;
    sw6_bis_Fq3 PY_twist_squared;
};

struct sw6_bis_affine_ate_coeffs {
    // TODO: trim (not all of them are needed)
    sw6_bis_Fq3 old_RX;
    sw6_bis_Fq3 old_RY;
    sw6_bis_Fq3 gamma;
    sw6_bis_Fq3 gamma_twist;
    sw6_bis_Fq3 gamma_X;
};

struct sw6_bis_affine_ate_G2_precomputation {
    sw6_bis_Fq3 QX;
    sw6_bis_Fq3 QY;
    std::vector<sw6_bis_affine_ate_coeffs> coeffs;
};

sw6_bis_affine_ate_G1_precomputation sw6_bis_affine_ate_precompute_G1(const sw6_bis_G1& P);
sw6_bis_affine_ate_G2_precomputation sw6_bis_affine_ate_precompute_G2(const sw6_bis_G2& Q);

sw6_bis_Fq6 sw6_bis_affine_ate_miller_loop(const sw6_bis_affine_ate_G1_precomputation &prec_P,
                                     const sw6_bis_affine_ate_G2_precomputation &prec_Q);

/* ate pairing */

struct sw6_bis_ate_G1_precomp {
    sw6_bis_Fq PX;
    sw6_bis_Fq PY;
    sw6_bis_Fq3 PX_twist;
    sw6_bis_Fq3 PY_twist;

    bool operator==(const sw6_bis_ate_G1_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const sw6_bis_ate_G1_precomp &prec_P);
    friend std::istream& operator>>(std::istream &in, sw6_bis_ate_G1_precomp &prec_P);
};

struct sw6_bis_ate_dbl_coeffs {
    sw6_bis_Fq3 c_H;
    sw6_bis_Fq3 c_4C;
    sw6_bis_Fq3 c_J;
    sw6_bis_Fq3 c_L;

    bool operator==(const sw6_bis_ate_dbl_coeffs &other) const;
    friend std::ostream& operator<<(std::ostream &out, const sw6_bis_ate_dbl_coeffs &dc);
    friend std::istream& operator>>(std::istream &in, sw6_bis_ate_dbl_coeffs &dc);
};

struct sw6_bis_ate_add_coeffs {
    sw6_bis_Fq3 c_L1;
    sw6_bis_Fq3 c_RZ;

    bool operator==(const sw6_bis_ate_add_coeffs &other) const;
    friend std::ostream& operator<<(std::ostream &out, const sw6_bis_ate_add_coeffs &dc);
    friend std::istream& operator>>(std::istream &in, sw6_bis_ate_add_coeffs &dc);
};

struct sw6_bis_ate_G2_precomp {
    sw6_bis_Fq3 QX;
    sw6_bis_Fq3 QY;
    sw6_bis_Fq3 QY2;
    sw6_bis_Fq3 QX_over_twist;
    sw6_bis_Fq3 QY_over_twist;
    std::vector<sw6_bis_ate_dbl_coeffs> dbl_coeffs;
    std::vector<sw6_bis_ate_add_coeffs> add_coeffs;

    bool operator==(const sw6_bis_ate_G2_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const sw6_bis_ate_G2_precomp &prec_Q);
    friend std::istream& operator>>(std::istream &in, sw6_bis_ate_G2_precomp &prec_Q);
};

sw6_bis_ate_G1_precomp sw6_bis_ate_precompute_G1(const sw6_bis_G1& P);
sw6_bis_ate_G2_precomp sw6_bis_ate_precompute_G2(const sw6_bis_G2& Q);

sw6_bis_Fq6 sw6_bis_ate_miller_loop(const sw6_bis_ate_G1_precomp &prec_P,
                              const sw6_bis_ate_G2_precomp &prec_Q);
sw6_bis_Fq6 sw6_bis_ate_double_miller_loop(const sw6_bis_ate_G1_precomp &prec_P1,
                                     const sw6_bis_ate_G2_precomp &prec_Q1,
                                     const sw6_bis_ate_G1_precomp &prec_P2,
                                     const sw6_bis_ate_G2_precomp &prec_Q2);

sw6_bis_Fq6 sw6_bis_ate_pairing(const sw6_bis_G1& P,
                          const sw6_bis_G2 &Q);
sw6_bis_GT sw6_bis_ate_reduced_pairing(const sw6_bis_G1 &P,
                                 const sw6_bis_G2 &Q);

/* choice of pairing */

typedef sw6_bis_ate_G1_precomp sw6_bis_G1_precomp;
typedef sw6_bis_ate_G2_precomp sw6_bis_G2_precomp;

sw6_bis_G1_precomp sw6_bis_precompute_G1(const sw6_bis_G1& P);

sw6_bis_G2_precomp sw6_bis_precompute_G2(const sw6_bis_G2& Q);

sw6_bis_Fq6 sw6_bis_miller_loop(const sw6_bis_G1_precomp &prec_P,
                          const sw6_bis_G2_precomp &prec_Q);

sw6_bis_Fq6 sw6_bis_double_miller_loop(const sw6_bis_G1_precomp &prec_P1,
                                 const sw6_bis_G2_precomp &prec_Q1,
                                 const sw6_bis_G1_precomp &prec_P2,
                                 const sw6_bis_G2_precomp &prec_Q2);

sw6_bis_Fq6 sw6_bis_pairing(const sw6_bis_G1& P,
                      const sw6_bis_G2 &Q);

sw6_bis_GT sw6_bis_reduced_pairing(const sw6_bis_G1 &P,
                             const sw6_bis_G2 &Q);

sw6_bis_GT sw6_bis_affine_reduced_pairing(const sw6_bis_G1 &P,
                                    const sw6_bis_G2 &Q);

} // libff

#endif // SW6_BIS_PAIRING_HPP_
