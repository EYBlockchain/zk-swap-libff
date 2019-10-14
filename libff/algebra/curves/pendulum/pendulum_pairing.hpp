#ifndef PENDULUM_PAIRING_HPP_
#define PENDULUM_PAIRING_HPP_

#include <vector>

#include <libff/algebra/curves/pendulum/pendulum_init.hpp>

namespace libff {

/* final exponentiation */

pendulum_Fq6 pendulum_final_exponentiation_last_chunk(const pendulum_Fq6 &elt,
                                              const pendulum_Fq6 &elt_inv);
pendulum_Fq6 pendulum_final_exponentiation_first_chunk(const pendulum_Fq6 &elt,
                                               const pendulum_Fq6 &elt_inv);
pendulum_GT pendulum_final_exponentiation(const pendulum_Fq6 &elt);

/* affine ate miller loop */

struct pendulum_affine_ate_G1_precomputation {
    pendulum_Fq PX;
    pendulum_Fq PY;
    pendulum_Fq3 PY_twist_squared;
};

struct pendulum_affine_ate_coeffs {
    // TODO: trim (not all of them are needed)
    pendulum_Fq3 old_RX;
    pendulum_Fq3 old_RY;
    pendulum_Fq3 gamma;
    pendulum_Fq3 gamma_twist;
    pendulum_Fq3 gamma_X;
};

struct pendulum_affine_ate_G2_precomputation {
    pendulum_Fq3 QX;
    pendulum_Fq3 QY;
    std::vector<pendulum_affine_ate_coeffs> coeffs;
};

pendulum_affine_ate_G1_precomputation pendulum_affine_ate_precompute_G1(const pendulum_G1& P);
pendulum_affine_ate_G2_precomputation pendulum_affine_ate_precompute_G2(const pendulum_G2& Q);

pendulum_Fq6 pendulum_affine_ate_miller_loop(const pendulum_affine_ate_G1_precomputation &prec_P,
                                     const pendulum_affine_ate_G2_precomputation &prec_Q);

/* ate pairing */

struct pendulum_ate_G1_precomp {
    pendulum_Fq PX;
    pendulum_Fq PY;
    pendulum_Fq3 PX_twist;
    pendulum_Fq3 PY_twist;

    bool operator==(const pendulum_ate_G1_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const pendulum_ate_G1_precomp &prec_P);
    friend std::istream& operator>>(std::istream &in, pendulum_ate_G1_precomp &prec_P);
};

struct pendulum_ate_dbl_coeffs {
    pendulum_Fq3 c_H;
    pendulum_Fq3 c_4C;
    pendulum_Fq3 c_J;
    pendulum_Fq3 c_L;

    bool operator==(const pendulum_ate_dbl_coeffs &other) const;
    friend std::ostream& operator<<(std::ostream &out, const pendulum_ate_dbl_coeffs &dc);
    friend std::istream& operator>>(std::istream &in, pendulum_ate_dbl_coeffs &dc);
};

struct pendulum_ate_add_coeffs {
    pendulum_Fq3 c_L1;
    pendulum_Fq3 c_RZ;

    bool operator==(const pendulum_ate_add_coeffs &other) const;
    friend std::ostream& operator<<(std::ostream &out, const pendulum_ate_add_coeffs &dc);
    friend std::istream& operator>>(std::istream &in, pendulum_ate_add_coeffs &dc);
};

struct pendulum_ate_G2_precomp {
    pendulum_Fq3 QX;
    pendulum_Fq3 QY;
    pendulum_Fq3 QY2;
    pendulum_Fq3 QX_over_twist;
    pendulum_Fq3 QY_over_twist;
    std::vector<pendulum_ate_dbl_coeffs> dbl_coeffs;
    std::vector<pendulum_ate_add_coeffs> add_coeffs;

    bool operator==(const pendulum_ate_G2_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const pendulum_ate_G2_precomp &prec_Q);
    friend std::istream& operator>>(std::istream &in, pendulum_ate_G2_precomp &prec_Q);
};

pendulum_ate_G1_precomp pendulum_ate_precompute_G1(const pendulum_G1& P);
pendulum_ate_G2_precomp pendulum_ate_precompute_G2(const pendulum_G2& Q);

pendulum_Fq6 pendulum_ate_miller_loop(const pendulum_ate_G1_precomp &prec_P,
                              const pendulum_ate_G2_precomp &prec_Q);
pendulum_Fq6 pendulum_ate_double_miller_loop(const pendulum_ate_G1_precomp &prec_P1,
                                     const pendulum_ate_G2_precomp &prec_Q1,
                                     const pendulum_ate_G1_precomp &prec_P2,
                                     const pendulum_ate_G2_precomp &prec_Q2);

pendulum_Fq6 pendulum_ate_pairing(const pendulum_G1& P,
                          const pendulum_G2 &Q);
pendulum_GT pendulum_ate_reduced_pairing(const pendulum_G1 &P,
                                 const pendulum_G2 &Q);

/* choice of pairing */

typedef pendulum_ate_G1_precomp pendulum_G1_precomp;
typedef pendulum_ate_G2_precomp pendulum_G2_precomp;

pendulum_G1_precomp pendulum_precompute_G1(const pendulum_G1& P);

pendulum_G2_precomp pendulum_precompute_G2(const pendulum_G2& Q);

pendulum_Fq6 pendulum_miller_loop(const pendulum_G1_precomp &prec_P,
                          const pendulum_G2_precomp &prec_Q);

pendulum_Fq6 pendulum_double_miller_loop(const pendulum_G1_precomp &prec_P1,
                                 const pendulum_G2_precomp &prec_Q1,
                                 const pendulum_G1_precomp &prec_P2,
                                 const pendulum_G2_precomp &prec_Q2);

pendulum_Fq6 pendulum_pairing(const pendulum_G1& P,
                      const pendulum_G2 &Q);

pendulum_GT pendulum_reduced_pairing(const pendulum_G1 &P,
                             const pendulum_G2 &Q);

pendulum_GT pendulum_affine_reduced_pairing(const pendulum_G1 &P,
                                    const pendulum_G2 &Q);

} // libff

#endif // PENDULUM_PAIRING_HPP_
