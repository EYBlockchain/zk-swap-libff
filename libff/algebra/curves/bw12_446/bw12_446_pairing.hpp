/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BW12_446_PAIRING_HPP_
#define BW12_446_PAIRING_HPP_
#include <vector>

#include <libff/algebra/curves/bw12_446/bw12_446_init.hpp>

namespace libff {

/* final exponentiation */

bw12_446_GT bw12_446_final_exponentiation(const bw12_446_Fq12 &elt);

/* ate pairing */

struct bw12_446_ate_G1_precomp {
    bw12_446_Fq PX;
    bw12_446_Fq PY;

    bool operator==(const bw12_446_ate_G1_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const bw12_446_ate_G1_precomp &prec_P);
    friend std::istream& operator>>(std::istream &in, bw12_446_ate_G1_precomp &prec_P);
};

struct bw12_446_ate_ell_coeffs {
    bw12_446_Fq2 ell_0;
    bw12_446_Fq2 ell_VW;
    bw12_446_Fq2 ell_VV;

    bool operator==(const bw12_446_ate_ell_coeffs &other) const;
    friend std::ostream& operator<<(std::ostream &out, const bw12_446_ate_ell_coeffs &dc);
    friend std::istream& operator>>(std::istream &in, bw12_446_ate_ell_coeffs &dc);
};

struct bw12_446_ate_G2_precomp {
    bw12_446_Fq2 QX;
    bw12_446_Fq2 QY;
    std::vector<bw12_446_ate_ell_coeffs> coeffs;

    bool operator==(const bw12_446_ate_G2_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const bw12_446_ate_G2_precomp &prec_Q);
    friend std::istream& operator>>(std::istream &in, bw12_446_ate_G2_precomp &prec_Q);
};

bw12_446_ate_G1_precomp bw12_446_ate_precompute_G1(const bw12_446_G1& P);
bw12_446_ate_G2_precomp bw12_446_ate_precompute_G2(const bw12_446_G2& Q);

bw12_446_Fq12 bw12_446_ate_miller_loop(const bw12_446_ate_G1_precomp &prec_P,
                              const bw12_446_ate_G2_precomp &prec_Q);
bw12_446_Fq12 bw12_446_ate_double_miller_loop(const bw12_446_ate_G1_precomp &prec_P1,
                                     const bw12_446_ate_G2_precomp &prec_Q1,
                                     const bw12_446_ate_G1_precomp &prec_P2,
                                     const bw12_446_ate_G2_precomp &prec_Q2);

bw12_446_Fq12 bw12_446_ate_pairing(const bw12_446_G1& P,
                          const bw12_446_G2 &Q);
bw12_446_GT bw12_446_ate_reduced_pairing(const bw12_446_G1 &P,
                                 const bw12_446_G2 &Q);

/* choice of pairing */

typedef bw12_446_ate_G1_precomp bw12_446_G1_precomp;
typedef bw12_446_ate_G2_precomp bw12_446_G2_precomp;

bw12_446_G1_precomp bw12_446_precompute_G1(const bw12_446_G1& P);

bw12_446_G2_precomp bw12_446_precompute_G2(const bw12_446_G2& Q);

bw12_446_Fq12 bw12_446_miller_loop(const bw12_446_G1_precomp &prec_P,
                          const bw12_446_G2_precomp &prec_Q);

bw12_446_Fq12 bw12_446_double_miller_loop(const bw12_446_G1_precomp &prec_P1,
                                 const bw12_446_G2_precomp &prec_Q1,
                                 const bw12_446_G1_precomp &prec_P2,
                                 const bw12_446_G2_precomp &prec_Q2);

bw12_446_Fq12 bw12_446_pairing(const bw12_446_G1& P,
                      const bw12_446_G2 &Q);

bw12_446_GT bw12_446_reduced_pairing(const bw12_446_G1 &P,
                             const bw12_446_G2 &Q);

bw12_446_GT bw12_446_affine_reduced_pairing(const bw12_446_G1 &P,
                                    const bw12_446_G2 &Q);

} // libff
#endif // BW12_446_PAIRING_HPP_
