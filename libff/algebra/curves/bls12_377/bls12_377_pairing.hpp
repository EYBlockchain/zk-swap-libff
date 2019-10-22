/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BLS12_377_PAIRING_HPP_
#define BLS12_377_PAIRING_HPP_
#include <vector>

#include <libff/algebra/curves/bls12_377/bls12_377_init.hpp>

namespace libff {

/* final exponentiation */

bls12_377_GT bls12_377_final_exponentiation(const bls12_377_Fq12 &elt);

/* ate pairing */

struct bls12_377_ate_G1_precomp {
    bls12_377_Fq PX;
    bls12_377_Fq PY;

    bool operator==(const bls12_377_ate_G1_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const bls12_377_ate_G1_precomp &prec_P);
    friend std::istream& operator>>(std::istream &in, bls12_377_ate_G1_precomp &prec_P);
};

struct bls12_377_ate_ell_coeffs {
    bls12_377_Fq2 ell_0;
    bls12_377_Fq2 ell_VW;
    bls12_377_Fq2 ell_VV;

    bool operator==(const bls12_377_ate_ell_coeffs &other) const;
    friend std::ostream& operator<<(std::ostream &out, const bls12_377_ate_ell_coeffs &dc);
    friend std::istream& operator>>(std::istream &in, bls12_377_ate_ell_coeffs &dc);
};

struct bls12_377_ate_G2_precomp {
    bls12_377_Fq2 QX;
    bls12_377_Fq2 QY;
    std::vector<bls12_377_ate_ell_coeffs> coeffs;

    bool operator==(const bls12_377_ate_G2_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const bls12_377_ate_G2_precomp &prec_Q);
    friend std::istream& operator>>(std::istream &in, bls12_377_ate_G2_precomp &prec_Q);
};

bls12_377_ate_G1_precomp bls12_377_ate_precompute_G1(const bls12_377_G1& P);
bls12_377_ate_G2_precomp bls12_377_ate_precompute_G2(const bls12_377_G2& Q);

bls12_377_Fq12 bls12_377_ate_miller_loop(const bls12_377_ate_G1_precomp &prec_P,
                              const bls12_377_ate_G2_precomp &prec_Q);
bls12_377_Fq12 bls12_377_ate_double_miller_loop(const bls12_377_ate_G1_precomp &prec_P1,
                                     const bls12_377_ate_G2_precomp &prec_Q1,
                                     const bls12_377_ate_G1_precomp &prec_P2,
                                     const bls12_377_ate_G2_precomp &prec_Q2);

bls12_377_Fq12 bls12_377_ate_pairing(const bls12_377_G1& P,
                          const bls12_377_G2 &Q);
bls12_377_GT bls12_377_ate_reduced_pairing(const bls12_377_G1 &P,
                                 const bls12_377_G2 &Q);

/* choice of pairing */

typedef bls12_377_ate_G1_precomp bls12_377_G1_precomp;
typedef bls12_377_ate_G2_precomp bls12_377_G2_precomp;

bls12_377_G1_precomp bls12_377_precompute_G1(const bls12_377_G1& P);

bls12_377_G2_precomp bls12_377_precompute_G2(const bls12_377_G2& Q);

bls12_377_Fq12 bls12_377_miller_loop(const bls12_377_G1_precomp &prec_P,
                          const bls12_377_G2_precomp &prec_Q);

bls12_377_Fq12 bls12_377_double_miller_loop(const bls12_377_G1_precomp &prec_P1,
                                 const bls12_377_G2_precomp &prec_Q1,
                                 const bls12_377_G1_precomp &prec_P2,
                                 const bls12_377_G2_precomp &prec_Q2);

bls12_377_Fq12 bls12_377_pairing(const bls12_377_G1& P,
                      const bls12_377_G2 &Q);

bls12_377_GT bls12_377_reduced_pairing(const bls12_377_G1 &P,
                             const bls12_377_G2 &Q);

bls12_377_GT bls12_377_affine_reduced_pairing(const bls12_377_G1 &P,
                                    const bls12_377_G2 &Q);

} // libff
#endif // BLS12_377_PAIRING_HPP_
