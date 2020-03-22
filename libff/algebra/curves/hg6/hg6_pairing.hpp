/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef SW6_BIS_PAIRING_HPP_
#define SW6_BIS_PAIRING_HPP_
#include <vector>

#include <libff/algebra/curves/hg6/hg6_init.hpp>

namespace libff {

/* final exponentiation */

hg6_Fq6 hg6_final_exponentiation_last_chunk(const hg6_Fq6 &elt);
hg6_Fq6 hg6_final_exponentiation_first_chunk(const hg6_Fq6 &elt,
                                                   const hg6_Fq6 &elt_inv);
hg6_GT hg6_final_exponentiation(const hg6_Fq6 &elt);

/* ate pairing */

struct hg6_ate_G1_precomp {
    hg6_Fq PX;
    hg6_Fq PY;

    bool operator==(const hg6_ate_G1_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const hg6_ate_G1_precomp &prec_P);
    friend std::istream& operator>>(std::istream &in, hg6_ate_G1_precomp &prec_P);
};

struct hg6_ate_ell_coeffs {
    hg6_Fq ell_0;
    hg6_Fq ell_VW;
    hg6_Fq ell_VV;

    bool operator==(const hg6_ate_ell_coeffs &other) const;
    friend std::ostream& operator<<(std::ostream &out, const hg6_ate_ell_coeffs &dc);
    friend std::istream& operator>>(std::istream &in, hg6_ate_ell_coeffs &dc);
};

struct hg6_ate_G2_precomp {
    hg6_Fq QX;
    hg6_Fq QY;
    std::vector<hg6_ate_ell_coeffs> coeffs;

    bool operator==(const hg6_ate_G2_precomp &other) const;
    friend std::ostream& operator<<(std::ostream &out, const hg6_ate_G2_precomp &prec_Q);
    friend std::istream& operator>>(std::istream &in, hg6_ate_G2_precomp &prec_Q);
};

hg6_ate_G1_precomp hg6_ate_precompute_G1(const hg6_G1& P);
hg6_ate_G2_precomp hg6_ate_precompute_G2(const hg6_G2& Q, const bigint<hg6_Fq::num_limbs> &loop_count);

hg6_Fq6 hg6_ate_miller_loop(const hg6_ate_G1_precomp &prec_P,
                              const hg6_ate_G2_precomp &prec_Q_1,
                              const hg6_ate_G2_precomp &prec_Q_2);
/*
hg6_Fq6 hg6_ate_double_miller_loop(const hg6_ate_G1_precomp &prec_P1,
                                     const hg6_ate_G2_precomp &prec_Q1,
                                     const hg6_ate_G1_precomp &prec_P2,
                                     const hg6_ate_G2_precomp &prec_Q2);
*/

hg6_Fq6 hg6_ate_pairing(const hg6_G1& P,
                          const hg6_G2 &Q);
hg6_GT hg6_ate_reduced_pairing(const hg6_G1 &P,
                                 const hg6_G2 &Q);

/* choice of pairing */

typedef hg6_ate_G1_precomp hg6_G1_precomp;
typedef hg6_ate_G2_precomp hg6_G2_precomp;

hg6_G1_precomp hg6_precompute_G1(const hg6_G1& P);

hg6_G2_precomp hg6_precompute_G2(const hg6_G2& Q, const bigint<hg6_Fq::num_limbs> &loop_count);

hg6_Fq6 hg6_miller_loop(const hg6_G1_precomp &prec_P,
                          const hg6_G2_precomp &prec_Q_1,
                          const hg6_G2_precomp &prec_Q_2);

/*
hg6_Fq6 hg6_double_miller_loop(const hg6_G1_precomp &prec_P1,
                                 const hg6_G2_precomp &prec_Q1,
                                 const hg6_G1_precomp &prec_P2,
                                 const hg6_G2_precomp &prec_Q2);
*/

hg6_Fq6 hg6_pairing(const hg6_G1& P,
                      const hg6_G2 &Q);

hg6_GT hg6_reduced_pairing(const hg6_G1 &P,
                             const hg6_G2 &Q);

} // libff
#endif // SW6_BIS_PAIRING_HPP_
