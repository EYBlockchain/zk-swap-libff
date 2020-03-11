/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <cassert>

#include <libff/algebra/curves/hg6/hg6_g1.hpp>
#include <libff/algebra/curves/hg6/hg6_g2.hpp>
#include <libff/algebra/curves/hg6/hg6_init.hpp>
#include <libff/algebra/curves/hg6/hg6_pairing.hpp>
#include <libff/common/profiling.hpp>

namespace libff {

bool hg6_ate_G1_precomp::operator==(const hg6_ate_G1_precomp &other) const
{
    return (this->PX == other.PX &&
            this->PY == other.PY);
}

std::ostream& operator<<(std::ostream &out, const hg6_ate_G1_precomp &prec_P)
{
    out << prec_P.PX << OUTPUT_SEPARATOR << prec_P.PY;

    return out;
}

std::istream& operator>>(std::istream &in, hg6_ate_G1_precomp &prec_P)
{
    in >> prec_P.PX;
    consume_OUTPUT_SEPARATOR(in);
    in >> prec_P.PY;

    return in;
}

bool  hg6_ate_ell_coeffs::operator==(const hg6_ate_ell_coeffs &other) const
{
    return (this->ell_0 == other.ell_0 &&
            this->ell_VW == other.ell_VW &&
            this->ell_VV == other.ell_VV);
}

std::ostream& operator<<(std::ostream &out, const hg6_ate_ell_coeffs &c)
{
    out << c.ell_0 << OUTPUT_SEPARATOR << c.ell_VW << OUTPUT_SEPARATOR << c.ell_VV;
    return out;
}

std::istream& operator>>(std::istream &in, hg6_ate_ell_coeffs &c)
{
    in >> c.ell_0;
    consume_OUTPUT_SEPARATOR(in);
    in >> c.ell_VW;
    consume_OUTPUT_SEPARATOR(in);
    in >> c.ell_VV;

    return in;
}

bool hg6_ate_G2_precomp::operator==(const hg6_ate_G2_precomp &other) const
{
    return (this->QX == other.QX &&
            this->QY == other.QY &&
            this->coeffs == other.coeffs);
}

std::ostream& operator<<(std::ostream& out, const hg6_ate_G2_precomp &prec_Q)
{
    out << prec_Q.QX << OUTPUT_SEPARATOR << prec_Q.QY << "\n";
    out << prec_Q.coeffs.size() << "\n";
    for (const hg6_ate_ell_coeffs &c : prec_Q.coeffs)
    {
        out << c << OUTPUT_NEWLINE;
    }
    return out;
}

std::istream& operator>>(std::istream& in, hg6_ate_G2_precomp &prec_Q)
{
    in >> prec_Q.QX;
    consume_OUTPUT_SEPARATOR(in);
    in >> prec_Q.QY;
    consume_newline(in);

    prec_Q.coeffs.clear();
    size_t s;
    in >> s;

    consume_newline(in);

    prec_Q.coeffs.reserve(s);

    for (size_t i = 0; i < s; ++i)
    {
        hg6_ate_ell_coeffs c;
        in >> c;
        consume_OUTPUT_NEWLINE(in);
        prec_Q.coeffs.emplace_back(c);
    }

    return in;
}

/* final exponentiations */

hg6_Fq6 hg6_final_exponentiation_first_chunk(const hg6_Fq6 &elt, const hg6_Fq6 &elt_inv)
{
    enter_block("Call to hg6_final_exponentiation_first_chunk");

    /* (q^3-1)*(q+1) */

    /* elt_q3 = elt^(q^3) */
    const hg6_Fq6 elt_q3 = elt.Frobenius_map(3);
    /* elt_q3_over_elt = elt^(q^3-1) */
    const hg6_Fq6 elt_q3_over_elt = elt_q3 * elt_inv;
    /* alpha = elt^((q^3-1) * q) */
    const hg6_Fq6 alpha = elt_q3_over_elt.Frobenius_map(1);
    /* beta = elt^((q^3-1)*(q+1) */
    const hg6_Fq6 beta = alpha * elt_q3_over_elt;
    leave_block("Call to hg6_final_exponentiation_first_chunk");
    return beta;
}

hg6_Fq6 hg6_exp_by_z(const hg6_Fq6 &elt)
{
  enter_block("Call to hg6_exp_by_z");

  hg6_Fq6 result = elt.cyclotomic_exp(hg6_final_exponent_z);
  if (hg6_final_exponent_is_z_neg)
  {
    result = result.unitary_inverse();
  }

  leave_block("Call to hg6_exp_by_z");

  return result;
}

hg6_Fq6 hg6_final_exponentiation_last_chunk(const hg6_Fq6 &elt)
{
    enter_block("Call to hg6_final_exponentiation_last_chunk");

    /*
     * R0(x) := (-103*x^7 + 70*x^6 + 269*x^5 - 197*x^4 - 314*x^3 - 73*x^2 - 263*x - 220)
     * R1(x) := (103*x^9 - 276*x^8 + 77*x^7 + 492*x^6 - 445*x^5 - 65*x^4 + 452*x^3 - 181*x^2 + 34*x + 229)
     *
     * La dernière partie de l'exponentiation finale est
     * f^R0(u)*(f^p)^R1(u) avec R0 et R1 ci-dessus et f^p est un Frobenius dans GF(p^6).
     */

    const hg6_Fq6 f0 = elt;
    const hg6_Fq6 f0p = f0.Frobenius_map(1);
    const hg6_Fq6 f1 = hg6_exp_by_z(f0);
    const hg6_Fq6 f1p = f1.Frobenius_map(1);
    const hg6_Fq6 f2 = hg6_exp_by_z(f1);
    const hg6_Fq6 f2p = f2.Frobenius_map(1);
    const hg6_Fq6 f3 = hg6_exp_by_z(f2);
    const hg6_Fq6 f3p = f3.Frobenius_map(1);
    const hg6_Fq6 f4 = hg6_exp_by_z(f3);
    const hg6_Fq6 f4p = f4.Frobenius_map(1);
    const hg6_Fq6 f5 = hg6_exp_by_z(f4);
    const hg6_Fq6 f5p = f5.Frobenius_map(1);
    const hg6_Fq6 f6 = hg6_exp_by_z(f5);
    const hg6_Fq6 f6p = f6.Frobenius_map(1);
    const hg6_Fq6 f7 = hg6_exp_by_z(f6);
    const hg6_Fq6 f7p = f7.Frobenius_map(1);

    // 4
    const hg6_Fq6 f8p = hg6_exp_by_z(f7p);
    const hg6_Fq6 f9p = hg6_exp_by_z(f8p);

    // 5
    const hg6_Fq6 result1 = f3p * f6p * f5p.Frobenius_map(3);

    // 6
    const hg6_Fq6 result2 = result1.squared();
    const hg6_Fq6 f4_2p = f4 * f2p;
    const hg6_Fq6 result3 = result2 * f5 * f0p * (f0 * f1 * f3 * f4_2p * f8p).Frobenius_map(3);

    // 7
    const hg6_Fq6 result4 = result3.squared();
    const hg6_Fq6 result5 = result4 * f9p * f7.Frobenius_map(3);

    // 8
    const hg6_Fq6 result6 = result5.squared();
    const hg6_Fq6 f2_4p = f2 * f4p;
    const hg6_Fq6 f4_2p_5p = f4_2p * f5p;
    const hg6_Fq6 result7 = result6 * f4_2p_5p * f6 * f7p * (f2_4p * f3 * f3p).Frobenius_map(3);

    // 9
    const hg6_Fq6 result8 = result7.squared();
    const hg6_Fq6 result9 = result8 * f0 * f7 * f1p * (f0p * f9p).Frobenius_map(3);

    // 10
    const hg6_Fq6 result10 = result9.squared();
    const hg6_Fq6 f6p_8p = f6p * f8p;
    const hg6_Fq6 f5_7p = f5 * f7p;
    const hg6_Fq6 result11 = result10 * f5_7p * f2p * (f6p_8p).Frobenius_map(3);

    // 11
    const hg6_Fq6 result12 = result11.squared();
    const hg6_Fq6 f3_6 = f3 * f6;
    const hg6_Fq6 f1_7 = f1 * f7;
    const hg6_Fq6 result13 = result12 * f3_6 * f9p * (f1_7 * f2).Frobenius_map(3);

    // 12
    const hg6_Fq6 result14 = result13.squared();
    const hg6_Fq6 result15 = result14 * f0 * f0p * f3p * f5p * (f4_2p * f5_7p * f6p_8p).Frobenius_map(3);

    // 13
    const hg6_Fq6 result16 = result15.squared();
    const hg6_Fq6 result17 = result16 * f1p * (f3_6).Frobenius_map(3);

    // 14
    const hg6_Fq6 result18 = result17.squared();
    const hg6_Fq6 result19 = result18 * f1_7 * f5_7p * f0p * (f2_4p * f4_2p_5p * f9p).Frobenius_map(3);

    leave_block("Call to hg6_final_exponentiation_last_chunk");

    return result19;
}

hg6_GT hg6_final_exponentiation(const hg6_Fq6 &elt)
{
    enter_block("Call to hg6_final_exponentiation");
    const hg6_Fq6 elt_inv = elt.inverse();
    hg6_Fq6 elt_to_first_chunk = hg6_final_exponentiation_first_chunk(elt, elt_inv);
    hg6_GT result = hg6_final_exponentiation_last_chunk(elt_to_first_chunk);
    // hg6_GT result = elt^hg6_final_exponent;
    leave_block("Call to hg6_final_exponentiation");

    return result;
}

/* ate pairing */

void doubling_step_for_miller_loop(const hg6_Fq two_inv,
                                           hg6_G2 &current,
                                           hg6_ate_ell_coeffs &c)
{
    const hg6_Fq X = current.X_, Y = current.Y_, Z = current.Z_;

    const hg6_Fq A = two_inv * (X * Y);                     // A = X1 * Y1 / 2
    const hg6_Fq B = Y.squared();                           // B = Y1^2
    const hg6_Fq C = Z.squared();                           // C = Z1^2
    const hg6_Fq D = C+C+C;                                 // D = 3 * C
    const hg6_Fq E = hg6_twist_coeff_b * D;           // E = twist_b * D
    const hg6_Fq F = E+E+E;                                 // F = 3 * E
    const hg6_Fq G = two_inv * (B+F);                       // G = (B+F)/2
    const hg6_Fq H = (Y+Z).squared() - (B+C);               // H = (Y1+Z1)^2-(B+C)
    const hg6_Fq I = E-B;                                   // I = E-B
    const hg6_Fq J = X.squared();                           // J = X1^2
    const hg6_Fq E_squared = E.squared();                   // E_squared = E^2

    current.X_ = A * (B-F);                                       // X3 = A * (B-F)
    current.Y_ = G.squared() - (E_squared+E_squared+E_squared);   // Y3 = G^2 - 3*E^2
    current.Z_ = B * H;                                           // Z3 = B * H
    c.ell_0 = hg6_twist * I;                               // ell_0 = xi * I
    c.ell_VW = -H;                                               // ell_VW = - H (later: * yP)
    c.ell_VV = J+J+J;                                            // ell_VV = 3*J (later: * xP)
}

void mixed_addition_step_for_miller_loop(const hg6_G2 base,
                                                 hg6_G2 &current,
                                                 hg6_ate_ell_coeffs &c)
{
    const hg6_Fq X1 = current.X_, Y1 = current.Y_, Z1 = current.Z_;
    const hg6_Fq &x2 = base.X_, &y2 = base.Y_;

    const hg6_Fq D = X1 - x2 * Z1;          // D = X1 - X2*Z1
    const hg6_Fq E = Y1 - y2 * Z1;          // E = Y1 - Y2*Z1
    const hg6_Fq F = D.squared();           // F = D^2
    const hg6_Fq G = E.squared();           // G = E^2
    const hg6_Fq H = D*F;                   // H = D*F
    const hg6_Fq I = X1 * F;                // I = X1 * F
    const hg6_Fq J = H + Z1*G - (I+I);      // J = H + Z1*G - (I+I)

    current.X_ = D * J;                           // X3 = D*J
    current.Y_ = E * (I-J)-(H * Y1);              // Y3 = E*(I-J)-(H*Y1)
    current.Z_ = Z1 * H;                          // Z3 = Z1*H
    c.ell_0 = hg6_twist * (E * x2 - D * y2); // ell_0 = xi * (E * X2 - D * Y2)
    c.ell_VV = - E;                              // ell_VV = - E (later: * xP)
    c.ell_VW = D;                                // ell_VW = D (later: * yP    )

}

hg6_ate_G1_precomp hg6_ate_precompute_G1(const hg6_G1& P)
{
    enter_block("Call to hg6_ate_precompute_G1");

    hg6_G1 Pcopy = P;
    Pcopy.to_affine_coordinates();

    hg6_ate_G1_precomp result;
    result.PX = Pcopy.X();
    result.PY = Pcopy.Y();

    leave_block("Call to hg6_ate_precompute_G1");
    return result;
}

hg6_ate_G2_precomp hg6_ate_precompute_G2(const hg6_G2& Q, const bigint<hg6_Fq::num_limbs> &loop_count)
{
    enter_block("Call to hg6_ate_precompute_G2");

    hg6_G2 Qcopy(Q);
    Qcopy.to_affine_coordinates();

    hg6_Fq two_inv = (hg6_Fq("2").inverse()); // could add to global params if needed

    hg6_ate_G2_precomp result;
    result.QX = Qcopy.X_;
    result.QY = Qcopy.Y_;

    hg6_G2 R;
    R.X_ = Qcopy.X_;
    R.Y_ = Qcopy.Y_;
    R.Z_ = hg6_Fq::one();

    bool found_nonzero = false;
    hg6_ate_ell_coeffs c;

    std::vector<long> NAF = find_wnaf(1, loop_count);
    for (long i = NAF.size() - 1; i >= 0; --i)
    {
        if (!found_nonzero)
        {
            /* this skips the MSB itself */
            found_nonzero |= (NAF[i] != 0);
            continue;
        }

        doubling_step_for_miller_loop(two_inv, R, c);
        result.coeffs.push_back(c);

        if (NAF[i] != 0)
        {
          if (NAF[i] > 0)
          {
            mixed_addition_step_for_miller_loop(Qcopy, R, c);
          }
          else
          {
            mixed_addition_step_for_miller_loop(-Qcopy, R, c);
          }
          result.coeffs.push_back(c);
        }
    }

    leave_block("Call to hg6_ate_precompute_G2");
    return result;
}

hg6_Fq6 hg6_ate_miller_loop(const hg6_ate_G1_precomp &prec_P,
                                     const hg6_ate_G2_precomp &prec_Q_1,
                                     const hg6_ate_G2_precomp &prec_Q_2)
{
    enter_block("Call to hg6_ate_miller_loop f_{u+1,Q}(P)");

    // f_{u+1,Q}(P)
    hg6_Fq6 f_1 = hg6_Fq6::one();

    bool found_nonzero_1 = false;
    size_t idx_1 = 0;

    const bigint<hg6_Fq::num_limbs> &loop_count_1 = hg6_ate_loop_count1;
    hg6_ate_ell_coeffs c_1;

    std::vector<long> NAF_1 = find_wnaf(1, loop_count_1);
    for (long i = NAF_1.size() - 1; i >= 0; --i)
    {
        if (!found_nonzero_1)
        {
            /* this skips the MSB itself */
            found_nonzero_1 |= (NAF_1[i] != 0);
            continue;
        }

        /* code below gets executed for all bits (EXCEPT the MSB itself) of
           hg6_param_p (skipping leading zeros) in MSB to LSB
           order */

        c_1 = prec_Q_1.coeffs[idx_1++];
        f_1 = f_1.squared();
        f_1 = f_1.mul_by_024(c_1.ell_0, prec_P.PY * c_1.ell_VW, prec_P.PX * c_1.ell_VV);

        if (NAF_1[i] != 0)
        {
            c_1 = prec_Q_1.coeffs[idx_1++];
            f_1 = f_1.mul_by_024(c_1.ell_0, prec_P.PY * c_1.ell_VW, prec_P.PX * c_1.ell_VV);
        }

    }

    if (hg6_ate_is_loop_count_neg)
    {
    	f_1 = f_1.inverse();
    }
    leave_block("Call to hg6_ate_miller_loop f_{u+1,Q}(P)");

    enter_block("Call to hg6_ate_miller_loop f_{u^3-u^2-u,Q}(P)");

    // f_{u^3-u^2-u,Q}(P)
    hg6_Fq6 f_2 = hg6_Fq6::one();

    bool found_nonzero_2 = false;
    size_t idx_2 = 0;

    const bigint<hg6_Fq::num_limbs> &loop_count_2 = hg6_ate_loop_count2;
    hg6_ate_ell_coeffs c_2;

    std::vector<long> NAF_2 = find_wnaf(1, loop_count_2);
    for (long i = NAF_2.size() - 1; i >= 0; --i)
    {
        if (!found_nonzero_2)
        {
            /* this skips the MSB itself */
            found_nonzero_2 |= (NAF_2[i] != 0);
            continue;
        }

        /* code below gets executed for all bits (EXCEPT the MSB itself) of
           hg6_param_p (skipping leading zeros) in MSB to LSB
           order */

        c_2 = prec_Q_2.coeffs[idx_2++];
        f_2 = f_2.squared();
        f_2 = f_2.mul_by_024(c_2.ell_0, prec_P.PY * c_2.ell_VW, prec_P.PX * c_2.ell_VV);

        if (NAF_2[i] != 0)
        {
            c_2 = prec_Q_2.coeffs[idx_2++];
            f_2 = f_2.mul_by_024(c_2.ell_0, prec_P.PY * c_2.ell_VW, prec_P.PX * c_2.ell_VV);
        }

    }

    if (hg6_ate_is_loop_count_neg)
    {
    	f_2 = f_2.inverse();
    }

    leave_block("Call to hg6_ate_miller_loop f_{u^3-u^2-u,Q}(P)");

    f_2 = f_2.Frobenius_map(1);

    return f_1 * f_2;
}

/*
hg6_Fq6 hg6_ate_double_miller_loop(const hg6_ate_G1_precomp &prec_P1,
                                     const hg6_ate_G2_precomp &prec_Q1,
                                     const hg6_ate_G1_precomp &prec_P2,
                                     const hg6_ate_G2_precomp &prec_Q2)
{
    enter_block("Call to hg6_ate_double_miller_loop");

    hg6_Fq6 f = hg6_Fq6::one();

    bool found_one = false;
    size_t idx = 0;

    const bigint<hg6_Fq::num_limbs> &loop_count = hg6_ate_loop_count;
    for (long i = loop_count.max_bits(); i >= 0; --i)
    {
        const bool bit = loop_count.test_bit(i);
        if (!found_one)
        {
            // this skips the MSB itself
            found_one |= bit;
            continue;
        }

        // code below gets executed for all bits (EXCEPT the MSB itself) of
        // hg6_param_p (skipping leading zeros) in MSB to LSB
        // order

        hg6_ate_ell_coeffs c1 = prec_Q1.coeffs[idx];
        hg6_ate_ell_coeffs c2 = prec_Q2.coeffs[idx];
        ++idx;

        f = f.squared();

        f = f.mul_by_024(c1.ell_0, prec_P1.PY * c1.ell_VW, prec_P1.PX * c1.ell_VV);
        f = f.mul_by_024(c2.ell_0, prec_P2.PY * c2.ell_VW, prec_P2.PX * c2.ell_VV);

        if (bit)
        {
            hg6_ate_ell_coeffs c1 = prec_Q1.coeffs[idx];
            hg6_ate_ell_coeffs c2 = prec_Q2.coeffs[idx];
            ++idx;

            f = f.mul_by_024(c1.ell_0, prec_P1.PY * c1.ell_VW, prec_P1.PX * c1.ell_VV);
            f = f.mul_by_024(c2.ell_0, prec_P2.PY * c2.ell_VW, prec_P2.PX * c2.ell_VV);
        }
    }

    if (hg6_ate_is_loop_count_neg)
    {
    	f = f.inverse();
    }

    leave_block("Call to hg6_ate_double_miller_loop");

    return f;
}
*/

hg6_Fq6 hg6_ate_pairing(const hg6_G1& P, const hg6_G2 &Q)
{
    enter_block("Call to hg6_ate_pairing");
    hg6_ate_G1_precomp prec_P = hg6_ate_precompute_G1(P);
    hg6_ate_G2_precomp prec_Q_1 = hg6_ate_precompute_G2(Q, hg6_ate_loop_count1);
    hg6_ate_G2_precomp prec_Q_2 = hg6_ate_precompute_G2(Q, hg6_ate_loop_count2);
    hg6_Fq6 result = hg6_ate_miller_loop(prec_P, prec_Q_1, prec_Q_2);
    leave_block("Call to hg6_ate_pairing");
    return result;
}

hg6_GT hg6_ate_reduced_pairing(const hg6_G1 &P, const hg6_G2 &Q)
{
    enter_block("Call to hg6_ate_reduced_pairing");
    const hg6_Fq6 f = hg6_ate_pairing(P, Q);
    const hg6_GT result = hg6_final_exponentiation(f);
    leave_block("Call to hg6_ate_reduced_pairing");
    return result;
}

/* choice of pairing */

hg6_G1_precomp hg6_precompute_G1(const hg6_G1& P)
{
    return hg6_ate_precompute_G1(P);
}

hg6_G2_precomp hg6_precompute_G2(const hg6_G2& Q, const bigint<hg6_Fq::num_limbs> &loop_count)
{
    return hg6_ate_precompute_G2(Q, loop_count);
}

hg6_Fq6 hg6_miller_loop(const hg6_G1_precomp &prec_P,
                          const hg6_G2_precomp &prec_Q_1,
                          const hg6_G2_precomp &prec_Q_2)
{
    return hg6_ate_miller_loop(prec_P, prec_Q_1, prec_Q_2);
}

/*
hg6_Fq6 hg6_double_miller_loop(const hg6_G1_precomp &prec_P1,
                                 const hg6_G2_precomp &prec_Q1,
                                 const hg6_G1_precomp &prec_P2,
                                 const hg6_G2_precomp &prec_Q2)
{
    return hg6_ate_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}
*/

hg6_Fq6 hg6_pairing(const hg6_G1& P,
                      const hg6_G2 &Q)
{
    return hg6_ate_pairing(P, Q);
}

hg6_GT hg6_reduced_pairing(const hg6_G1 &P,
                             const hg6_G2 &Q)
{
    return hg6_ate_reduced_pairing(P, Q);
}
} // libff
