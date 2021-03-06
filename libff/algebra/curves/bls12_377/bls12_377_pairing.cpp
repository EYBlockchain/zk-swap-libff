/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <cassert>

#include <libff/algebra/curves/bls12_377/bls12_377_g1.hpp>
#include <libff/algebra/curves/bls12_377/bls12_377_g2.hpp>
#include <libff/algebra/curves/bls12_377/bls12_377_init.hpp>
#include <libff/algebra/curves/bls12_377/bls12_377_pairing.hpp>
#include <libff/common/profiling.hpp>

namespace libff {

bool bls12_377_ate_G1_precomp::operator==(const bls12_377_ate_G1_precomp &other) const
{
    return (this->PX == other.PX &&
            this->PY == other.PY);
}

std::ostream& operator<<(std::ostream &out, const bls12_377_ate_G1_precomp &prec_P)
{
    out << prec_P.PX << OUTPUT_SEPARATOR << prec_P.PY;

    return out;
}

std::istream& operator>>(std::istream &in, bls12_377_ate_G1_precomp &prec_P)
{
    in >> prec_P.PX;
    consume_OUTPUT_SEPARATOR(in);
    in >> prec_P.PY;

    return in;
}

bool  bls12_377_ate_ell_coeffs::operator==(const bls12_377_ate_ell_coeffs &other) const
{
    return (this->ell_0 == other.ell_0 &&
            this->ell_VW == other.ell_VW &&
            this->ell_VV == other.ell_VV);
}

std::ostream& operator<<(std::ostream &out, const bls12_377_ate_ell_coeffs &c)
{
    out << c.ell_0 << OUTPUT_SEPARATOR << c.ell_VW << OUTPUT_SEPARATOR << c.ell_VV;
    return out;
}

std::istream& operator>>(std::istream &in, bls12_377_ate_ell_coeffs &c)
{
    in >> c.ell_0;
    consume_OUTPUT_SEPARATOR(in);
    in >> c.ell_VW;
    consume_OUTPUT_SEPARATOR(in);
    in >> c.ell_VV;

    return in;
}

bool bls12_377_ate_G2_precomp::operator==(const bls12_377_ate_G2_precomp &other) const
{
    return (this->QX == other.QX &&
            this->QY == other.QY &&
            this->coeffs == other.coeffs);
}

std::ostream& operator<<(std::ostream& out, const bls12_377_ate_G2_precomp &prec_Q)
{
    out << prec_Q.QX << OUTPUT_SEPARATOR << prec_Q.QY << "\n";
    out << prec_Q.coeffs.size() << "\n";
    for (const bls12_377_ate_ell_coeffs &c : prec_Q.coeffs)
    {
        out << c << OUTPUT_NEWLINE;
    }
    return out;
}

std::istream& operator>>(std::istream& in, bls12_377_ate_G2_precomp &prec_Q)
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
        bls12_377_ate_ell_coeffs c;
        in >> c;
        consume_OUTPUT_NEWLINE(in);
        prec_Q.coeffs.emplace_back(c);
    }

    return in;
}

/* final exponentiations */

bls12_377_Fq12 bls12_377_final_exponentiation_first_chunk(const bls12_377_Fq12 &elt)
{
    enter_block("Call to bls12_377_final_exponentiation_first_chunk");

    /*
      Computes result = elt^((q^6-1)*(q^2+1)).
      Follows, e.g., Beuchat et al page 9, by computing result as follows:
         elt^((q^6-1)*(q^2+1)) = (conj(elt) * elt^(-1))^(q^2+1)
      More precisely:
      A = conj(elt)
      B = elt.inverse()
      C = A * B
      D = C.Frobenius_map(2)
      result = D * C
    */

    const bls12_377_Fq12 A = bls12_377_Fq12(elt.c0,-elt.c1);
    const bls12_377_Fq12 B = elt.inverse();
    const bls12_377_Fq12 C = A * B;
    const bls12_377_Fq12 D = C.Frobenius_map(2);
    const bls12_377_Fq12 result = D * C;

    leave_block("Call to bls12_377_final_exponentiation_first_chunk");

    return result;
}

bls12_377_Fq12 bls12_377_exp_by_z(const bls12_377_Fq12 &elt)
{
    enter_block("Call to bls12_377_exp_by_z");

    bls12_377_Fq12 result = elt.cyclotomic_exp(bls12_377_final_exponent_z);
    if (bls12_377_final_exponent_is_z_neg)
    {
        result = result.unitary_inverse();
    }

    leave_block("Call to bls12_377_exp_by_z");

    return result;
}

bls12_377_Fq12 bls12_377_final_exponentiation_last_chunk(const bls12_377_Fq12 &elt)
{
    enter_block("Call to bls12_377_final_exponentiation_last_chunk");

    const bls12_377_Fq12 A = elt.cyclotomic_squared();   // elt^2
    const bls12_377_Fq12 B = A.unitary_inverse();        // elt^(-2)
    const bls12_377_Fq12 C = bls12_377_exp_by_z(elt);    // elt^z
    const bls12_377_Fq12 D = C.cyclotomic_squared();     // elt^(2z)
    const bls12_377_Fq12 E = B * C;                      // elt^(z-2)
    const bls12_377_Fq12 F = bls12_377_exp_by_z(E);      // elt^(z^2-2z)
    const bls12_377_Fq12 G = bls12_377_exp_by_z(F);      // elt^(z^3-2z^2)
    const bls12_377_Fq12 H = bls12_377_exp_by_z(G);      // elt^(z^4-2z^3)
    const bls12_377_Fq12 I = H * D;                      // elt^(z^4-2z^3+2z)
    const bls12_377_Fq12 J = bls12_377_exp_by_z(I);      // elt^(z^5-2z^4+2z^2)
    const bls12_377_Fq12 K = E.unitary_inverse();        // elt^(-z+2)
    const bls12_377_Fq12 L = K * J;                      // elt^(z^5-2z^4+2z^2) * elt^(-z+2)
    const bls12_377_Fq12 M = elt * L;                    // elt^(z^5-2z^4+2z^2) * elt^(-z+2) * elt
    const bls12_377_Fq12 N = elt.unitary_inverse();      // elt^(-1)
    const bls12_377_Fq12 O = F * elt;                    // elt^(z^2-2z) * elt
    const bls12_377_Fq12 P = O.Frobenius_map(3);         // (elt^(z^2-2z) * elt)^(q^3)
    const bls12_377_Fq12 Q = I * N;                      // elt^(z^4-2z^3+2z) * elt^(-1)
    const bls12_377_Fq12 R = Q.Frobenius_map(1);         // (elt^(z^4-2z^3+2z) * elt^(-1))^q
    const bls12_377_Fq12 S = C * G;                      // elt^(z^3-2z^2) * elt^z
    const bls12_377_Fq12 T = S.Frobenius_map(2);         // (elt^(z^3-2z^2) * elt^z)^(q^2)
    const bls12_377_Fq12 U = T * P;                      // (elt^(z^2-2z) * elt)^(q^3) * (elt^(z^3-2z^2) * elt^z)^(q^2)
    const bls12_377_Fq12 V = U * R;                      // (elt^(z^2-2z) * elt)^(q^3) * (elt^(z^3-2z^2) * elt^z)^(q^2) * (elt^(z^4-2z^3+2z) * elt^(-1))^q
    const bls12_377_Fq12 W = V * M;                      // (elt^(z^2-2z) * elt)^(q^3) * (elt^(z^3-2z^2) * elt^z)^(q^2) * (elt^(z^4-2z^3+2z) * elt^(-1))^q * elt^(z^5-2z^4+2z^2) * elt^(-z+2) * elt

    const bls12_377_Fq12 result = W;

    leave_block("Call to bls12_377_final_exponentiation_last_chunk");

    return result;
}

bls12_377_GT bls12_377_final_exponentiation(const bls12_377_Fq12 &elt)
{
    enter_block("Call to bls12_377_final_exponentiation");
    /* OLD naive version:
        bls12_377_GT result = elt^bls12_377_final_exponent;
    */
    bls12_377_Fq12 A = bls12_377_final_exponentiation_first_chunk(elt);
    bls12_377_GT result = bls12_377_final_exponentiation_last_chunk(A);

    leave_block("Call to bls12_377_final_exponentiation");
    return result;
}

/* ate pairing */

void doubling_step_for_miller_loop(const bls12_377_Fq two_inv,
                                           bls12_377_G2 &current,
                                           bls12_377_ate_ell_coeffs &c)
{
    const bls12_377_Fq2 X = current.X, Y = current.Y, Z = current.Z;

    const bls12_377_Fq2 A = two_inv * (X * Y);                     // A = X1 * Y1 / 2
    const bls12_377_Fq2 B = Y.squared();                           // B = Y1^2
    const bls12_377_Fq2 C = Z.squared();                           // C = Z1^2
    const bls12_377_Fq2 D = C+C+C;                                 // D = 3 * C
    const bls12_377_Fq2 E = bls12_377_twist_coeff_b * D;           // E = twist_b * D
    const bls12_377_Fq2 F = E+E+E;                                 // F = 3 * E
    const bls12_377_Fq2 G = two_inv * (B+F);                       // G = (B+F)/2
    const bls12_377_Fq2 H = (Y+Z).squared() - (B+C);               // H = (Y1+Z1)^2-(B+C)
    const bls12_377_Fq2 I = E-B;                                   // I = E-B
    const bls12_377_Fq2 J = X.squared();                           // J = X1^2
    const bls12_377_Fq2 E_squared = E.squared();                   // E_squared = E^2

    current.X = A * (B-F);                                       // X3 = A * (B-F)
    current.Y = G.squared() - (E_squared+E_squared+E_squared);   // Y3 = G^2 - 3*E^2
    current.Z = B * H;                                           // Z3 = B * H
    c.ell_0 = bls12_377_twist * I;                               // ell_0 = xi * I
    c.ell_VW = -H;                                               // ell_VW = - H (later: * yP)
    c.ell_VV = J+J+J;                                            // ell_VV = 3*J (later: * xP)
}

void mixed_addition_step_for_miller_loop(const bls12_377_G2 base,
                                                 bls12_377_G2 &current,
                                                 bls12_377_ate_ell_coeffs &c)
{
    const bls12_377_Fq2 X1 = current.X, Y1 = current.Y, Z1 = current.Z;
    const bls12_377_Fq2 &x2 = base.X, &y2 = base.Y;

    const bls12_377_Fq2 D = X1 - x2 * Z1;          // D = X1 - X2*Z1
    const bls12_377_Fq2 E = Y1 - y2 * Z1;          // E = Y1 - Y2*Z1
    const bls12_377_Fq2 F = D.squared();           // F = D^2
    const bls12_377_Fq2 G = E.squared();           // G = E^2
    const bls12_377_Fq2 H = D*F;                   // H = D*F
    const bls12_377_Fq2 I = X1 * F;                // I = X1 * F
    const bls12_377_Fq2 J = H + Z1*G - (I+I);      // J = H + Z1*G - (I+I)

    current.X = D * J;                           // X3 = D*J
    current.Y = E * (I-J)-(H * Y1);              // Y3 = E*(I-J)-(H*Y1)
    current.Z = Z1 * H;                          // Z3 = Z1*H
    c.ell_0 = bls12_377_twist * (E * x2 - D * y2); // ell_0 = xi * (E * X2 - D * Y2)
    c.ell_VV = - E;                              // ell_VV = - E (later: * xP)
    c.ell_VW = D;                                // ell_VW = D (later: * yP    )
}

bls12_377_ate_G1_precomp bls12_377_ate_precompute_G1(const bls12_377_G1& P)
{
    enter_block("Call to bls12_377_ate_precompute_G1");

    bls12_377_G1 Pcopy = P;
    Pcopy.to_affine_coordinates();

    bls12_377_ate_G1_precomp result;
    result.PX = Pcopy.X;
    result.PY = Pcopy.Y;

    leave_block("Call to bls12_377_ate_precompute_G1");
    return result;
}

bls12_377_ate_G2_precomp bls12_377_ate_precompute_G2(const bls12_377_G2& Q)
{
    enter_block("Call to bls12_377_ate_precompute_G2");

    bls12_377_G2 Qcopy(Q);
    Qcopy.to_affine_coordinates();

    bls12_377_Fq two_inv = (bls12_377_Fq("2").inverse()); // could add to global params if needed

    bls12_377_ate_G2_precomp result;
    result.QX = Qcopy.X;
    result.QY = Qcopy.Y;

    bls12_377_G2 R;
    R.X = Qcopy.X;
    R.Y = Qcopy.Y;
    R.Z = bls12_377_Fq2::one();

    const bigint<bls12_377_Fq::num_limbs> &loop_count = bls12_377_ate_loop_count;
    bool found_one = false;
    bls12_377_ate_ell_coeffs c;

    for (long i = loop_count.max_bits(); i >= 0; --i)
    {
        const bool bit = loop_count.test_bit(i);
        if (!found_one)
        {
            /* this skips the MSB itself */
            found_one |= bit;
            continue;
        }

        doubling_step_for_miller_loop(two_inv, R, c);
        result.coeffs.push_back(c);

        if (bit)
        {
            mixed_addition_step_for_miller_loop(Qcopy, R, c);
            result.coeffs.push_back(c);
        }
    }

    leave_block("Call to bls12_377_ate_precompute_G2");
    return result;
}

bls12_377_Fq12 bls12_377_ate_miller_loop(const bls12_377_ate_G1_precomp &prec_P,
                                     const bls12_377_ate_G2_precomp &prec_Q)
{
    enter_block("Call to bls12_377_ate_miller_loop");

    bls12_377_Fq12 f = bls12_377_Fq12::one();

    bool found_one = false;
    size_t idx = 0;

    const bigint<bls12_377_Fq::num_limbs> &loop_count = bls12_377_ate_loop_count;
    bls12_377_ate_ell_coeffs c;

    for (long i = loop_count.max_bits(); i >= 0; --i)
    {
        const bool bit = loop_count.test_bit(i);
        if (!found_one)
        {
            /* this skips the MSB itself */
            found_one |= bit;
            continue;
        }

        /* code below gets executed for all bits (EXCEPT the MSB itself) of
           bls12_377_param_p (skipping leading zeros) in MSB to LSB
           order */

        c = prec_Q.coeffs[idx++];
        f = f.squared();
        f = f.mul_by_024(c.ell_0, prec_P.PY * c.ell_VW, prec_P.PX * c.ell_VV);

        if (bit)
        {
            c = prec_Q.coeffs[idx++];
            f = f.mul_by_024(c.ell_0, prec_P.PY * c.ell_VW, prec_P.PX * c.ell_VV);
        }

    }

    if (bls12_377_ate_is_loop_count_neg)
    {
    	f = f.inverse();
    }

    /*
    c = prec_Q.coeffs[idx++];
    f = f.mul_by_024(c.ell_0,prec_P.PY * c.ell_VW,prec_P.PX * c.ell_VV);

    c = prec_Q.coeffs[idx++];
    f = f.mul_by_024(c.ell_0,prec_P.PY * c.ell_VW,prec_P.PX * c.ell_VV);
    */

    leave_block("Call to bls12_377_ate_miller_loop");
    return f;
}

bls12_377_Fq12 bls12_377_ate_double_miller_loop(const bls12_377_ate_G1_precomp &prec_P1,
                                     const bls12_377_ate_G2_precomp &prec_Q1,
                                     const bls12_377_ate_G1_precomp &prec_P2,
                                     const bls12_377_ate_G2_precomp &prec_Q2)
{
    enter_block("Call to bls12_377_ate_double_miller_loop");

    bls12_377_Fq12 f = bls12_377_Fq12::one();

    bool found_one = false;
    size_t idx = 0;

    const bigint<bls12_377_Fq::num_limbs> &loop_count = bls12_377_ate_loop_count;
    for (long i = loop_count.max_bits(); i >= 0; --i)
    {
        const bool bit = loop_count.test_bit(i);
        if (!found_one)
        {
            /* this skips the MSB itself */
            found_one |= bit;
            continue;
        }

        /* code below gets executed for all bits (EXCEPT the MSB itself) of
           bls12_377_param_p (skipping leading zeros) in MSB to LSB
           order */

        bls12_377_ate_ell_coeffs c1 = prec_Q1.coeffs[idx];
        bls12_377_ate_ell_coeffs c2 = prec_Q2.coeffs[idx];
        ++idx;

        f = f.squared();

        f = f.mul_by_024(c1.ell_0, prec_P1.PY * c1.ell_VW, prec_P1.PX * c1.ell_VV);
        f = f.mul_by_024(c2.ell_0, prec_P2.PY * c2.ell_VW, prec_P2.PX * c2.ell_VV);

        if (bit)
        {
            bls12_377_ate_ell_coeffs c1 = prec_Q1.coeffs[idx];
            bls12_377_ate_ell_coeffs c2 = prec_Q2.coeffs[idx];
            ++idx;

            f = f.mul_by_024(c1.ell_0, prec_P1.PY * c1.ell_VW, prec_P1.PX * c1.ell_VV);
            f = f.mul_by_024(c2.ell_0, prec_P2.PY * c2.ell_VW, prec_P2.PX * c2.ell_VV);
        }
    }

    if (bls12_377_ate_is_loop_count_neg)
    {
    	f = f.inverse();
    }

    /*
    bls12_377_ate_ell_coeffs c1 = prec_Q1.coeffs[idx];
    bls12_377_ate_ell_coeffs c2 = prec_Q2.coeffs[idx];
    ++idx;
    f = f.mul_by_024(c1.ell_0, prec_P1.PY * c1.ell_VW, prec_P1.PX * c1.ell_VV);
    f = f.mul_by_024(c2.ell_0, prec_P2.PY * c2.ell_VW, prec_P2.PX * c2.ell_VV);

    c1 = prec_Q1.coeffs[idx];
    c2 = prec_Q2.coeffs[idx];
    ++idx;
    f = f.mul_by_024(c1.ell_0, prec_P1.PY * c1.ell_VW, prec_P1.PX * c1.ell_VV);
    f = f.mul_by_024(c2.ell_0, prec_P2.PY * c2.ell_VW, prec_P2.PX * c2.ell_VV);
    */

    leave_block("Call to bls12_377_ate_double_miller_loop");

    return f;
}

bls12_377_Fq12 bls12_377_ate_pairing(const bls12_377_G1& P, const bls12_377_G2 &Q)
{
    enter_block("Call to bls12_377_ate_pairing");
    bls12_377_ate_G1_precomp prec_P = bls12_377_ate_precompute_G1(P);
    bls12_377_ate_G2_precomp prec_Q = bls12_377_ate_precompute_G2(Q);
    bls12_377_Fq12 result = bls12_377_ate_miller_loop(prec_P, prec_Q);
    leave_block("Call to bls12_377_ate_pairing");
    return result;
}

bls12_377_GT bls12_377_ate_reduced_pairing(const bls12_377_G1 &P, const bls12_377_G2 &Q)
{
    enter_block("Call to bls12_377_ate_reduced_pairing");
    const bls12_377_Fq12 f = bls12_377_ate_pairing(P, Q);
    const bls12_377_GT result = bls12_377_final_exponentiation(f);
    leave_block("Call to bls12_377_ate_reduced_pairing");
    return result;
}

/* choice of pairing */

bls12_377_G1_precomp bls12_377_precompute_G1(const bls12_377_G1& P)
{
    return bls12_377_ate_precompute_G1(P);
}

bls12_377_G2_precomp bls12_377_precompute_G2(const bls12_377_G2& Q)
{
    return bls12_377_ate_precompute_G2(Q);
}

bls12_377_Fq12 bls12_377_miller_loop(const bls12_377_G1_precomp &prec_P,
                          const bls12_377_G2_precomp &prec_Q)
{
    return bls12_377_ate_miller_loop(prec_P, prec_Q);
}

bls12_377_Fq12 bls12_377_double_miller_loop(const bls12_377_G1_precomp &prec_P1,
                                 const bls12_377_G2_precomp &prec_Q1,
                                 const bls12_377_G1_precomp &prec_P2,
                                 const bls12_377_G2_precomp &prec_Q2)
{
    return bls12_377_ate_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

bls12_377_Fq12 bls12_377_pairing(const bls12_377_G1& P,
                      const bls12_377_G2 &Q)
{
    return bls12_377_ate_pairing(P, Q);
}

bls12_377_GT bls12_377_reduced_pairing(const bls12_377_G1 &P,
                             const bls12_377_G2 &Q)
{
    return bls12_377_ate_reduced_pairing(P, Q);
}
} // libff
