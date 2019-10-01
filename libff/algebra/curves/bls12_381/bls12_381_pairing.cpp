#include <cassert>

#include <libff/algebra/curves/bls12_381/bls12_381_g1.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_g2.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_init.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_pairing.hpp>
#include <libff/common/profiling.hpp>

namespace libff {

bool bls12_381_ate_G1_precomp::operator==(const bls12_381_ate_G1_precomp &other) const
{
    return (this->PX == other.PX &&
            this->PY == other.PY);
}

std::ostream& operator<<(std::ostream &out, const bls12_381_ate_G1_precomp &prec_P)
{
    out << prec_P.PX << OUTPUT_SEPARATOR << prec_P.PY;

    return out;
}

std::istream& operator>>(std::istream &in, bls12_381_ate_G1_precomp &prec_P)
{
    in >> prec_P.PX;
    consume_OUTPUT_SEPARATOR(in);
    in >> prec_P.PY;

    return in;
}

bool  bls12_381_ate_ell_coeffs::operator==(const bls12_381_ate_ell_coeffs &other) const
{
    return (this->ell_0 == other.ell_0 &&
            this->ell_VW == other.ell_VW &&
            this->ell_VV == other.ell_VV);
}

std::ostream& operator<<(std::ostream &out, const bls12_381_ate_ell_coeffs &c)
{
    out << c.ell_0 << OUTPUT_SEPARATOR << c.ell_VW << OUTPUT_SEPARATOR << c.ell_VV;
    return out;
}

std::istream& operator>>(std::istream &in, bls12_381_ate_ell_coeffs &c)
{
    in >> c.ell_0;
    consume_OUTPUT_SEPARATOR(in);
    in >> c.ell_VW;
    consume_OUTPUT_SEPARATOR(in);
    in >> c.ell_VV;

    return in;
}

bool bls12_381_ate_G2_precomp::operator==(const bls12_381_ate_G2_precomp &other) const
{
    return (this->QX == other.QX &&
            this->QY == other.QY &&
            this->coeffs == other.coeffs);
}

std::ostream& operator<<(std::ostream& out, const bls12_381_ate_G2_precomp &prec_Q)
{
    out << prec_Q.QX << OUTPUT_SEPARATOR << prec_Q.QY << "\n";
    out << prec_Q.coeffs.size() << "\n";
    for (const bls12_381_ate_ell_coeffs &c : prec_Q.coeffs)
    {
        out << c << OUTPUT_NEWLINE;
    }
    return out;
}

std::istream& operator>>(std::istream& in, bls12_381_ate_G2_precomp &prec_Q)
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
        bls12_381_ate_ell_coeffs c;
        in >> c;
        consume_OUTPUT_NEWLINE(in);
        prec_Q.coeffs.emplace_back(c);
    }

    return in;
}

/* final exponentiations */

bls12_381_Fq12 bls12_381_final_exponentiation_first_chunk(const bls12_381_Fq12 &elt)
{
    enter_block("Call to bls12_381_final_exponentiation_first_chunk");

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

    const bls12_381_Fq12 A = bls12_381_Fq12(elt.c0,-elt.c1);
    const bls12_381_Fq12 B = elt.inverse();
    const bls12_381_Fq12 C = A * B;
    const bls12_381_Fq12 D = C.Frobenius_map(2);
    const bls12_381_Fq12 result = D * C;

    leave_block("Call to bls12_381_final_exponentiation_first_chunk");

    return result;
}

bls12_381_Fq12 bls12_381_exp_by_neg_z(const bls12_381_Fq12 &elt)
{
    enter_block("Call to bls12_381_exp_by_neg_z");

    bls12_381_Fq12 result = elt.cyclotomic_exp(bls12_381_final_exponent_z);
    if (!bls12_381_final_exponent_is_z_neg)
    {
        result = result.unitary_inverse();
    }

    leave_block("Call to bls12_381_exp_by_neg_z");

    return result;
}

bls12_381_Fq12 bls12_381_final_exponentiation_last_chunk(const bls12_381_Fq12 &elt)
{
    enter_block("Call to bls12_381_final_exponentiation_last_chunk");

    /*
      Follows Laura Fuentes-Castaneda et al. "Faster hashing to G2"
      by computing:

      result = elt^(q^3 * (12*z^3 + 6z^2 + 4z - 1) +
                    q^2 * (12*z^3 + 6z^2 + 6z) +
                    q   * (12*z^3 + 6z^2 + 4z) +
                    1   * (12*z^3 + 12z^2 + 6z + 1))
      which equals

      result = elt^( 2z * ( 6z^2 + 3z + 1 ) * (q^4 - q^2 + 1)/r ).

      Using the following addition chain:

      A = exp_by_neg_z(elt)  // = elt^(-z)
      B = A^2                // = elt^(-2*z)
      C = B^2                // = elt^(-4*z)
      D = C * B              // = elt^(-6*z)
      E = exp_by_neg_z(D)    // = elt^(6*z^2)
      F = E^2                // = elt^(12*z^2)
      G = epx_by_neg_z(F)    // = elt^(-12*z^3)
      H = conj(D)            // = elt^(6*z)
      I = conj(G)            // = elt^(12*z^3)
      J = I * E              // = elt^(12*z^3 + 6*z^2)
      K = J * H              // = elt^(12*z^3 + 6*z^2 + 6*z)
      L = K * B              // = elt^(12*z^3 + 6*z^2 + 4*z)
      M = K * E              // = elt^(12*z^3 + 12*z^2 + 6*z)
      N = M * elt            // = elt^(12*z^3 + 12*z^2 + 6*z + 1)
      O = L.Frobenius_map(1) // = elt^(q*(12*z^3 + 6*z^2 + 4*z))
      P = O * N              // = elt^(q*(12*z^3 + 6*z^2 + 4*z) * (12*z^3 + 12*z^2 + 6*z + 1))
      Q = K.Frobenius_map(2) // = elt^(q^2 * (12*z^3 + 6*z^2 + 6*z))
      R = Q * P              // = elt^(q^2 * (12*z^3 + 6*z^2 + 6*z) + q*(12*z^3 + 6*z^2 + 4*z) * (12*z^3 + 12*z^2 + 6*z + 1))
      S = conj(elt)          // = elt^(-1)
      T = S * L              // = elt^(12*z^3 + 6*z^2 + 4*z - 1)
      U = T.Frobenius_map(3) // = elt^(q^3(12*z^3 + 6*z^2 + 4*z - 1))
      V = U * R              // = elt^(q^3(12*z^3 + 6*z^2 + 4*z - 1) + q^2 * (12*z^3 + 6*z^2 + 6*z) + q*(12*z^3 + 6*z^2 + 4*z) * (12*z^3 + 12*z^2 + 6*z + 1))
      result = V

    */

    const bls12_381_Fq12 A = bls12_381_exp_by_neg_z(elt);
    const bls12_381_Fq12 B = A.cyclotomic_squared();
    const bls12_381_Fq12 C = B.cyclotomic_squared();
    const bls12_381_Fq12 D = C * B;
    const bls12_381_Fq12 E = bls12_381_exp_by_neg_z(D);
    const bls12_381_Fq12 F = E.cyclotomic_squared();
    const bls12_381_Fq12 G = bls12_381_exp_by_neg_z(F);
    const bls12_381_Fq12 H = D.unitary_inverse();
    const bls12_381_Fq12 I = G.unitary_inverse();
    const bls12_381_Fq12 J = I * E;
    const bls12_381_Fq12 K = J * H;
    const bls12_381_Fq12 L = K * B;
    const bls12_381_Fq12 M = K * E;
    const bls12_381_Fq12 N = M * elt;
    const bls12_381_Fq12 O = L.Frobenius_map(1);
    const bls12_381_Fq12 P = O * N;
    const bls12_381_Fq12 Q = K.Frobenius_map(2);
    const bls12_381_Fq12 R = Q * P;
    const bls12_381_Fq12 S = elt.unitary_inverse();
    const bls12_381_Fq12 T = S * L;
    const bls12_381_Fq12 U = T.Frobenius_map(3);
    const bls12_381_Fq12 V = U * R;

    const bls12_381_Fq12 result = V;

    leave_block("Call to bls12_381_final_exponentiation_last_chunk");

    return result;
}

bls12_381_GT bls12_381_final_exponentiation(const bls12_381_Fq12 &elt)
{
    enter_block("Call to bls12_381_final_exponentiation");
    /* OLD naive version:
        bls12_381_GT result = elt^bls12_381_final_exponent;
    */
    bls12_381_Fq12 A = bls12_381_final_exponentiation_first_chunk(elt);
    bls12_381_GT result = bls12_381_final_exponentiation_last_chunk(A);

    leave_block("Call to bls12_381_final_exponentiation");
    return result;
}

/* ate pairing */

void doubling_step_for_flipped_miller_loop(const bls12_381_Fq two_inv,
                                           bls12_381_G2 &current,
                                           bls12_381_ate_ell_coeffs &c)
{
    const bls12_381_Fq2 X = current.X, Y = current.Y, Z = current.Z;

    const bls12_381_Fq2 A = two_inv * (X * Y);                     // A = X1 * Y1 / 2
    const bls12_381_Fq2 B = Y.squared();                           // B = Y1^2
    const bls12_381_Fq2 C = Z.squared();                           // C = Z1^2
    const bls12_381_Fq2 D = C+C+C;                                 // D = 3 * C
    const bls12_381_Fq2 E = bls12_381_twist_coeff_b * D;             // E = twist_b * D
    const bls12_381_Fq2 F = E+E+E;                                 // F = 3 * E
    const bls12_381_Fq2 G = two_inv * (B+F);                       // G = (B+F)/2
    const bls12_381_Fq2 H = (Y+Z).squared() - (B+C);               // H = (Y1+Z1)^2-(B+C)
    const bls12_381_Fq2 I = E-B;                                   // I = E-B
    const bls12_381_Fq2 J = X.squared();                           // J = X1^2
    const bls12_381_Fq2 E_squared = E.squared();                   // E_squared = E^2

    current.X = A * (B-F);                                       // X3 = A * (B-F)
    current.Y = G.squared() - (E_squared+E_squared+E_squared);   // Y3 = G^2 - 3*E^2
    current.Z = B * H;                                           // Z3 = B * H
    c.ell_0 = bls12_381_twist * I;                                 // ell_0 = xi * I
    c.ell_VW = -H;                                               // ell_VW = - H (later: * yP)
    c.ell_VV = J+J+J;                                            // ell_VV = 3*J (later: * xP)
}

void mixed_addition_step_for_flipped_miller_loop(const bls12_381_G2 base,
                                                 bls12_381_G2 &current,
                                                 bls12_381_ate_ell_coeffs &c)
{
    const bls12_381_Fq2 X1 = current.X, Y1 = current.Y, Z1 = current.Z;
    const bls12_381_Fq2 &x2 = base.X, &y2 = base.Y;

    const bls12_381_Fq2 D = X1 - x2 * Z1;          // D = X1 - X2*Z1
    const bls12_381_Fq2 E = Y1 - y2 * Z1;          // E = Y1 - Y2*Z1
    const bls12_381_Fq2 F = D.squared();           // F = D^2
    const bls12_381_Fq2 G = E.squared();           // G = E^2
    const bls12_381_Fq2 H = D*F;                   // H = D*F
    const bls12_381_Fq2 I = X1 * F;                // I = X1 * F
    const bls12_381_Fq2 J = H + Z1*G - (I+I);      // J = H + Z1*G - (I+I)

    current.X = D * J;                           // X3 = D*J
    current.Y = E * (I-J)-(H * Y1);              // Y3 = E*(I-J)-(H*Y1)
    current.Z = Z1 * H;                          // Z3 = Z1*H
    c.ell_0 = bls12_381_twist * (E * x2 - D * y2); // ell_0 = xi * (E * X2 - D * Y2)
    c.ell_VV = - E;                              // ell_VV = - E (later: * xP)
    c.ell_VW = D;                                // ell_VW = D (later: * yP    )
}

bls12_381_ate_G1_precomp bls12_381_ate_precompute_G1(const bls12_381_G1& P)
{
    enter_block("Call to bls12_381_ate_precompute_G1");

    bls12_381_G1 Pcopy = P;
    Pcopy.to_affine_coordinates();

    bls12_381_ate_G1_precomp result;
    result.PX = Pcopy.X;
    result.PY = Pcopy.Y;

    leave_block("Call to bls12_381_ate_precompute_G1");
    return result;
}

bls12_381_ate_G2_precomp bls12_381_ate_precompute_G2(const bls12_381_G2& Q)
{
    enter_block("Call to bls12_381_ate_precompute_G2");

    bls12_381_G2 Qcopy(Q);
    Qcopy.to_affine_coordinates();

    bls12_381_Fq two_inv = (bls12_381_Fq("2").inverse()); // could add to global params if needed

    bls12_381_ate_G2_precomp result;
    result.QX = Qcopy.X;
    result.QY = Qcopy.Y;

    bls12_381_G2 R;
    R.X = Qcopy.X;
    R.Y = Qcopy.Y;
    R.Z = bls12_381_Fq2::one();

    const bigint<bls12_381_Fq::num_limbs> &loop_count = bls12_381_ate_loop_count;
    bool found_one = false;
    bls12_381_ate_ell_coeffs c;

    for (long i = loop_count.max_bits(); i >= 0; --i)
    {
        const bool bit = loop_count.test_bit(i);
        if (!found_one)
        {
            /* this skips the MSB itself */
            found_one |= bit;
            continue;
        }

        doubling_step_for_flipped_miller_loop(two_inv, R, c);
        result.coeffs.push_back(c);

        if (bit)
        {
            mixed_addition_step_for_flipped_miller_loop(Qcopy, R, c);
            result.coeffs.push_back(c);
        }
    }

    bls12_381_G2 Q1 = Qcopy.mul_by_q();
    assert(Q1.Z == bls12_381_Fq2::one());
    bls12_381_G2 Q2 = Q1.mul_by_q();
    assert(Q2.Z == bls12_381_Fq2::one());

    if (bls12_381_ate_is_loop_count_neg)
    {
        R.Y = - R.Y;
    }
    Q2.Y = - Q2.Y;

    mixed_addition_step_for_flipped_miller_loop(Q1, R, c);
    result.coeffs.push_back(c);

    mixed_addition_step_for_flipped_miller_loop(Q2, R, c);
    result.coeffs.push_back(c);

    leave_block("Call to bls12_381_ate_precompute_G2");
    return result;
}

bls12_381_Fq12 bls12_381_ate_miller_loop(const bls12_381_ate_G1_precomp &prec_P,
                                     const bls12_381_ate_G2_precomp &prec_Q)
{
    enter_block("Call to bls12_381_ate_miller_loop");

    bls12_381_Fq12 f = bls12_381_Fq12::one();

    bool found_one = false;
    size_t idx = 0;

    const bigint<bls12_381_Fq::num_limbs> &loop_count = bls12_381_ate_loop_count;
    bls12_381_ate_ell_coeffs c;

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
           bls12_381_param_p (skipping leading zeros) in MSB to LSB
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

    if (bls12_381_ate_is_loop_count_neg)
    {
    	f = f.inverse();
    }

    c = prec_Q.coeffs[idx++];
    f = f.mul_by_024(c.ell_0,prec_P.PY * c.ell_VW,prec_P.PX * c.ell_VV);

    c = prec_Q.coeffs[idx++];
    f = f.mul_by_024(c.ell_0,prec_P.PY * c.ell_VW,prec_P.PX * c.ell_VV);

    leave_block("Call to bls12_381_ate_miller_loop");
    return f;
}

bls12_381_Fq12 bls12_381_ate_double_miller_loop(const bls12_381_ate_G1_precomp &prec_P1,
                                     const bls12_381_ate_G2_precomp &prec_Q1,
                                     const bls12_381_ate_G1_precomp &prec_P2,
                                     const bls12_381_ate_G2_precomp &prec_Q2)
{
    enter_block("Call to bls12_381_ate_double_miller_loop");

    bls12_381_Fq12 f = bls12_381_Fq12::one();

    bool found_one = false;
    size_t idx = 0;

    const bigint<bls12_381_Fq::num_limbs> &loop_count = bls12_381_ate_loop_count;
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
           bls12_381_param_p (skipping leading zeros) in MSB to LSB
           order */

        bls12_381_ate_ell_coeffs c1 = prec_Q1.coeffs[idx];
        bls12_381_ate_ell_coeffs c2 = prec_Q2.coeffs[idx];
        ++idx;

        f = f.squared();

        f = f.mul_by_024(c1.ell_0, prec_P1.PY * c1.ell_VW, prec_P1.PX * c1.ell_VV);
        f = f.mul_by_024(c2.ell_0, prec_P2.PY * c2.ell_VW, prec_P2.PX * c2.ell_VV);

        if (bit)
        {
            bls12_381_ate_ell_coeffs c1 = prec_Q1.coeffs[idx];
            bls12_381_ate_ell_coeffs c2 = prec_Q2.coeffs[idx];
            ++idx;

            f = f.mul_by_024(c1.ell_0, prec_P1.PY * c1.ell_VW, prec_P1.PX * c1.ell_VV);
            f = f.mul_by_024(c2.ell_0, prec_P2.PY * c2.ell_VW, prec_P2.PX * c2.ell_VV);
        }
    }

    if (bls12_381_ate_is_loop_count_neg)
    {
    	f = f.inverse();
    }

    bls12_381_ate_ell_coeffs c1 = prec_Q1.coeffs[idx];
    bls12_381_ate_ell_coeffs c2 = prec_Q2.coeffs[idx];
    ++idx;
    f = f.mul_by_024(c1.ell_0, prec_P1.PY * c1.ell_VW, prec_P1.PX * c1.ell_VV);
    f = f.mul_by_024(c2.ell_0, prec_P2.PY * c2.ell_VW, prec_P2.PX * c2.ell_VV);

    c1 = prec_Q1.coeffs[idx];
    c2 = prec_Q2.coeffs[idx];
    ++idx;
    f = f.mul_by_024(c1.ell_0, prec_P1.PY * c1.ell_VW, prec_P1.PX * c1.ell_VV);
    f = f.mul_by_024(c2.ell_0, prec_P2.PY * c2.ell_VW, prec_P2.PX * c2.ell_VV);

    leave_block("Call to bls12_381_ate_double_miller_loop");

    return f;
}

bls12_381_Fq12 bls12_381_ate_pairing(const bls12_381_G1& P, const bls12_381_G2 &Q)
{
    enter_block("Call to bls12_381_ate_pairing");
    bls12_381_ate_G1_precomp prec_P = bls12_381_ate_precompute_G1(P);
    bls12_381_ate_G2_precomp prec_Q = bls12_381_ate_precompute_G2(Q);
    bls12_381_Fq12 result = bls12_381_ate_miller_loop(prec_P, prec_Q);
    leave_block("Call to bls12_381_ate_pairing");
    return result;
}

bls12_381_GT bls12_381_ate_reduced_pairing(const bls12_381_G1 &P, const bls12_381_G2 &Q)
{
    enter_block("Call to bls12_381_ate_reduced_pairing");
    const bls12_381_Fq12 f = bls12_381_ate_pairing(P, Q);
    f.print();
    const bls12_381_GT result = bls12_381_final_exponentiation(f);
    leave_block("Call to bls12_381_ate_reduced_pairing");
    return result;
}

/* choice of pairing */

bls12_381_G1_precomp bls12_381_precompute_G1(const bls12_381_G1& P)
{
    return bls12_381_ate_precompute_G1(P);
}

bls12_381_G2_precomp bls12_381_precompute_G2(const bls12_381_G2& Q)
{
    return bls12_381_ate_precompute_G2(Q);
}

bls12_381_Fq12 bls12_381_miller_loop(const bls12_381_G1_precomp &prec_P,
                          const bls12_381_G2_precomp &prec_Q)
{
    return bls12_381_ate_miller_loop(prec_P, prec_Q);
}

bls12_381_Fq12 bls12_381_double_miller_loop(const bls12_381_G1_precomp &prec_P1,
                                 const bls12_381_G2_precomp &prec_Q1,
                                 const bls12_381_G1_precomp &prec_P2,
                                 const bls12_381_G2_precomp &prec_Q2)
{
    return bls12_381_ate_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

bls12_381_Fq12 bls12_381_pairing(const bls12_381_G1& P,
                      const bls12_381_G2 &Q)
{
    return bls12_381_ate_pairing(P, Q);
}

bls12_381_GT bls12_381_reduced_pairing(const bls12_381_G1 &P,
                             const bls12_381_G2 &Q)
{
    return bls12_381_ate_reduced_pairing(P, Q);
}
} // libff
