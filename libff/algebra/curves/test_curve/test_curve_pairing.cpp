#include <cassert>

#include <libff/algebra/curves/test_curve/test_curve_g1.hpp>
#include <libff/algebra/curves/test_curve/test_curve_g2.hpp>
#include <libff/algebra/curves/test_curve/test_curve_init.hpp>
#include <libff/algebra/curves/test_curve/test_curve_pairing.hpp>
#include <libff/common/profiling.hpp>

namespace libff {

bool test_curve_ate_G1_precomp::operator==(const test_curve_ate_G1_precomp &other) const
{
    return (this->PX == other.PX &&
            this->PY == other.PY);
}

std::ostream& operator<<(std::ostream &out, const test_curve_ate_G1_precomp &prec_P)
{
    out << prec_P.PX << OUTPUT_SEPARATOR << prec_P.PY;

    return out;
}

std::istream& operator>>(std::istream &in, test_curve_ate_G1_precomp &prec_P)
{
    in >> prec_P.PX;
    consume_OUTPUT_SEPARATOR(in);
    in >> prec_P.PY;

    return in;
}

bool  test_curve_ate_ell_coeffs::operator==(const test_curve_ate_ell_coeffs &other) const
{
    return (this->ell_0 == other.ell_0 &&
            this->ell_VW == other.ell_VW &&
            this->ell_VV == other.ell_VV);
}

std::ostream& operator<<(std::ostream &out, const test_curve_ate_ell_coeffs &c)
{
    out << c.ell_0 << OUTPUT_SEPARATOR << c.ell_VW << OUTPUT_SEPARATOR << c.ell_VV;
    return out;
}

std::istream& operator>>(std::istream &in, test_curve_ate_ell_coeffs &c)
{
    in >> c.ell_0;
    consume_OUTPUT_SEPARATOR(in);
    in >> c.ell_VW;
    consume_OUTPUT_SEPARATOR(in);
    in >> c.ell_VV;

    return in;
}

bool test_curve_ate_G2_precomp::operator==(const test_curve_ate_G2_precomp &other) const
{
    return (this->QX == other.QX &&
            this->QY == other.QY &&
            this->coeffs == other.coeffs);
}

std::ostream& operator<<(std::ostream& out, const test_curve_ate_G2_precomp &prec_Q)
{
    out << prec_Q.QX << OUTPUT_SEPARATOR << prec_Q.QY << "\n";
    out << prec_Q.coeffs.size() << "\n";
    for (const test_curve_ate_ell_coeffs &c : prec_Q.coeffs)
    {
        out << c << OUTPUT_NEWLINE;
    }
    return out;
}

std::istream& operator>>(std::istream& in, test_curve_ate_G2_precomp &prec_Q)
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
        test_curve_ate_ell_coeffs c;
        in >> c;
        consume_OUTPUT_NEWLINE(in);
        prec_Q.coeffs.emplace_back(c);
    }

    return in;
}

/* final exponentiations */

test_curve_Fq12 test_curve_final_exponentiation_first_chunk(const test_curve_Fq12 &elt)
{
    enter_block("Call to test_curve_final_exponentiation_first_chunk");

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

    const test_curve_Fq12 A = test_curve_Fq12(elt.c0,-elt.c1);
    const test_curve_Fq12 B = elt.inverse();
    const test_curve_Fq12 C = A * B;
    const test_curve_Fq12 D = C.Frobenius_map(2);
    const test_curve_Fq12 result = D * C;

    leave_block("Call to test_curve_final_exponentiation_first_chunk");

    return result;
}

test_curve_Fq12 test_curve_exp_by_z(const test_curve_Fq12 &elt)
{
    enter_block("Call to test_curve_exp_by_z");

    test_curve_Fq12 result = elt.cyclotomic_exp(test_curve_final_exponent_z);
    if (test_curve_final_exponent_is_z_neg)
    {
        result = result.unitary_inverse();
    }

    leave_block("Call to test_curve_exp_by_z");

    return result;
}

test_curve_Fq12 test_curve_final_exponentiation_last_chunk(const test_curve_Fq12 &elt)
{
    enter_block("Call to test_curve_final_exponentiation_last_chunk");

    const test_curve_Fq12 A = elt.cyclotomic_squared();   // elt^2
    const test_curve_Fq12 B = A.unitary_inverse();        // elt^(-2)
    const test_curve_Fq12 C = test_curve_exp_by_z(elt);    // elt^z
    const test_curve_Fq12 D = C.cyclotomic_squared();     // elt^(2z)
    const test_curve_Fq12 E = B * C;                      // elt^(z-2)
    const test_curve_Fq12 F = test_curve_exp_by_z(E);      // elt^(z^2-2z)
    const test_curve_Fq12 G = test_curve_exp_by_z(F);      // elt^(z^3-2z^2)
    const test_curve_Fq12 H = test_curve_exp_by_z(G);      // elt^(z^4-2z^3)
    const test_curve_Fq12 I = H * D;                      // elt^(z^4-2z^3+2z)
    const test_curve_Fq12 J = test_curve_exp_by_z(I);      // elt^(z^5-2z^4+2z^2)
    const test_curve_Fq12 K = E.unitary_inverse();        // elt^(-z+2)
    const test_curve_Fq12 L = K * J;                      // elt^(z^5-2z^4+2z^2) * elt^(-z+2)
    const test_curve_Fq12 M = elt * L;                    // elt^(z^5-2z^4+2z^2) * elt^(-z+2) * elt
    const test_curve_Fq12 N = elt.unitary_inverse();      // elt^(-1)
    const test_curve_Fq12 O = F * elt;                    // elt^(z^2-2z) * elt
    const test_curve_Fq12 P = O.Frobenius_map(3);         // (elt^(z^2-2z) * elt)^(q^3)
    const test_curve_Fq12 Q = I * N;                      // elt^(z^4-2z^3+2z) * elt^(-1)
    const test_curve_Fq12 R = Q.Frobenius_map(1);         // (elt^(z^4-2z^3+2z) * elt^(-1))^q
    const test_curve_Fq12 S = C * G;                      // elt^(z^3-2z^2) * elt^z
    const test_curve_Fq12 T = S.Frobenius_map(2);         // (elt^(z^3-2z^2) * elt^z)^(q^2)
    const test_curve_Fq12 U = T * P;                      // (elt^(z^2-2z) * elt)^(q^3) * (elt^(z^3-2z^2) * elt^z)^(q^2)
    const test_curve_Fq12 V = U * R;                      // (elt^(z^2-2z) * elt)^(q^3) * (elt^(z^3-2z^2) * elt^z)^(q^2) * (elt^(z^4-2z^3+2z) * elt^(-1))^q
    const test_curve_Fq12 W = V * M;                      // (elt^(z^2-2z) * elt)^(q^3) * (elt^(z^3-2z^2) * elt^z)^(q^2) * (elt^(z^4-2z^3+2z) * elt^(-1))^q * elt^(z^5-2z^4+2z^2) * elt^(-z+2) * elt

    const test_curve_Fq12 result = W;

    leave_block("Call to test_curve_final_exponentiation_last_chunk");

    return result;
}

test_curve_GT test_curve_final_exponentiation(const test_curve_Fq12 &elt)
{
    enter_block("Call to test_curve_final_exponentiation");
    /* OLD naive version:
        test_curve_GT result = elt^test_curve_final_exponent;
    */
    test_curve_Fq12 A = test_curve_final_exponentiation_first_chunk(elt);
    test_curve_GT result = test_curve_final_exponentiation_last_chunk(A);

    leave_block("Call to test_curve_final_exponentiation");
    return result;
}

/* ate pairing */

void doubling_step_for_miller_loop(const test_curve_Fq two_inv,
                                           test_curve_G2 &current,
                                           test_curve_ate_ell_coeffs &c)
{
    const test_curve_Fq2 X = current.X, Y = current.Y, Z = current.Z;

    const test_curve_Fq2 A = two_inv * (X * Y);                     // A = X1 * Y1 / 2
    const test_curve_Fq2 B = Y.squared();                           // B = Y1^2
    const test_curve_Fq2 C = Z.squared();                           // C = Z1^2
    const test_curve_Fq2 D = C+C+C;                                 // D = 3 * C
    const test_curve_Fq2 E = test_curve_twist_coeff_b * D;           // E = twist_b * D
    const test_curve_Fq2 F = E+E+E;                                 // F = 3 * E
    const test_curve_Fq2 G = two_inv * (B+F);                       // G = (B+F)/2
    const test_curve_Fq2 H = (Y+Z).squared() - (B+C);               // H = (Y1+Z1)^2-(B+C)
    const test_curve_Fq2 I = E-B;                                   // I = E-B
    const test_curve_Fq2 J = X.squared();                           // J = X1^2
    const test_curve_Fq2 E_squared = E.squared();                   // E_squared = E^2

    current.X = A * (B-F);                                       // X3 = A * (B-F)
    current.Y = G.squared() - (E_squared+E_squared+E_squared);   // Y3 = G^2 - 3*E^2
    current.Z = B * H;                                           // Z3 = B * H
    c.ell_0 = test_curve_twist * I;                               // ell_0 = xi * I
    c.ell_VW = -H;                                               // ell_VW = - H (later: * yP)
    c.ell_VV = J+J+J;                                            // ell_VV = 3*J (later: * xP)
}

void mixed_addition_step_for_miller_loop(const test_curve_G2 base,
                                                 test_curve_G2 &current,
                                                 test_curve_ate_ell_coeffs &c)
{
    const test_curve_Fq2 X1 = current.X, Y1 = current.Y, Z1 = current.Z;
    const test_curve_Fq2 &x2 = base.X, &y2 = base.Y;

    const test_curve_Fq2 D = X1 - x2 * Z1;          // D = X1 - X2*Z1
    const test_curve_Fq2 E = Y1 - y2 * Z1;          // E = Y1 - Y2*Z1
    const test_curve_Fq2 F = D.squared();           // F = D^2
    const test_curve_Fq2 G = E.squared();           // G = E^2
    const test_curve_Fq2 H = D*F;                   // H = D*F
    const test_curve_Fq2 I = X1 * F;                // I = X1 * F
    const test_curve_Fq2 J = H + Z1*G - (I+I);      // J = H + Z1*G - (I+I)

    current.X = D * J;                           // X3 = D*J
    current.Y = E * (I-J)-(H * Y1);              // Y3 = E*(I-J)-(H*Y1)
    current.Z = Z1 * H;                          // Z3 = Z1*H
    c.ell_0 = test_curve_twist * (E * x2 - D * y2); // ell_0 = xi * (E * X2 - D * Y2)
    c.ell_VV = - E;                              // ell_VV = - E (later: * xP)
    c.ell_VW = D;                                // ell_VW = D (later: * yP    )
}

test_curve_ate_G1_precomp test_curve_ate_precompute_G1(const test_curve_G1& P)
{
    enter_block("Call to test_curve_ate_precompute_G1");

    test_curve_G1 Pcopy = P;
    Pcopy.to_affine_coordinates();

    test_curve_ate_G1_precomp result;
    result.PX = Pcopy.X;
    result.PY = Pcopy.Y;

    leave_block("Call to test_curve_ate_precompute_G1");
    return result;
}

test_curve_ate_G2_precomp test_curve_ate_precompute_G2(const test_curve_G2& Q)
{
    enter_block("Call to test_curve_ate_precompute_G2");

    test_curve_G2 Qcopy(Q);
    Qcopy.to_affine_coordinates();

    test_curve_Fq two_inv = (test_curve_Fq("2").inverse()); // could add to global params if needed

    test_curve_ate_G2_precomp result;
    result.QX = Qcopy.X;
    result.QY = Qcopy.Y;

    test_curve_G2 R;
    R.X = Qcopy.X;
    R.Y = Qcopy.Y;
    R.Z = test_curve_Fq2::one();

    const bigint<test_curve_Fq::num_limbs> &loop_count = test_curve_ate_loop_count;
    bool found_one = false;
    test_curve_ate_ell_coeffs c;

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

    leave_block("Call to test_curve_ate_precompute_G2");
    return result;
}

test_curve_Fq12 test_curve_ate_miller_loop(const test_curve_ate_G1_precomp &prec_P,
                                     const test_curve_ate_G2_precomp &prec_Q)
{
    enter_block("Call to test_curve_ate_miller_loop");

    test_curve_Fq12 f = test_curve_Fq12::one();

    bool found_one = false;
    size_t idx = 0;

    const bigint<test_curve_Fq::num_limbs> &loop_count = test_curve_ate_loop_count;
    test_curve_ate_ell_coeffs c;

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
           test_curve_param_p (skipping leading zeros) in MSB to LSB
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

    if (test_curve_ate_is_loop_count_neg)
    {
    	f = f.inverse();
    }

    /*
    c = prec_Q.coeffs[idx++];
    f = f.mul_by_024(c.ell_0,prec_P.PY * c.ell_VW,prec_P.PX * c.ell_VV);

    c = prec_Q.coeffs[idx++];
    f = f.mul_by_024(c.ell_0,prec_P.PY * c.ell_VW,prec_P.PX * c.ell_VV);
    */

    leave_block("Call to test_curve_ate_miller_loop");
    return f;
}

test_curve_Fq12 test_curve_ate_double_miller_loop(const test_curve_ate_G1_precomp &prec_P1,
                                     const test_curve_ate_G2_precomp &prec_Q1,
                                     const test_curve_ate_G1_precomp &prec_P2,
                                     const test_curve_ate_G2_precomp &prec_Q2)
{
    enter_block("Call to test_curve_ate_double_miller_loop");

    test_curve_Fq12 f = test_curve_Fq12::one();

    bool found_one = false;
    size_t idx = 0;

    const bigint<test_curve_Fq::num_limbs> &loop_count = test_curve_ate_loop_count;
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
           test_curve_param_p (skipping leading zeros) in MSB to LSB
           order */

        test_curve_ate_ell_coeffs c1 = prec_Q1.coeffs[idx];
        test_curve_ate_ell_coeffs c2 = prec_Q2.coeffs[idx];
        ++idx;

        f = f.squared();

        f = f.mul_by_024(c1.ell_0, prec_P1.PY * c1.ell_VW, prec_P1.PX * c1.ell_VV);
        f = f.mul_by_024(c2.ell_0, prec_P2.PY * c2.ell_VW, prec_P2.PX * c2.ell_VV);

        if (bit)
        {
            test_curve_ate_ell_coeffs c1 = prec_Q1.coeffs[idx];
            test_curve_ate_ell_coeffs c2 = prec_Q2.coeffs[idx];
            ++idx;

            f = f.mul_by_024(c1.ell_0, prec_P1.PY * c1.ell_VW, prec_P1.PX * c1.ell_VV);
            f = f.mul_by_024(c2.ell_0, prec_P2.PY * c2.ell_VW, prec_P2.PX * c2.ell_VV);
        }
    }

    if (test_curve_ate_is_loop_count_neg)
    {
    	f = f.inverse();
    }

    /*
    test_curve_ate_ell_coeffs c1 = prec_Q1.coeffs[idx];
    test_curve_ate_ell_coeffs c2 = prec_Q2.coeffs[idx];
    ++idx;
    f = f.mul_by_024(c1.ell_0, prec_P1.PY * c1.ell_VW, prec_P1.PX * c1.ell_VV);
    f = f.mul_by_024(c2.ell_0, prec_P2.PY * c2.ell_VW, prec_P2.PX * c2.ell_VV);

    c1 = prec_Q1.coeffs[idx];
    c2 = prec_Q2.coeffs[idx];
    ++idx;
    f = f.mul_by_024(c1.ell_0, prec_P1.PY * c1.ell_VW, prec_P1.PX * c1.ell_VV);
    f = f.mul_by_024(c2.ell_0, prec_P2.PY * c2.ell_VW, prec_P2.PX * c2.ell_VV);
    */

    leave_block("Call to test_curve_ate_double_miller_loop");

    return f;
}

test_curve_Fq12 test_curve_ate_pairing(const test_curve_G1& P, const test_curve_G2 &Q)
{
    enter_block("Call to test_curve_ate_pairing");
    test_curve_ate_G1_precomp prec_P = test_curve_ate_precompute_G1(P);
    test_curve_ate_G2_precomp prec_Q = test_curve_ate_precompute_G2(Q);
    test_curve_Fq12 result = test_curve_ate_miller_loop(prec_P, prec_Q);
    leave_block("Call to test_curve_ate_pairing");
    return result;
}

test_curve_GT test_curve_ate_reduced_pairing(const test_curve_G1 &P, const test_curve_G2 &Q)
{
    enter_block("Call to test_curve_ate_reduced_pairing");
    const test_curve_Fq12 f = test_curve_ate_pairing(P, Q);
    const test_curve_GT result = test_curve_final_exponentiation(f);
    leave_block("Call to test_curve_ate_reduced_pairing");
    return result;
}

/* choice of pairing */

test_curve_G1_precomp test_curve_precompute_G1(const test_curve_G1& P)
{
    return test_curve_ate_precompute_G1(P);
}

test_curve_G2_precomp test_curve_precompute_G2(const test_curve_G2& Q)
{
    return test_curve_ate_precompute_G2(Q);
}

test_curve_Fq12 test_curve_miller_loop(const test_curve_G1_precomp &prec_P,
                          const test_curve_G2_precomp &prec_Q)
{
    return test_curve_ate_miller_loop(prec_P, prec_Q);
}

test_curve_Fq12 test_curve_double_miller_loop(const test_curve_G1_precomp &prec_P1,
                                 const test_curve_G2_precomp &prec_Q1,
                                 const test_curve_G1_precomp &prec_P2,
                                 const test_curve_G2_precomp &prec_Q2)
{
    return test_curve_ate_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

test_curve_Fq12 test_curve_pairing(const test_curve_G1& P,
                      const test_curve_G2 &Q)
{
    return test_curve_ate_pairing(P, Q);
}

test_curve_GT test_curve_reduced_pairing(const test_curve_G1 &P,
                             const test_curve_G2 &Q)
{
    return test_curve_ate_reduced_pairing(P, Q);
}
} // libff
