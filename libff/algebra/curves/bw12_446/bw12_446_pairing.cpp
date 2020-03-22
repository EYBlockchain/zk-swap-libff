#include <cassert>

#include <libff/algebra/curves/bw12_446/bw12_446_g1.hpp>
#include <libff/algebra/curves/bw12_446/bw12_446_g2.hpp>
#include <libff/algebra/curves/bw12_446/bw12_446_init.hpp>
#include <libff/algebra/curves/bw12_446/bw12_446_pairing.hpp>
#include <libff/common/profiling.hpp>

namespace libff {

  bool bw12_446_ate_G1_precomp::operator==(const bw12_446_ate_G1_precomp &other) const
  {
    return (this->PX == other.PX &&
        this->PY == other.PY);
  }

  std::ostream& operator<<(std::ostream &out, const bw12_446_ate_G1_precomp &prec_P)
  {
    out << prec_P.PX << OUTPUT_SEPARATOR << prec_P.PY;

    return out;
  }

  std::istream& operator>>(std::istream &in, bw12_446_ate_G1_precomp &prec_P)
  {
    in >> prec_P.PX;
    consume_OUTPUT_SEPARATOR(in);
    in >> prec_P.PY;

    return in;
  }

  bool  bw12_446_ate_ell_coeffs::operator==(const bw12_446_ate_ell_coeffs &other) const
  {
    return (this->ell_0 == other.ell_0 &&
        this->ell_VW == other.ell_VW &&
        this->ell_VV == other.ell_VV);
  }

  std::ostream& operator<<(std::ostream &out, const bw12_446_ate_ell_coeffs &c)
  {
    out << c.ell_0 << OUTPUT_SEPARATOR << c.ell_VW << OUTPUT_SEPARATOR << c.ell_VV;
    return out;
  }

  std::istream& operator>>(std::istream &in, bw12_446_ate_ell_coeffs &c)
  {
    in >> c.ell_0;
    consume_OUTPUT_SEPARATOR(in);
    in >> c.ell_VW;
    consume_OUTPUT_SEPARATOR(in);
    in >> c.ell_VV;

    return in;
  }

  bool bw12_446_ate_G2_precomp::operator==(const bw12_446_ate_G2_precomp &other) const
  {
    return (this->QX == other.QX &&
        this->QY == other.QY &&
        this->coeffs == other.coeffs);
  }

  std::ostream& operator<<(std::ostream& out, const bw12_446_ate_G2_precomp &prec_Q)
  {
    out << prec_Q.QX << OUTPUT_SEPARATOR << prec_Q.QY << "\n";
    out << prec_Q.coeffs.size() << "\n";
    for (const bw12_446_ate_ell_coeffs &c : prec_Q.coeffs)
    {
      out << c << OUTPUT_NEWLINE;
    }
    return out;
  }

  std::istream& operator>>(std::istream& in, bw12_446_ate_G2_precomp &prec_Q)
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
      bw12_446_ate_ell_coeffs c;
      in >> c;
      consume_OUTPUT_NEWLINE(in);
      prec_Q.coeffs.emplace_back(c);
    }

    return in;
  }

  /* final exponentiations */

  bw12_446_Fq12 bw12_446_final_exponentiation_first_chunk(const bw12_446_Fq12 &elt)
  {
    enter_block("Call to bw12_446_final_exponentiation_first_chunk");

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

    const bw12_446_Fq12 A = bw12_446_Fq12(elt.c0,-elt.c1);
    const bw12_446_Fq12 B = elt.inverse();
    const bw12_446_Fq12 C = A * B;
    const bw12_446_Fq12 D = C.Frobenius_map(2);
    const bw12_446_Fq12 result = D * C;

    leave_block("Call to bw12_446_final_exponentiation_first_chunk");

    return result;
  }

  bw12_446_Fq12 bw12_446_exp_by_neg_z(const bw12_446_Fq12 &elt)
  {
    enter_block("Call to bw12_446_exp_by_neg_z");

    bw12_446_Fq12 result = elt.cyclotomic_exp(bw12_446_final_exponent_z);
    if (!bw12_446_final_exponent_is_z_neg)
    {
      result = result.unitary_inverse();
    }

    leave_block("Call to bw12_446_exp_by_neg_z");

    return result;
  }

  bw12_446_Fq12 bw12_446_final_exponentiation_last_chunk(const bw12_446_Fq12 &elt)
  {
    enter_block("Call to bw12_446_final_exponentiation_last_chunk");

    /*
     * ++ powers of f
     * A = f^2
     * B = A * f (3)
     * C = B^2 (6)
     * D = C * f (7)
     *
     * ++ powers of f^z
     * E = exp_by_neg_z(f)
     * F = E^2 (2z)
     * G = F^2 (4z)
     * H = G * E (5z)
     * I = H * E (6z)
     * J = G^2 (8z)
     * K = J * H (13z)
     * L = K^2 (26z)
     * M = K * L (39z)
     *
     * ++ powers of f^(z^2)
     * N = exp_by_neg_z(I) (6z^2)
     * O = N^2 (12z^2)
     * P = O^2 (24z^2)
     * Q = P^2 (48z^2)
     * R = Q * N (54z^2)
     * S = R * O (66z^2)
     * T = S^2 (132z^2)
     * U = T * N (138z^2)
     *
     * ++ powers of f^(z^3)
     * V = exp_by_neg_z(N) (6z^3)
     * W = V^2 (12z^3)
     * X = W * V (18z^3)
     * Y = X^2 (36z^3)
     * Z = Y^2 (72z^3)
     * AA = Z^2 (144z^3)
     * BB = AA * X (162z^3)
     * CC = BB * Y * V (204z^3)
     * DD = CC * V (210z^3)
     * EE = AA^2 (288z^3)
     * FF = EE * X (306z^3)
     *
     * ++ powers of f^(z^4)
     * GG = exp_by_neg_z(Z) (72z^4)
     * HH = GG^2 (144z^4)
     * II = HH^2 (288z^4)
     * JJ = II * GG (360z^4)
     * KK = JJ * GG (432z^4)
     *
     * ++ powers of f^(z^5)
     * LL = exp_by_neg_z(II) (288z^5)
     * MM = LL^2 (576z^5)
     *
     * ++ list of positive powers of f
     * NN = R * JJ
     * OO = S * JJ
     * PP = D * U * JJ
     * QQ = H * CC * MM
     * RR = NN * OO.Frobenius_map(1) * PP.Frobeinus_map(2) * QQ.Frobenius_map(3)
     *
     * ++ list of negative powers of f
     * SS = I * BB * LL
     * TT = K * DD * LL
     * UU = M * FF * LL
     * VV = f * R * KK
     * WW = SS * TT.Frobenius_map(1) * UU.Frobeinus_map(2) * VV.Frobenius_map(3)
     * XX = WW.inverse()
     * YY = XX * RR
     */
    const bw12_446_Fq12 A = elt.squared();
    const bw12_446_Fq12 B = A * elt;
    const bw12_446_Fq12 C = B.squared();
    const bw12_446_Fq12 D = C * elt;
    const bw12_446_Fq12 E = bw12_446_exp_by_neg_z(elt);
    const bw12_446_Fq12 F = E.squared();
    const bw12_446_Fq12 G = F.squared();
    const bw12_446_Fq12 H = G * E;
    const bw12_446_Fq12 I = H * E;
    const bw12_446_Fq12 J = G.squared();
    const bw12_446_Fq12 K = J * H;
    const bw12_446_Fq12 L = K.squared();
    const bw12_446_Fq12 M = K * L;
    const bw12_446_Fq12 N = bw12_446_exp_by_neg_z(I);
    const bw12_446_Fq12 O = N.squared();
    const bw12_446_Fq12 P = O.squared();
    const bw12_446_Fq12 Q = P.squared();
    const bw12_446_Fq12 R = Q * N;
    const bw12_446_Fq12 S = R * O;
    const bw12_446_Fq12 T = S.squared();
    const bw12_446_Fq12 U = T * N;
    const bw12_446_Fq12 V = bw12_446_exp_by_neg_z(N);
    const bw12_446_Fq12 W = V.squared();
    const bw12_446_Fq12 X = W * V;
    const bw12_446_Fq12 Y = X.squared();
    const bw12_446_Fq12 Z = Y.squared();
    const bw12_446_Fq12 AA = Z.squared();
    const bw12_446_Fq12 BB = AA * X;
    const bw12_446_Fq12 CC = BB * Y * V;
    const bw12_446_Fq12 DD = CC * V;
    const bw12_446_Fq12 EE = AA.squared();
    const bw12_446_Fq12 FF = EE * X;
    const bw12_446_Fq12 GG = bw12_446_exp_by_neg_z(Z);
    const bw12_446_Fq12 HH = GG.squared();
    const bw12_446_Fq12 II = HH.squared();
    const bw12_446_Fq12 JJ = II * GG;
    const bw12_446_Fq12 KK = JJ * GG;
    const bw12_446_Fq12 LL = bw12_446_exp_by_neg_z(II);
    const bw12_446_Fq12 MM = LL.squared();
    const bw12_446_Fq12 NN = R * JJ;
    const bw12_446_Fq12 OO = S * JJ;
    const bw12_446_Fq12 PP = D * U * JJ;
    const bw12_446_Fq12 QQ = H * CC * MM;
    const bw12_446_Fq12 RR = NN * OO.Frobenius_map(1) * PP.Frobenius_map(2) * QQ.Frobenius_map(3);
    const bw12_446_Fq12 SS = I * BB * LL;
    const bw12_446_Fq12 TT = K * DD * LL;
    const bw12_446_Fq12 UU = M * FF * LL;
    const bw12_446_Fq12 VV = elt * R * KK;
    const bw12_446_Fq12 WW = SS * TT.Frobenius_map(1) * UU.Frobenius_map(2) * VV.Frobenius_map(3);
    const bw12_446_Fq12 XX = WW.unitary_inverse();
    const bw12_446_Fq12 YY = XX * RR;

    const bw12_446_Fq12 result = YY;

    leave_block("Call to bw12_446_final_exponentiation_last_chunk");

    return result;
}

bw12_446_GT bw12_446_final_exponentiation(const bw12_446_Fq12 &elt)
{
    enter_block("Call to bw12_446_final_exponentiation");
    /* OLD naive version:
        bw12_446_GT result = elt^bw12_446_final_exponent;
    */
    bw12_446_Fq12 A = bw12_446_final_exponentiation_first_chunk(elt);
    bw12_446_GT result = bw12_446_final_exponentiation_last_chunk(A);

    leave_block("Call to bw12_446_final_exponentiation");
    return result;
}

/* ate pairing */

void doubling_step_for_flipped_miller_loop(const bw12_446_Fq two_inv,
                                           bw12_446_G2 &current,
                                           bw12_446_ate_ell_coeffs &c)
{
    const bw12_446_Fq2 X = current.X, Y = current.Y, Z = current.Z;

    const bw12_446_Fq2 A = two_inv * (X * Y);                     // A = X1 * Y1 / 2
    const bw12_446_Fq2 B = Y.squared();                           // B = Y1^2
    const bw12_446_Fq2 C = Z.squared();                           // C = Z1^2
    const bw12_446_Fq2 D = C+C+C;                                 // D = 3 * C
    const bw12_446_Fq2 E = bw12_446_twist_coeff_b * D;             // E = twist_b * D
    const bw12_446_Fq2 F = E+E+E;                                 // F = 3 * E
    const bw12_446_Fq2 G = two_inv * (B+F);                       // G = (B+F)/2
    const bw12_446_Fq2 H = (Y+Z).squared() - (B+C);               // H = (Y1+Z1)^2-(B+C)
    const bw12_446_Fq2 I = E-B;                                   // I = E-B
    const bw12_446_Fq2 J = X.squared();                           // J = X1^2
    const bw12_446_Fq2 E_squared = E.squared();                   // E_squared = E^2

    current.X = A * (B-F);                                       // X3 = A * (B-F)
    current.Y = G.squared() - (E_squared+E_squared+E_squared);   // Y3 = G^2 - 3*E^2
    current.Z = B * H;                                           // Z3 = B * H
    c.ell_0 = bw12_446_twist * I;                                 // ell_0 = xi * I
    c.ell_VW = -H;                                               // ell_VW = - H (later: * yP)
    c.ell_VV = J+J+J;                                            // ell_VV = 3*J (later: * xP)
}

void mixed_addition_step_for_flipped_miller_loop(const bw12_446_G2 base,
                                                 bw12_446_G2 &current,
                                                 bw12_446_ate_ell_coeffs &c)
{
    const bw12_446_Fq2 X1 = current.X, Y1 = current.Y, Z1 = current.Z;
    const bw12_446_Fq2 &x2 = base.X, &y2 = base.Y;

    const bw12_446_Fq2 D = X1 - x2 * Z1;          // D = X1 - X2*Z1
    const bw12_446_Fq2 E = Y1 - y2 * Z1;          // E = Y1 - Y2*Z1
    const bw12_446_Fq2 F = D.squared();           // F = D^2
    const bw12_446_Fq2 G = E.squared();           // G = E^2
    const bw12_446_Fq2 H = D*F;                   // H = D*F
    const bw12_446_Fq2 I = X1 * F;                // I = X1 * F
    const bw12_446_Fq2 J = H + Z1*G - (I+I);      // J = H + Z1*G - (I+I)

    current.X = D * J;                           // X3 = D*J
    current.Y = E * (I-J)-(H * Y1);              // Y3 = E*(I-J)-(H*Y1)
    current.Z = Z1 * H;                          // Z3 = Z1*H
    c.ell_0 = bw12_446_twist * (E * x2 - D * y2); // ell_0 = xi * (E * X2 - D * Y2)
    c.ell_VV = - E;                              // ell_VV = - E (later: * xP)
    c.ell_VW = D;                                // ell_VW = D (later: * yP    )
}

bw12_446_ate_G1_precomp bw12_446_ate_precompute_G1(const bw12_446_G1& P)
{
    enter_block("Call to bw12_446_ate_precompute_G1");

    bw12_446_G1 Pcopy = P;
    Pcopy.to_affine_coordinates();

    bw12_446_ate_G1_precomp result;
    result.PX = Pcopy.X;
    result.PY = Pcopy.Y;

    leave_block("Call to bw12_446_ate_precompute_G1");
    return result;
}

bw12_446_ate_G2_precomp bw12_446_ate_precompute_G2(const bw12_446_G2& Q)
{
    enter_block("Call to bw12_446_ate_precompute_G2");

    bw12_446_G2 Qcopy(Q);
    Qcopy.to_affine_coordinates();

    bw12_446_Fq two_inv = (bw12_446_Fq("2").inverse()); // could add to global params if needed

    bw12_446_ate_G2_precomp result;
    result.QX = Qcopy.X;
    result.QY = Qcopy.Y;

    bw12_446_G2 R;
    R.X = Qcopy.X;
    R.Y = Qcopy.Y;
    R.Z = bw12_446_Fq2::one();

    const bigint<bw12_446_Fq::num_limbs> &loop_count = bw12_446_ate_loop_count;
    bool found_one = false;
    bw12_446_ate_ell_coeffs c;

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

    bw12_446_G2 Q1 = Qcopy.mul_by_q();
    Q1.Y = - Q1.Y;
    assert(Q1.Z == bw12_446_Fq2::one());
    bw12_446_G2 Q2 = Q1.mul_by_q();
    assert(Q2.Z == bw12_446_Fq2::one());

    if (bw12_446_ate_is_loop_count_neg)
    {
        R.Y = - R.Y;
    }

    mixed_addition_step_for_flipped_miller_loop(Q1, R, c);
    result.coeffs.push_back(c);

    mixed_addition_step_for_flipped_miller_loop(Q2, R, c);
    result.coeffs.push_back(c);

    leave_block("Call to bw12_446_ate_precompute_G2");
    return result;
}

bw12_446_Fq12 bw12_446_ate_miller_loop(const bw12_446_ate_G1_precomp &prec_P,
                                     const bw12_446_ate_G2_precomp &prec_Q)
{
    enter_block("Call to bw12_446_ate_miller_loop");

    bw12_446_Fq12 f = bw12_446_Fq12::one();

    bool found_one = false;
    size_t idx = 0;

    const bigint<bw12_446_Fq::num_limbs> &loop_count = bw12_446_ate_loop_count;
    bw12_446_ate_ell_coeffs c;

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
           bw12_446_param_p (skipping leading zeros) in MSB to LSB
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

    if (bw12_446_ate_is_loop_count_neg)
    {
    	f = f.inverse();
    }

    c = prec_Q.coeffs[idx++];
    f = f.mul_by_024(c.ell_0,prec_P.PY * c.ell_VW,prec_P.PX * c.ell_VV);

    c = prec_Q.coeffs[idx++];
    f = f.mul_by_024(c.ell_0,prec_P.PY * c.ell_VW,prec_P.PX * c.ell_VV);

    leave_block("Call to bw12_446_ate_miller_loop");
    return f;
}

bw12_446_Fq12 bw12_446_ate_double_miller_loop(const bw12_446_ate_G1_precomp &prec_P1,
                                     const bw12_446_ate_G2_precomp &prec_Q1,
                                     const bw12_446_ate_G1_precomp &prec_P2,
                                     const bw12_446_ate_G2_precomp &prec_Q2)
{
    enter_block("Call to bw12_446_ate_double_miller_loop");

    bw12_446_Fq12 f = bw12_446_Fq12::one();

    bool found_one = false;
    size_t idx = 0;

    const bigint<bw12_446_Fq::num_limbs> &loop_count = bw12_446_ate_loop_count;
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
           bw12_446_param_p (skipping leading zeros) in MSB to LSB
           order */

        bw12_446_ate_ell_coeffs c1 = prec_Q1.coeffs[idx];
        bw12_446_ate_ell_coeffs c2 = prec_Q2.coeffs[idx];
        ++idx;

        f = f.squared();

        f = f.mul_by_024(c1.ell_0, prec_P1.PY * c1.ell_VW, prec_P1.PX * c1.ell_VV);
        f = f.mul_by_024(c2.ell_0, prec_P2.PY * c2.ell_VW, prec_P2.PX * c2.ell_VV);

        if (bit)
        {
            bw12_446_ate_ell_coeffs c1 = prec_Q1.coeffs[idx];
            bw12_446_ate_ell_coeffs c2 = prec_Q2.coeffs[idx];
            ++idx;

            f = f.mul_by_024(c1.ell_0, prec_P1.PY * c1.ell_VW, prec_P1.PX * c1.ell_VV);
            f = f.mul_by_024(c2.ell_0, prec_P2.PY * c2.ell_VW, prec_P2.PX * c2.ell_VV);
        }
    }

    if (bw12_446_ate_is_loop_count_neg)
    {
    	f = f.inverse();
    }

    bw12_446_ate_ell_coeffs c1 = prec_Q1.coeffs[idx];
    bw12_446_ate_ell_coeffs c2 = prec_Q2.coeffs[idx];
    ++idx;
    f = f.mul_by_024(c1.ell_0, prec_P1.PY * c1.ell_VW, prec_P1.PX * c1.ell_VV);
    f = f.mul_by_024(c2.ell_0, prec_P2.PY * c2.ell_VW, prec_P2.PX * c2.ell_VV);

    c1 = prec_Q1.coeffs[idx];
    c2 = prec_Q2.coeffs[idx];
    ++idx;
    f = f.mul_by_024(c1.ell_0, prec_P1.PY * c1.ell_VW, prec_P1.PX * c1.ell_VV);
    f = f.mul_by_024(c2.ell_0, prec_P2.PY * c2.ell_VW, prec_P2.PX * c2.ell_VV);

    leave_block("Call to bw12_446_ate_double_miller_loop");

    return f;
}

bw12_446_Fq12 bw12_446_ate_pairing(const bw12_446_G1& P, const bw12_446_G2 &Q)
{
    enter_block("Call to bw12_446_ate_pairing");
    bw12_446_ate_G1_precomp prec_P = bw12_446_ate_precompute_G1(P);
    bw12_446_ate_G2_precomp prec_Q = bw12_446_ate_precompute_G2(Q);
    bw12_446_Fq12 result = bw12_446_ate_miller_loop(prec_P, prec_Q);
    leave_block("Call to bw12_446_ate_pairing");
    return result;
}

bw12_446_GT bw12_446_ate_reduced_pairing(const bw12_446_G1 &P, const bw12_446_G2 &Q)
{
    enter_block("Call to bw12_446_ate_reduced_pairing");
    const bw12_446_Fq12 f = bw12_446_ate_pairing(P, Q);
    const bw12_446_GT result = bw12_446_final_exponentiation(f);
    leave_block("Call to bw12_446_ate_reduced_pairing");
    return result;
}

/* choice of pairing */

bw12_446_G1_precomp bw12_446_precompute_G1(const bw12_446_G1& P)
{
    return bw12_446_ate_precompute_G1(P);
}

bw12_446_G2_precomp bw12_446_precompute_G2(const bw12_446_G2& Q)
{
    return bw12_446_ate_precompute_G2(Q);
}

bw12_446_Fq12 bw12_446_miller_loop(const bw12_446_G1_precomp &prec_P,
                          const bw12_446_G2_precomp &prec_Q)
{
    return bw12_446_ate_miller_loop(prec_P, prec_Q);
}

bw12_446_Fq12 bw12_446_double_miller_loop(const bw12_446_G1_precomp &prec_P1,
                                 const bw12_446_G2_precomp &prec_Q1,
                                 const bw12_446_G1_precomp &prec_P2,
                                 const bw12_446_G2_precomp &prec_Q2)
{
    return bw12_446_ate_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

bw12_446_Fq12 bw12_446_pairing(const bw12_446_G1& P,
                      const bw12_446_G2 &Q)
{
    return bw12_446_ate_pairing(P, Q);
}

bw12_446_GT bw12_446_reduced_pairing(const bw12_446_G1 &P,
                             const bw12_446_G2 &Q)
{
    return bw12_446_ate_reduced_pairing(P, Q);
}
} // libff
