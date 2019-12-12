#include <cassert>

#include <libff/algebra/curves/sw6_bis/sw6_bis_g1.hpp>
#include <libff/algebra/curves/sw6_bis/sw6_bis_g2.hpp>
#include <libff/algebra/curves/sw6_bis/sw6_bis_init.hpp>
#include <libff/algebra/curves/sw6_bis/sw6_bis_pairing.hpp>
#include <libff/algebra/scalar_multiplication/wnaf.hpp>
#include <libff/common/profiling.hpp>

namespace libff {

bool sw6_bis_ate_G1_precomp::operator==(const sw6_bis_ate_G1_precomp &other) const
{
    return (this->PX == other.PX &&
            this->PY == other.PY &&
            this->PX_twist == other.PX_twist &&
            this->PY_twist == other.PY_twist);
}

std::ostream& operator<<(std::ostream &out, const sw6_bis_ate_G1_precomp &prec_P)
{
    out << prec_P.PX << OUTPUT_SEPARATOR << prec_P.PY << OUTPUT_SEPARATOR << prec_P.PX_twist << OUTPUT_SEPARATOR << prec_P.PY_twist;

    return out;
}

std::istream& operator>>(std::istream &in, sw6_bis_ate_G1_precomp &prec_P)
{
    in >> prec_P.PX;
    consume_OUTPUT_SEPARATOR(in);
    in >> prec_P.PY;
    consume_OUTPUT_SEPARATOR(in);
    in >> prec_P.PX_twist;
    consume_OUTPUT_SEPARATOR(in);
    in >> prec_P.PY_twist;

    return in;
}

bool sw6_bis_ate_dbl_coeffs::operator==(const sw6_bis_ate_dbl_coeffs &other) const
{
    return (this->c_H == other.c_H &&
            this->c_4C == other.c_4C &&
            this->c_J == other.c_J &&
            this->c_L == other.c_L);
}

std::ostream& operator<<(std::ostream &out, const sw6_bis_ate_dbl_coeffs &dc)
{
    out << dc.c_H << OUTPUT_SEPARATOR << dc.c_4C << OUTPUT_SEPARATOR << dc.c_J << OUTPUT_SEPARATOR << dc.c_L;
    return out;
}

std::istream& operator>>(std::istream &in, sw6_bis_ate_dbl_coeffs &dc)
{
    in >> dc.c_H;
    consume_OUTPUT_SEPARATOR(in);
    in >> dc.c_4C;
    consume_OUTPUT_SEPARATOR(in);
    in >> dc.c_J;
    consume_OUTPUT_SEPARATOR(in);
    in >> dc.c_L;

    return in;
}

bool sw6_bis_ate_add_coeffs::operator==(const sw6_bis_ate_add_coeffs &other) const
{
    return (this->c_L1 == other.c_L1 &&
            this->c_RZ == other.c_RZ);
}

std::ostream& operator<<(std::ostream &out, const sw6_bis_ate_add_coeffs &ac)
{
    out << ac.c_L1 << OUTPUT_SEPARATOR << ac.c_RZ;
    return out;
}

std::istream& operator>>(std::istream &in, sw6_bis_ate_add_coeffs &ac)
{
    in >> ac.c_L1;
    consume_OUTPUT_SEPARATOR(in);
    in >> ac.c_RZ;

    return in;
}


bool sw6_bis_ate_G2_precomp::operator==(const sw6_bis_ate_G2_precomp &other) const
{
    return (this->QX == other.QX &&
            this->QY == other.QY &&
            this->QY2 == other.QY2 &&
            this->QX_over_twist == other.QX_over_twist &&
            this->QY_over_twist == other.QY_over_twist &&
            this->dbl_coeffs == other.dbl_coeffs &&
            this->add_coeffs == other.add_coeffs);
}

std::ostream& operator<<(std::ostream& out, const sw6_bis_ate_G2_precomp &prec_Q)
{
    out << prec_Q.QX << OUTPUT_SEPARATOR
        << prec_Q.QY << OUTPUT_SEPARATOR
        << prec_Q.QY2  << OUTPUT_SEPARATOR
        << prec_Q.QX_over_twist << OUTPUT_SEPARATOR
        << prec_Q.QY_over_twist << "\n";
    out << prec_Q.dbl_coeffs.size() << "\n";
    for (const sw6_bis_ate_dbl_coeffs &dc : prec_Q.dbl_coeffs)
    {
        out << dc << OUTPUT_NEWLINE;
    }
    out << prec_Q.add_coeffs.size() << "\n";
    for (const sw6_bis_ate_add_coeffs &ac : prec_Q.add_coeffs)
    {
        out << ac << OUTPUT_NEWLINE;
    }

    return out;
}

std::istream& operator>>(std::istream& in, sw6_bis_ate_G2_precomp &prec_Q)
{
    in >> prec_Q.QX;
    consume_OUTPUT_SEPARATOR(in);
    in >> prec_Q.QY;
    consume_OUTPUT_SEPARATOR(in);
    in >> prec_Q.QY2;
    consume_OUTPUT_SEPARATOR(in);
    in >> prec_Q.QX_over_twist;
    consume_OUTPUT_SEPARATOR(in);
    in >> prec_Q.QY_over_twist;
    consume_newline(in);

    prec_Q.dbl_coeffs.clear();
    size_t dbl_s;
    in >> dbl_s;
    consume_newline(in);

    prec_Q.dbl_coeffs.reserve(dbl_s);

    for (size_t i = 0; i < dbl_s; ++i)
    {
        sw6_bis_ate_dbl_coeffs dc;
        in >> dc;
        consume_OUTPUT_NEWLINE(in);
        prec_Q.dbl_coeffs.emplace_back(dc);
    }

    prec_Q.add_coeffs.clear();
    size_t add_s;
    in >> add_s;
    consume_newline(in);

    prec_Q.add_coeffs.reserve(add_s);

    for (size_t i = 0; i < add_s; ++i)
    {
        sw6_bis_ate_add_coeffs ac;
        in >> ac;
        consume_OUTPUT_NEWLINE(in);
        prec_Q.add_coeffs.emplace_back(ac);
    }

    return in;
}

sw6_bis_Fq6 sw6_bis_exp_by_z(const sw6_bis_Fq6 &elt)
{
  enter_block("Call to sw6_bis_exp_by_z");

  sw6_bis_Fq6 result = elt.cyclotomic_exp(sw6_bis_final_exponent_z);
  if (sw6_bis_final_exponent_is_z_neg)
  {
    result = result.unitary_inverse();
  }

  leave_block("Call to sw6_bis_exp_by_z");

  return result;
}

/* final exponentiations */

sw6_bis_Fq6 old_sw6_bis_final_exponentiation_last_chunk(const sw6_bis_Fq6 &elt, const sw6_bis_Fq6 &elt_inv)
{
    enter_block("Call to old_sw6_bis_final_exponentiation_last_chunk");
    const sw6_bis_Fq6 elt_q = elt.Frobenius_map(1);
    sw6_bis_Fq6 w1_part = elt_q.cyclotomic_exp(sw6_bis_final_exponent_last_chunk_w1);
    sw6_bis_Fq6 w0_part;
    if (sw6_bis_final_exponent_last_chunk_is_w0_neg)
    {
    	w0_part = elt_inv.cyclotomic_exp(sw6_bis_final_exponent_last_chunk_abs_of_w0);
    } else {
    	w0_part = elt.cyclotomic_exp(sw6_bis_final_exponent_last_chunk_abs_of_w0);
    }
    sw6_bis_Fq6 result = w1_part * w0_part;
    leave_block("Call to old_sw6_bis_final_exponentiation_last_chunk");

    return result;
}

sw6_bis_Fq6 sw6_bis_final_exponentiation_last_chunk(const sw6_bis_Fq6 &elt)
{
    enter_block("Call to sw6_bis_final_exponentiation_last_chunk");
    /*
     * R0(x) := (-103*x^7 + 70*x^6 + 269*x^5 - 197*x^4 - 314*x^3 - 73*x^2 - 263*x - 220)
     * R1(x) := (103*x^9 - 276*x^8 + 77*x^7 + 492*x^6 - 445*x^5 - 65*x^4 + 452*x^3 - 181*x^2 + 34*x + 229)
     *
     * La dernière partie de l'exponentiation finale est
     * f^R0(u)*(f^p)^R1(u) avec R0 et R1 ci-dessus et f^p est un Frobenius dans GF(p^6).
     *
     * ++ powers of f
     * A0 = f^2 (2)
     * B0 = A0^2 (4)
     * C0 = B0^2 (8)
     * D0 = C0^2 (16)
     * E0 = D0^2 (32)
     * F0 = E0^2 (64)
     * G0 = F0^2 (128)
     * H0 = G0*F0 (192)
     * I0 = H0*D0 (208)
     * J0 = I0*C0 (216)
     * (R0) K0 = J0*B0 (220)
     * L0 = K0*C0 (228)
     * (R1) M0 = L0*f (229)
     *
     * ++ powers of f^z
     * A1 = exp_by_z(A0) (2z)
     * B1 = exp_by_z(E0) (32z)
     * (R1) C1 = A1*B1        (34z)
     * D1 = exp_by_z(M0) (229z)
     * (R0) E1 = D1*C1        (263z)
     * F1 = A1^2         (4z)
     * G1 = exp_by_z(f)  (z)
     * H1 = F1*G1        (5z)
     *
     * ++ powers of f^(z^2)
     * A2 = exp_by_z(C1) (34z^2)
     * B2 = C1^2         (68z^2)
     * C2 = exp_by_z(H1) (5z^2)
     * (R0) D2 = C2*B2        (73z^2)
     * E2 = D2^2         (146z^2)
     * F2 = E2*A2        (180z^2)
     * G2 = exp_by_z(G1) (z^2)
     * H2 = G2*F2        (181z^2)
     * (R1) H2_inv = H2.inverse()   (-181z^2)
     * I2 = H2*D2        (254z^2)
     *
     * ++ powers of f^(z^3)
     * A3 = exp_by_z(I2) (254z^3)
     * B3 = exp_by_z(A2) (34z^3)
     * C3 = exp_by_z(C2) (5z^3)
     * D3 = C3^2         (10z^3)
     * E3 = D3*C3        (15z^3)
     * F3 = exp_by_z(G2) (z^3)
     * G3 = F3*E3        (16z^3)
     * H3 = G3*A3        (270z^3)
     * error H30 = H3*B3       (304z^3)
     * error (R0) H31 = H30*D3      (314z^3)
     * I3 = B3^2         (68z^3)
     * J3 = I3^2         (136z^3)
     * K3 = F3^2         (2z^3)
     * L3 = K3*J3        (138z^3)
     * error (R1) M3 = H31*L3        (452z^3)
     * N3 = E3^2         (30z^3)
     * O3 = N3*B3        (64z^3)
     *
     * ++ powers of f^(z^4)
     * A4 = exp_by_z(F3) (z^4)
     * B4 = A4^2         (2z^4)
     * C4 = exp_by_z(O3) (64z^4)
     * D4 = C4*A4        (65z^4)
     * (R1) D4_inv = D4.inverse()        (-65z^4)
     * E4 = D4^2         (130z^4)
     * F4 = E4*D4        (195z^4)
     * (R0) G4 = F4*B4        (197z^4)
     * error H4 = E4^2         (260z^4)
     * I4 = B4^2         (4z^4)
     * J4 = exp_by_z(C3) (5z^4)
     * K4 = J4*I4        (9z^4)
     * L4 = K4*H4        (269z^4)
     *
     * ++ powers of f^(z^5)
     * A5 = exp_by_z(L4) (269z^5)
     * (R0) A5_inv = A5.inverse() (-269z^5)
     * B5 = exp_by_z(F4) (195z^5)
     * C5 = B5^2         (390z^5)
     * E5 = exp_by_z(K4) (9z^5)
     * F5 = E5^2         (18z^5)
     * G5 = F5^2         (36z^5)
     * H5 = G5*F5        (54z^5)
     * I5 = exp_by_z(A4) (z^5)
     * J5 = I5*H5        (55z^5)
     * K5 = J5*C5        (445z^5)
     * (R1) K5_inv = K5.inverse()        (-445z^5)
     * error L5 = J5*E5        (64z^5)
     * M5 = exp_by_z(J4) (5z^5)
     * error N5 = M5*L5        (69z^5)
     *
     * ++ powers of f^(z^6)
     * error A6 = exp_by_z(N5) (69z^6)
     * B6 = exp_by_z(I5) (z^6)
     * C6 = B5^2         (2z^6)
     * error D6 = B6*A6        (70z^6)
     * (R0) D6_inv = D6.inverse() (70z^6)
     * E6 = exp_by_z(K5) (445z^6)
     * F6 = exp_by_z(F5) (18z^6)
     * G6 = F6^2         (36z^6)
     * H6 = exp_by_z(E5) (9z^6)
     * I6 = H6*G6        (45z^6)
     * J6 = I6*C6        (47z^6)
     * (R1) K6 = J6*E6        (492z^6)
     * L6 = J6*F6        (65z^6)
     * M6 = C6^2         (4z^6)
     * N6 = M6*B6        (5z^6)
     *
     * ++ powers of f^(z^7)
     * A7 = exp_by_z(J6) (47z^7)
     * B7 = A7^2         (94z^7)
     * C7 = exp_by_z(H6) (9z^7)
     * (R0) D7 = C7*B7        (103z^7)
     * E7 = exp_by_z(L6) (65z^7)
     * F7 = exp_by_z(N6) (5z^7)
     * G7 = E7*F7        (70z^7)
     * H7 = exp_by_z(C6) (2z^7)
     * I7 = H7*F7        (7z^7)
     * (R1) J7 = I7*G7        (77z^7)
     *
     * ++ powers of f^(z^8)
     * A8 = exp_by_z(D7) (103z^8)
     * B8 = A8^2         (206z^8)
     * C8 = exp_by_z(G7) (70z^8)
     * D8 = B8*C8        (276z^8)
     * (R1) D8_inv = D8.inverse()        (-276z^8)
     *
     * ++ powers of f^(z^9)
     * (R1) A9 = exp_by_z(A8) (103z^9)
     *
     * ++ f^R0(u)
     * R01 = K0*E1
     * R02 = D2*R01
     * error R03 = H31*R02
     * R04 = G4*R03
     * R05 = A5_inv*R04
     * R06 = D6_inv*R05
     * R07 = D7*R06
     * R0 = R07.inverse()  (f^R0(u))
     *
     * ++ (f^p)^R1(u)
     * R11 = M0*C1
     * R12 = H2_inv*R11
     * R13 = M3*R12
     * R14 = D4_inv*R13
     * R15 = K5_inv*R14
     * R16 = K6*R15
     * R17 = J7*R16
     * R18 = D8_inv*R17
     * R19 = A9*R18
     * R1 = R19.Frobenius_map(1)
     *
     * result = R0*R1
     */

    const sw6_bis_Fq6 A0 = elt.squared();
    const sw6_bis_Fq6 B0 = A0.squared();
    const sw6_bis_Fq6 C0 = B0.squared();
    const sw6_bis_Fq6 D0 = C0.squared();
    const sw6_bis_Fq6 E0 = D0.squared();
    const sw6_bis_Fq6 F0 = E0.squared();
    const sw6_bis_Fq6 G0 = F0.squared();
    const sw6_bis_Fq6 H0 = G0*F0;
    const sw6_bis_Fq6 I0 = H0*D0;
    const sw6_bis_Fq6 J0 = I0*C0;
    const sw6_bis_Fq6 K0 = J0*B0;
    const sw6_bis_Fq6 L0 = K0*C0;
    const sw6_bis_Fq6 M0 = L0*elt;

    const sw6_bis_Fq6 A1 = sw6_bis_exp_by_z(A0);
    const sw6_bis_Fq6 B1 = sw6_bis_exp_by_z(E0);
    const sw6_bis_Fq6 C1 = A1*B1;
    const sw6_bis_Fq6 D1 = sw6_bis_exp_by_z(M0);
    const sw6_bis_Fq6 E1 = D1*C1;
    const sw6_bis_Fq6 F1 = A1.squared();
    const sw6_bis_Fq6 G1 = sw6_bis_exp_by_z(elt);
    const sw6_bis_Fq6 H1 = F1*G1;

    const sw6_bis_Fq6 A2 = sw6_bis_exp_by_z(C1);
    const sw6_bis_Fq6 B2 = C1.squared();
    const sw6_bis_Fq6 C2 = sw6_bis_exp_by_z(H1);
    const sw6_bis_Fq6 D2 = C2*B2;
    const sw6_bis_Fq6 E2 = D2.squared();
    const sw6_bis_Fq6 F2 = E2*A2;
    const sw6_bis_Fq6 G2 = sw6_bis_exp_by_z(G1);
    const sw6_bis_Fq6 H2 = G2*F2;
    const sw6_bis_Fq6 H2_inv = H2.unitary_inverse();
    const sw6_bis_Fq6 I2 = H2*D2;

    const sw6_bis_Fq6 A3 = sw6_bis_exp_by_z(I2);
    const sw6_bis_Fq6 B3 = sw6_bis_exp_by_z(A2);
    const sw6_bis_Fq6 C3 = sw6_bis_exp_by_z(C2);
    const sw6_bis_Fq6 D3 = C3.squared();
    const sw6_bis_Fq6 E3 = D3*C3;
    const sw6_bis_Fq6 F3 = sw6_bis_exp_by_z(G2);
    const sw6_bis_Fq6 G3 = F3*E3;
    const sw6_bis_Fq6 H3 = G3*A3;
    const sw6_bis_Fq6 H30 = H3*B3;
    const sw6_bis_Fq6 H31 = H30*D3;
    const sw6_bis_Fq6 I3 = B3.squared();
    const sw6_bis_Fq6 J3 = I3.squared();
    const sw6_bis_Fq6 K3 = F3.squared();
    const sw6_bis_Fq6 L3 = K3*J3;
    const sw6_bis_Fq6 M3 = H31*L3;
    const sw6_bis_Fq6 N3 = E3.squared();
    const sw6_bis_Fq6 O3 = N3*B3;

    const sw6_bis_Fq6 A4 = sw6_bis_exp_by_z(F3);
    const sw6_bis_Fq6 B4 = A4.squared();
    const sw6_bis_Fq6 C4 = sw6_bis_exp_by_z(O3);
    const sw6_bis_Fq6 D4 = C4*A4;
    const sw6_bis_Fq6 D4_inv = D4.unitary_inverse();
    const sw6_bis_Fq6 E4 = D4.squared();
    const sw6_bis_Fq6 F4 = E4*D4;
    const sw6_bis_Fq6 G4 = F4*B4;
    const sw6_bis_Fq6 H4 = E4.squared();
    const sw6_bis_Fq6 I4 = B4.squared();
    const sw6_bis_Fq6 J4 = sw6_bis_exp_by_z(C3);
    const sw6_bis_Fq6 K4 = J4*I4;
    const sw6_bis_Fq6 L4 = K4*H4;

    const sw6_bis_Fq6 A5 = sw6_bis_exp_by_z(L4);
    const sw6_bis_Fq6 A5_inv = A5.unitary_inverse();
    const sw6_bis_Fq6 B5 = sw6_bis_exp_by_z(F4);
    const sw6_bis_Fq6 C5 = B5.squared();
    const sw6_bis_Fq6 E5 = sw6_bis_exp_by_z(K4);
    const sw6_bis_Fq6 F5 = E5.squared();
    const sw6_bis_Fq6 G5 = F5.squared();
    const sw6_bis_Fq6 H5 = G5*F5;
    const sw6_bis_Fq6 I5 = sw6_bis_exp_by_z(A4);
    const sw6_bis_Fq6 J5 = I5*H5;
    const sw6_bis_Fq6 K5 = J5*C5;
    const sw6_bis_Fq6 K5_inv = K5.unitary_inverse();
    const sw6_bis_Fq6 L5 = J5*E5;
    const sw6_bis_Fq6 M5 = sw6_bis_exp_by_z(J4);
    const sw6_bis_Fq6 N5 = M5*L5;

    const sw6_bis_Fq6 A6 = sw6_bis_exp_by_z(N5);
    const sw6_bis_Fq6 B6 = sw6_bis_exp_by_z(I5);
    const sw6_bis_Fq6 C6 = B5.squared();
    const sw6_bis_Fq6 D6 = B6*A6;
    const sw6_bis_Fq6 D6_inv = D6.unitary_inverse();
    const sw6_bis_Fq6 E6 = sw6_bis_exp_by_z(K5);
    const sw6_bis_Fq6 F6 = sw6_bis_exp_by_z(F5);
    const sw6_bis_Fq6 G6 = F6.squared();
    const sw6_bis_Fq6 H6 = sw6_bis_exp_by_z(E5);
    const sw6_bis_Fq6 I6 = H6*G6;
    const sw6_bis_Fq6 J6 = I6*C6;
    const sw6_bis_Fq6 K6 = J6*E6;
    const sw6_bis_Fq6 L6 = J6*F6;
    const sw6_bis_Fq6 M6 = C6.squared();
    const sw6_bis_Fq6 N6 = M6*B6;

    const sw6_bis_Fq6 A7 = sw6_bis_exp_by_z(J6);
    const sw6_bis_Fq6 B7 = A7.squared();
    const sw6_bis_Fq6 C7 = sw6_bis_exp_by_z(H6);
    const sw6_bis_Fq6 D7 = C7*B7;
    const sw6_bis_Fq6 E7 = sw6_bis_exp_by_z(L6);
    const sw6_bis_Fq6 F7 = sw6_bis_exp_by_z(N6);
    const sw6_bis_Fq6 G7 = E7*F7;
    const sw6_bis_Fq6 H7 = sw6_bis_exp_by_z(C6);
    const sw6_bis_Fq6 I7 = H7*F7;
    const sw6_bis_Fq6 J7 = I7*G7;

    const sw6_bis_Fq6 A8 = sw6_bis_exp_by_z(D7);
    const sw6_bis_Fq6 B8 = A8.squared();
    const sw6_bis_Fq6 C8 = sw6_bis_exp_by_z(G7);
    const sw6_bis_Fq6 D8 = B8*C8;
    const sw6_bis_Fq6 D8_inv = D8.unitary_inverse();

    const sw6_bis_Fq6 A9 = sw6_bis_exp_by_z(A8);

    const sw6_bis_Fq6 R01 = K0*E1;
    const sw6_bis_Fq6 R02 = D2*R01;
    const sw6_bis_Fq6 R03 = H31*R02;
    const sw6_bis_Fq6 R04 = G4*R03;
    const sw6_bis_Fq6 R05 = A5_inv*R04;
    const sw6_bis_Fq6 R06 = D6_inv*R05;
    const sw6_bis_Fq6 R07 = D7*R06;
    const sw6_bis_Fq6 R0 = R07.unitary_inverse();

    const sw6_bis_Fq6 R11 = M0*C1;
    const sw6_bis_Fq6 R12 = H2_inv*R11;
    const sw6_bis_Fq6 R13 = M3*R12;
    const sw6_bis_Fq6 R14 = D4_inv*R13;
    const sw6_bis_Fq6 R15 = K5_inv*R14;
    const sw6_bis_Fq6 R16 = K6*R15;
    const sw6_bis_Fq6 R17 = J7*R16;
    const sw6_bis_Fq6 R18 = D8_inv*R17;
    const sw6_bis_Fq6 R19 = A9*R18;
    const sw6_bis_Fq6 R1 = R19.Frobenius_map(1);

    const sw6_bis_Fq6 result = R0*R1;

    leave_block("Call to sw6_bis_final_exponentiation_last_chunk");

    return result;
}

sw6_bis_Fq6 sw6_bis_final_exponentiation_first_chunk(const sw6_bis_Fq6 &elt, const sw6_bis_Fq6 &elt_inv)
{
    enter_block("Call to sw6_bis_final_exponentiation_first_chunk");

    /* (q^3-1)*(q+1) */

    /* elt_q3 = elt^(q^3) */
    const sw6_bis_Fq6 elt_q3 = elt.Frobenius_map(3);
    /* elt_q3_over_elt = elt^(q^3-1) */
    const sw6_bis_Fq6 elt_q3_over_elt = elt_q3 * elt_inv;
    /* alpha = elt^((q^3-1) * q) */
    const sw6_bis_Fq6 alpha = elt_q3_over_elt.Frobenius_map(1);
    /* beta = elt^((q^3-1)*(q+1) */
    const sw6_bis_Fq6 beta = alpha * elt_q3_over_elt;
    leave_block("Call to sw6_bis_final_exponentiation_first_chunk");
    return beta;
}

sw6_bis_GT sw6_bis_final_exponentiation(const sw6_bis_Fq6 &elt)
{
    enter_block("Call to sw6_bis_final_exponentiation");
    const sw6_bis_Fq6 elt_inv = elt.inverse();
    sw6_bis_Fq6 elt_to_first_chunk = sw6_bis_final_exponentiation_first_chunk(elt, elt_inv);
    const sw6_bis_Fq6 elt_inv_to_first_chunk = sw6_bis_final_exponentiation_first_chunk(elt_inv, elt);
    sw6_bis_GT result = old_sw6_bis_final_exponentiation_last_chunk(elt_to_first_chunk, elt_inv_to_first_chunk);
    // sw6_bis_GT result = sw6_bis_final_exponentiation_last_chunk(elt_to_first_chunk);

    leave_block("Call to sw6_bis_final_exponentiation");

    return result;
}

/* affine ate miller loop */

sw6_bis_affine_ate_G1_precomputation sw6_bis_affine_ate_precompute_G1(const sw6_bis_G1& P)
{
    enter_block("Call to sw6_bis_affine_ate_precompute_G1");

    sw6_bis_G1 Pcopy = P;
    Pcopy.to_affine_coordinates();

    sw6_bis_affine_ate_G1_precomputation result;
    result.PX = Pcopy.X();
    result.PY = Pcopy.Y();
    result.PY_twist_squared = Pcopy.Y() * sw6_bis_twist.squared();

    leave_block("Call to sw6_bis_affine_ate_precompute_G1");
    return result;
}

sw6_bis_affine_ate_G2_precomputation sw6_bis_affine_ate_precompute_G2(const sw6_bis_G2& Q)
{
    enter_block("Call to sw6_bis_affine_ate_precompute_G2");

    sw6_bis_G2 Qcopy(Q);
    Qcopy.to_affine_coordinates();

    sw6_bis_affine_ate_G2_precomputation result;
    result.QX = Qcopy.X();
    result.QY = Qcopy.Y();

    sw6_bis_Fq3 RX = Qcopy.X();
    sw6_bis_Fq3 RY = Qcopy.Y();

    const bigint<sw6_bis_Fq::num_limbs> &loop_count = sw6_bis_ate_loop_count;
    bool found_nonzero = false;

    std::vector<long> NAF = find_wnaf(1, loop_count);
    for (long i = NAF.size() - 1; i >= 0; --i)
    {
        if (!found_nonzero)
        {
            /* this skips the MSB itself */
            found_nonzero |= (NAF[i] != 0);
            continue;
        }

        sw6_bis_affine_ate_coeffs c;
        c.old_RX = RX;
        c.old_RY = RY;
        sw6_bis_Fq3 old_RX_2 = c.old_RX.squared();
        c.gamma = (old_RX_2 + old_RX_2 + old_RX_2 + sw6_bis_twist_coeff_a) * (c.old_RY + c.old_RY).inverse();
        c.gamma_twist = c.gamma * sw6_bis_twist;
        c.gamma_X = c.gamma * c.old_RX;
        result.coeffs.push_back(c);

        RX = c.gamma.squared() - (c.old_RX+c.old_RX);
        RY = c.gamma * (c.old_RX - RX) - c.old_RY;

        if (NAF[i] != 0)
        {
            sw6_bis_affine_ate_coeffs c;
            c.old_RX = RX;
            c.old_RY = RY;
            if (NAF[i] > 0)
            {
                c.gamma = (c.old_RY - result.QY) * (c.old_RX - result.QX).inverse();
            }
            else
            {
                c.gamma = (c.old_RY + result.QY) * (c.old_RX - result.QX).inverse();
            }
            c.gamma_twist = c.gamma * sw6_bis_twist;
            c.gamma_X = c.gamma * result.QX;
            result.coeffs.push_back(c);

            RX = c.gamma.squared() - (c.old_RX+result.QX);
            RY = c.gamma * (c.old_RX - RX) - c.old_RY;
        }
    }

    /* TODO: maybe handle neg
    if (sw6_bis_ate_is_loop_count_neg)
    {
    	sw6_bis_ate_add_coeffs ac;
		sw6_bis_affine_ate_dbl_coeffs c;
		c.old_RX = RX;
		c.old_RY = -RY;
		old_RX_2 = c.old_RY.squared();
		c.gamma = (old_RX_2 + old_RX_2 + old_RX_2 + sw6_bis_coeff_a) * (c.old_RY + c.old_RY).inverse();
		c.gamma_twist = c.gamma * sw6_bis_twist;
		c.gamma_X = c.gamma * c.old_RX;
		result.coeffs.push_back(c);
    }
    */

    leave_block("Call to sw6_bis_affine_ate_precompute_G2");
    return result;
}

sw6_bis_Fq6 sw6_bis_affine_ate_miller_loop(const sw6_bis_affine_ate_G1_precomputation &prec_P,
                                     const sw6_bis_affine_ate_G2_precomputation &prec_Q)
{
    enter_block("Call to sw6_bis_affine_ate_miller_loop");

    sw6_bis_Fq6 f = sw6_bis_Fq6::one();

    const bigint<sw6_bis_Fq::num_limbs> &loop_count = sw6_bis_ate_loop_count;
    bool found_nonzero = false;
    size_t idx = 0;

    std::vector<long> NAF = find_wnaf(1, loop_count);
    for (long i = NAF.size() - 1; i >= 0; --i)
    {
        if (!found_nonzero)
        {
            /* this skips the MSB itself */
            found_nonzero |= (NAF[i] != 0);
            continue;
        }

        /* code below gets executed for all bits (EXCEPT the MSB itself) of
           sw6_bis_param_p (skipping leading zeros) in MSB to LSB
           order */
        sw6_bis_affine_ate_coeffs c = prec_Q.coeffs[idx++];

        sw6_bis_Fq6 g_RR_at_P = sw6_bis_Fq6(prec_P.PY_twist_squared,
                                      - prec_P.PX * c.gamma_twist + c.gamma_X - c.old_RY);
        f = f.squared().mul_by_2345(g_RR_at_P);

        if (NAF[i] != 0)
        {
            sw6_bis_affine_ate_coeffs c = prec_Q.coeffs[idx++];
            sw6_bis_Fq6 g_RQ_at_P;
            if (NAF[i] > 0)
            {
                g_RQ_at_P = sw6_bis_Fq6(prec_P.PY_twist_squared,
                                     - prec_P.PX * c.gamma_twist + c.gamma_X - prec_Q.QY);
            }
            else
            {
                g_RQ_at_P = sw6_bis_Fq6(prec_P.PY_twist_squared,
                                     - prec_P.PX * c.gamma_twist + c.gamma_X + prec_Q.QY);
            }
            f = f.mul_by_2345(g_RQ_at_P);
        }

    }

    /* TODO: maybe handle neg
    if (sw6_bis_ate_is_loop_count_neg)
    {
    	// TODO:
    	sw6_bis_affine_ate_coeffs ac = prec_Q.coeffs[idx++];
    	sw6_bis_Fq6 g_RnegR_at_P = sw6_bis_Fq6(prec_P.PY_twist_squared,
                                          - prec_P.PX * c.gamma_twist + c.gamma_X - c.old_RY);
    	f = (f * g_RnegR_at_P).inverse();
    }
    */

    leave_block("Call to sw6_bis_affine_ate_miller_loop");

    return f;
}

/* ate pairing */

struct extended_sw6_bis_G2_projective {
    sw6_bis_Fq3 X;
    sw6_bis_Fq3 Y;
    sw6_bis_Fq3 Z;
    sw6_bis_Fq3 T;

    void print() const
    {
        printf("extended sw6_bis_G2 projective X/Y/Z/T:\n");
        X.print();
        Y.print();
        Z.print();
        T.print();
    }

    void test_invariant() const
    {
        assert(T == Z.squared());
    }
};

void doubling_step_for_flipped_miller_loop(extended_sw6_bis_G2_projective &current,
                                           sw6_bis_ate_dbl_coeffs &dc)
{
    const sw6_bis_Fq3 X = current.X, Y = current.Y, Z = current.Z, T = current.T;

    const sw6_bis_Fq3 A = T.squared(); // A = T1^2
    const sw6_bis_Fq3 B = X.squared(); // B = X1^2
    const sw6_bis_Fq3 C = Y.squared(); // C = Y1^2
    const sw6_bis_Fq3 D = C.squared(); // D = C^2
    const sw6_bis_Fq3 E = (X+C).squared() - B - D; // E = (X1+C)^2-B-D
    const sw6_bis_Fq3 F = (B+B+B) + sw6_bis_twist_coeff_a * A; // F = 3*B +  a  *A
    const sw6_bis_Fq3 G = F.squared(); // G = F^2

    current.X = -(E+E+E+E) + G; // X3 = -4*E+G
    current.Y = -sw6_bis_Fq("8")*D + F*(E+E-current.X); // Y3 = -8*D+F*(2*E-X3)
    current.Z = (Y+Z).squared() - C - Z.squared(); // Z3 = (Y1+Z1)^2-C-Z1^2
    current.T = current.Z.squared(); // T3 = Z3^2

    dc.c_H = (current.Z + T).squared() - current.T - A; // H = (Z3+T1)^2-T3-A
    dc.c_4C = C+C+C+C; // fourC = 4*C
    dc.c_J = (F+T).squared() - G - A; // J = (F+T1)^2-G-A
    dc.c_L = (F+X).squared() - G - B; // L = (F+X1)^2-G-B

#ifdef DEBUG
    current.test_invariant();
#endif
}

void mixed_addition_step_for_flipped_miller_loop(const sw6_bis_Fq3 base_X, const sw6_bis_Fq3 base_Y, const sw6_bis_Fq3 base_Y_squared,
                                                 extended_sw6_bis_G2_projective &current,
                                                 sw6_bis_ate_add_coeffs &ac)
{
    const sw6_bis_Fq3 X1 = current.X, Y1 = current.Y, Z1 = current.Z, T1 = current.T;
    const sw6_bis_Fq3 &x2 = base_X,    &y2 =  base_Y, &y2_squared = base_Y_squared;

    const sw6_bis_Fq3 B = x2 * T1; // B = x2 * T1
    const sw6_bis_Fq3 D = ((y2 + Z1).squared() - y2_squared - T1) * T1; // D = ((y2 + Z1)^2 - y2squared - T1) * T1
    const sw6_bis_Fq3 H = B - X1; // H = B - X1
    const sw6_bis_Fq3 I = H.squared(); // I = H^2
    const sw6_bis_Fq3 E = I + I + I + I; // E = 4*I
    const sw6_bis_Fq3 J = H * E; // J = H * E
    const sw6_bis_Fq3 V = X1 * E; // V = X1 * E
    const sw6_bis_Fq3 L1 = D - (Y1 + Y1); // L1 = D - 2 * Y1

    current.X = L1.squared() - J - (V+V); // X3 = L1^2 - J - 2*V
    current.Y = L1 * (V-current.X) - (Y1+Y1) * J; // Y3 = L1 * (V-X3) - 2*Y1 * J
    current.Z = (Z1+H).squared() - T1 - I; // Z3 = (Z1 + H)^2 - T1 - I
    current.T = current.Z.squared(); // T3 = Z3^2

    ac.c_L1 = L1;
    ac.c_RZ = current.Z;
#ifdef DEBUG
    current.test_invariant();
#endif
}

sw6_bis_ate_G1_precomp sw6_bis_ate_precompute_G1(const sw6_bis_G1& P)
{
    enter_block("Call to sw6_bis_ate_precompute_G1");

    sw6_bis_G1 Pcopy = P;
    Pcopy.to_affine_coordinates();

    sw6_bis_ate_G1_precomp result;
    result.PX = Pcopy.X();
    result.PY = Pcopy.Y();
    result.PX_twist = Pcopy.X() * sw6_bis_twist;
    result.PY_twist = Pcopy.Y() * sw6_bis_twist;

    leave_block("Call to sw6_bis_ate_precompute_G1");
    return result;
}

sw6_bis_ate_G2_precomp sw6_bis_ate_precompute_G2(const sw6_bis_G2& Q)
{
    enter_block("Call to sw6_bis_ate_precompute_G2");

    sw6_bis_G2 Qcopy(Q);
    Qcopy.to_affine_coordinates();

    sw6_bis_Fq3 sw6_bis_twist_inv = sw6_bis_twist.inverse(); // could add to global params if needed

    sw6_bis_ate_G2_precomp result;
    result.QX = Qcopy.X();
    result.QY = Qcopy.Y();
    result.QY2 = Qcopy.Y().squared();
    result.QX_over_twist = Qcopy.X() * sw6_bis_twist_inv;
    result.QY_over_twist = Qcopy.Y() * sw6_bis_twist_inv;

    extended_sw6_bis_G2_projective R;
    R.X = Qcopy.X();
    R.Y = Qcopy.Y();
    R.Z = sw6_bis_Fq3::one();
    R.T = sw6_bis_Fq3::one();

    const bigint<sw6_bis_Fq::num_limbs> &loop_count = sw6_bis_ate_loop_count;
    bool found_one = false;
    for (long i = loop_count.max_bits() - 1; i >= 0; --i)
    {
        const bool bit = loop_count.test_bit(i);

        if (!found_one)
        {
            /* this skips the MSB itself */
            found_one |= bit;
            continue;
        }

        sw6_bis_ate_dbl_coeffs dc;
        doubling_step_for_flipped_miller_loop(R, dc);
        result.dbl_coeffs.push_back(dc);

        if (bit)
        {
            sw6_bis_ate_add_coeffs ac;
            mixed_addition_step_for_flipped_miller_loop(result.QX, result.QY, result.QY2, R, ac);
            result.add_coeffs.push_back(ac);
        }
    }

    if (sw6_bis_ate_is_loop_count_neg)
    {
    	sw6_bis_Fq3 RZ_inv = R.Z.inverse();
    	sw6_bis_Fq3 RZ2_inv = RZ_inv.squared();
    	sw6_bis_Fq3 RZ3_inv = RZ2_inv * RZ_inv;
    	sw6_bis_Fq3 minus_R_affine_X = R.X * RZ2_inv;
    	sw6_bis_Fq3 minus_R_affine_Y = - R.Y * RZ3_inv;
    	sw6_bis_Fq3 minus_R_affine_Y2 = minus_R_affine_Y.squared();
    	sw6_bis_ate_add_coeffs ac;
        mixed_addition_step_for_flipped_miller_loop(minus_R_affine_X, minus_R_affine_Y, minus_R_affine_Y2, R, ac);
        result.add_coeffs.push_back(ac);
    }

    leave_block("Call to sw6_bis_ate_precompute_G2");
    return result;
}

sw6_bis_Fq6 sw6_bis_ate_miller_loop(const sw6_bis_ate_G1_precomp &prec_P,
                              const sw6_bis_ate_G2_precomp &prec_Q)
{
    enter_block("Call to sw6_bis_ate_miller_loop");

    sw6_bis_Fq3 L1_coeff = sw6_bis_Fq3(prec_P.PX, sw6_bis_Fq::zero(), sw6_bis_Fq::zero()) - prec_Q.QX_over_twist;

    sw6_bis_Fq6 f = sw6_bis_Fq6::one();

    bool found_one = false;
    size_t dbl_idx = 0;
    size_t add_idx = 0;

    const bigint<sw6_bis_Fq::num_limbs> &loop_count = sw6_bis_ate_loop_count;

    for (long i = loop_count.max_bits() - 1; i >= 0; --i)
    {
        const bool bit = loop_count.test_bit(i);

        if (!found_one)
        {
            /* this skips the MSB itself */
            found_one |= bit;
            continue;
        }

        /* code below gets executed for all bits (EXCEPT the MSB itself) of
           sw6_bis_param_p (skipping leading zeros) in MSB to LSB
           order */
        sw6_bis_ate_dbl_coeffs dc = prec_Q.dbl_coeffs[dbl_idx++];

        sw6_bis_Fq6 g_RR_at_P = sw6_bis_Fq6(- dc.c_4C - dc.c_J * prec_P.PX_twist + dc.c_L,
                                      dc.c_H * prec_P.PY_twist);
        f = f.squared() * g_RR_at_P;

        if (bit)
        {
            sw6_bis_ate_add_coeffs ac = prec_Q.add_coeffs[add_idx++];
            sw6_bis_Fq6 g_RQ_at_P = sw6_bis_Fq6(ac.c_RZ * prec_P.PY_twist,
                                          -(prec_Q.QY_over_twist * ac.c_RZ + L1_coeff * ac.c_L1));
            f = f * g_RQ_at_P;
        }

    }

    if (sw6_bis_ate_is_loop_count_neg)
    {
    	sw6_bis_ate_add_coeffs ac = prec_Q.add_coeffs[add_idx++];
    	sw6_bis_Fq6 g_RnegR_at_P = sw6_bis_Fq6(ac.c_RZ * prec_P.PY_twist,
                                         -(prec_Q.QY_over_twist * ac.c_RZ + L1_coeff * ac.c_L1));
    	f = (f * g_RnegR_at_P).inverse();
    }

    leave_block("Call to sw6_bis_ate_miller_loop");

    return f;
}

sw6_bis_Fq6 sw6_bis_ate_double_miller_loop(const sw6_bis_ate_G1_precomp &prec_P1,
                                     const sw6_bis_ate_G2_precomp &prec_Q1,
                                     const sw6_bis_ate_G1_precomp &prec_P2,
                                     const sw6_bis_ate_G2_precomp &prec_Q2)
{
    enter_block("Call to sw6_bis_ate_double_miller_loop");

    sw6_bis_Fq3 L1_coeff1 = sw6_bis_Fq3(prec_P1.PX, sw6_bis_Fq::zero(), sw6_bis_Fq::zero()) - prec_Q1.QX_over_twist;
    sw6_bis_Fq3 L1_coeff2 = sw6_bis_Fq3(prec_P2.PX, sw6_bis_Fq::zero(), sw6_bis_Fq::zero()) - prec_Q2.QX_over_twist;

    sw6_bis_Fq6 f = sw6_bis_Fq6::one();

    bool found_one = false;
    size_t dbl_idx = 0;
    size_t add_idx = 0;

    const bigint<sw6_bis_Fq::num_limbs> &loop_count = sw6_bis_ate_loop_count;

    for (long i = loop_count.max_bits() - 1; i >= 0; --i)
    {
        const bool bit = loop_count.test_bit(i);

        if (!found_one)
        {
            /* this skips the MSB itself */
            found_one |= bit;
            continue;
        }

        /* code below gets executed for all bits (EXCEPT the MSB itself) of
           sw6_bis_param_p (skipping leading zeros) in MSB to LSB
           order */
        sw6_bis_ate_dbl_coeffs dc1 = prec_Q1.dbl_coeffs[dbl_idx];
        sw6_bis_ate_dbl_coeffs dc2 = prec_Q2.dbl_coeffs[dbl_idx];
        ++dbl_idx;

        sw6_bis_Fq6 g_RR_at_P1 = sw6_bis_Fq6(- dc1.c_4C - dc1.c_J * prec_P1.PX_twist + dc1.c_L,
                                       dc1.c_H * prec_P1.PY_twist);

        sw6_bis_Fq6 g_RR_at_P2 = sw6_bis_Fq6(- dc2.c_4C - dc2.c_J * prec_P2.PX_twist + dc2.c_L,
                                       dc2.c_H * prec_P2.PY_twist);

        f = f.squared() * g_RR_at_P1 * g_RR_at_P2;

        if (bit)
        {
            sw6_bis_ate_add_coeffs ac1 = prec_Q1.add_coeffs[add_idx];
            sw6_bis_ate_add_coeffs ac2 = prec_Q2.add_coeffs[add_idx];
            ++add_idx;

            sw6_bis_Fq6 g_RQ_at_P1 = sw6_bis_Fq6(ac1.c_RZ * prec_P1.PY_twist,
                                           -(prec_Q1.QY_over_twist * ac1.c_RZ + L1_coeff1 * ac1.c_L1));
            sw6_bis_Fq6 g_RQ_at_P2 = sw6_bis_Fq6(ac2.c_RZ * prec_P2.PY_twist,
                                           -(prec_Q2.QY_over_twist * ac2.c_RZ + L1_coeff2 * ac2.c_L1));

            f = f * g_RQ_at_P1 * g_RQ_at_P2;
        }
    }

    if (sw6_bis_ate_is_loop_count_neg)
    {
    	sw6_bis_ate_add_coeffs ac1 = prec_Q1.add_coeffs[add_idx];
        sw6_bis_ate_add_coeffs ac2 = prec_Q2.add_coeffs[add_idx];
    	++add_idx;
    	sw6_bis_Fq6 g_RnegR_at_P1 = sw6_bis_Fq6(ac1.c_RZ * prec_P1.PY_twist,
                                          -(prec_Q1.QY_over_twist * ac1.c_RZ + L1_coeff1 * ac1.c_L1));
    	sw6_bis_Fq6 g_RnegR_at_P2 = sw6_bis_Fq6(ac2.c_RZ * prec_P2.PY_twist,
                                          -(prec_Q2.QY_over_twist * ac2.c_RZ + L1_coeff2 * ac2.c_L1));

    	f = (f * g_RnegR_at_P1 * g_RnegR_at_P2).inverse();
    }

    leave_block("Call to sw6_bis_ate_double_miller_loop");

    return f;
}

sw6_bis_Fq6 sw6_bis_ate_pairing(const sw6_bis_G1& P, const sw6_bis_G2 &Q)
{
    enter_block("Call to sw6_bis_ate_pairing");
    sw6_bis_ate_G1_precomp prec_P = sw6_bis_ate_precompute_G1(P);
    sw6_bis_ate_G2_precomp prec_Q = sw6_bis_ate_precompute_G2(Q);
    sw6_bis_Fq6 result = sw6_bis_ate_miller_loop(prec_P, prec_Q);
    leave_block("Call to sw6_bis_ate_pairing");
    return result;
}

sw6_bis_GT sw6_bis_ate_reduced_pairing(const sw6_bis_G1 &P, const sw6_bis_G2 &Q)
{
    enter_block("Call to sw6_bis_ate_reduced_pairing");
    const sw6_bis_Fq6 f = sw6_bis_ate_pairing(P, Q);
    const sw6_bis_GT result = sw6_bis_final_exponentiation(f);
    leave_block("Call to sw6_bis_ate_reduced_pairing");
    return result;
}

sw6_bis_G1_precomp sw6_bis_precompute_G1(const sw6_bis_G1& P)
{
    return sw6_bis_ate_precompute_G1(P);
}

sw6_bis_G2_precomp sw6_bis_precompute_G2(const sw6_bis_G2& Q)
{
    return sw6_bis_ate_precompute_G2(Q);
}

sw6_bis_Fq6 sw6_bis_miller_loop(const sw6_bis_G1_precomp &prec_P,
                          const sw6_bis_G2_precomp &prec_Q)
{
    return sw6_bis_ate_miller_loop(prec_P, prec_Q);
}

sw6_bis_Fq6 sw6_bis_double_miller_loop(const sw6_bis_G1_precomp &prec_P1,
                                 const sw6_bis_G2_precomp &prec_Q1,
                                 const sw6_bis_G1_precomp &prec_P2,
                                 const sw6_bis_G2_precomp &prec_Q2)
{
    return sw6_bis_ate_double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

sw6_bis_Fq6 sw6_bis_pairing(const sw6_bis_G1& P,
                      const sw6_bis_G2 &Q)
{
    return sw6_bis_ate_pairing(P, Q);
}

sw6_bis_GT sw6_bis_reduced_pairing(const sw6_bis_G1 &P,
                             const sw6_bis_G2 &Q)
{
    return sw6_bis_ate_reduced_pairing(P, Q);
}

sw6_bis_GT sw6_bis_affine_reduced_pairing(const sw6_bis_G1 &P,
                                    const sw6_bis_G2 &Q)
{
    const sw6_bis_affine_ate_G1_precomputation prec_P = sw6_bis_affine_ate_precompute_G1(P);
    const sw6_bis_affine_ate_G2_precomputation prec_Q = sw6_bis_affine_ate_precompute_G2(Q);
    const sw6_bis_Fq6 f = sw6_bis_affine_ate_miller_loop(prec_P, prec_Q);
    const sw6_bis_GT result = sw6_bis_final_exponentiation(f);
    return result;
}

} // libff
