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

/*
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
*/

sw6_bis_Fq6 sw6_bis_final_exponentiation_last_chunk(const sw6_bis_Fq6 &elt)
{
    enter_block("Call to sw6_bis_final_exponentiation_last_chunk");
    /*
     * R0(x) := (-103*x^7 + 70*x^6 + 269*x^5 - 197*x^4 - 314*x^3 - 73*x^2 - 263*x - 220)
     * R1(x) := (103*x^9 - 276*x^8 + 77*x^7 + 492*x^6 - 445*x^5 - 65*x^4 + 452*x^3 - 181*x^2 + 34*x + 229)
     *
     * La dernière partie de l'exponentiation finale est
     * f^R0(u)*(f^p)^R1(u) avec R0 et R1 ci-dessus et f^p est un Frobenius dans GF(p^6).
     */

    const sw6_bis_Fq6 f0 = elt;
    const sw6_bis_Fq6 f0p = f0.Frobenius_map(1);
    const sw6_bis_Fq6 f1 = sw6_bis_exp_by_z(f0);
    const sw6_bis_Fq6 f1p = f1.Frobenius_map(1);
    const sw6_bis_Fq6 f2 = sw6_bis_exp_by_z(f1);
    const sw6_bis_Fq6 f2p = f2.Frobenius_map(1);
    const sw6_bis_Fq6 f3 = sw6_bis_exp_by_z(f2);
    const sw6_bis_Fq6 f3p = f3.Frobenius_map(1);
    const sw6_bis_Fq6 f4 = sw6_bis_exp_by_z(f3);
    const sw6_bis_Fq6 f4p = f4.Frobenius_map(1);
    const sw6_bis_Fq6 f5 = sw6_bis_exp_by_z(f4);
    const sw6_bis_Fq6 f5p = f5.Frobenius_map(1);
    const sw6_bis_Fq6 f6 = sw6_bis_exp_by_z(f5);
    const sw6_bis_Fq6 f6p = f6.Frobenius_map(1);
    const sw6_bis_Fq6 f7 = sw6_bis_exp_by_z(f6);
    const sw6_bis_Fq6 f7p = f7.Frobenius_map(1);

    // 4
    const sw6_bis_Fq6 f8p = sw6_bis_exp_by_z(f7p);
    const sw6_bis_Fq6 f9p = sw6_bis_exp_by_z(f8p);

    // 5
    const sw6_bis_Fq6 result1 = f3p * f6p * f5p.Frobenius_map(3);

    // 6
    const sw6_bis_Fq6 result2 = result1.squared();
    const sw6_bis_Fq6 f4_2p = f4 * f2p;
    const sw6_bis_Fq6 result3 = result2 * f5 * f0p * (f0 * f1 * f3 * f4_2p * f8p).Frobenius_map(3);

    // 7
    const sw6_bis_Fq6 result4 = result3.squared();
    const sw6_bis_Fq6 result5 = result4 * f9p * f7.Frobenius_map(3);

    // 8
    const sw6_bis_Fq6 result6 = result5.squared();
    const sw6_bis_Fq6 f2_4p = f2 * f4p;
    const sw6_bis_Fq6 f4_2p_5p = f4_2p * f5p;
    const sw6_bis_Fq6 result7 = result6 * f4_2p_5p * f6 * f7p * (f2_4p * f3 * f3p).Frobenius_map(3);

    // 9
    const sw6_bis_Fq6 result8 = result7.squared();
    const sw6_bis_Fq6 result9 = result8 * f0 * f7 * f1p * (f0p * f9p).Frobenius_map(3);

    // 10
    const sw6_bis_Fq6 result10 = result9.squared();
    const sw6_bis_Fq6 f6p_8p = f6p * f8p;
    const sw6_bis_Fq6 f5_7p = f5 * f7p;
    const sw6_bis_Fq6 result11 = result10 * f5_7p * f2p * (f6p_8p).Frobenius_map(3);

    // 11
    const sw6_bis_Fq6 result12 = result11.squared();
    const sw6_bis_Fq6 f3_6 = f3 * f6;
    const sw6_bis_Fq6 f1_7 = f1 * f7;
    const sw6_bis_Fq6 result13 = result12 * f3_6 * f9p * (f1_7 * f2).Frobenius_map(3);

    // 12
    const sw6_bis_Fq6 result14 = result13.squared();
    const sw6_bis_Fq6 result15 = result14 * f0 * f0p * f3p * f5p * (f4_2p * f5_7p * f6p_8p).Frobenius_map(3);

    // 13
    const sw6_bis_Fq6 result16 = result15.squared();
    const sw6_bis_Fq6 result17 = result16 * f1p * (f3_6).Frobenius_map(3);

    // 14
    const sw6_bis_Fq6 result18 = result17.squared();
    const sw6_bis_Fq6 result19 = result18 * f1_7 * f5_7p * f0p * (f2_4p * f4_2p_5p * f9p).Frobenius_map(3);

    leave_block("Call to sw6_bis_final_exponentiation_last_chunk");

    return result19;
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
    // sw6_bis_GT result = old_sw6_bis_final_exponentiation_last_chunk(elt_to_first_chunk, elt_inv_to_first_chunk);
    sw6_bis_GT result = sw6_bis_final_exponentiation_last_chunk(elt_to_first_chunk);

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

sw6_bis_ate_G2_precomp sw6_bis_ate_precompute_G2(const sw6_bis_G2& Q, const bigint<sw6_bis_Fq::num_limbs> &loop_count)
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
                              const sw6_bis_ate_G2_precomp &prec_Q_1, const sw6_bis_ate_G2_precomp &prec_Q_2)
{
    // f_{u+1,Q}(P)
    enter_block("Call to sw6_bis_ate_miller_loop 1");

    sw6_bis_Fq3 L1_coeff_1 = sw6_bis_Fq3(prec_P.PX, sw6_bis_Fq::zero(), sw6_bis_Fq::zero()) - prec_Q_1.QX_over_twist;

    sw6_bis_Fq6 f_1 = sw6_bis_Fq6::one();

    bool found_one_1 = false;
    size_t dbl_idx_1 = 0;
    size_t add_idx_1 = 0;

    const bigint<sw6_bis_Fq::num_limbs> &loop_count_1 = sw6_bis_ate_loop_count1;

    for (long i = loop_count_1.max_bits() - 1; i >= 0; --i)
    {
        const bool bit = loop_count_1.test_bit(i);

        if (!found_one_1)
        {
            /* this skips the MSB itself */
            found_one_1 |= bit;
            continue;
        }

        /* code below gets executed for all bits (EXCEPT the MSB itself) of
           sw6_bis_param_p (skipping leading zeros) in MSB to LSB
           order */
        sw6_bis_ate_dbl_coeffs dc = prec_Q_1.dbl_coeffs[dbl_idx_1++];

        sw6_bis_Fq6 g_RR_at_P = sw6_bis_Fq6(- dc.c_4C - dc.c_J * prec_P.PX_twist + dc.c_L,
                                      dc.c_H * prec_P.PY_twist);
        f_1 = f_1.squared() * g_RR_at_P;

        if (bit)
        {
            sw6_bis_ate_add_coeffs ac = prec_Q_1.add_coeffs[add_idx_1++];
            sw6_bis_Fq6 g_RQ_at_P = sw6_bis_Fq6(ac.c_RZ * prec_P.PY_twist,
                                          -(prec_Q_1.QY_over_twist * ac.c_RZ + L1_coeff_1 * ac.c_L1));
            f_1 = f_1 * g_RQ_at_P;
        }

    }

    if (sw6_bis_ate_is_loop_count_neg)
    {
    	sw6_bis_ate_add_coeffs ac = prec_Q_1.add_coeffs[add_idx_1++];
    	sw6_bis_Fq6 g_RnegR_at_P = sw6_bis_Fq6(ac.c_RZ * prec_P.PY_twist,
                                         -(prec_Q_1.QY_over_twist * ac.c_RZ + L1_coeff_1 * ac.c_L1));
    	f_1 = (f_1 * g_RnegR_at_P).inverse();
    }

    leave_block("Call to sw6_bis_ate_miller_loop 1");

    // f_{u^3-u^2-u,Q}(P)^q
    enter_block("Call to sw6_bis_ate_miller_loop 2");

    sw6_bis_Fq3 L1_coeff_2 = sw6_bis_Fq3(prec_P.PX, sw6_bis_Fq::zero(), sw6_bis_Fq::zero()) - prec_Q_2.QX_over_twist;

    sw6_bis_Fq6 f_2 = sw6_bis_Fq6::one();

    bool found_one_2 = false;
    size_t dbl_idx_2 = 0;
    size_t add_idx_2 = 0;

    const bigint<sw6_bis_Fq::num_limbs> &loop_count_2 = sw6_bis_ate_loop_count2;

    for (long i = loop_count_2.max_bits() - 1; i >= 0; --i)
    {
        const bool bit = loop_count_2.test_bit(i);

        if (!found_one_2)
        {
            /* this skips the MSB itself */
            found_one_2 |= bit;
            continue;
        }

        /* code below gets executed for all bits (EXCEPT the MSB itself) of
           sw6_bis_param_p (skipping leading zeros) in MSB to LSB
           order */
        sw6_bis_ate_dbl_coeffs dc = prec_Q_2.dbl_coeffs[dbl_idx_2++];

        sw6_bis_Fq6 g_RR_at_P = sw6_bis_Fq6(- dc.c_4C - dc.c_J * prec_P.PX_twist + dc.c_L,
                                      dc.c_H * prec_P.PY_twist);
        f_2 = f_2.squared() * g_RR_at_P;

        if (bit)
        {
            sw6_bis_ate_add_coeffs ac = prec_Q_2.add_coeffs[add_idx_2++];
            sw6_bis_Fq6 g_RQ_at_P = sw6_bis_Fq6(ac.c_RZ * prec_P.PY_twist,
                                          -(prec_Q_2.QY_over_twist * ac.c_RZ + L1_coeff_2 * ac.c_L1));
            f_2 = f_2 * g_RQ_at_P;
        }

    }

    if (sw6_bis_ate_is_loop_count_neg)
    {
    	sw6_bis_ate_add_coeffs ac = prec_Q_2.add_coeffs[add_idx_2++];
    	sw6_bis_Fq6 g_RnegR_at_P = sw6_bis_Fq6(ac.c_RZ * prec_P.PY_twist,
                                         -(prec_Q_2.QY_over_twist * ac.c_RZ + L1_coeff_2 * ac.c_L1));
    	f_2 = (f_2 * g_RnegR_at_P).inverse();
    }
    f_2 = f_2.Frobenius_map(1);

    leave_block("Call to sw6_bis_ate_miller_loop 2");

    return f_1 * f_2;
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
    // f_{u+1,Q}
    sw6_bis_ate_G2_precomp prec_Q_1 = sw6_bis_ate_precompute_G2(Q, sw6_bis_ate_loop_count1);
    // f_{u^3-u^2-u,Q}
    sw6_bis_ate_G2_precomp prec_Q_2 = sw6_bis_ate_precompute_G2(Q, sw6_bis_ate_loop_count2);
    sw6_bis_Fq6 result = sw6_bis_ate_miller_loop(prec_P, prec_Q_1, prec_Q_2);
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

sw6_bis_G2_precomp sw6_bis_precompute_G2(const sw6_bis_G2& Q, const bigint<sw6_bis_Fq::num_limbs> &loop_count)
{
    return sw6_bis_ate_precompute_G2(Q, loop_count);
}

sw6_bis_Fq6 sw6_bis_miller_loop(const sw6_bis_G1_precomp &prec_P,
                          const sw6_bis_G2_precomp &prec_Q_1, const sw6_bis_G2_precomp &prec_Q_2)
{
    return sw6_bis_ate_miller_loop(prec_P, prec_Q_1, prec_Q_2);
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
