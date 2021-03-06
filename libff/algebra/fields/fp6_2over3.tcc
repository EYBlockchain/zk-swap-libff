/** @file
 *****************************************************************************
 Implementation of arithmetic in the finite field F[(p^3)^2].
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef FP6_2OVER3_TCC_
#define FP6_2OVER3_TCC_
#include <libff/algebra/fields/field_utils.hpp>
#include <libff/algebra/scalar_multiplication/wnaf.hpp>

namespace libff {

template<mp_size_t n, const bigint<n>& modulus>
Fp3_model<n,modulus> Fp6_2over3_model<n, modulus>::mul_by_non_residue(const Fp3_model<n,modulus> &elem)
{
    return Fp3_model<n, modulus>(non_residue * elem.c2, elem.c0, elem.c1);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp6_2over3_model<n, modulus> Fp6_2over3_model<n, modulus>::zero()
{
    return Fp6_2over3_model<n,modulus>(my_Fp3::zero(),
                                my_Fp3::zero());
}

template<mp_size_t n, const bigint<n>& modulus>
Fp6_2over3_model<n, modulus> Fp6_2over3_model<n, modulus>::one()
{
    return Fp6_2over3_model<n,modulus>(my_Fp3::one(),
                                my_Fp3::zero());
}

template<mp_size_t n, const bigint<n>& modulus>
Fp6_2over3_model<n,modulus> Fp6_2over3_model<n,modulus>::random_element()
{
    Fp6_2over3_model<n, modulus> r;
    r.c0 = my_Fp3::random_element();
    r.c1 = my_Fp3::random_element();

    return r;
}

template<mp_size_t n, const bigint<n>& modulus>
bool Fp6_2over3_model<n,modulus>::operator==(const Fp6_2over3_model<n,modulus> &other) const
{
    return (this->c0 == other.c0 && this->c1 == other.c1);
}

template<mp_size_t n, const bigint<n>& modulus>
bool Fp6_2over3_model<n,modulus>::operator!=(const Fp6_2over3_model<n,modulus> &other) const
{
    return !(operator==(other));
}

template<mp_size_t n, const bigint<n>& modulus>
Fp6_2over3_model<n,modulus> Fp6_2over3_model<n,modulus>::operator+(const Fp6_2over3_model<n,modulus> &other) const
{
    return Fp6_2over3_model<n,modulus>(this->c0 + other.c0,
                                this->c1 + other.c1);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp6_2over3_model<n,modulus> Fp6_2over3_model<n,modulus>::operator-(const Fp6_2over3_model<n,modulus> &other) const
{
    return Fp6_2over3_model<n,modulus>(this->c0 - other.c0,
                                this->c1 - other.c1);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp6_2over3_model<n, modulus> operator*(const Fp_model<n, modulus> &lhs, const Fp6_2over3_model<n, modulus> &rhs)
{
    return Fp6_2over3_model<n,modulus>(lhs*rhs.c0,
                                lhs*rhs.c1);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp6_2over3_model<n,modulus> Fp6_2over3_model<n,modulus>::operator*(const Fp6_2over3_model<n,modulus> &other) const
{
    /* Devegili OhEig Scott Dahab --- Multiplication and Squaring on Pairing-Friendly Fields.pdf; Section 3 (Karatsuba) */

    const my_Fp3 &B = other.c1, &A = other.c0,
                 &b = this->c1, &a = this->c0;
    const my_Fp3 aA = a*A;
    const my_Fp3 bB = b*B;
    const my_Fp3 beta_bB = Fp6_2over3_model<n,modulus>::mul_by_non_residue(bB);

    return Fp6_2over3_model<n,modulus>(aA + beta_bB,
                                       (a+b)*(A+B) - aA  - bB);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp6_2over3_model<n,modulus> Fp6_2over3_model<n,modulus>::mul_by_045(const Fp_model<n, modulus> &ell_0,
                                                                    const Fp_model<n, modulus> &ell_VW,
                                                                    const Fp_model<n, modulus> &ell_VV) const
{
  /*
    // OLD
    Fp6_2over3_model<n,modulus> a(my_Fp3(ell_VW, my_Fp::zero(), my_Fp::zero()),
                                  my_Fp3(my_Fp::zero(), ell_0, ell_VV));

    return (*this) * a;
  */

    my_Fp z0 = this->c0.c0;
    my_Fp z1 = this->c0.c1;
    my_Fp z2 = this->c0.c2;
    my_Fp z3 = this->c1.c0;
    my_Fp z4 = this->c1.c1;
    my_Fp z5 = this->c1.c2;

    my_Fp x0 = ell_VW;
    my_Fp x4 = ell_0;
    my_Fp x5 = ell_VV;

    my_Fp t0, t1, t2, t3, t4, t5;
    my_Fp tmp1, tmp2;

    tmp1 = my_Fp3::non_residue * x4;
    tmp2 = my_Fp3::non_residue * x5;

    t0 = x0 * z0 + tmp1 * z4 + tmp2 * z3;
    t1 = x0 * z1 + tmp1 * z5 + tmp2 * z4;
    t2 = x0 * z2 + x4 * z3 + tmp2 * z5;
    t3 = x0 * z3 + tmp1 * z2 + tmp2 * z1;
    t4 = x0 * z4 + x4 * z0 + tmp2 * z2;
    t5 = x0 * z5 + x4 * z1 + x5 * z0;

    return Fp6_2over3_model<n,modulus>(my_Fp3(t0,t1,t2),my_Fp3(t3,t4,t5));
}

template<mp_size_t n, const bigint<n>& modulus>
Fp6_2over3_model<n,modulus> Fp6_2over3_model<n,modulus>::mul_by_024(const Fp_model<n, modulus> &ell_0,
                                                                                const Fp_model<n, modulus> &ell_VW,
                                                                                const Fp_model<n, modulus> &ell_VV) const
{
  /*
       Fp6_2over3_model<n,modulus> a(my_Fp3(ell_0, my_Fp::zero(), ell_VV),
       my_Fp3(my_Fp::zero(), ell_VW, my_Fp::zero()));

       return (*this) * a;
  */
    my_Fp z0 = this->c0.c0;
    my_Fp z1 = this->c0.c1;
    my_Fp z2 = this->c0.c2;
    my_Fp z3 = this->c1.c0;
    my_Fp z4 = this->c1.c1;
    my_Fp z5 = this->c1.c2;

    my_Fp x0 = ell_0;
    my_Fp x2 = ell_VV;
    my_Fp x4 = ell_VW;

    my_Fp t0, t1, t2, s0, T3, T4, D0, D2, D4, S1;

    D0 = z0 * x0;
    D2 = z2 * x2;
    D4 = z4 * x4;
    t2 = z0 + z4;
    t1 = z0 + z2;
    s0 = z1 + z3 + z5;

    // For z.a_.a_ = z0.
    S1 = z1 * x2;
    T3 = S1 + D4;
    T4 = my_Fp3::non_residue * T3 + D0;
    z0 = T4;

    // For z.a_.b_ = z1
    T3 = z5 * x4;
    S1 = S1 + T3;
    T3 = T3 + D2;
    T4 = my_Fp3::non_residue * T3;
    T3 = z1 * x0;
    S1 = S1 + T3;
    T4 = T4 + T3;
    z1 = T4;

    // For z.a_.c_ = z2
    t0 = x0 + x2;
    T3 = t1 * t0 - D0 - D2;
    T4 = z3 * x4;
    S1 = S1 + T4;
    T3 = T3 + T4;

    // For z.b_.a_ = z3 (z3 needs z2)
    t0 = z2 + z4;
    z2 = T3;
    t1 = x2 + x4;
    T3 = t0 * t1 - D2 - D4;
    T4 = my_Fp3::non_residue * T3;
    T3 = z3 * x0;
    S1 = S1 + T3;
    T4 = T4 + T3;
    z3 = T4;

    // For z.b_.b_ = z4
    T3 = z5 * x2;
    S1 = S1 + T3;
    T4 = my_Fp3::non_residue * T3;
    t0 = x0 + x4;
    T3 = t2 * t0 - D0 - D4;
    T4 = T4 + T3;
    z4 = T4;

    // For z.b_.c_ = z5.
    t0 = x0 + x2 + x4;
    T3 = s0 * t0 - S1;
    z5 = T3;

    return Fp6_2over3_model<n,modulus>(my_Fp3(z0,z1,z2),my_Fp3(z3,z4,z5));
}

template<mp_size_t n, const bigint<n>& modulus>
Fp6_2over3_model<n,modulus> Fp6_2over3_model<n,modulus>::mul_by_2345(const Fp6_2over3_model<n,modulus> &other) const
{
    /* Devegili OhEig Scott Dahab --- Multiplication and Squaring on Pairing-Friendly Fields.pdf; Section 3 (Karatsuba) */
    assert(other.c0.c0.is_zero());
    assert(other.c0.c1.is_zero());

    const my_Fp3 &B = other.c1, &A = other.c0,
                 &b = this->c1, &a = this->c0;
    const my_Fp3 aA = my_Fp3(a.c1 * A.c2 * non_residue, a.c2 * A.c2 * non_residue, a.c0 * A.c2);
    const my_Fp3 bB = b*B;
    const my_Fp3 beta_bB = Fp6_2over3_model<n,modulus>::mul_by_non_residue(bB);

    return Fp6_2over3_model<n,modulus>(aA + beta_bB,
                                       (a+b)*(A+B) - aA  - bB);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp6_2over3_model<n,modulus> Fp6_2over3_model<n,modulus>::operator-() const
{
    return Fp6_2over3_model<n,modulus>(-this->c0,
                                -this->c1);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp6_2over3_model<n,modulus> Fp6_2over3_model<n,modulus>::squared() const
{
    /* Devegili OhEig Scott Dahab --- Multiplication and Squaring on Pairing-Friendly Fields.pdf; Section 3 (Complex) */
    const my_Fp3 &b = this->c1, &a = this->c0;
    const my_Fp3 ab = a * b;

    return Fp6_2over3_model<n,modulus>((a+b)*(a+Fp6_2over3_model<n,modulus>::mul_by_non_residue(b))-ab-Fp6_2over3_model<n,modulus>::mul_by_non_residue(ab),
                                ab + ab);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp6_2over3_model<n,modulus> Fp6_2over3_model<n,modulus>::inverse() const
{
    /* From "High-Speed Software Implementation of the Optimal Ate Pairing over Barreto-Naehrig Curves"; Algorithm 8 */

    const my_Fp3 &b = this->c1, &a = this->c0;
    const my_Fp3 t1 = b.squared();
    const my_Fp3 t0 = a.squared() - Fp6_2over3_model<n,modulus>::mul_by_non_residue(t1);
    const my_Fp3 new_t1 = t0.inverse();

    return Fp6_2over3_model<n,modulus>(a * new_t1,
                                       - (b * new_t1));
}

template<mp_size_t n, const bigint<n>& modulus>
Fp6_2over3_model<n,modulus> Fp6_2over3_model<n,modulus>::Frobenius_map(unsigned long power) const
{
    return Fp6_2over3_model<n,modulus>(c0.Frobenius_map(power),
                                       Frobenius_coeffs_c1[power % 6] * c1.Frobenius_map(power));
}

template<mp_size_t n, const bigint<n>& modulus>
Fp6_2over3_model<n,modulus> Fp6_2over3_model<n,modulus>::unitary_inverse() const
{
    return Fp6_2over3_model<n,modulus>(this->c0,
                                -this->c1);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp6_2over3_model<n,modulus> Fp6_2over3_model<n,modulus>::cyclotomic_squared() const
{
    my_Fp2 a = my_Fp2(c0.c0, c1.c1);
    //my_Fp a_a = c0.c0; // a = Fp2([c0[0],c1[1]])
    //my_Fp a_b = c1.c1;

    my_Fp2 b = my_Fp2(c1.c0, c0.c2);
    //my_Fp b_a = c1.c0; // b = Fp2([c1[0],c0[2]])
    //my_Fp b_b = c0.c2;

    my_Fp2 c = my_Fp2(c0.c1, c1.c2);
    //my_Fp c_a = c0.c1; // c = Fp2([c0[1],c1[2]])
    //my_Fp c_b = c1.c2;

    my_Fp2 asq = a.squared();
    my_Fp2 bsq = b.squared();
    my_Fp2 csq = c.squared();

    // A = vector(3*a^2 - 2*Fp2([vector(a)[0],-vector(a)[1]]))
    //my_Fp A_a = my_Fp(3l) * asq_a - my_Fp(2l) * a_a;
    my_Fp A_a = asq.c0 - a.c0;
    A_a = A_a + A_a + asq.c0;
    //my_Fp A_b = my_Fp(3l) * asq_b + my_Fp(2l) * a_b;
    my_Fp A_b = asq.c1 + a.c1;
    A_b = A_b + A_b + asq.c1;

    // B = vector(3*Fp2([non_residue*c2[1],c2[0]]) + 2*Fp2([vector(b)[0],-vector(b)[1]]))
    //my_Fp B_a = my_Fp(3l) * my_Fp3::non_residue * csq_b + my_Fp(2l) * b_a;
    my_Fp B_tmp = my_Fp3::non_residue * csq.c1;
    my_Fp B_a = B_tmp + b.c0;
    B_a = B_a + B_a + B_tmp;

    //my_Fp B_b = my_Fp(3l) * csq_a - my_Fp(2l) * b_b;
    my_Fp B_b = csq.c0 - b.c1;
    B_b = B_b + B_b + csq.c0;

    // C = vector(3*b^2 - 2*Fp2([vector(c)[0],-vector(c)[1]]))
    //my_Fp C_a = my_Fp(3l) * bsq_a - my_Fp(2l) * c_a;
    my_Fp C_a = bsq.c0 - c.c0;
    C_a = C_a + C_a + bsq.c0;
    // my_Fp C_b = my_Fp(3l) * bsq_b + my_Fp(2l) * c_b;
    my_Fp C_b = bsq.c1 + c.c1;
    C_b = C_b + C_b + bsq.c1;

    // e0 = Fp3([A[0],C[0],B[1]])
    // e1 = Fp3([B[0],A[1],C[1]])
    // fin = Fp6e([e0,e1])
    // return fin

    return Fp6_2over3_model<n, modulus>(my_Fp3(A_a, C_a, B_b),
                                        my_Fp3(B_a, A_b, C_b));
}

template<mp_size_t n, const bigint<n>& modulus>
template<mp_size_t m>
Fp6_2over3_model<n, modulus> Fp6_2over3_model<n,modulus>::cyclotomic_exp(const bigint<m> &exponent) const
{
    Fp6_2over3_model<n,modulus> res = Fp6_2over3_model<n,modulus>::one();
    Fp6_2over3_model<n,modulus> this_inverse = this->unitary_inverse();

    bool found_nonzero = false;
    std::vector<long> NAF = find_wnaf(1, exponent);

    for (long i = static_cast<long>(NAF.size() - 1); i >= 0; --i)
    {
        if (found_nonzero)
        {
            res = res.cyclotomic_squared();
        }

        if (NAF[i] != 0)
        {
            found_nonzero = true;

            if (NAF[i] > 0)
            {
                res = res * (*this);
            }
            else
            {
                res = res * this_inverse;
            }
        }
    }

    return res;
}

template<mp_size_t n, const bigint<n>& modulus>
std::ostream& operator<<(std::ostream &out, const Fp6_2over3_model<n, modulus> &el)
{
    out << el.c0 << OUTPUT_SEPARATOR << el.c1;
    return out;
}

template<mp_size_t n, const bigint<n>& modulus>
std::istream& operator>>(std::istream &in, Fp6_2over3_model<n, modulus> &el)
{
    in >> el.c0 >> el.c1;
    return in;
}

template<mp_size_t n, const bigint<n>& modulus, mp_size_t m>
Fp6_2over3_model<n, modulus> operator^(const Fp6_2over3_model<n, modulus> &self, const bigint<m> &exponent)
{
    return power<Fp6_2over3_model<n, modulus>, m>(self, exponent);
}

template<mp_size_t n, const bigint<n>& modulus, mp_size_t m, const bigint<m>& exp_modulus>
Fp6_2over3_model<n, modulus> operator^(const Fp6_2over3_model<n, modulus> &self, const Fp_model<m, exp_modulus> &exponent)
{
    return self^(exponent.as_bigint());
}

} // libff
#endif // FP6_2OVER3_TCC_
