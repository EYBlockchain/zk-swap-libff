/** @file
 *****************************************************************************
 Declaration of arithmetic in the finite field F[((p^2)^3)^2].
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef FP12_2OVER3OVER2_HPP_
#define FP12_2OVER3OVER2_HPP_
#include <vector>

#include <libff/algebra/fields/fp.hpp>
#include <libff/algebra/fields/fp2.hpp>
#include <libff/algebra/fields/fp6_3over2.hpp>

namespace libff {

template<mp_size_t n, const bigint<n>& modulus>
class Fp12_2over3over2_model;

template<mp_size_t n, const bigint<n>& modulus>
std::ostream& operator<<(std::ostream &, const Fp12_2over3over2_model<n, modulus> &);

template<mp_size_t n, const bigint<n>& modulus>
std::istream& operator>>(std::istream &, Fp12_2over3over2_model<n, modulus> &);

/**
 * Arithmetic in the finite field F[((p^2)^3)^2].
 *
 * Let p := modulus. This interface provides arithmetic for the extension field
 * Fp12 = Fp6[W]/(W^2-V) where Fp6 = Fp2[V]/(V^3-non_residue) and non_residue is in Fp2
 *
 * ASSUMPTION: p = 1 (mod 6)
 */
template<mp_size_t n, const bigint<n>& modulus>
class Fp12_2over3over2_model {
public:
    typedef Fp_model<n, modulus> my_Fp;
    typedef Fp2_model<n, modulus> my_Fp2;
    typedef Fp6_3over2_model<n, modulus> my_Fp6;
    typedef Fp12_2over3over2_model<n, modulus> my_Fp12;

    static Fp2_model<n, modulus> non_residue;
    static Fp2_model<n, modulus> Frobenius_coeffs_c1[12]; // non_residue^((modulus^i-1)/6) for i=0,...,11

    my_Fp6 c0, c1;
    Fp12_2over3over2_model() {};
    Fp12_2over3over2_model(const my_Fp6& c0, const my_Fp6& c1) : c0(c0), c1(c1) {};

    void clear() { c0.clear(); c1.clear(); }
    void print(int indent = 0, const bool as_hex=false) const {
        std::cout << std::string(indent, '\t') << "Fp12_2over3over2.c0:\n";
        c0.print(indent + 1, as_hex);
        std::cout << std::string(indent, '\t') << "Fp12_2over3over2.c1:\n";
        c1.print(indent + 1, as_hex);
    }

    static Fp12_2over3over2_model<n, modulus> zero();
    static Fp12_2over3over2_model<n, modulus> one();
    static Fp12_2over3over2_model<n, modulus> random_element();

    bool is_zero() const { return c0.is_zero() && c1.is_zero(); }
    bool operator==(const Fp12_2over3over2_model &other) const;
    bool operator!=(const Fp12_2over3over2_model &other) const;

    Fp12_2over3over2_model operator+(const Fp12_2over3over2_model &other) const;
    Fp12_2over3over2_model operator-(const Fp12_2over3over2_model &other) const;
    Fp12_2over3over2_model operator*(const Fp12_2over3over2_model &other) const;
    Fp12_2over3over2_model operator-() const;
    void add(const my_Fp12 &x, const my_Fp12 &y);
    void subtract(const my_Fp12 &x, const my_Fp12 &y);
    void negate(const my_Fp12 &x);
    void copy(const my_Fp12 &x);
    void square(const my_Fp12 &x);
    Fp12_2over3over2_model squared() const; // default is squared_complex
    Fp12_2over3over2_model squared_karatsuba() const;
    Fp12_2over3over2_model squared_complex() const;
    void inverse(const my_Fp12 &x);
    Fp12_2over3over2_model inverse() const;
    void multiply(const my_Fp12 &x, const my_Fp12 &y);
    void frobenius_map(const my_Fp12 &x, unsigned long power);
    Fp12_2over3over2_model Frobenius_map(unsigned long power) const;
    Fp12_2over3over2_model unitary_inverse() const;
    void conjugate(const my_Fp12 &x);
    Fp12_2over3over2_model cyclotomic_squared() const;

    void multiply_by_c034(const my_Fp12 &a, const my_Fp2 &c0, const my_Fp2 &c1, const my_Fp2 &c4);
    void multiply_by_c014(const my_Fp12 &a, const my_Fp2 &c0, const my_Fp2 &c1, const my_Fp2 &c4);
    Fp12_2over3over2_model mul_by_024(const my_Fp2 &ell_0, const my_Fp2 &ell_VW, const my_Fp2 &ell_VV) const;

    static my_Fp6 mul_by_non_residue(const my_Fp6 &elt);

    template<mp_size_t m>
    Fp12_2over3over2_model cyclotomic_exp(const bigint<m> &exponent) const;

    static bigint<n> base_field_char() { return modulus; }
    static size_t extension_degree() { return 12; }

    friend std::ostream& operator<< <n, modulus>(std::ostream &out, const Fp12_2over3over2_model<n, modulus> &el);
    friend std::istream& operator>> <n, modulus>(std::istream &in, Fp12_2over3over2_model<n, modulus> &el);
};

template<mp_size_t n, const bigint<n>& modulus>
std::ostream& operator<<(std::ostream& out, const std::vector<Fp12_2over3over2_model<n, modulus> > &v);

template<mp_size_t n, const bigint<n>& modulus>
std::istream& operator>>(std::istream& in, std::vector<Fp12_2over3over2_model<n, modulus> > &v);

template<mp_size_t n, const bigint<n>& modulus>
Fp12_2over3over2_model<n, modulus> operator*(const Fp_model<n, modulus> &lhs, const Fp12_2over3over2_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n>& modulus>
Fp12_2over3over2_model<n, modulus> operator*(const Fp2_model<n, modulus> &lhs, const Fp12_2over3over2_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n>& modulus>
Fp12_2over3over2_model<n, modulus> operator*(const Fp6_3over2_model<n, modulus> &lhs, const Fp12_2over3over2_model<n, modulus> &rhs);

template<mp_size_t n, const bigint<n>& modulus, mp_size_t m>
Fp12_2over3over2_model<n, modulus> operator^(const Fp12_2over3over2_model<n, modulus> &self, const bigint<m> &exponent);

template<mp_size_t n, const bigint<n>& modulus, mp_size_t m, const bigint<m>& exp_modulus>
Fp12_2over3over2_model<n, modulus> operator^(const Fp12_2over3over2_model<n, modulus> &self, const Fp_model<m, exp_modulus> &exponent);

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n, modulus> Fp12_2over3over2_model<n, modulus>::non_residue;

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n, modulus> Fp12_2over3over2_model<n, modulus>::Frobenius_coeffs_c1[12];

} // libff
#include <libff/algebra/fields/fp12_2over3over2.tcc>
#endif // FP12_2OVER3OVER2_HPP_
