/** @file
 *****************************************************************************

 Declaration of interfaces for the SW6 G1 group.

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef SW6_G1_HPP_
#define SW6_G1_HPP_

#include <vector>

#include <libff/algebra/curves/curve_utils.hpp>
#include <libff/algebra/curves/sw6/sw6_init.hpp>

namespace libff {

class sw6_G1;
std::ostream& operator<<(std::ostream &, const sw6_G1&);
std::istream& operator>>(std::istream &, sw6_G1&);

class sw6_G1 {
public:
    sw6_Fq X, Y, Z;
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static sw6_G1 G1_zero;
    static sw6_G1 G1_one;
    static sw6_Fq coeff_a;
    static sw6_Fq coeff_b;

    typedef sw6_Fq base_field;
    typedef sw6_Fr scalar_field;

    // using projective coordinates
    sw6_G1();
    sw6_G1(const sw6_Fq& X, const sw6_Fq& Y) : X(X), Y(Y), Z(base_field::one()) {}
    sw6_G1(const sw6_Fq& X, const sw6_Fq& Y, const sw6_Fq& Z) : X(X), Y(Y), Z(Z) {}

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const sw6_G1 &other) const;
    bool operator!=(const sw6_G1 &other) const;

    sw6_G1 operator+(const sw6_G1 &other) const;
    sw6_G1 operator-() const;
    sw6_G1 operator-(const sw6_G1 &other) const;

    sw6_G1 add(const sw6_G1 &other) const;
    sw6_G1 mixed_add(const sw6_G1 &other) const;
    sw6_G1 dbl() const;

    bool is_well_formed() const;

    static sw6_G1 zero();
    static sw6_G1 one();
    static sw6_G1 random_element();

    static size_t size_in_bits() { return base_field::size_in_bits() + 1; }
    static bigint<base_field::num_limbs> base_field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    friend std::ostream& operator<<(std::ostream &out, const sw6_G1 &g);
    friend std::istream& operator>>(std::istream &in, sw6_G1 &g);

    static void batch_to_special_all_non_zeros(std::vector<sw6_G1> &vec);
};

template<mp_size_t m>
sw6_G1 operator*(const bigint<m> &lhs, const sw6_G1 &rhs)
{
    return scalar_mul<sw6_G1, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
sw6_G1 operator*(const Fp_model<m,modulus_p> &lhs, const sw6_G1 &rhs)
{
    return scalar_mul<sw6_G1, m>(rhs, lhs.as_bigint());
}

std::ostream& operator<<(std::ostream& out, const std::vector<sw6_G1> &v);
std::istream& operator>>(std::istream& in, std::vector<sw6_G1> &v);

} // libff

#endif // SW6_G1_HPP_
