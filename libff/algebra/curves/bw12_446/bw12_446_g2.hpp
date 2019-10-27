/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BW12_446_G2_HPP_
#define BW12_446_G2_HPP_
#include <vector>

#include <libff/algebra/curves/bw12_446/bw12_446_init.hpp>
#include <libff/algebra/curves/curve_utils.hpp>

namespace libff {

class bw12_446_G2;
std::ostream& operator<<(std::ostream &, const bw12_446_G2&);
std::istream& operator>>(std::istream &, bw12_446_G2&);

class bw12_446_G2 {
public:
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static bw12_446_G2 G2_zero;
    static bw12_446_G2 G2_one;

    typedef bw12_446_Fq base_field;
    typedef bw12_446_Fq2 twist_field;
    typedef bw12_446_Fr scalar_field;

    bw12_446_Fq2 X, Y, Z;

    // using Jacobian coordinates
    bw12_446_G2();
    bw12_446_G2(const bw12_446_Fq2& X, const bw12_446_Fq2& Y, const bw12_446_Fq2& Z) : X(X), Y(Y), Z(Z) {};

    static bw12_446_Fq2 mul_by_b(const bw12_446_Fq2 &elt);

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const bw12_446_G2 &other) const;
    bool operator!=(const bw12_446_G2 &other) const;

    bw12_446_G2 operator+(const bw12_446_G2 &other) const;
    bw12_446_G2 operator-() const;
    bw12_446_G2 operator-(const bw12_446_G2 &other) const;

    bw12_446_G2 add(const bw12_446_G2 &other) const;
    bw12_446_G2 mixed_add(const bw12_446_G2 &other) const;
    bw12_446_G2 dbl() const;
    bw12_446_G2 mul_by_q() const;

    bool is_well_formed() const;

    static bw12_446_G2 zero();
    static bw12_446_G2 one();
    static bw12_446_G2 random_element();

    static size_t size_in_bits() { return twist_field::size_in_bits() + 1; }
    static bigint<base_field::num_limbs> base_field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    friend std::ostream& operator<<(std::ostream &out, const bw12_446_G2 &g);
    friend std::istream& operator>>(std::istream &in, bw12_446_G2 &g);

    static void batch_to_special_all_non_zeros(std::vector<bw12_446_G2> &vec);
};

template<mp_size_t m>
bw12_446_G2 operator*(const bigint<m> &lhs, const bw12_446_G2 &rhs)
{
    return scalar_mul<bw12_446_G2, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
bw12_446_G2 operator*(const Fp_model<m,modulus_p> &lhs, const bw12_446_G2 &rhs)
{
    return scalar_mul<bw12_446_G2, m>(rhs, lhs.as_bigint());
}


} // libff
#endif // BW12_446_G2_HPP_
