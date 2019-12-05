/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BW12_446_G1_HPP_
#define BW12_446_G1_HPP_
#include <vector>

#include <libff/algebra/curves/bw12_446/bw12_446_init.hpp>
#include <libff/algebra/curves/curve_utils.hpp>

namespace libff {

class bw12_446_G1;
std::ostream& operator<<(std::ostream &, const bw12_446_G1&);
std::istream& operator>>(std::istream &, bw12_446_G1&);

class bw12_446_G1 {
public:
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static bw12_446_G1 G1_zero;
    static bw12_446_G1 G1_one;

    typedef bw12_446_Fq base_field;
    typedef bw12_446_Fr scalar_field;

    bw12_446_Fq X, Y, Z;

    // using Jacobian coordinates
    bw12_446_G1();
    bw12_446_G1(const bw12_446_Fq& X, const bw12_446_Fq& Y, const bw12_446_Fq& Z) : X(X), Y(Y), Z(Z) {};

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const bw12_446_G1 &other) const;
    bool operator!=(const bw12_446_G1 &other) const;

    bw12_446_G1 operator+(const bw12_446_G1 &other) const;
    bw12_446_G1 operator-() const;
    bw12_446_G1 operator-(const bw12_446_G1 &other) const;

    bw12_446_G1 add(const bw12_446_G1 &other) const;
    bw12_446_G1 mixed_add(const bw12_446_G1 &other) const;
    bw12_446_G1 dbl() const;

    bool is_well_formed() const;

    static bw12_446_G1 zero();
    static bw12_446_G1 one();
    static bw12_446_G1 random_element();

    static size_t size_in_bits() { return base_field::size_in_bits() + 1; }
    static bigint<base_field::num_limbs> base_field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    friend std::ostream& operator<<(std::ostream &out, const bw12_446_G1 &g);
    friend std::istream& operator>>(std::istream &in, bw12_446_G1 &g);

    static void batch_to_special_all_non_zeros(std::vector<bw12_446_G1> &vec);
};

template<mp_size_t m>
bw12_446_G1 operator*(const bigint<m> &lhs, const bw12_446_G1 &rhs)
{
    return scalar_mul<bw12_446_G1, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
bw12_446_G1 operator*(const Fp_model<m,modulus_p> &lhs, const bw12_446_G1 &rhs)
{
    return scalar_mul<bw12_446_G1, m>(rhs, lhs.as_bigint());
}

std::ostream& operator<<(std::ostream& out, const std::vector<bw12_446_G1> &v);
std::istream& operator>>(std::istream& in, std::vector<bw12_446_G1> &v);

} // libff
#endif // BW12_446_G1_HPP_
