/** @file
 *****************************************************************************

 Declaration of interfaces for the SW6_BIS G1 group.

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef SW6_BIS_G1_HPP_
#define SW6_BIS_G1_HPP_

#include <vector>

#include <libff/algebra/curves/curve_utils.hpp>
#include <libff/algebra/curves/hg6/hg6_init.hpp>

namespace libff {

class hg6_G1;
std::ostream& operator<<(std::ostream &, const hg6_G1&);
std::istream& operator>>(std::istream &, hg6_G1&);

class hg6_G1 {
private:
    hg6_Fq X_, Y_, Z_;
public:
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static hg6_G1 G1_zero;
    static hg6_G1 G1_one;
    static hg6_Fq coeff_a;
    static hg6_Fq coeff_b;

    typedef hg6_Fq base_field;
    typedef hg6_Fr scalar_field;

    // using projective coordinates
    hg6_G1();
    hg6_G1(const hg6_Fq& X, const hg6_Fq& Y) : X_(X), Y_(Y), Z_(base_field::one()) {}
    hg6_G1(const hg6_Fq& X, const hg6_Fq& Y, const hg6_Fq& Z) : X_(X), Y_(Y), Z_(Z) {}

    hg6_Fq X() const { return X_; }
    hg6_Fq Y() const { return Y_; }
    hg6_Fq Z() const { return Z_; }

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const hg6_G1 &other) const;
    bool operator!=(const hg6_G1 &other) const;

    hg6_G1 operator+(const hg6_G1 &other) const;
    hg6_G1 operator-() const;
    hg6_G1 operator-(const hg6_G1 &other) const;

    hg6_G1 add(const hg6_G1 &other) const;
    hg6_G1 mixed_add(const hg6_G1 &other) const;
    hg6_G1 dbl() const;

    bool is_well_formed() const;

    static hg6_G1 zero();
    static hg6_G1 one();
    static hg6_G1 random_element();

    static size_t size_in_bits() { return base_field::size_in_bits() + 1; }
    static bigint<base_field::num_limbs> base_field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    friend std::ostream& operator<<(std::ostream &out, const hg6_G1 &g);
    friend std::istream& operator>>(std::istream &in, hg6_G1 &g);

    static void batch_to_special_all_non_zeros(std::vector<hg6_G1> &vec);
};

template<mp_size_t m>
hg6_G1 operator*(const bigint<m> &lhs, const hg6_G1 &rhs)
{
    return scalar_mul<hg6_G1, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
hg6_G1 operator*(const Fp_model<m,modulus_p> &lhs, const hg6_G1 &rhs)
{
    return scalar_mul<hg6_G1, m>(rhs, lhs.as_bigint());
}

std::ostream& operator<<(std::ostream& out, const std::vector<hg6_G1> &v);
std::istream& operator>>(std::istream& in, std::vector<hg6_G1> &v);

} // libff

#endif // SW6_BIS_G1_HPP_
