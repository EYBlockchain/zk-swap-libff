/** @file
 *****************************************************************************

 Declaration of interfaces for the PENDULUM G1 group.

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef PENDULUM_G1_HPP_
#define PENDULUM_G1_HPP_

#include <vector>

#include <libff/algebra/curves/curve_utils.hpp>
#include <libff/algebra/curves/pendulum/pendulum_init.hpp>

namespace libff {

class pendulum_G1;
std::ostream& operator<<(std::ostream &, const pendulum_G1&);
std::istream& operator>>(std::istream &, pendulum_G1&);

class pendulum_G1 {
private:
    pendulum_Fq X_, Y_, Z_;
public:
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static pendulum_G1 G1_zero;
    static pendulum_G1 G1_one;
    static pendulum_Fq coeff_a;
    static pendulum_Fq coeff_b;

    typedef pendulum_Fq base_field;
    typedef pendulum_Fr scalar_field;

    // using projective coordinates
    pendulum_G1();
    pendulum_G1(const pendulum_Fq& X, const pendulum_Fq& Y) : X_(X), Y_(Y), Z_(base_field::one()) {}
    pendulum_G1(const pendulum_Fq& X, const pendulum_Fq& Y, const pendulum_Fq& Z) : X_(X), Y_(Y), Z_(Z) {}

    pendulum_Fq X() const { return X_; }
    pendulum_Fq Y() const { return Y_; }
    pendulum_Fq Z() const { return Z_; }

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const pendulum_G1 &other) const;
    bool operator!=(const pendulum_G1 &other) const;

    pendulum_G1 operator+(const pendulum_G1 &other) const;
    pendulum_G1 operator-() const;
    pendulum_G1 operator-(const pendulum_G1 &other) const;

    pendulum_G1 add(const pendulum_G1 &other) const;
    pendulum_G1 mixed_add(const pendulum_G1 &other) const;
    pendulum_G1 dbl() const;

    bool is_well_formed() const;

    static pendulum_G1 zero();
    static pendulum_G1 one();
    static pendulum_G1 random_element();

    static size_t size_in_bits() { return base_field::size_in_bits() + 1; }
    static bigint<base_field::num_limbs> base_field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    friend std::ostream& operator<<(std::ostream &out, const pendulum_G1 &g);
    friend std::istream& operator>>(std::istream &in, pendulum_G1 &g);

    static void batch_to_special_all_non_zeros(std::vector<pendulum_G1> &vec);
};

template<mp_size_t m>
pendulum_G1 operator*(const bigint<m> &lhs, const pendulum_G1 &rhs)
{
    return scalar_mul<pendulum_G1, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
pendulum_G1 operator*(const Fp_model<m,modulus_p> &lhs, const pendulum_G1 &rhs)
{
    return scalar_mul<pendulum_G1, m>(rhs, lhs.as_bigint());
}

std::ostream& operator<<(std::ostream& out, const std::vector<pendulum_G1> &v);
std::istream& operator>>(std::istream& in, std::vector<pendulum_G1> &v);

} // libff

#endif // PENDULUM_G1_HPP_
