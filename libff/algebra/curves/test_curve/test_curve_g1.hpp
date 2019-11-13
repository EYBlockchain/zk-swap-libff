#ifndef TEST_CURVE_G1_HPP_
#define TEST_CURVE_G1_HPP_
#include <vector>

#include <libff/algebra/curves/test_curve/test_curve_init.hpp>
#include <libff/algebra/curves/curve_utils.hpp>

namespace libff {

class test_curve_G1;
std::ostream& operator<<(std::ostream &, const test_curve_G1&);
std::istream& operator>>(std::istream &, test_curve_G1&);

class test_curve_G1 {
public:
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static test_curve_G1 G1_zero;
    static test_curve_G1 G1_one;

    typedef test_curve_Fq base_field;
    typedef test_curve_Fr scalar_field;

    test_curve_Fq X, Y, Z;

    // using Jacobian coordinates
    test_curve_G1();
    test_curve_G1(const test_curve_Fq& X, const test_curve_Fq& Y, const test_curve_Fq& Z) : X(X), Y(Y), Z(Z) {};

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const test_curve_G1 &other) const;
    bool operator!=(const test_curve_G1 &other) const;

    test_curve_G1 operator+(const test_curve_G1 &other) const;
    test_curve_G1 operator-() const;
    test_curve_G1 operator-(const test_curve_G1 &other) const;

    test_curve_G1 add(const test_curve_G1 &other) const;
    test_curve_G1 mixed_add(const test_curve_G1 &other) const;
    test_curve_G1 dbl() const;

    bool is_well_formed() const;

    static test_curve_G1 zero();
    static test_curve_G1 one();
    static test_curve_G1 random_element();

    static size_t size_in_bits() { return base_field::size_in_bits() + 1; }
    static bigint<base_field::num_limbs> base_field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    friend std::ostream& operator<<(std::ostream &out, const test_curve_G1 &g);
    friend std::istream& operator>>(std::istream &in, test_curve_G1 &g);

    static void batch_to_special_all_non_zeros(std::vector<test_curve_G1> &vec);
};

template<mp_size_t m>
test_curve_G1 operator*(const bigint<m> &lhs, const test_curve_G1 &rhs)
{
    return scalar_mul<test_curve_G1, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
test_curve_G1 operator*(const Fp_model<m,modulus_p> &lhs, const test_curve_G1 &rhs)
{
    return scalar_mul<test_curve_G1, m>(rhs, lhs.as_bigint());
}

std::ostream& operator<<(std::ostream& out, const std::vector<test_curve_G1> &v);
std::istream& operator>>(std::istream& in, std::vector<test_curve_G1> &v);

} // libff
#endif // TEST_CURVE_G1_HPP_
