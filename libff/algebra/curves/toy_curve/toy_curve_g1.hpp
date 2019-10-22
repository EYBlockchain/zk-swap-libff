#ifndef TOY_CURVE_G1_HPP_
#define TOY_CURVE_G1_HPP_
#include <vector>

#include <libff/algebra/curves/toy_curve/toy_curve_init.hpp>
#include <libff/algebra/curves/curve_utils.hpp>

namespace libff {

class toy_curve_G1;
std::ostream& operator<<(std::ostream &, const toy_curve_G1&);
std::istream& operator>>(std::istream &, toy_curve_G1&);

class toy_curve_G1 {
public:
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static toy_curve_G1 G1_zero;
    static toy_curve_G1 G1_one;

    typedef toy_curve_Fq base_field;
    typedef toy_curve_Fr scalar_field;

    toy_curve_Fq X, Y, Z;

    // using Jacobian coordinates
    toy_curve_G1();
    toy_curve_G1(const toy_curve_Fq& X, const toy_curve_Fq& Y, const toy_curve_Fq& Z) : X(X), Y(Y), Z(Z) {};

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const toy_curve_G1 &other) const;
    bool operator!=(const toy_curve_G1 &other) const;

    toy_curve_G1 operator+(const toy_curve_G1 &other) const;
    toy_curve_G1 operator-() const;
    toy_curve_G1 operator-(const toy_curve_G1 &other) const;

    toy_curve_G1 add(const toy_curve_G1 &other) const;
    toy_curve_G1 mixed_add(const toy_curve_G1 &other) const;
    toy_curve_G1 dbl() const;

    bool is_well_formed() const;

    static toy_curve_G1 zero();
    static toy_curve_G1 one();
    static toy_curve_G1 random_element();

    static size_t size_in_bits() { return base_field::size_in_bits() + 1; }
    static bigint<base_field::num_limbs> base_field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    friend std::ostream& operator<<(std::ostream &out, const toy_curve_G1 &g);
    friend std::istream& operator>>(std::istream &in, toy_curve_G1 &g);

    static void batch_to_special_all_non_zeros(std::vector<toy_curve_G1> &vec);
};

template<mp_size_t m>
toy_curve_G1 operator*(const bigint<m> &lhs, const toy_curve_G1 &rhs)
{
    return scalar_mul<toy_curve_G1, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
toy_curve_G1 operator*(const Fp_model<m,modulus_p> &lhs, const toy_curve_G1 &rhs)
{
    return scalar_mul<toy_curve_G1, m>(rhs, lhs.as_bigint());
}

std::ostream& operator<<(std::ostream& out, const std::vector<toy_curve_G1> &v);
std::istream& operator>>(std::istream& in, std::vector<toy_curve_G1> &v);

} // libff
#endif // TOY_CURVE_G1_HPP_
