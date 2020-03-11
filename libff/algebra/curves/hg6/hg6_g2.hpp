#ifndef SW6_BIS_G2_HPP_
#define SW6_BIS_G2_HPP_

#include <vector>

#include <libff/algebra/curves/curve_utils.hpp>
#include <libff/algebra/curves/hg6/hg6_init.hpp>

namespace libff {

class hg6_G2;
std::ostream& operator<<(std::ostream &, const hg6_G2&);
std::istream& operator>>(std::istream &, hg6_G2&);

class hg6_G2 {
public:
    hg6_Fq X_, Y_, Z_;
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static hg6_G2 G2_zero;
    static hg6_G2 G2_one;
    static hg6_Fq twist;
    static hg6_Fq coeff_a;
    static hg6_Fq coeff_b;

    typedef hg6_Fq base_field;
    typedef hg6_Fq twist_field;
    typedef hg6_Fr scalar_field;

    // using projective coordinates
    hg6_G2();
    hg6_G2(const hg6_Fq& X, const hg6_Fq& Y, const hg6_Fq& Z) : X_(X), Y_(Y), Z_(Z) {}

    hg6_Fq X() const { return X_; }
    hg6_Fq Y() const { return Y_; }
    hg6_Fq Z() const { return Z_; }

    static hg6_Fq mul_by_a(const hg6_Fq &elt);
    static hg6_Fq mul_by_b(const hg6_Fq &elt);

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const hg6_G2 &other) const;
    bool operator!=(const hg6_G2 &other) const;

    hg6_G2 operator+(const hg6_G2 &other) const;
    hg6_G2 operator-() const;
    hg6_G2 operator-(const hg6_G2 &other) const;

    hg6_G2 add(const hg6_G2 &other) const;
    hg6_G2 mixed_add(const hg6_G2 &other) const;
    hg6_G2 dbl() const;
    hg6_G2 mul_by_q() const;

    bool is_well_formed() const;

    static hg6_G2 zero();
    static hg6_G2 one();
    static hg6_G2 random_element();

    static size_t size_in_bits() { return twist_field::size_in_bits() + 1; }
    static bigint<base_field::num_limbs> base_field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    friend std::ostream& operator<<(std::ostream &out, const hg6_G2 &g);
    friend std::istream& operator>>(std::istream &in, hg6_G2 &g);

    static void batch_to_special_all_non_zeros(std::vector<hg6_G2> &vec);
};

template<mp_size_t m>
hg6_G2 operator*(const bigint<m> &lhs, const hg6_G2 &rhs)
{
    return scalar_mul<hg6_G2, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
hg6_G2 operator*(const Fp_model<m,modulus_p> &lhs, const hg6_G2 &rhs)
{
    return scalar_mul<hg6_G2, m>(rhs, lhs.as_bigint());
}

} // libff

#endif // SW6_BIS_G2_HPP_
