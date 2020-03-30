#ifndef SW6_BIS_G2_HPP_
#define SW6_BIS_G2_HPP_

#include <vector>

#include <libff/algebra/curves/curve_utils.hpp>
#include <libff/algebra/curves/bw6_761/bw6_761_init.hpp>

namespace libff {

class bw6_761_G2;
std::ostream& operator<<(std::ostream &, const bw6_761_G2&);
std::istream& operator>>(std::istream &, bw6_761_G2&);

class bw6_761_G2 {
public:
    bw6_761_Fq X_, Y_, Z_;
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static bw6_761_G2 G2_zero;
    static bw6_761_G2 G2_one;
    static bw6_761_Fq twist;
    static bw6_761_Fq coeff_b;

    // endomorphism
    static bw6_761_Fq cube_root_of_unity; // w1
    static bw6_761_Fr eigen_value; // lambda1

    typedef bw6_761_Fq base_field;
    typedef bw6_761_Fq twist_field;
    typedef bw6_761_Fr scalar_field;

    // using projective coordinates
    bw6_761_G2();
    bw6_761_G2(const bw6_761_Fq& X, const bw6_761_Fq& Y, const bw6_761_Fq& Z) : X_(X), Y_(Y), Z_(Z) {}

    bw6_761_Fq X() const { return X_; }
    bw6_761_Fq Y() const { return Y_; }
    bw6_761_Fq Z() const { return Z_; }

    static bw6_761_Fq mul_by_b(const bw6_761_Fq &elt);

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const bw6_761_G2 &other) const;
    bool operator!=(const bw6_761_G2 &other) const;

    bw6_761_G2 operator+(const bw6_761_G2 &other) const;
    bw6_761_G2 operator-() const;
    bw6_761_G2 operator-(const bw6_761_G2 &other) const;

    bw6_761_G2 add(const bw6_761_G2 &other) const;
    bw6_761_G2 mixed_add(const bw6_761_G2 &other) const;
    bw6_761_G2 dbl() const;
    bw6_761_G2 mul_by_q() const;
    bw6_761_G2 endomorphism() const;
    bw6_761_G2 clear_cofactor() const;
    bw6_761_G2 mul_by_r() const;

    bool is_on_subgroup() const;
    bool is_well_formed() const;

    static bw6_761_G2 zero();
    static bw6_761_G2 one();
    static bw6_761_G2 random_element();

    static size_t size_in_bits() { return twist_field::size_in_bits() + 1; }
    static bigint<base_field::num_limbs> base_field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    friend std::ostream& operator<<(std::ostream &out, const bw6_761_G2 &g);
    friend std::istream& operator>>(std::istream &in, bw6_761_G2 &g);

    static void batch_to_special_all_non_zeros(std::vector<bw6_761_G2> &vec);
};

template<mp_size_t m>
bw6_761_G2 operator*(const bigint<m> &lhs, const bw6_761_G2 &rhs)
{
    return scalar_mul<bw6_761_G2, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
bw6_761_G2 operator*(const Fp_model<m,modulus_p> &lhs, const bw6_761_G2 &rhs)
{
    return scalar_mul<bw6_761_G2, m>(rhs, lhs.as_bigint());
}

} // libff

#endif // SW6_BIS_G2_HPP_
