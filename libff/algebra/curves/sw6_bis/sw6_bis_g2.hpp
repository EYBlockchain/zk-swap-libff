#ifndef SW6_BIS_G2_HPP_
#define SW6_BIS_G2_HPP_

#include <vector>

#include <libff/algebra/curves/curve_utils.hpp>
#include <libff/algebra/curves/sw6_bis/sw6_bis_init.hpp>

namespace libff {

class sw6_bis_G2;
std::ostream& operator<<(std::ostream &, const sw6_bis_G2&);
std::istream& operator>>(std::istream &, sw6_bis_G2&);

class sw6_bis_G2 {
private:
    sw6_bis_Fq3 X_, Y_, Z_;
public:
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static sw6_bis_G2 G2_zero;
    static sw6_bis_G2 G2_one;
    static sw6_bis_Fq3 twist;
    static sw6_bis_Fq3 coeff_a;
    static sw6_bis_Fq3 coeff_b;

    typedef sw6_bis_Fq base_field;
    typedef sw6_bis_Fq3 twist_field;
    typedef sw6_bis_Fr scalar_field;

    // using projective coordinates
    sw6_bis_G2();
    sw6_bis_G2(const sw6_bis_Fq3& X, const sw6_bis_Fq3& Y, const sw6_bis_Fq3& Z) : X_(X), Y_(Y), Z_(Z) {}

    sw6_bis_Fq3 X() const { return X_; }
    sw6_bis_Fq3 Y() const { return Y_; }
    sw6_bis_Fq3 Z() const { return Z_; }

    static sw6_bis_Fq3 mul_by_a(const sw6_bis_Fq3 &elt);
    static sw6_bis_Fq3 mul_by_b(const sw6_bis_Fq3 &elt);

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const sw6_bis_G2 &other) const;
    bool operator!=(const sw6_bis_G2 &other) const;

    sw6_bis_G2 operator+(const sw6_bis_G2 &other) const;
    sw6_bis_G2 operator-() const;
    sw6_bis_G2 operator-(const sw6_bis_G2 &other) const;

    sw6_bis_G2 add(const sw6_bis_G2 &other) const;
    sw6_bis_G2 mixed_add(const sw6_bis_G2 &other) const;
    sw6_bis_G2 dbl() const;
    sw6_bis_G2 mul_by_q() const;

    bool is_well_formed() const;

    static sw6_bis_G2 zero();
    static sw6_bis_G2 one();
    static sw6_bis_G2 random_element();

    static size_t size_in_bits() { return twist_field::size_in_bits() + 1; }
    static bigint<base_field::num_limbs> base_field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    friend std::ostream& operator<<(std::ostream &out, const sw6_bis_G2 &g);
    friend std::istream& operator>>(std::istream &in, sw6_bis_G2 &g);

    static void batch_to_special_all_non_zeros(std::vector<sw6_bis_G2> &vec);
};

template<mp_size_t m>
sw6_bis_G2 operator*(const bigint<m> &lhs, const sw6_bis_G2 &rhs)
{
    return scalar_mul<sw6_bis_G2, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p>
sw6_bis_G2 operator*(const Fp_model<m,modulus_p> &lhs, const sw6_bis_G2 &rhs)
{
    return scalar_mul<sw6_bis_G2, m>(rhs, lhs.as_bigint());
}

} // libff

#endif // SW6_BIS_G2_HPP_
