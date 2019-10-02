#include <cassert>

#include <libff/algebra/curves/bls12_381/bls12_381_g1.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_g2.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_init.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_pairing.hpp>
#include <libff/common/profiling.hpp>


typedef libff::bls12_381_G1 G1;
typedef libff::bls12_381_G2 G2;
typedef libff::bls12_381_Fq2 Fq2;
typedef libff::bls12_381_Fq12 Fq12;


struct MillerTriple {
    Fq2 a;
    Fq2 b;
    Fq2 c;
};


/*
 * There is no benefit to being constant time in the exponent, since
 * "bls_x" is assumed to be publicly known.
 */
template <unsigned int right_shift, bool square_at_end>
static inline void exp_by_x_restrict(Fq12& result, const Fq12& a) {
    result = Fq12::one();

    for (int i = (int)libff::bls12_381_x_highest_set_bit; i != (right_shift - 1); i--) {
        std::cout << "..." << i << " " << right_shift << " " << (right_shift - 1) << std::endl;
        result.square(result);
        if (libff::bls12_381_x.test_bit(i)) {
            result.multiply(result, a);
        }
    }

    if (square_at_end) {
        result.square(result);
    }

    if (libff::bls12_381_x_is_negative) {
        result.conjugate(result);
    }
}


static void final_exponentiation(Fq12& result, const Fq12& a)
{
    Fq12 f1;
    f1.conjugate(a);

    Fq12 f2;
    f2.inverse(a);
    Fq12 r;
    r.multiply(f1, f2);
    f2.copy(r);
    r.frobenius_map(r, 2);
    r.multiply(r, f2);

    Fq12 y0;
    y0.square(r);

    Fq12& y1 = result;
    exp_by_x_restrict<0, false>(y1, y0);

    Fq12 y2;
    exp_by_x_restrict<1, false>(y2, y1);

    Fq12 y3;
    y3.conjugate(r);
    y1.multiply(y1, y3);
    y1.conjugate(y1);
    y1.multiply(y1, y2);
    exp_by_x_restrict<1, true>(y2, y1);
    exp_by_x_restrict<1, true>(y3, y2);
    y1.conjugate(y1);
    y3.multiply(y3, y1);
    y1.conjugate(y1);
    y1.frobenius_map(y1, 3);
    y2.frobenius_map(y2, 2);
    y1.multiply(y1, y2);
    exp_by_x_restrict<1, true>(y2, y3);
    y2.multiply(y2, y0);
    y2.multiply(y2, r);
    y1.multiply(y1, y2);
    y2.frobenius_map(y3, 1);
    y1.multiply(y1, y2);
}


static void miller_doubling_step(MillerTriple& result, G2& r)
{
    Fq2& tmp0 = result.a;
    tmp0.square(r.X);

    Fq2 tmp1;
    tmp1.square(r.Y);

    Fq2 tmp2;
    tmp2.square(tmp1);

    Fq2& tmp3 = result.b;
    tmp3.add(tmp1, r.X);
    tmp3.square(tmp3);
    tmp3.subtract(tmp3, tmp0);
    tmp3.subtract(tmp3, tmp2);
    tmp3.multiply2(tmp3);

    Fq2 tmp4;
    tmp4.multiply2(tmp0);
    tmp4.add(tmp4, tmp0);

    Fq2& tmp6 = result.c;
    tmp6.add(r.X, tmp4);

    Fq2 tmp5;
    tmp5.square(tmp4);

    Fq2 zsquared;
    zsquared.square(r.Z);

    r.X.subtract(tmp5, tmp3);
    r.X.subtract(r.X, tmp3);

    r.Z.add(r.Z, r.Y);
    r.Z.square(r.Z);
    r.Z.subtract(r.Z, tmp1);
    r.Z.subtract(r.Z, zsquared);

    r.Y.subtract(tmp3, r.X);
    r.Y.multiply(r.Y, tmp4);

    tmp2.multiply2(tmp2);
    tmp2.multiply2(tmp2);
    tmp2.multiply2(tmp2);

    r.Y.subtract(r.Y, tmp2);

    // Calculate result.b
    tmp3.multiply(tmp4, zsquared);
    tmp3.multiply2(tmp3);
    tmp3.negate(tmp3);

    // Calculate result.c
    tmp6.square(tmp6);
    tmp6.subtract(tmp6, tmp0);
    tmp6.subtract(tmp6, tmp5);

    tmp1.multiply2(tmp1);
    tmp1.multiply2(tmp1);

    tmp6.subtract(tmp6, tmp1);

    // Calculate result.a
    tmp0.multiply(r.Z, zsquared);
    tmp0.multiply2(tmp0);
}


static void miller_addition_step(MillerTriple& result, G2& r, const G2& g2) {
    Fq2 zsquared;
    zsquared.square(r.Z);

    Fq2 ysquared;
    ysquared.square(g2.Y);

    Fq2 t0;
    t0.multiply(zsquared, g2.X);

    Fq2& t1 = result.b;
    t1.add(g2.Y, r.Z);
    t1.square(t1);
    t1.subtract(t1, ysquared);
    t1.subtract(t1, zsquared);
    t1.multiply(t1, zsquared);

    Fq2 t2;
    t2.subtract(t0, r.X);

    Fq2 t3;
    t3.square(t2);

    Fq2 t4;
    t4.multiply2(t3);
    t4.multiply2(t4);

    Fq2 t5;
    t5.multiply(t4, t2);

    Fq2 t6;
    t6.subtract(t1, r.Y);
    t6.subtract(t6, r.Y);

    Fq2& t9 = result.c;
    t9.multiply(t6, g2.X);

    Fq2 t7;
    t7.multiply(t4, r.X);

    r.X.square(t6);
    r.X.subtract(r.X, t5);
    r.X.subtract(r.X, t7);
    r.X.subtract(r.X, t7);

    r.Z.add(r.Z, t2);
    r.Z.square(r.Z);
    r.Z.subtract(r.Z, zsquared);
    r.Z.subtract(r.Z, t3);

    Fq2& t10 = result.a;
    t10.add(g2.Y, r.Z);

    Fq2 t8;
    t8.subtract(t7, r.X);
    t8.multiply(t8, t6);

    t0.multiply(r.Y, t5);
    t0.multiply2(t0);

    r.Y.subtract(t8, t0);

    t10.square(t10);
    t10.subtract(t10, ysquared);

    Fq2 ztsquared;
    ztsquared.square(r.Z);

    t10.subtract(t10, ztsquared);

    t9.multiply2(t9);
    t9.subtract(t9, t10);

    t10.multiply2(r.Z);

    t6.negate(t6);

    t1.multiply2(t6);
}


/*
 * NOTE: This structure is approximately 20 KiB in size, so use with
 * caution.
 */
struct G2Prepared {
    static constexpr unsigned int num_coeffs = libff::bls12_381_x_highest_set_bit + libff::bls12_381_x_num_set_bits - 1;

    MillerTriple coeffs[num_coeffs];
    bool infinity;

    bool is_zero() const {
        return this->infinity;
    }

    G2Prepared( const G2 &g2 )
    {
        prepare(g2);
    }

    void prepare(const G2& g2) {
        G2 g2_copy = g2;
        g2_copy.to_affine_coordinates();

        G2 r = g2_copy;
        int coeff_idx = 0;

        /* Skips the least significant bit and most significant set bit. */
        for (unsigned int i = libff::bls12_381_x_highest_set_bit - 1; i != 0; i--) {
            miller_doubling_step(this->coeffs[coeff_idx++], r);
            if (libff::bls12_381_x.test_bit(i)) {
                miller_addition_step(this->coeffs[coeff_idx++], r, g2_copy);
            }
        }

        miller_doubling_step(this->coeffs[coeff_idx], r);
        this->infinity = g2_copy.is_zero();
    }
};


struct AffinePair {
    G1 g1;
    G2 g2;
    G2 r;   // private

    AffinePair( const G1 &a, const G2 &b )
    :
        g1(a),
        g2(b)
    {
        g1.to_affine_coordinates();
        g2.to_affine_coordinates();
    }
};


struct PreparedPair {
    G1 g1;
    G2Prepared g2;

    PreparedPair( const G1 &a, const G2Prepared &b )
    :
        g1(a),
        g2(b)
    {
        g1.to_affine_coordinates();
    }

    size_t coeff_idx;   // private
};



static void ell(Fq12& f, const MillerTriple& coeffs, const G1& g1)
{
    Fq2 c0;
    Fq2 c1;

    c0.c0 = coeffs.a.c0 * g1.Y;
    c0.c1 = coeffs.a.c1 * g1.Y;

    c1.c0 = coeffs.b.c0 * g1.X;
    c1.c1 = coeffs.b.c1 * g1.X;

    f.multiply_by_c014(f, coeffs.c, c1, c0);
}



static void miller_loop(Fq12& result, AffinePair* affine_pairs, size_t num_affine_pairs, PreparedPair* prepared_pairs, size_t num_prepared_pairs) {
    MillerTriple coeffs;

    result.copy(Fq12::one());

    for (size_t j = 0; j != num_affine_pairs; j++) {
        AffinePair& pair = affine_pairs[j];
        pair.r = pair.g2;
        pair.r.to_affine_coordinates();
    }
    for (size_t j = 0; j != num_prepared_pairs; j++) {
        PreparedPair& pair = prepared_pairs[j];
        pair.coeff_idx = 0;
    }

    /* Skips the least significant bit and most significant set bit. */
    for (unsigned int i = libff::bls12_381_x_highest_set_bit - 1; i != 0; i--) {
        for (size_t j = 0; j != num_affine_pairs; j++) {
            AffinePair& pair = affine_pairs[j];
            if (!pair.g1.is_zero() && !pair.g2.is_zero()) {
                miller_doubling_step(coeffs, pair.r);
                ell(result, coeffs, pair.g1);
            }
        }
        for (size_t j = 0; j != num_prepared_pairs; j++) {
            PreparedPair& pair = prepared_pairs[j];
            if (!pair.g1.is_zero() && !pair.g2.is_zero()) {
                ell(result, pair.g2.coeffs[pair.coeff_idx++], pair.g1);
            }
        }

        if (libff::bls12_381_x.test_bit(i)) {
            for (size_t j = 0; j != num_affine_pairs; j++) {
                AffinePair& pair = affine_pairs[j];
                if (!pair.g1.is_zero() && !pair.g2.is_zero()) {
                    miller_addition_step(coeffs, pair.r, pair.g2);
                    ell(result, coeffs, pair.g1);
                }
            }
            for (size_t j = 0; j != num_prepared_pairs; j++) {
                PreparedPair& pair = prepared_pairs[j];
                if (!pair.g1.is_zero() && !pair.g2.is_zero()) {
                    ell(result, pair.g2.coeffs[pair.coeff_idx++], pair.g1);
                }
            }
        }

        result.square(result);
    }

    for (size_t j = 0; j != num_affine_pairs; j++) {
        AffinePair& pair = affine_pairs[j];
        if (!pair.g1.is_zero() && !pair.g2.is_zero()) {
            miller_doubling_step(coeffs, pair.r);
            ell(result, coeffs, pair.g1);
        }
    }
    for (size_t j = 0; j != num_prepared_pairs; j++) {
        PreparedPair& pair = prepared_pairs[j];
        if (!pair.g1.is_zero() && !pair.g2.is_zero()) {
            ell(result, pair.g2.coeffs[pair.coeff_idx++], pair.g1);
        }
    }

    if (libff::bls12_381_x_is_negative) {
        result.conjugate(result);
    }
}


void miller_loop(Fq12 &result, const G1 &g1, const G2 &g2) {
    AffinePair pair(g1, g2);
    miller_loop(result, &pair, 1, nullptr, 0);
}


void miller_loop(Fq12 &result, const G1 &g1, const G2Prepared &g2) {
    PreparedPair pair(g1, g2);
    miller_loop(result, nullptr, 0, &pair, 1);
}


namespace libff {


bls12_381_Fq12 bls12_381_pairing(const bls12_381_G1& P, const bls12_381_G2 &Q)
{
    bls12_381_Fq12 result;
    miller_loop(result, P, Q);
    return result;
}


bls12_381_Fq12 bls12_381_reduced_pairing(const bls12_381_G1& P, const bls12_381_G2 &Q)
{
    bls12_381_Fq12 result = bls12_381_pairing(P, Q);
    final_exponentiation(result, result);
    return result;
}


} // libff
