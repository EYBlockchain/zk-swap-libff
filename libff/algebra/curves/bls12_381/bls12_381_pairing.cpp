#include <cassert>

#include <libff/algebra/curves/bls12_381/bls12_381_g1.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_g2.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_init.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_pp.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_pairing.hpp>
#include <libff/common/profiling.hpp>


using libff::BlsTwistType;


template<typename ppT>
struct MillerTriple {
    typename ppT::Fqe_type a;
    typename ppT::Fqe_type b;
    typename ppT::Fqe_type c;
};


template<typename ppT, typename Fqk_type=typename ppT::Fqk_type>
static void exp_by_x(Fqk_type &output, const Fqk_type &a)
{
    auto res = Fqk_type::one();
    bool found_one = false;
    for( int i = 63; i >= 0; i-- )
    {
        if( found_one ) {
            res = res.squared();
        }

        if( ppT::X & (1ul<<i) ) {
            found_one = true;
            res = res * a;
        }
    }

    if (ppT::X_IS_NEG) {
        res.conjugate(res);
    }

    output = res;
}


template<typename ppT, typename Fqk_type=typename ppT::Fqk_type>
static void final_exponentiation(Fqk_type &result, const Fqk_type &f)
{
    // f1 = r.conjugate() = f^(p^6)
    Fqk_type f1;                        // let mut f1 = *f;
    f1.frobenius_map(f, 6);             // f1.frobenius_map(6);

    Fqk_type f2;
    f2.inverse(f);

    // f2 = f^(-1);
    // r = f^(p^6 - 1)
    Fqk_type r = f1 * f2;               // let mut r = f1 * &f2;

    // f2 = f^(p^6 - 1)
    f2 = r;                             // f2 = r;
    // r = f^((p^6 - 1)(p^2))
    r.frobenius_map(r, 2);              // r.frobenius_map(2);

    // r = f^((p^6 - 1)(p^2) + (p^6 - 1))
    // r = f^((p^6 - 1)(p^2 + 1))
    r = r * f2;                         // r *= &f2;

    // Hard part of the final exponentation is below:
    // From https://eprint.iacr.org/2016/130.pdf, Table 1
                                        // let mut y0 = r.cyclotomic_square();
                                        // y0.conjugate();
    Fqk_type y0 = r.squared().unitary_inverse();

    Fqk_type y5;
    exp_by_x<ppT>(y5, r);               // let mut y5 = Self::exp_by_x(r);
    Fqk_type y1 = y5.squared();         // let mut y1 = y5.cyclotomic_square();
    Fqk_type y3 = y0 * y5;              // let mut y3 = y0 * &y5;
    exp_by_x<ppT>(y0, y3);              // y0 = Self::exp_by_x(y3);
    Fqk_type y2;
    exp_by_x<ppT>(y2, y0);              // let y2 = Self::exp_by_x(y0);
    Fqk_type y4;
    exp_by_x<ppT>(y4, y2);              // let mut y4 = Self::exp_by_x(y2);
    y4 = y4 * y1;                       // y4 *= &y1;

    exp_by_x<ppT>(y1, y4);              // y1 = Self::exp_by_x(y4);
    y3.conjugate(y3);                   // y3.conjugate();
    y1 = y1 * y3;                       // y1 *= &y3;
    y1 = y1 * r;                        // y1 *= &r;
    y3 = r;                             // y3 = r;
    y3.conjugate(y3);                   // y3.conjugate();

    y0 = y0 * r;                        // y0 *= &r;
    y0.frobenius_map(y0, 3);            // y0.frobenius_map(3);
    y4 = y4 * y3;                       // y4 *= &y3;
    y4.frobenius_map(y4, 1);            // y4.frobenius_map(1);
    y5 = y5 * y2;                       // y5 *= &y2;
    y5.frobenius_map(y5, 2);            // y5.frobenius_map(2);
    y5 = y5 * y0;                       // y5 *= &y0;
    y5 = y5 * y4;                       // y5 *= &y4;
    y5 = y5 * y1;                       // y5 *= &y1;

    result = y5;
}


template<typename ppT>
static void miller_doubling_step(MillerTriple<ppT> &result, typename ppT::G2_type &r, const typename ppT::Fq_type &two_inv)
{
    // Formula for line function when working with homogeneous projective coordinates.

    const auto a = (two_inv * r.X) * r.Y;           // let mut a = r.x * &r.y;
                                                    // a.mul_by_fp(two_inv);    
    const auto b = r.Y.squared();                   // let b = r.y.square();
    const auto c = r.Z.squared();                   // let c = r.z.square();
    const auto e = ppT::TWIST_COEFF_B * (c.multiply2() + c);   // let e = B::G2Parameters::COEFF_B * &(c.double() + &c);
    const auto f = e.multiply2() + e;               // let f = e.double() + &e;
    const auto g = two_inv * (b + f);               // let mut g = b + &f;
                                                    // g.mul_by_fp(two_inv);
    const auto h = (r.Y + r.Z).squared() - (b + c); // let h = (r.y + &r.z).square() - &(b + &c);
    const auto i = e - b;                           // let i = e - &b;
    const auto j = r.X.squared();                   // let j = r.x.square();
    const auto e_square = e.squared();              // let e_square = e.square();

    r.X = a * (b - f);                              // r.x = a * &(b - &f);
    r.Y = g.squared() - (e_square.multiply2() + e_square); //r.y = g.square() - &(e_square.double() + &e_square);
    r.Z = b * h;                                    // r.z = b * &h;

    // TODO: compile-time, remove branch
    if( ppT::TWIST_TYPE == BlsTwistType::M ) {
        result.a = i;
        result.b = j.multiply2() + j;
        result.c = -h;
    }
    else {
        result.a = -h;
        result.b = j.multiply2() + j;
        result.c = i;
    }
}


template<typename ppT>
static void miller_addition_step(MillerTriple<ppT> &result, typename ppT::G2_type &r, const typename ppT::G2_type &q)
{
    // Formula for line function when working with homogeneous projective coordinates.

    const auto theta = r.Y - (q.Y * r.Z);           // let theta = r.y - &(q.y * &r.z);
    const auto lambda = r.X - (q.X * r.Z);          // let lambda = r.x - &(q.x * &r.z);
    const auto c = theta.squared();                 // let c = theta.square();
    const auto d = lambda.squared();                // let d = lambda.square();
    const auto e = lambda * d;                      // let e = lambda * &d;
    const auto f = r.Z * c;                         // let f = r.z * &c;
    const auto g = r.X * d;                         // let g = r.x * &d;
    const auto h = e + f - g.multiply2();           // let h = e + &f - &g.double();
    r.X = lambda * h;                               // r.x = lambda * &h;
    r.Y = theta * (g - h) - (e * r.Y);              // r.y = theta * &(g - &h) - &(e * &r.y);
    r.Z = r.Z * e;                                  // r.z *= &e;
    const auto j = (theta * q.X) - (lambda * q.Y);  // let j = theta * &q.x - &(lambda * &q.y);

    // TODO: compile-time, remove branch
    if( ppT::TWIST_TYPE == BlsTwistType::M ) {
        result.a = j;
        result.b = -theta;
        result.c = lambda;
    }
    else {
        result.a = lambda;
        result.b = -theta;
        result.c = j;
    }
}


/* NOTE: This structure is approximately 20 KiB in size, so use with caution. */
template<typename ppT>
struct G2Prepared
{
    static constexpr unsigned int num_coeffs = ppT::X_HIGHEST_BIT + ppT::X_NUM_BITS - 1;
    using G2_type = typename ppT::G2_type;

    std::array<MillerTriple<ppT>, num_coeffs> coeffs;
    bool infinity;

    bool is_zero() const {
        return this->infinity;
    }

    G2Prepared( const G2_type &g2 )
    {
        prepare(g2);
    }

    void prepare(const G2_type &input_point)
    {
        if( (this->infinity = input_point.is_zero()) ) {
            return;
        }

        G2_type q = input_point;
        q.to_affine_coordinates();

        G2_type r = q;
        int coeff_idx = 0;

        // TODO: pre-compute two_inv... rather than every time it's prepared
        const auto two_inv = ppT::Fq_type::one().multiply2().inverse();

        // Skip the 1st bit
        for (int i = 62; i >= 0; i--)
        {
            miller_doubling_step(this->coeffs[coeff_idx++], r, two_inv);
     
            if ( ppT::X & (1ul<<i) )
            {
                miller_addition_step(this->coeffs[coeff_idx++], r, q);
            }
        }
    }
};


template<typename ppT>
struct PreparedPair {
    typename ppT::G1_type g1;
    G2Prepared<ppT> g2;

    PreparedPair( const decltype(g1) &a, const decltype(g2) &b )
    :
        g1(a), g2(b)
    {
        g1.to_affine_coordinates();
    }
};


// Twisting isomorphism from E to E'
template<typename ppT>
static inline void ell(typename ppT::Fqk_type &f, const MillerTriple<ppT> &coeffs, const typename ppT::G1_type &g1)
{
    if( ppT::TWIST_TYPE == BlsTwistType::M ) {
        // Multiply type twist
        f.multiply_by_c014(f, coeffs.a, g1.X * coeffs.b, g1.Y * coeffs.c);       
    }
    else {
        // Divide type twist
        f.multiply_by_c034(f, g1.Y * coeffs.a, g1.X * coeffs.b, coeffs.c);       
    }
}


template<typename ppT>
static typename ppT::Fqk_type miller_loop( const typename ppT::G1_type &P, const typename ppT::G2_type &Q )
{
    const PreparedPair<ppT> pair(P, Q);
    int coeff_idx = 0;
    auto f = ppT::Fqk_type::one();

    for( int i = 62; i >= 0; i-- )
    {        
        f.square(f);

        ell<ppT>(f, pair.g2.coeffs[coeff_idx++], pair.g1);

        if ( ppT::X & (1ul<<i) ) {
            ell<ppT>(f, pair.g2.coeffs[coeff_idx++], pair.g1);
        }
    }

    if( ppT::X_IS_NEG ) {
        f.conjugate(f);
    }

    return f;
}


namespace libff {


bls12_381_Fq12 bls12_381_pairing(const bls12_381_G1& P, const bls12_381_G2 &Q)
{
    return miller_loop<libff::bls12_381_pp>(P, Q);
}


bls12_381_Fq12 bls12_381_reduced_pairing(const bls12_381_G1& P, const bls12_381_G2 &Q)
{
    bls12_381_Fq12 result = bls12_381_pairing(P, Q);
    final_exponentiation<libff::bls12_381_pp>(result, result);
    return result;
}


} // libff
