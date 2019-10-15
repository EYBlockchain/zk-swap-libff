#include <libff/algebra/curves/bls12.hpp>

namespace libff {

namespace bls12 {


template<typename ppT, typename Fqk_type>
static inline Fqk_type exp_by_x(const Fqk_type &a)
{
    auto res = Fqk_type::one();
    for( int i = ppT::X_HIGHEST_BIT; i >= 0; i-- )
    {
        res = res.cyclotomic_squared();

        if( ppT::X & (1ul<<i) ) {
            res = res * a;
        }
    }

    if constexpr(ppT::X_IS_NEG) {
        return res.unitary_inverse();
    }

    return res;
}


template<typename ppT, typename Fqk_type>
static Fqk_type final_exponentiation(const Fqk_type &f)
{
    Fqk_type r = f.Frobenius_map(6) * f.inverse();
    r = r.Frobenius_map(2) * r;
    // Hard part of the final exponentation is below:
    // From https://eprint.iacr.org/2016/130.pdf, Table 1
    auto y0 = r.cyclotomic_squared().unitary_inverse();
    const auto y5 = exp_by_x<ppT>(r);
    const auto y3 = y0 * y5;
    y0 = exp_by_x<ppT>(y3);
    const auto y2 = exp_by_x<ppT>(y0);
    const auto y4 = exp_by_x<ppT>(y2) * y5.cyclotomic_squared();
    return (  (y0 * r).Frobenius_map(3)
            * (y5 * y2).Frobenius_map(2)
            * (y4 * r.unitary_inverse()).Frobenius_map(1)
            * exp_by_x<ppT>(y4)
            * y3.unitary_inverse()
            * r);
}


template<typename ppT>
static void miller_doubling_step(MillerTriple<ppT> &result, typename ppT::G2_type &r, const typename ppT::Fq_type &two_inv)
{
    // Formula for line function when working with homogeneous projective coordinates.
    const auto b = r.Y.squared();
    const auto c = r.Z.squared();
    const auto e = ppT::G2_type::coeff_b * c.multiply3();
    const auto f = e.multiply3();
    const auto h = (r.Y + r.Z).squared() - (b + c);

    if constexpr( ppT::TWIST_TYPE == TwistType::M ) {
        result.a = e - b;
        result.b = r.X.squared().multiply3();
        result.c = -h;
    }
    else {
        result.a = -h;
        result.b = r.X.squared().multiply3();
        result.c = e - b;
    }

    r.X = ((two_inv * r.X) * r.Y) * (b - f);
    r.Y = (two_inv * (b + f)).squared() - e.squared().multiply3();
    r.Z = b * h;
}


template<typename ppT>
static void miller_addition_step(MillerTriple<ppT> &result, typename ppT::G2_type &r, const typename ppT::G2_type &q)
{
    // Formula for line function when working with homogeneous projective coordinates.
    const auto theta = r.Y - (q.Y * r.Z);
    const auto lambda = r.X - (q.X * r.Z);
    const auto d = lambda.squared();
    const auto e = lambda * d;
    const auto g = r.X * d;
    const auto h = e + (r.Z * theta.squared()) - g.multiply2();
    const auto j = (theta * q.X) - (lambda * q.Y);

    r.X = lambda * h;
    r.Y = theta * (g - h) - (e * r.Y);
    r.Z = r.Z * e;

    if constexpr( ppT::TWIST_TYPE == TwistType::M ) {
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


template<typename ppT>
G2Prepared<ppT>::G2Prepared( const G2_type &g2 )
:
    infinity(g2.is_zero())
{
    if( ! infinity ) {
        _prepare(g2);
    }
}


template<typename ppT>
void G2Prepared<ppT>::_prepare(const G2_type &input_point)
{
    G2_type q = input_point;
    q.to_affine_coordinates();

    G2_type r = q;
    int coeff_idx = 0;

    // TODO: pre-compute two_inv... rather than every time it's prepared
    const auto two_inv = ppT::Fq_type::one().multiply2().inverse();

    // Skip the 1st bit
    for (int i = (ppT::X_HIGHEST_BIT-1); i >= 0; i--)
    {
        miller_doubling_step(this->coeffs[coeff_idx++], r, two_inv);

        if ( ppT::X & (1ul<<i) )
        {
            miller_addition_step(this->coeffs[coeff_idx++], r, q);
        }
    }
}


template<typename ppT>
struct PreparedPair {
    typename ppT::G1_type g1;
    const G2Prepared<ppT> g2;

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
    if constexpr( ppT::TWIST_TYPE == TwistType::M ) {
        f.multiply_by_c014(f, coeffs.a, g1.X * coeffs.b, g1.Y * coeffs.c);       
    }
    else {
        f.multiply_by_c034(f, g1.Y * coeffs.a, g1.X * coeffs.b, coeffs.c);       
    }
}


template<typename ppT>
static typename ppT::Fqk_type miller_loop( const std::vector<PreparedPair<ppT>> &pairs )
{
    int coeff_idx = 0;
    auto f = ppT::Fqk_type::one();

    for( int i = (ppT::X_HIGHEST_BIT-1); i >= 0; i-- )
    {
        f.square(f);

        for( const auto &p : pairs ) {
            ell<ppT>(f, p.g2.coeffs[coeff_idx], p.g1);
        }
        coeff_idx += 1;

        if ( ppT::X & (1ul<<i) )
        {
            for( const auto &p : pairs ) {
                ell<ppT>(f, p.g2.coeffs[coeff_idx], p.g1);
            }
            coeff_idx += 1;
        }
    }

    if constexpr( ppT::X_IS_NEG ) {
        f.conjugate(f);
    }

    return f;
}


template<typename ppT>
static inline typename ppT::Fqk_type miller_loop( const typename ppT::G1_type &P, const typename ppT::G2_type &Q )
{
    return miller_loop<ppT>({PreparedPair<ppT>(P, Q)});
}


// namespace bls12
}

// namespace libff
}
