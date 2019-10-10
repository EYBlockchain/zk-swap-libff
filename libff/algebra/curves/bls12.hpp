#pragma once

#include <cstdint>
#include <limits>
#include <array>

namespace libff {

namespace bls12 {


template<typename ppT>
struct MillerTriple {
    typename ppT::Fqe_type a;
    typename ppT::Fqe_type b;
    typename ppT::Fqe_type c;
};

enum TwistType {
    M,  // Multiply at twist
    D   // Divide at twist
};

// For X_NUM_ONES
template <uint64_t x>
struct CountOnes {
    static constexpr unsigned int n = (x & 1) + CountOnes<(x>>1)>::n;
};
template <>
struct CountOnes<0> {
    static constexpr unsigned int n = 0;
};

// for X_HIGHEST_BIT
template<uint64_t N, uint64_t B=std::numeric_limits<uint64_t>::digits-1>
struct FindMSB {
    static constexpr unsigned int MSB = (N&(1ul<<B)) ? B: FindMSB<N,B-1>::MSB;
};
template<uint64_t N>
struct FindMSB<N,0> {
    static constexpr unsigned int MSB = (N==0) ? -1 : 0;
};


/* NOTE: This structure is approximately 20 KiB in size, so use with caution. */
template<typename ppT>
struct G2Prepared
{
    static constexpr unsigned int num_coeffs = ppT::X_HIGHEST_BIT + ppT::X_NUM_ONES - 1;
    using G2_type = typename ppT::G2_type;

    std::array<MillerTriple<ppT>, num_coeffs> coeffs;
    const bool infinity;

    bool is_zero() const;

    G2Prepared( const G2_type &g2 );

    void _prepare(const G2_type &input_point);
};



// namespace bls12
}

// namespace libff
}