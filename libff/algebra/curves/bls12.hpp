#pragma once

#include <cstdint>
#include <limits>

namespace libff {

namespace bls12 {


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


// namespace bls12
}

// namespace libff
}