#pragma once
#include <vector>

#include <libff/algebra/curves/sw.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_pp.hpp>

namespace libff {

class bls12_381_G2;

using bls12_381_G2_base = sw::SWJacobianPoint<bls12_381_pp, bls12_381_pp::Fqe_type, bls12_381_pp::Fp_type, bls12_381_twist_coeff_b, bls12_381_G2>;

class bls12_381_G2 : public bls12_381_G2_base
{
	using bls12_381_G2_base::SWJacobianPoint;

public:
	// Pass-thru Copy constructor from base type
	bls12_381_G2( const bls12_381_G2_base &base )
	: bls12_381_G2_base(base) {}

	int sign_bit() const {
		return Y.c0.as_bigint().data[0] & 1;
	}

	// Define storage in our own class, otherwise it gets messy
	static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static bls12_381_G2 _zero;
    static bls12_381_G2 _one;
};

// namespace libff
}
