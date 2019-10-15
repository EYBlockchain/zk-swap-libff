#include <libff/algebra/curves/bls12_381/bls12_381_g1.hpp>

namespace libff {

bls12_381_G1 bls12_381_G1::_one;
bls12_381_G1 bls12_381_G1::_zero;
std::vector<size_t> bls12_381_G1::wnaf_window_table;
std::vector<size_t> bls12_381_G1::fixed_base_exp_window_table;

}