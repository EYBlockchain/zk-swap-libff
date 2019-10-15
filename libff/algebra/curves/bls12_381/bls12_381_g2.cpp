#include <libff/algebra/curves/bls12_381/bls12_381_g2.hpp>

namespace libff {

bls12_381_G2 bls12_381_G2::_one;
bls12_381_G2 bls12_381_G2::_zero;
std::vector<size_t> bls12_381_G2::wnaf_window_table;
std::vector<size_t> bls12_381_G2::fixed_base_exp_window_table;

}