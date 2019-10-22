# Curve TOY_CURVE (BN40)
A toy example of a Barreto-Naehrig pairing friendly-elliptic curve (embedding degree `k=12`) with a field size of 81-bits.

### Progress:
- [x] TOY_CURVE_g1.cpp
- [x] TOY_CURVE_g1.hpp
- [x] TOY_CURVE_g2.cpp
- [x] TOY_CURVE_g2.hpp
- [x] TOY_CURVE_init.cpp
- [x] TOY_CURVE_init.hpp
- [x] TOY_CURVE_pairing.cpp
- [x] TOY_CURVE_pairing.hpp
- [x] TOY_CURVE_pp.cpp
- [x] TOY_CURVE_pp.hpp

## TODO:
In `TOY_CURVE_init.hpp`:

* recompute optimal `wnaf_window_table` for G1 and G2

* recompute optimal `fixed_base_exp_window_table` in G1 and G2

## Tests:
- [x] algebra_field_test
- [x] algebra_group_test
- [x] algebra_bilinearity_test

