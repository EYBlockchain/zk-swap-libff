# Curve BLS12_381
Zcash pairing-friendly curve (embedding degree `k=12`).

### Progress:
<<<<<<< HEAD
- [x] BLS12_377_g1.cpp
- [x] BLS12_377_g1.hpp
- [x] BLS12_377_g2.cpp
- [x] BLS12_377_g2.hpp
- [x] BLS12_377_init.cpp
- [x] BLS12_377_init.hpp
- [ ] BLS12_377_pairing.cpp
- [x] BLS12_377_pairing.hpp
- [x] BLS12_377_pp.cpp
=======
- [x] BLS12_377_g1.cpp
- [x] BLS12_377_g1.hpp
- [x] BLS12_377_g2.cpp
- [x] BLS12_377_g2.hpp
- [x] BLS12_377_init.cpp
- [x] BLS12_377_init.hpp
- [x] BLS12_377_pairing.cpp
- [x] BLS12_377_pairing.hpp
- [x] BLS12_377_pp.cpp
>>>>>>> origin/new_curves
- [x] BLS12_377_pp.hpp

## TODO:
In `BLS12_381_init.hpp`:

* fill in `wnaf_window_table` for G1 and G2

* fill in `fixed_base_exp_window_table` in G1 and G2

* implement pairing for M-type curves (need to add `mul_by_014()` in `libff/algebra/fields/fp12_2over3over2.tcc` and change `ell_0`, `ell_VV` and `ell_VW` in Miller loop doubling/addition)

## Tests:
- [x] algebra_field_test
- [x] algebra_group_test
- [ ] algebra_bilinearity_test

