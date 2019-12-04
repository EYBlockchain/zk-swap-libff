# Curve SW6_bis
An alternative curve to Zexe's `SW6`. A pairing-friendly curve (embedding degree `k=6`), constructed over `BLS12_377` via modified Cocks-Pinch method.

It has 3 enhacements over `SW6`:
- Smaller field size (`761` bits vs. `782` bits, or in Montgemery domain `768` bits vs. `832` bits)
- Smaller `G2` representation (coordinates over `Fq` vs. coordinates over `Fq^3`, thanks to the use of sextic twist vs. use of quadratic twist)
- Optimized final exponentiation (polynomial representation of `q` vs. N/A)
- use of GLV fast multiplication (polynomial representation of `q` vs. N/A)

### Progress:
- [x] sw6_bis_g1.cpp
- [x] sw6_bis_g1.hpp
- [x] sw6_bis_g2.cpp
- [x] sw6_bis_g2.hpp
- [x] sw6_bis_init.cpp
- [x] sw6_bis_init.hpp
- [x] sw6_bis_pairing.cpp
- [x] sw6_bis_pairing.hpp
- [x] sw6_bis_pp.cpp
- [x] sw6_bis_pp.hpp

## Tests:
- [x] algebra_field_test
- [x] algebra_group_test
- [x] algebra_bilinearity_test

## TODO:
In `sw6_bis_init.hpp`:

* [ ] fill in `wnaf_window_table` for G1 and G2

* [ ] fill in `fixed_base_exp_window_table` in G1 and G2

* [ ] implement and use the sextic twist

* [ ] implement optimized final exponentiation

* [ ] implement GLV fast multiplication
