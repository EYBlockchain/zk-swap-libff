# Curve sw6_bis
An alternative curve to Zexe's SW6 (slightly faster). A pairing-friendly curve (embedding degree `k=6`), constructed over BLS12_377 via Cocks-Pinch method.

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

* recompute optimal `wnaf_window_table` for G1 and G2

* recompute optimal `fixed_base_exp_window_table` in G1 and G2

* use a sextic twist to sw6_bis so that elements in `G2` will be in `Fq`

* optimize FE (wip with Aurore)
