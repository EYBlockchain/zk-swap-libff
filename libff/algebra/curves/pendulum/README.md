# Pendulum curve
A pairing-friendly curve (embedding degree `k=6`), constructed over MNT6 via Cocks-Pinch method.

### Progress:
- [x] pendulum_g1.cpp
- [x] pendulum_g1.hpp
- [x] pendulum_g2.cpp
- [x] pendulum_g2.hpp
- [x] pendulum_init.cpp
- [x] pendulum_init.hpp
- [x] pendulum_pairing.cpp
- [x] pendulum_pairing.hpp
- [x] pendulum_pp.cpp
- [x] pendulum_pp.hpp

## Tests:
- [x] algebra_field_test
- [x] algebra_group_test
- [x] algebra_bilinearity_test

## TODO:
In `pendulum_init.hpp`:

* recompute optimal `wnaf_window_table` for G1 and G2

* recompute optimal `fixed_base_exp_window_table` in G1 and G2
