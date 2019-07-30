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
- [ ] algebra_field_test
- [ ] algebra_group_test
- [ ] algebra_bilinearity_test

## TODO:
In `pendulum_init.hpp`:

* fill in `wnaf_window_table` for G1 and G2 

* fill in `fixed_base_exp_window_table` in G1 and G2

* double chek pairing parameters 

