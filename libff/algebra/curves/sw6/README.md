# Curve SW6 
Zexe second pairing-friendly curve (embedding degree `k=6`), constructed over BLS12_377 via Cocks-Pinch method.

### Progress:
- [x] SW6_g1.cpp  
- [x] SW6_g1.hpp  
- [x] SW6_g2.cpp  
- [x] SW6_g2.hpp  
- [x] SW6_init.cpp  
- [x] SW6_init.hpp  
- [x] SW6_pairing.cpp  
- [x] SW6_pairing.hpp  
- [x] SW6_pp.cpp  
- [x] SW6_pp.hpp

## Tests:
- [x] algebra_field_test
- [x] algebra_group_test
- [x] algebra_bilinearity_test

## TODO:
In `sw6_init.hpp`:

* fill in `wnaf_window_table` for G1 and G2 

* fill in `fixed_base_exp_window_table` in G1 and G2

