# Curve MNT4753 
Coda MNT4753 pairing-friendly curve (embedding degree `k=4`). It forms a cycle with MNT6753.

### Progress:
- [x] mnt4753_g1.cpp  
- [x] mnt4753_g1.hpp  
- [x] mnt4753_g2.cpp  
- [x] mnt4753_g2.hpp  
- [x] mnt4753_init.cpp  
- [x] mnt4753_init.hpp  
- [x] mnt4753_pairing.cpp  
- [x] mnt4753_pairing.hpp  
- [x] mnt4753_pp.cpp  
- [x] mnt4753_pp.hpp

## Tests:
- [x] algebra_field_test
- [x] algebra_group_test
- [x] algebra_bilinearity_test

## TODO:
In `mnt4753_init.hpp`:

* recompute optimal `wnaf_window_table` for G1 and G2 

* recompute optimal `fixed_base_exp_window_table` in G1 and G2

