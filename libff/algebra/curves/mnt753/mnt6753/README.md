# Curve MNT6753 
Coda MNT6753 pairing-friendly curve (embedding degree `k=6`). I form a cycle with MNT4753.

### Progress:
- [x] mnt6753_g1.cpp  
- [x] mnt6753_g1.hpp  
- [x] mnt6753_g2.cpp  
- [x] mnt6753_g2.hpp  
- [x] mnt6753_init.cpp  
- [x] mnt6753_init.hpp  
- [x] mnt6753_pairing.cpp  
- [x] mnt6753_pairing.hpp  
- [x] mnt6753_pp.cpp  
- [x] mnt6753_pp.hpp

## Tests:
- [x] algebra_field_test
- [x] algebra_group_test
- [x] algebra_bilinearity_test

## TODO:
In `mnt6753_init.hpp`:

* fill in `wnaf_window_table` for G1 and G2 

* fill in `fixed_base_exp_window_table` in G1 and G2

