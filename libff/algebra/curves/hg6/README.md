# Curve HG6
An alternative curve to Zexe's SW6 which is faster. A pairing-friendly curve (embedding degree `k=6`), constructed over BLS12-377 via *modified* Brezing-Weng method.
Beside having a field size of 761-bit and G2 elements on Fq, parameters are expressed as polynomials allowing very efficient pairing computation.

### Progress:
- [x] hg6_g1.cpp
- [x] hg6_g1.hpp
- [x] hg6_g2.cpp
- [x] hg6_g2.hpp
- [x] hg6_init.cpp
- [x] hg6_init.hpp
- [x] hg6_pairing.cpp
- [x] hg6_pairing.hpp
- [x] hg6_pp.cpp
- [x] hg6_pp.hpp

## Tests:
- [x] algebra_field_test
- [x] algebra_group_test
- [x] algebra_bilinearity_test

## TODO:
- [x] change curve coefficient to `b=-1`
- [x] use a sextic twist to hg6 so that elements in `G2` are be in `Fq`
- [x] optimal ate pairing (Alg.5)
- [ ] optimized Alg.5
  - [x] NAF
  - [ ] f_{u^2-u-1,Q}
- [x] Alg.5 in affine coordinates
- [x] optimize FE hard part
- [x] optimize mul_by_024 in Fq6
- [ ] factor a square in mul pairing (double_pairing)

- [ ] Fast multiplication with GLV endomorphisms (patent until 09/20)
- [ ] Faster hash into curve with endomorphisms (patent until 09/20)
- [ ] Faster points checks with endomorphisms (patent until 09/20)

- [ ] recompute optimal `wnaf_window_table` for G1 and G2
- [ ] recompute optimal `fixed_base_exp_window_table` in G1 and G2
