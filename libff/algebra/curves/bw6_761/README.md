# Curve BW6_761
An alternative curve to Zexe's SW6 which is faster. A pairing-friendly curve (embedding degree `k=6`), constructed over BLS12-377 via *modified* Brezing-Weng method.
Beside having a field size of 761-bit and G2 elements on Fq, parameters are expressed as polynomials allowing very efficient pairing computation.

### Progress:
- [x] bw6_761_g1.cpp
- [x] bw6_761_g1.hpp
- [x] bw6_761_g2.cpp
- [x] bw6_761_g2.hpp
- [x] bw6_761_init.cpp
- [x] bw6_761_init.hpp
- [x] bw6_761_pairing.cpp
- [x] bw6_761_pairing.hpp
- [x] bw6_761_pp.cpp
- [x] bw6_761_pp.hpp

## Tests:
- [x] algebra_field_test
- [x] algebra_group_test
- [x] algebra_bilinearity_test

## TODO:
- [x] change curve coefficient to `b=-1`
- [x] use a sextic twist to bw6_761 so that elements in `G2` are in `Fq`
- [x] optimal ate pairing (Alg.5)
- [ ] optimized Alg.5
  - [x] NAF
  - [ ] f_{u^2-u-1,Q}
- [x] optimize FE hard part
- [x] optimize `mul_by_024` in Fq6
- [x] implement the D-twist
  - [x] G2 on `y^2=x^3-1/2` over `Fq`
  - [x] pairing with `mul_by_045`
- [ ] factor a square in mul pairing (double_pairing)

- [ ] Fast multiplication with GLV endomorphisms (patent until 09/20)
- [ ] Faster hash into curve with endomorphisms (patent until 09/20)
- [ ] Faster points checks with endomorphisms (patent until 09/20)

- [ ] recompute optimal `wnaf_window_table` for G1 and G2
- [ ] recompute optimal `fixed_base_exp_window_table` in G1 and G2
