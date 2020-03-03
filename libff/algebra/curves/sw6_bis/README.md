# Curve sw6_bis
An alternative curve to Zexe's SW6 (faster). A pairing-friendly curve (embedding degree `k=6`), constructed over BLS12_377 via *modified* Brezing-Weng method.

### Progress:
- [x] sw6_bis_g1.cpp
- [x] sw6_bis_g1.hpp
- [x] sw6_bis_g2.cpp
- [x] sw6_bis_g2.hpp
- [x] sw6_bis_init.cpp
- [x] sw6_bis_init.hpp
- [ ] sw6_bis_pairing.cpp(*)
- [x] sw6_bis_pairing.hpp
- [x] sw6_bis_pp.cpp
- [x] sw6_bis_pp.hpp

## Tests:
- [x] algebra_field_test
- [x] algebra_group_test
- [ ] algebra_bilinearity_test(*)

## TODO:
- [x] change curve coefficient to `b=-1`
- [ ] use a sextic twist to sw6_bis so that elements in `G2` will be in `Fq`
- [x] optimal ate pairing (Alg.5)
- [ ] optimized Alg.5
  - [x] NAF
  - [ ] f_{u^2-u-1,Q}
- [x] Alg.5 in affine coordinates
- [x] optimize FE hard part

- [ ] Fast multiplication with GLV endomorphisms (patent until 09/20)
- [ ] Faster hash into curve with endomorphisms (patent until 09/20)
- [ ] Faster points checks with endomorphisms (patent until 09/20)

- [ ] recompute optimal `wnaf_window_table` for G1 and G2
- [ ] recompute optimal `fixed_base_exp_window_table` in G1 and G2

## Note
(*) pairing works with the non-optimized miller loop versions
