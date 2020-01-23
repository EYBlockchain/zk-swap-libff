# Curve sw6_bis
An alternative curve to Zexe's SW6 (faster). A pairing-friendly curve (embedding degree `k=6`), constructed over BLS12_377 via *modified* Cocks-Pinch method.

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
- [ ] optimize FE hard part
- [ ] Fast multiplication with GLV endomorphisms
- [ ] Faster hash into curve with endomorphisms
- [ ] Faster points checks with endomorphisms
- [ ] recompute optimal `wnaf_window_table` for G1 and G2
- [ ] recompute optimal `fixed_base_exp_window_table` in G1 and G2

## Note
(*) pairing works with naive final exponentiation:
  - [x] first chunk: `(q^3-1)(q+1)`
  - [x] last chunk: `(q^2-q+1)/r = w0+q*w1`
but we want to implement the last chunk as `R0(u)+p*R1(u)` with `R0` et `R1` polynomials and `p` a Frobenius in `GF(p^6)`
