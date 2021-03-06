# Curve BW12_446
A new optimized Brezing-Weng pairing-friendly curve of embedding degree `k=12` that conservatively targets 128-bit security level. The subgroup is highly 2-adic (`s=37`).
*source*: Family 17 choice (b) in [this paper](https://eprint.iacr.org/2019/555.pdf).

### Progress:
- [x] BW12_446_g1.cpp
- [x] BW12_446_g1.hpp
- [x] BW12_446_g2.cpp
- [x] BW12_446_g2.hpp
- [x] BW12_446_init.cpp
- [x] BW12_446_init.hpp
- [x] BW12_446_pairing.cpp
- [x] BW12_446_pairing.hpp
- [x] BW12_446_pp.cpp
- [x] BW12_446_pp.hpp

## TODO:

## Tests:
- [x] algebra_field_test
- [x] algebra_group_test
- [x] algebra_bilinearity_test
