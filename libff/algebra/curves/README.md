# Original curves supported by libff
- BN128
- ALT_BN128
- EDWARDS
- MNT4
- MNT6

# New curves added to libff
- [x] BLS12_381: ZCash SNARK pairing-friendly elliptic curve.
- [x] BW12_446: A new optimized Brezing-Weng pairing-friendly curve of embedding degree `k=12` that conservatively targets 128-bit security level. The subgroup is highly 2-adic (`s=37`). *source*: Family 17 choice (b) in [this paper](https://eprint.iacr.org/2019/555.pdf).
- [x] BLS12_377: ZEXE inner curve (highly 2-adic w.r.t. both the subgroup order and the field characteristic)
- [x] SW6: ZEXE outter curve (constructed via Cocks-Pinch method)
- [x] HG6: A new curve to replace SW6, that is much faster.
- [x] PENDULUM: An outter curve for MNT6 (constructed via Cocks-Pinch method)
- [x] MNT4753: Coda cycle
- [x] MNT6753: Coda cycle
- [x] TOY_CURVE: A BN toy example
