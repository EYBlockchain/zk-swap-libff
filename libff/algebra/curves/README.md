# Original curves supported by libff
- BN128
- ALT_BN128
- EDWARDS
- MNT4
- MNT6

# New curves added to libff
- [x] BLS12_377: ZEXE inner curve (highly 2-adic w.r.t. both the subgroup order and the field characteristic)
- [x] SW6: ZEXE outter curve (constructed via Cocks-Pinch method)
- [x] SW6_BIS: A new curve to replace SW6
- [x] PENDULUM: An outter curve for MNT6 (constructed via Cocks-Pinch method)
- [x] MNT4753: Coda cycle
- [x] MNT6753: Coda cycle
- [x] TOY_CURVE: A BN toy example

## TODO:
- [ ] BLS12_381: ZCash curve
- [ ] JUBJUB: Edwards curve on top of BLS12_381
- [ ] BABY_JUBJUB: Edwards curve on top of ALT_BN128
