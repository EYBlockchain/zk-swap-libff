# libff curves
## Implemented curves 
Originally in `libff`, the authors implemented 5 curves:
* BN128
This is a Barreto-Naehrig elliptic curve that targets 128-bit security level, with a group of points and a field characteristic of order 254 bits and an embedding degree `k=12`. The code is dynamically generated. 
* ALT_BN128
This is the same Barreto-Naehrig elliptic curve but without dynamically generated code (a bit slower that BN128). It is the one implemented in Ethereum Virtual Machine.
* MNT6 and MNT4
These two curves form a pairing-friendly cycle introduced in this paper. They are used for recursive SNARKs but do not target 128-bit security level as the field characteristic is of 298 bits.
* EDWARDS 
This is an Edwards curve built on top of ALT_BN128 to permit efficient implementation of pedersen hashes in `libsnark`. It goes by the name of Baby_Jubjub within Ethereum community.

## New curves
We implemented 7 new curves in `libff`:
- [x] BLS12_381
This a Barreto-Lynn-Scott curve implemented in ZCash that targets 128-bit security level, with a group of points of order 254 bits, a field characteristic of 381 bits and an embedding degree `k=12`.
- [x] MNT4753 and MNT6753
This is a similar pairing-friendly cycle as MNT4/MNT6 that targets 128-bit security using fields of 753 bits.
- [x] BLS12_377
This is the first ZEXE curve that enables bounded recursive SNARKs. It is a 128-bit BLS curve that is highly 2-adic w.r.t. both the subgroup order and field size and an embedding degree `k=12`.
- [x] SW6
This is the second ZEXE curve that enables bounded recursive SNARKs. It is constructed via Cocks-Pinch method with an embedding degree `k=6`.
- [x] SW6_BIS
This is a new curve built on top of BLS12_377 and is a bit more efficient than SW6.
- [x] PENDULUM
This is a new curve built on top of MNT6 and enable efficient implementation of ZKP aggregation.
- [ ] JUBJUB 
This is an Edwards curve built on top of BLS12_381 (implemented in ZCash) to enable efficient implementation of Pedersen hashes in `libsnark`.
