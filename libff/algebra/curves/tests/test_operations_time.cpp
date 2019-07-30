/* 
 *************************************************************************************
 * Curves comparison on the basis of time execution of ECADD, ECDBL and ECMUL in both 
 * G1 and G2 groups and of ECPAIR. Where:
 * ECADD is the addition of two points on the curve,
 * ECDBL is the doubling of a point on the curve,
 * ECMUL is the multiplication of a point by a random scalar in the finite fields, and 
 * ECPAIR is the pairing computation of a point of G1 with a point of G2.
 *
 * It is interesting to compare the following curves against each other: 
 * BN128 vs. ALT_BN128
 * ALT_BN128 vs. BLS12_381
 * ALT_BN128 vs. BLS12_377
 * BLS12_381 vs. BLS12_377
 * SW6 vs. SW6_BIS
 * MNT4 vs. MNT6
 * MNT4752 vs. MNT6753
 * MNT4752 vs. MNT4
 * MNT6752 vs. MNT6
 *************************************************************************************
 */

#include <time.h>

#include <libff/algebra/curves/edwards/edwards_pp.hpp>
#include <libff/common/profiling.hpp>
#ifdef CURVE_BN128
#include <libff/algebra/curves/bn128/bn128_pp.hpp>
#endif
#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>
#include <libff/algebra/curves/mnt/mnt4/mnt4_pp.hpp>
#include <libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp>

#include <libff/algebra/curves/bls12_377/bls12_377_pp.hpp>
#include <libff/algebra/curves/sw6/sw6_pp.hpp>
#include <libff/algebra/curves/sw6_bis/sw6_bis_pp.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_pp.hpp>
#include <libff/algebra/curves/mnt753/mnt6753/mnt6753_pp.hpp>
#include <libff/algebra/curves/mnt753/mnt4753/mnt4753_pp.hpp>

using namespace libff;

template<typename ppT>
void time_test()
{
    G1<ppT> G1_one = G1<ppT>::one();
    G2<ppT> G2_one = G2<ppT>::one();
    Fr<ppT> rand1 = Fr<ppT>::random_element();
    Fr<ppT> rand2 = Fr<ppT>::random_element();
    G1<ppT> P1 = rand1 * G1_one;
    G1<ppT> P2 = rand2 * G1_one;
    G2<ppT> Q1 = rand1 * G2_one;
    G2<ppT> Q2 = rand2 * G2_one;

    printf("ECADD, ECDBL and ECMUL in G1:\n");
    clock_t tStart = clock();
    G1<ppT> P_add = P1 + P2;
    printf("time taken for ECADD in G1 is %fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    tStart = clock();
    G1<ppT> P_mul = rand1 * P1;
    printf("time taken for ECMUL in G1 is %fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    tStart = clock();
    G1<ppT> P_dbl = P1.dbl();
    printf("time taken for ECDBL in G1 is %fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    printf("ECADD, ECDBL and ECMUL in G2:\n");
    tStart = clock();
    G2<ppT> Q_add = Q1 + Q2;
    printf("time taken for ECADD in G2 is %fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    tStart = clock();
    G2<ppT> Q_mul = rand1 * Q1;
    printf("time taken for ECMUL in G2 is %fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    tStart = clock();
    G2<ppT> Q_dbl = Q1.dbl();
    printf("time taken for ECDBL in G2 is %fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    printf("ECPAIR:\n");
    tStart = clock();
    GT<ppT> pair = ppT::reduced_pairing(P1,Q1);
    printf("time taken for ECPAIR is %fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
}


int main(void)
{
    printf("EDWARDS curve:\n");
    edwards_pp::init_public_params();
    time_test<edwards_pp>();
    printf("*************************\n\n");

    printf("MNT6 curve:\n");
    mnt6_pp::init_public_params();
    time_test<mnt6_pp>();
    printf("*************************\n\n");
    
    printf("MNT4 curve:\n");
    mnt4_pp::init_public_params();
    time_test<mnt4_pp>();
    printf("*************************\n\n");
    
    printf("MNT6753 curve:\n");
    mnt6753_pp::init_public_params();
    time_test<mnt6753_pp>();
    printf("*************************\n\n");
    
    printf("MNT4753 curve:\n");
    mnt4753_pp::init_public_params();
    time_test<mnt4753_pp>();
    printf("*************************\n\n");
    
    printf("ALT_BN128 curve:\n");
    alt_bn128_pp::init_public_params();
    time_test<alt_bn128_pp>();
    printf("*************************\n\n");
    
    printf("BLS12_381 curve:\n");
    bls12_381_pp::init_public_params();
    time_test<bls12_381_pp>();
    printf("*************************\n\n");
    
    printf("BLS12_377 curve:\n");
    bls12_377_pp::init_public_params();
    time_test<bls12_377_pp>();
    printf("*************************\n\n");
    
    printf("SW6 curve:\n");
    sw6_pp::init_public_params();
    time_test<sw6_pp>();
    printf("*************************\n\n");
    
    printf("SW6_BIS curve:\n");
    sw6_bis_pp::init_public_params();
    time_test<sw6_bis_pp>();
    printf("*************************\n\n");

#ifdef CURVE_BN128
    printf("BN128 curve:\n");
    bn128_pp::init_public_params();
    time_test<bn128_pp>();
    printf("*************************\n\n");
#endif 
}
