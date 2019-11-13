/**
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#include <libff/algebra/curves/test_curve/test_curve_pp.hpp>
#include <libff/algebra/curves/bw12_446/bw12_446_pp.hpp>
#include <libff/algebra/curves/sw6/sw6_pp.hpp>
#include <libff/algebra/curves/sw6_bis/sw6_bis_pp.hpp>
#include <libff/algebra/curves/pendulum/pendulum_pp.hpp>
#include <libff/algebra/curves/edwards/edwards_pp.hpp>
#include <libff/common/profiling.hpp>
#ifdef CURVE_BN128
#include <libff/algebra/curves/bn128/bn128_pp.hpp>
#endif
#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>
#include <libff/algebra/curves/bls12_377/bls12_377_pp.hpp>
#include <libff/algebra/curves/toy_curve/toy_curve_pp.hpp>
#include <libff/algebra/curves/mnt/mnt4/mnt4_pp.hpp>
#include <libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp>
#include <libff/algebra/curves/mnt753/mnt6753/mnt6753_pp.hpp>
#include <libff/algebra/curves/mnt753/mnt4753/mnt4753_pp.hpp>

using namespace libff;

template<typename ppT>
void pairing_batching_test()
{
    // instantiate... 
    Fr<ppT> VKx_poly = (Fr<ppT>::random_element());
    Fr<ppT> VKy_poly = (Fr<ppT>::random_element());
    Fr<ppT> VKz_poly = (Fr<ppT>::random_element());

    Fr<ppT> A1_poly = Fr<ppT>::random_element();
    Fr<ppT> B1_poly = Fr<ppT>::random_element();
    Fr<ppT> C1_poly = (A1_poly * B1_poly - VKx_poly * VKy_poly) * VKz_poly.inverse();
    Fr<ppT> A2_poly = Fr<ppT>::random_element();
    Fr<ppT> B2_poly = Fr<ppT>::random_element();
    Fr<ppT> C2_poly = (A2_poly * B2_poly - VKx_poly * VKy_poly) * VKz_poly.inverse();

    // Proof
    G1<ppT> A1 = A1_poly * G1<ppT>::one(); 
    G2<ppT> B1 = B1_poly * G2<ppT>::one();   
    G1<ppT> C1 = C1_poly * G1<ppT>::one();
    G1<ppT> A2 = A2_poly * G1<ppT>::one(); 
    G2<ppT> B2 = B2_poly * G2<ppT>::one();   
    G1<ppT> C2 = C2_poly * G1<ppT>::one();
    
    // Verification key
    G1<ppT> VKx = VKx_poly * G1<ppT>::one(); 
    G2<ppT> VKy = VKy_poly * G2<ppT>::one();   
    G2<ppT> VKz = VKz_poly * G2<ppT>::one();

    // Verifier equation
    printf("\nverifier checks proof#1: e(A1,B1)=e(VK1,VK2)*e(C1,VK3)\n");
    assert(ppT::reduced_pairing(A1, B1) == ppT::reduced_pairing(VKx, VKy) * ppT::reduced_pairing(C1, VKz));
    printf("\nverifier checks proof#2: e(A2,B2)=e(VK1,VK2)*e(C2,VK3)\n");
    assert(ppT::reduced_pairing(A2, B2) == ppT::reduced_pairing(VKx, VKy) * ppT::reduced_pairing(C2, VKz));
    // TODO: need some randomness?
    printf("\nverifier checks a batching of proof#1 and proof#2: e(A1,B1)*e(A2,B2)=e(2*VK1,VK2)*e(C1+C2,VK3)\n");
    assert(ppT::reduced_pairing(A1, B1) * ppT::reduced_pairing(A2, B2) == ppT::reduced_pairing(Fr<ppT>("2")*VKx, VKy) * ppT::reduced_pairing(C1+C2, VKz));
}

template<typename ppT>
void pairing_test()
{
    GT<ppT> GT_one = GT<ppT>::one();

    printf("Running bilinearity tests:\n");
    G1<ppT> P = (Fr<ppT>::random_element()) * G1<ppT>::one();
    G2<ppT> Q = (Fr<ppT>::random_element()) * G2<ppT>::one();

    printf("P:\n");
    P.print();
    P.print_coordinates();
    printf("Q:\n");
    Q.print();
    Q.print_coordinates();
    printf("\n\n");

    Fr<ppT> s = Fr<ppT>::random_element();
    G1<ppT> sP = s * P;
    G2<ppT> sQ = s * Q;

    printf("Pairing bilinearity tests (three must match):\n");
    GT<ppT> ans1 = ppT::reduced_pairing(sP, Q);
    GT<ppT> ans2 = ppT::reduced_pairing(P, sQ);
    GT<ppT> ans3 = ppT::reduced_pairing(P, Q)^s;
    ans1.print();
    ans2.print();
    ans3.print();
    assert(ans1 == ans2);
    assert(ans2 == ans3);

    assert(ans1 != GT_one);
    assert((ans1^Fr<ppT>::field_char()) == GT_one);
    printf("\n\n");
}

template<typename ppT>
void double_miller_loop_test()
{
    const G1<ppT> P1 = (Fr<ppT>::random_element()) * G1<ppT>::one();
    const G1<ppT> P2 = (Fr<ppT>::random_element()) * G1<ppT>::one();
    const G2<ppT> Q1 = (Fr<ppT>::random_element()) * G2<ppT>::one();
    const G2<ppT> Q2 = (Fr<ppT>::random_element()) * G2<ppT>::one();

    const G1_precomp<ppT> prec_P1 = ppT::precompute_G1(P1);
    const G1_precomp<ppT> prec_P2 = ppT::precompute_G1(P2);
    const G2_precomp<ppT> prec_Q1 = ppT::precompute_G2(Q1);
    const G2_precomp<ppT> prec_Q2 = ppT::precompute_G2(Q2);

    const Fqk<ppT> ans_1 = ppT::miller_loop(prec_P1, prec_Q1);
    const Fqk<ppT> ans_2 = ppT::miller_loop(prec_P2, prec_Q2);
    const Fqk<ppT> ans_12 = ppT::double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
    assert(ans_1 * ans_2 == ans_12);
}

template<typename ppT>
void affine_pairing_test()
{
    GT<ppT> GT_one = GT<ppT>::one();

    printf("Running bilinearity tests:\n");
    G1<ppT> P = (Fr<ppT>::random_element()) * G1<ppT>::one();
    G2<ppT> Q = (Fr<ppT>::random_element()) * G2<ppT>::one();

    printf("P:\n");
    P.print();
    printf("Q:\n");
    Q.print();
    printf("\n\n");

    Fr<ppT> s = Fr<ppT>::random_element();
    G1<ppT> sP = s * P;
    G2<ppT> sQ = s * Q;

    printf("Pairing bilinearity tests (three must match):\n");
    GT<ppT> ans1 = ppT::affine_reduced_pairing(sP, Q);
    GT<ppT> ans2 = ppT::affine_reduced_pairing(P, sQ);
    GT<ppT> ans3 = ppT::affine_reduced_pairing(P, Q)^s;
    ans1.print();
    ans2.print();
    ans3.print();
    assert(ans1 == ans2);
    assert(ans2 == ans3);

    assert(ans1 != GT_one);
    assert((ans1^Fr<ppT>::field_char()) == GT_one);
    printf("\n\n");
}

int main(void)
{
    start_profiling();

    printf("test_curve:\n");
    test_curve_pp::init_public_params();
    pairing_test<test_curve_pp>();
    double_miller_loop_test<test_curve_pp>();

    /*
    printf("bw12_446:\n");
    bw12_446_pp::init_public_params();
    pairing_test<bw12_446_pp>();
    double_miller_loop_test<bw12_446_pp>();

    printf("pendulum:\n");
    pendulum_pp::init_public_params();
    pairing_test<pendulum_pp>();
    double_miller_loop_test<pendulum_pp>();

    printf("sw6:\n");
    sw6_pp::init_public_params();
    pairing_test<sw6_pp>();
    double_miller_loop_test<sw6_pp>();

    printf("sw6_bis:\n");
    sw6_bis_pp::init_public_params();
    pairing_test<sw6_bis_pp>();
    double_miller_loop_test<sw6_bis_pp>();

    printf("edwards:\n");
    edwards_pp::init_public_params();
    pairing_test<edwards_pp>();
    double_miller_loop_test<edwards_pp>();

    printf("mnt6:\n");
    mnt6_pp::init_public_params();
    pairing_test<mnt6_pp>();
    double_miller_loop_test<mnt6_pp>();
    affine_pairing_test<mnt6_pp>();

    printf("mnt4:\n");
    mnt4_pp::init_public_params();
    pairing_test<mnt4_pp>();
    double_miller_loop_test<mnt4_pp>();
    affine_pairing_test<mnt4_pp>();

    printf("mnt4753:\n");
    mnt4753_pp::init_public_params();
    pairing_test<mnt4753_pp>();
    double_miller_loop_test<mnt4753_pp>();
    affine_pairing_test<mnt4753_pp>();

    printf("mnt6753:\n");
    mnt6753_pp::init_public_params();
    pairing_test<mnt6753_pp>();
    double_miller_loop_test<mnt6753_pp>();
    affine_pairing_test<mnt6753_pp>();

    printf("alt_bn128:\n");
    alt_bn128_pp::init_public_params();
    pairing_test<alt_bn128_pp>();
    double_miller_loop_test<alt_bn128_pp>();

    printf("toy_curve:\n");
    toy_curve_pp::init_public_params();
    pairing_test<toy_curve_pp>();
    double_miller_loop_test<toy_curve_pp>();

    printf("bls12_377:\n");
    bls12_377_pp::init_public_params();
    pairing_test<bls12_377_pp>();
    double_miller_loop_test<bls12_377_pp>();
    pairing_batching_test<bls12_377_pp>();

#ifdef CURVE_BN128       // BN128 has fancy dependencies so it may be disabled
    bn128_pp::init_public_params();
    pairing_test<bn128_pp>();
    double_miller_loop_test<bn128_pp>();
#endif
    */
}
