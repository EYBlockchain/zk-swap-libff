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
#include <libff/algebra/curves/pendulum/pendulum_pp.hpp>

using namespace libff;

template<typename ppT>
void pairing_test()
{
    GT<ppT> GT_one = GT<ppT>::one();
    GT<ppT> GT_random = GT<ppT>::random_element();
    printf("GT_one:\n");
    GT_one.print();
    // printf("GT_random:\n");
    // GT_random.print();

    G1<ppT> G1_one = G1<ppT>::one();
    printf("G1_one:\n");
    G1_one.print();
    G2<ppT> G2_one = G2<ppT>::one();
    printf("G2_one:\n");
    G2_one.print();

    printf("Running bilinearity tests:\n");
    //G1<ppT> P = (Fr<ppT>::random_element()) * G1_one;
    G1<ppT> P = Fr<ppT>("2") * G1<ppT>::one();
    //G2<ppT> Q = (Fr<ppT>::random_element()) * G2_one;
    G2<ppT> Q = Fr<ppT>("3") * G2<ppT>::one();

    P.to_affine_coordinates();
    Q.to_affine_coordinates();

    printf("P:\n");
    P.print();
    P.print_coordinates();
    printf("Q:\n");
    Q.print();
    Q.print_coordinates();
    printf("\n\n");

    //Fr<ppT> s = Fr<ppT>::random_element();
    Fr<ppT> s = Fr<ppT>("2");
    G1<ppT> sP = s * P;
    G2<ppT> sQ = s * Q;

    printf("Pairing bilinearity tests (three must match):\n");
    GT<ppT> ans1 = ppT::reduced_pairing(sP, Q);
    GT<ppT> ans2 = ppT::reduced_pairing(P, sQ);
    GT<ppT> ans3 = ppT::reduced_pairing(P, Q)^s;
    printf("e(sP,Q):\n");
    ans1.print();
    printf("e(P,sQ):\n");
    ans2.print();
    printf("e(P,Q)^s:\n");
    ans3.print();
    assert(ans1 == ans2);
    assert(ans2 == ans3);
    printf("e(sP,Q)=e(P,sQ)=e(sP,Q)^s checked!\n");
    assert(ans1 != GT_one);
    printf("non-degenracy e(sP,Q) != 1 checked!\n");
    assert((ans1^Fr<ppT>::field_char()) == GT_one);
    printf("Fq12 order e(sP,Q)^r = 1 checked!\n");
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

    /*
    // EDWARDS
    edwards_pp::init_public_params();
    pairing_test<edwards_pp>();
    double_miller_loop_test<edwards_pp>();

    // MNT6
    mnt6_pp::init_public_params();
    pairing_test<mnt6_pp>();
    double_miller_loop_test<mnt6_pp>();
    affine_pairing_test<mnt6_pp>();

    // MNT4
    mnt4_pp::init_public_params();
    pairing_test<mnt4_pp>();
    double_miller_loop_test<mnt4_pp>();
    affine_pairing_test<mnt4_pp>();

    // ALT_BN128 
    alt_bn128_pp::init_public_params();
    pairing_test<alt_bn128_pp>();
    double_miller_loop_test<alt_bn128_pp>();
    */
    
    // new curve: BLS12_381
    bls12_381_pp::init_public_params();
    pairing_test<bls12_381_pp>();
    //double_miller_loop_test<bls12_381_pp>();

    /*
    // new curve: BLS12_377
    bls12_377_pp::init_public_params();
    pairing_test<bls12_377_pp>();
    double_miller_loop_test<bls12_377_pp>();
    
    // new curve: SW6
    sw6_pp::init_public_params();
    pairing_test<sw6_pp>();
    double_miller_loop_test<sw6_pp>();

    
    // new curve: MNT6753 
    mnt6753_pp::init_public_params();
    pairing_test<mnt6753_pp>();
    double_miller_loop_test<mnt6753_pp>();
    
    // new curve: MNT4753 
    mnt4753_pp::init_public_params();
    pairing_test<mnt4753_pp>();
    double_miller_loop_test<mnt4753_pp>();
    
    // new curve: SW6_BIS
    sw6_bis_pp::init_public_params();
    pairing_test<sw6_bis_pp>();
    double_miller_loop_test<sw6_bis_pp>();
    
    // new curve: PENDULUM 
    pendulum_pp::init_public_params();
    pairing_test<pendulum_pp>();
    double_miller_loop_test<pendulum_pp>();
    */


#ifdef CURVE_BN128       // BN128 has fancy dependencies so it may be disabled
    // bn128_pp::init_public_params();
    // pairing_test<bn128_pp>();
    // double_miller_loop_test<bn128_pp>();
#endif
}
