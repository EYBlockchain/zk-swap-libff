/**
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#include <iostream>
#include <sstream>
#include <libff/algebra/curves/edwards/edwards_pp.hpp>
#include <libff/algebra/curves/mnt/mnt4/mnt4_pp.hpp>
#include <libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp>
#include <libff/common/profiling.hpp>
#ifdef CURVE_BN128
#include <libff/algebra/curves/bn128/bn128_pp.hpp>
#endif
#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>

#include <libff/algebra/curves/bls12_377/bls12_377_pp.hpp>
#include <libff/algebra/curves/sw6/sw6_pp.hpp>
#include <libff/algebra/curves/sw6_bis/sw6_bis_pp.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_pp.hpp>
#include <libff/algebra/curves/mnt753/mnt6753/mnt6753_pp.hpp>
#include <libff/algebra/curves/mnt753/mnt4753/mnt4753_pp.hpp>
#include <libff/algebra/curves/pendulum/pendulum_pp.hpp>
//#include <libff/algebra/curves/jubjub/jubjub_pp.hpp>

using namespace libff;



template<typename GroupT>
void test_mixed_add()
{
    GroupT base, el, result;

    base = GroupT::zero();
    el = GroupT::zero();
    el.to_special();
    result = base.mixed_add(el);
    assert(result == base + el);

    base = GroupT::zero();
    el = GroupT::random_element();
    el.to_special();
    result = base.mixed_add(el);
    assert(result == base + el);

    base = GroupT::random_element();
    el = GroupT::zero();
    el.to_special();
    result = base.mixed_add(el);
    assert(result == base + el);

    base = GroupT::random_element();
    el = GroupT::random_element();
    el.to_special();
    result = base.mixed_add(el);
    assert(result == base + el);

    base = GroupT::random_element();
    el = base;
    el.to_special();
    result = base.mixed_add(el);
    assert(result == base.dbl());
}

template<typename GroupT>
void test_group()
{
    bigint<1> rand1 = bigint<1>("76749407");
    bigint<1> rand2 = bigint<1>("44410867");
    bigint<1> randsum = bigint<1>("121160274");

    GroupT zero = GroupT::zero();
    assert(zero == zero);
    GroupT one = GroupT::one();
    assert(one == one);
    GroupT two = bigint<1>(2l) * GroupT::one();
    assert(two == two);
    GroupT five = bigint<1>(5l) * GroupT::one();

    GroupT three = bigint<1>(3l) * GroupT::one();
    GroupT four = bigint<1>(4l) * GroupT::one();
    GroupT fourtytwo = bigint<1>(42l) * GroupT::one();
    GroupT fifty = bigint<1>(50l) * GroupT::one();

    assert(two+five == three+four);
    assert(fifty-three == five+fourtytwo);

    GroupT a = GroupT::random_element();
    // std::cout << "random point 1:\n";
    // a.print();
    GroupT b = GroupT::random_element();
    // std::cout << "random point 2:\n";
    // b.print();

    assert(one != zero);
    assert(a != zero);
    assert(a != one);

    assert(b != zero);
    assert(b != one);

    assert(a.dbl() == a + a);
    assert(b.dbl() == b + b);
    assert(one + two == three);
    assert(two + one == three);
    assert(a + b == b + a);
    assert(a - a == zero);
    assert(a - b == a + (-b));
    assert(a - b == (-b) + a);

    // handle special cases
    assert(zero + (-a) == -a);
    assert(zero - a == -a);
    assert(a - zero == a);
    assert(a + zero == a);
    assert(zero + a == a);

    assert((a + b).dbl() == (a + b) + (b + a));
    assert(bigint<1>("2") * (a + b) == (a + b) + (b + a));

    assert((rand1 * a) + (rand2 * a) == (randsum * a));

    assert(GroupT::order() * a == zero);
    assert(GroupT::order() * one == zero);
    assert((GroupT::order() * a) - a != zero);
    assert((GroupT::order() * one) - one != zero);

    test_mixed_add<GroupT>();
}

template<typename GroupT>
void test_mul_by_q()
{
    GroupT a = GroupT::random_element();
    GroupT b = GroupT::base_field_char()*a;
    GroupT c = a.mul_by_q();
    // GroupT tmp = a.mul_by_q();
    // std::cout << "get X.Frobenius_map(1) and Y.Frobenius_map(1)\n";
    // tmp.print();
    // printf("random point a:\n");
    // a.print();
    // printf("q*a:\n");
    // b.print();
    // printf("a.mul_by_q():\n");
    // c.print();
    std::cout << "mul_by_q b";
    b.print();

    std::cout << "mul_by_q c";
    c.print();

    //assert(b == c);
}

template<typename GroupT>
void test_output()
{
    GroupT g = GroupT::one();
    // printf("one()\n");
    // g.print_coordinates();

    for (size_t i = 0; i < 100; ++i)
    {
        std::stringstream ss;
        ss << g;
        GroupT gg;
        // printf("here\n");
        ss >> gg;
        assert(g == gg);
        /* use a random point in next iteration */
        g = GroupT::random_element();
    }
}


template<typename ppT>
void test_known_points()
{
    using G1T = G1<ppT>;
    using G2T = G2<ppT>;
    using FpT = typename ppT::Fp_type;

    const FpT fourtytwo {"42"};

    const auto a = G1T::one();
    const auto b = G2T::one();
    auto c = fourtytwo * a;
    c.to_affine_coordinates();
    auto d = fourtytwo * b;
    d.to_affine_coordinates();
    auto e = a + a;
    e.to_affine_coordinates();
    auto f = b + b;
    f.to_affine_coordinates();

    const auto result = ppT::reduced_pairing(c, d);

    std::cout << "G1*42 = ";
    c.print();

    std::cout << "G2*42 = ";
    d.print();

    std::cout << "G1+G1 = ";
    e.print();

    std::cout << "G2+G2 = ";
    f.print();

    std::cout << "pairing(42,42) = ";
    result.print();
}


template<typename ppT>
void standard_tests(const char *name) {
    std::cout << name << std::endl;
    ppT::init_public_params();
    test_known_points<ppT>();
    test_group<G1<ppT> >();
    test_output<G1<ppT> >();
    test_group<G2<ppT> >();
    test_output<G2<ppT> >();
    //test_mul_by_q<G2<ppT> >();
    std::cout << std::endl << std::endl;
}


int main(void)
{
    standard_tests<bls12_381_pp>("bls12_381");
    standard_tests<bls12_377_pp>("bls12_377");
    standard_tests<edwards_pp>("edwards");
    standard_tests<mnt4_pp>("mnt4");
    standard_tests<mnt6_pp>("mnt6");
    standard_tests<alt_bn128_pp>("alt_bn128");
    standard_tests<mnt6753_pp>("mnt6753");
    standard_tests<mnt4753_pp>("mnt4753");
    /*
    standard_tests<sw6_pp>("sw6");
    standard_tests<sw6_bis_pp>("sw6_bis");
    standard_tests<pendulum_pp>("pendulum");
    */
    
    // new curve: JUBJUB 
    //standard_tests<jubjub_pp>("pendulum");

//#ifdef CURVE_BN128       // BN128 has fancy dependencies so it may be disabled
//    standard_tests<bn128_pp>("bn128");
//#endif
}
