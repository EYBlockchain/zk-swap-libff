#include <libff/algebra/curves/toy_curve/toy_curve_g1.hpp>
#include <libff/algebra/curves/toy_curve/toy_curve_g2.hpp>
#include <libff/algebra/curves/toy_curve/toy_curve_init.hpp>

namespace libff {

bigint<toy_curve_r_limbs> toy_curve_modulus_r;
bigint<toy_curve_q_limbs> toy_curve_modulus_q;

toy_curve_Fq toy_curve_coeff_b;
toy_curve_Fq2 toy_curve_twist;
toy_curve_Fq2 toy_curve_twist_coeff_b;
toy_curve_Fq toy_curve_twist_mul_by_b_c0;
toy_curve_Fq toy_curve_twist_mul_by_b_c1;
toy_curve_Fq2 toy_curve_twist_mul_by_q_X;
toy_curve_Fq2 toy_curve_twist_mul_by_q_Y;

bigint<toy_curve_q_limbs> toy_curve_ate_loop_count;
bool toy_curve_ate_is_loop_count_neg;
bigint<12*toy_curve_q_limbs> toy_curve_final_exponent;
bigint<toy_curve_q_limbs> toy_curve_final_exponent_z;
bool toy_curve_final_exponent_is_z_neg;

void init_toy_curve_params()
{
    typedef bigint<toy_curve_r_limbs> bigint_r;
    typedef bigint<toy_curve_q_limbs> bigint_q;

    assert(sizeof(mp_limb_t) == 8 || sizeof(mp_limb_t) == 4); // Montgomery assumes this

    /* parameters for scalar field Fr */

    toy_curve_modulus_r = bigint_r("1211208045514572150807217");
    assert(toy_curve_Fr::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        toy_curve_Fr::Rsquared = bigint_r("79671369965451331016344");
        toy_curve_Fr::Rcubed = bigint_r("530347193141001599194170");
        toy_curve_Fr::inv = 0xfdf45c86db4f11af;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        toy_curve_Fr::Rsquared = bigint_r("1100134687786334797816816");
        toy_curve_Fr::Rcubed = bigint_r("87100873514965819468348");
        toy_curve_Fr::inv = 0xdb4f11af;
    }
    toy_curve_Fr::num_bits = 81;
    toy_curve_Fr::euler = bigint_r("605604022757286075403608");
    toy_curve_Fr::s = 4;
    toy_curve_Fr::t = bigint_r("75700502844660759425451");
    toy_curve_Fr::t_minus_1_over_2 = bigint_r("37850251422330379712725");
    toy_curve_Fr::multiplicative_generator = toy_curve_Fr("5");
    toy_curve_Fr::root_of_unity = toy_curve_Fr("253350081684880917638143");
    toy_curve_Fr::nqr = toy_curve_Fr("5");
    toy_curve_Fr::nqr_to_t = toy_curve_Fr("253350081684880917638143");

    /* parameters for base field Fq */

    toy_curve_modulus_q = bigint_q("1211208045515672698496983");
    assert(toy_curve_Fq::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        toy_curve_Fq::Rsquared = bigint_q("650507409816592346205806");
        toy_curve_Fq::Rcubed = bigint_q("517440680038978689559956");
        toy_curve_Fq::inv = 0x3beb29520f794419;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        toy_curve_Fq::Rsquared = bigint_q("706344792498391613488943");
        toy_curve_Fq::Rcubed = bigint_q("242532178488217303219420");
        toy_curve_Fq::inv = 0xf794419;
    }
    toy_curve_Fq::num_bits = 81;
    toy_curve_Fq::euler = bigint_q("605604022757836349248491");
    toy_curve_Fq::s = 1;
    toy_curve_Fq::t = bigint_q("605604022757836349248491");
    toy_curve_Fq::t_minus_1_over_2 = bigint_q("302802011378918174624245");
    toy_curve_Fq::multiplicative_generator = toy_curve_Fq("3");
    toy_curve_Fq::root_of_unity = toy_curve_Fq("1211208045515672698496982");
    toy_curve_Fq::nqr = toy_curve_Fq("3");
    toy_curve_Fq::nqr_to_t = toy_curve_Fq("1211208045515672698496982");

    /* parameters for twist field Fq2 */
    toy_curve_Fq2::euler = bigint<2*toy_curve_q_limbs>("733512464760947933639364383815815166593630051144");
    toy_curve_Fq2::s = 4;
    toy_curve_Fq2::t = bigint<2*toy_curve_q_limbs>("91689058095118491704920547976976895824203756393");
    toy_curve_Fq2::t_minus_1_over_2 = bigint<2*toy_curve_q_limbs>("45844529047559245852460273988488447912101878196");
    toy_curve_Fq2::non_residue = toy_curve_Fq("1211208045515672698496982");
    toy_curve_Fq2::nqr = toy_curve_Fq2(toy_curve_Fq("2"),toy_curve_Fq("1"));
    toy_curve_Fq2::nqr_to_t = toy_curve_Fq2(toy_curve_Fq("1048851813973016762019162"),toy_curve_Fq("178987081698553272687971"));
    toy_curve_Fq2::Frobenius_coeffs_c1[0] = toy_curve_Fq("1");
    toy_curve_Fq2::Frobenius_coeffs_c1[1] = toy_curve_Fq("1211208045515672698496982");

    /* parameters for Fq6 */

    toy_curve_Fq6::non_residue = toy_curve_Fq2(toy_curve_Fq("9"),toy_curve_Fq("1"));
    toy_curve_Fq6::Frobenius_coeffs_c1[0] = toy_curve_Fq2(toy_curve_Fq("1"),toy_curve_Fq("0"));
    toy_curve_Fq6::Frobenius_coeffs_c1[1] = toy_curve_Fq2(toy_curve_Fq("920235461909360806387090"),toy_curve_Fq("557781887453241790070377"));
    toy_curve_Fq6::Frobenius_coeffs_c1[2] = toy_curve_Fq2(toy_curve_Fq("1414034297008940566"),toy_curve_Fq("0"));
    toy_curve_Fq6::Frobenius_coeffs_c1[3] = toy_curve_Fq2(toy_curve_Fq("562915637257563914341675"),toy_curve_Fq("408288931399064860138241"));
    toy_curve_Fq6::Frobenius_coeffs_c1[4] = toy_curve_Fq2(toy_curve_Fq("1211206631481375689556416"),toy_curve_Fq("0"));
    toy_curve_Fq6::Frobenius_coeffs_c1[5] = toy_curve_Fq2(toy_curve_Fq("939264991864420676265201"),toy_curve_Fq("245137226663366048288365"));
    toy_curve_Fq6::Frobenius_coeffs_c2[0] = toy_curve_Fq2(toy_curve_Fq("1"),toy_curve_Fq("0"));
    toy_curve_Fq6::Frobenius_coeffs_c2[1] = toy_curve_Fq2(toy_curve_Fq("478503903485433842330536"),toy_curve_Fq("943119613315175336864312"));
    toy_curve_Fq6::Frobenius_coeffs_c2[2] = toy_curve_Fq2(toy_curve_Fq("1211206631481375689556416"),toy_curve_Fq("0"));
    toy_curve_Fq6::Frobenius_coeffs_c2[3] = toy_curve_Fq2(toy_curve_Fq("902686482693806422192136"),toy_curve_Fq("1043810449293902033798424"));
    toy_curve_Fq6::Frobenius_coeffs_c2[4] = toy_curve_Fq2(toy_curve_Fq("1414034297008940566"),toy_curve_Fq("0"));
    toy_curve_Fq6::Frobenius_coeffs_c2[5] = toy_curve_Fq2(toy_curve_Fq("1041225704852105132471294"),toy_curve_Fq("435486028422268026331230"));

    /* parameters for Fq12 */

    toy_curve_Fq12::non_residue = toy_curve_Fq2(toy_curve_Fq("9"),toy_curve_Fq("1"));
    toy_curve_Fq12::Frobenius_coeffs_c1[0]  = toy_curve_Fq2(toy_curve_Fq("1"),toy_curve_Fq("0"));
    toy_curve_Fq12::Frobenius_coeffs_c1[1]  = toy_curve_Fq2(toy_curve_Fq("792583545297129969483201"),toy_curve_Fq("284873871376666073012322"));
    toy_curve_Fq12::Frobenius_coeffs_c1[2]  = toy_curve_Fq2(toy_curve_Fq("1211206631481375689556416"),toy_curve_Fq("0"));
    toy_curve_Fq12::Frobenius_coeffs_c1[3]  = toy_curve_Fq2(toy_curve_Fq("536479378566875297433672"),toy_curve_Fq("421052740272693079504280"));
    toy_curve_Fq12::Frobenius_coeffs_c1[4]  = toy_curve_Fq2(toy_curve_Fq("1414034297008940566"),toy_curve_Fq("0"));
    toy_curve_Fq12::Frobenius_coeffs_c1[5]  = toy_curve_Fq2(toy_curve_Fq("1093353167167340130077093"),toy_curve_Fq("505281433866313545980381"));
    toy_curve_Fq12::Frobenius_coeffs_c1[6]  = toy_curve_Fq2(toy_curve_Fq("1"),toy_curve_Fq("0"));
    toy_curve_Fq12::Frobenius_coeffs_c1[7]  = toy_curve_Fq2(toy_curve_Fq("792583545297129969483201"),toy_curve_Fq("284873871376666073012322"));
    toy_curve_Fq12::Frobenius_coeffs_c1[8]  = toy_curve_Fq2(toy_curve_Fq("1211206631481375689556416"),toy_curve_Fq("0"));
    toy_curve_Fq12::Frobenius_coeffs_c1[9]  = toy_curve_Fq2(toy_curve_Fq("536479378566875297433672"),toy_curve_Fq("421052740272693079504280"));
    toy_curve_Fq12::Frobenius_coeffs_c1[10] = toy_curve_Fq2(toy_curve_Fq("1414034297008940566"),toy_curve_Fq("0"));
    toy_curve_Fq12::Frobenius_coeffs_c1[11] = toy_curve_Fq2(toy_curve_Fq("1093353167167340130077093"),toy_curve_Fq("505281433866313545980381"));

    /* choice of short Weierstrass curve and its twist */

    toy_curve_coeff_b = toy_curve_Fq("11337408");
    toy_curve_twist = toy_curve_Fq2(toy_curve_Fq("2"), toy_curve_Fq("1"));
    toy_curve_twist_coeff_b = toy_curve_coeff_b * toy_curve_twist.inverse();
    // toy_curve_twist_coeff_b = toy_curve_Fq2(toy_curve_Fq("726724827309403623633153"), toy_curve_Fq("242241609103134537431915"));
    
    // TODO
    toy_curve_twist_mul_by_b_c0 = toy_curve_coeff_b * toy_curve_Fq2::non_residue;
    toy_curve_twist_mul_by_b_c1 = toy_curve_coeff_b * toy_curve_Fq2::non_residue;
    
    // (238034680536228445107023*z + 1199867752729316448016774 , 328515184787529717232823*z + 769861615151601207864903)
    toy_curve_twist_mul_by_q_X = toy_curve_Fq2(toy_curve_Fq("1199867752729316448016774"),
                                               toy_curve_Fq("238034680536228445107023"));
    toy_curve_twist_mul_by_q_Y = toy_curve_Fq2(toy_curve_Fq("769861615151601207864903"),
                                               toy_curve_Fq("328515184787529717232823"));

    /* choice of group G1 */
    toy_curve_G1::G1_zero = toy_curve_G1(toy_curve_Fq::zero(),
                                     toy_curve_Fq::one(),
                                     toy_curve_Fq::zero());
    toy_curve_G1::G1_one = toy_curve_G1(toy_curve_Fq("681251598544902089201317"),
                                    toy_curve_Fq("599134107450332280136757"),
                                    toy_curve_Fq::one());

    // TODO
    // toy_curve_G1::wnaf_window_table.resize(0);
    // toy_curve_G1::wnaf_window_table.push_back(11);
    // toy_curve_G1::wnaf_window_table.push_back(24);
    // toy_curve_G1::wnaf_window_table.push_back(60);
    // toy_curve_G1::wnaf_window_table.push_back(127);

    // TODO
    // toy_curve_G1::fixed_base_exp_window_table.resize(0);
    // // window 1 is unbeaten in [-inf, 4.99]
    // toy_curve_G1::fixed_base_exp_window_table.push_back(1);
    // // window 2 is unbeaten in [4.99, 10.99]
    // toy_curve_G1::fixed_base_exp_window_table.push_back(5);
    // // window 3 is unbeaten in [10.99, 32.29]
    // toy_curve_G1::fixed_base_exp_window_table.push_back(11);
    // // window 4 is unbeaten in [32.29, 55.23]
    // toy_curve_G1::fixed_base_exp_window_table.push_back(32);
    // // window 5 is unbeaten in [55.23, 162.03]
    // toy_curve_G1::fixed_base_exp_window_table.push_back(55);
    // // window 6 is unbeaten in [162.03, 360.15]
    // toy_curve_G1::fixed_base_exp_window_table.push_back(162);
    // // window 7 is unbeaten in [360.15, 815.44]
    // toy_curve_G1::fixed_base_exp_window_table.push_back(360);
    // // window 8 is unbeaten in [815.44, 2373.07]
    // toy_curve_G1::fixed_base_exp_window_table.push_back(815);
    // // window 9 is unbeaten in [2373.07, 6977.75]
    // toy_curve_G1::fixed_base_exp_window_table.push_back(2373);
    // // window 10 is unbeaten in [6977.75, 7122.23]
    // toy_curve_G1::fixed_base_exp_window_table.push_back(6978);
    // // window 11 is unbeaten in [7122.23, 57818.46]
    // toy_curve_G1::fixed_base_exp_window_table.push_back(7122);
    // // window 12 is never the best
    // toy_curve_G1::fixed_base_exp_window_table.push_back(0);
    // // window 13 is unbeaten in [57818.46, 169679.14]
    // toy_curve_G1::fixed_base_exp_window_table.push_back(57818);
    // // window 14 is never the best
    // toy_curve_G1::fixed_base_exp_window_table.push_back(0);
    // // window 15 is unbeaten in [169679.14, 439758.91]
    // toy_curve_G1::fixed_base_exp_window_table.push_back(169679);
    // // window 16 is unbeaten in [439758.91, 936073.41]
    // toy_curve_G1::fixed_base_exp_window_table.push_back(439759);
    // // window 17 is unbeaten in [936073.41, 4666554.74]
    // toy_curve_G1::fixed_base_exp_window_table.push_back(936073);
    // // window 18 is never the best
    // toy_curve_G1::fixed_base_exp_window_table.push_back(0);
    // // window 19 is unbeaten in [4666554.74, 7580404.42]
    // toy_curve_G1::fixed_base_exp_window_table.push_back(4666555);
    // // window 20 is unbeaten in [7580404.42, 34552892.20]
    // toy_curve_G1::fixed_base_exp_window_table.push_back(7580404);
    // // window 21 is never the best
    // toy_curve_G1::fixed_base_exp_window_table.push_back(0);
    // // window 22 is unbeaten in [34552892.20, inf]
    // toy_curve_G1::fixed_base_exp_window_table.push_back(34552892);

    /* choice of group G2 */

    toy_curve_G2::G2_zero = toy_curve_G2(toy_curve_Fq2::zero(),
                                     toy_curve_Fq2::one(),
                                     toy_curve_Fq2::zero());

    toy_curve_G2::G2_one = toy_curve_G2(toy_curve_Fq2(toy_curve_Fq("721310978809192310510685"),
                                                      toy_curve_Fq("1139677188336733367133892")),
                                        toy_curve_Fq2(toy_curve_Fq("1117722825184865630397924"),
                                                      toy_curve_Fq("794050534252954482295013")),
                                        toy_curve_Fq2::one());

    // TODO
    // toy_curve_G2::wnaf_window_table.resize(0);
    // toy_curve_G2::wnaf_window_table.push_back(5);
    // toy_curve_G2::wnaf_window_table.push_back(15);
    // toy_curve_G2::wnaf_window_table.push_back(39);
    // toy_curve_G2::wnaf_window_table.push_back(109);

    // TODO
    // toy_curve_G2::fixed_base_exp_window_table.resize(0);
    // // window 1 is unbeaten in [-inf, 5.10]
    // toy_curve_G2::fixed_base_exp_window_table.push_back(1);
    // // window 2 is unbeaten in [5.10, 10.43]
    // toy_curve_G2::fixed_base_exp_window_table.push_back(5);
    // // window 3 is unbeaten in [10.43, 25.28]
    // toy_curve_G2::fixed_base_exp_window_table.push_back(10);
    // // window 4 is unbeaten in [25.28, 59.00]
    // toy_curve_G2::fixed_base_exp_window_table.push_back(25);
    // // window 5 is unbeaten in [59.00, 154.03]
    // toy_curve_G2::fixed_base_exp_window_table.push_back(59);
    // // window 6 is unbeaten in [154.03, 334.25]
    // toy_curve_G2::fixed_base_exp_window_table.push_back(154);
    // // window 7 is unbeaten in [334.25, 742.58]
    // toy_curve_G2::fixed_base_exp_window_table.push_back(334);
    // // window 8 is unbeaten in [742.58, 2034.40]
    // toy_curve_G2::fixed_base_exp_window_table.push_back(743);
    // // window 9 is unbeaten in [2034.40, 4987.56]
    // toy_curve_G2::fixed_base_exp_window_table.push_back(2034);
    // // window 10 is unbeaten in [4987.56, 8888.27]
    // toy_curve_G2::fixed_base_exp_window_table.push_back(4988);
    // // window 11 is unbeaten in [8888.27, 26271.13]
    // toy_curve_G2::fixed_base_exp_window_table.push_back(8888);
    // // window 12 is unbeaten in [26271.13, 39768.20]
    // toy_curve_G2::fixed_base_exp_window_table.push_back(26271);
    // // window 13 is unbeaten in [39768.20, 106275.75]
    // toy_curve_G2::fixed_base_exp_window_table.push_back(39768);
    // // window 14 is unbeaten in [106275.75, 141703.40]
    // toy_curve_G2::fixed_base_exp_window_table.push_back(106276);
    // // window 15 is unbeaten in [141703.40, 462422.97]
    // toy_curve_G2::fixed_base_exp_window_table.push_back(141703);
    // // window 16 is unbeaten in [462422.97, 926871.84]
    // toy_curve_G2::fixed_base_exp_window_table.push_back(462423);
    // // window 17 is unbeaten in [926871.84, 4873049.17]
    // toy_curve_G2::fixed_base_exp_window_table.push_back(926872);
    // // window 18 is never the best
    // toy_curve_G2::fixed_base_exp_window_table.push_back(0);
    // // window 19 is unbeaten in [4873049.17, 5706707.88]
    // toy_curve_G2::fixed_base_exp_window_table.push_back(4873049);
    // // window 20 is unbeaten in [5706707.88, 31673814.95]
    // toy_curve_G2::fixed_base_exp_window_table.push_back(5706708);
    // // window 21 is never the best
    // toy_curve_G2::fixed_base_exp_window_table.push_back(0);
    // // window 22 is unbeaten in [31673814.95, inf]
    // toy_curve_G2::fixed_base_exp_window_table.push_back(31673815);

    /* pairing parameters */

    toy_curve_ate_loop_count = bigint_q("2569688"); // 6x+2
    toy_curve_ate_is_loop_count_neg = false;
    toy_curve_final_exponent = bigint<12*toy_curve_q_limbs>("8230120935534082392488788604567970523486776024120818076811068198858268218988401383942947953904346117128229154366615933535994788726253110512101997828009013749721407671124731898055629375957137719841356872465220670708850071238945381194646598264826646428911691925341280");
    toy_curve_final_exponent_z = bigint_q("428281"); // x
    toy_curve_final_exponent_is_z_neg = false;

}
} // libff
