#include <libff/algebra/curves/bw12_446/bw12_446_g1.hpp>
#include <libff/algebra/curves/bw12_446/bw12_446_g2.hpp>
#include <libff/algebra/curves/bw12_446/bw12_446_init.hpp>

namespace libff {

bigint<bw12_446_r_limbs> bw12_446_modulus_r;
bigint<bw12_446_q_limbs> bw12_446_modulus_q;

bw12_446_Fq bw12_446_coeff_b;
bw12_446_Fq2 bw12_446_twist;
bw12_446_Fq2 bw12_446_twist_coeff_b;
bw12_446_Fq bw12_446_twist_mul_by_b_c0;
bw12_446_Fq bw12_446_twist_mul_by_b_c1;
bw12_446_Fq2 bw12_446_twist_mul_by_q_X;
bw12_446_Fq2 bw12_446_twist_mul_by_q_Y;

bigint<bw12_446_q_limbs> bw12_446_ate_loop_count;
bool bw12_446_ate_is_loop_count_neg;
bigint<12*bw12_446_q_limbs> bw12_446_final_exponent;
bigint<bw12_446_q_limbs> bw12_446_final_exponent_z;
bool bw12_446_final_exponent_is_z_neg;

void init_bw12_446_params()
{
    typedef bigint<bw12_446_r_limbs> bigint_r;
    typedef bigint<bw12_446_q_limbs> bigint_q;

    assert(sizeof(mp_limb_t) == 8 || sizeof(mp_limb_t) == 4); // Montgomery assumes this

    /* bw12_446 Fr parameters */

    bw12_446_modulus_r = bigint_r("90637159839200800243440280491062214495118814519237505016426951489992654958877780496875521");
    assert(bw12_446_Fr::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
            bw12_446_Fr::Rsquared = bigint_r("50862560107253004884149423848909871830875171837624363981390021923094034434418442324228825");
            bw12_446_Fr::Rcubed = bigint_r("50675323606149040807823350422840507188867444679633985479770931554234217690988076251041885");
            bw12_446_Fr::inv = 0xffffff9fffffffff;
    }
    if (sizeof(mp_limb_t) == 4)
    {
            bw12_446_Fr::Rsquared = bigint_r("50862560107253004884149423848909871830875171837624363981390021923094034434418442324228825");
            bw12_446_Fr::Rcubed = bigint_r("50675323606149040807823350422840507188867444679633985479770931554234217690988076251041885");
            bw12_446_Fr::inv = 0xffffffff;
    }
    bw12_446_Fr::num_bits = 296;
    bw12_446_Fr::euler = bigint_r("45318579919600400121720140245531107247559407259618752508213475744996327479438890248437760");
    bw12_446_Fr::s = 37; 
    bw12_446_Fr::t = bigint_r("659472133259993281124045319237209438107084239303640097324583268273219109388285");
    bw12_446_Fr::t_minus_1_over_2 = bigint_r("329736066629996640562022659618604719053542119651820048662291634136609554694142");
    bw12_446_Fr::multiplicative_generator = bw12_446_Fr("2");
    bw12_446_Fr::root_of_unity = bw12_446_Fr("25401307393381059377016791523743214306167571729903391344128117975778465653379228184528783");
    bw12_446_Fr::nqr = bw12_446_Fr("7"); 
    bw12_446_Fr::nqr_to_t = bw12_446_Fr("17557157426661147254776552357054583200036963560273902201178260883255059622477573823170198");

     /* bw12_446 Fq parameters */

    bw12_446_modulus_q = bigint_q("218297830370226601612193514776502382704221475011035725792624279635059316315225569375126760584265176077347911285544462735949991646330887");
    assert(bw12_446_Fq::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
            bw12_446_Fq::Rsquared = bigint_q("214581801153929916038975092469364979435093898507770163846706714227498986532565755035532659537001899755210179971956409038736628463491791");
            bw12_446_Fq::Rcubed = bigint_q("28686717986377017394462670120833580730649520450112924097741077920950688866593634642126838918815964753569582926385973453664307485859093");
            bw12_446_Fq::inv = 0x687d633249249249;
    }
    if (sizeof(mp_limb_t) == 4)
    {
            bw12_446_Fq::Rsquared = bigint_q("214581801153929916038975092469364979435093898507770163846706714227498986532565755035532659537001899755210179971956409038736628463491791");
            bw12_446_Fq::Rcubed = bigint_q("28686717986377017394462670120833580730649520450112924097741077920950688866593634642126838918815964753569582926385973453664307485859093");
            bw12_446_Fq::inv = 0x49249249;
    }
    bw12_446_Fq::num_bits = 447;
    bw12_446_Fq::euler = bigint_q("109148915185113300806096757388251191352110737505517862896312139817529658157612784687563380292132588038673955642772231367974995823165443");
    bw12_446_Fq::s = 1; 
    bw12_446_Fq::t = bigint_q("109148915185113300806096757388251191352110737505517862896312139817529658157612784687563380292132588038673955642772231367974995823165443");
    bw12_446_Fq::t_minus_1_over_2 = bigint_q("54574457592556650403048378694125595676055368752758931448156069908764829078806392343781690146066294019336977821386115683987497911582721");
    bw12_446_Fq::multiplicative_generator = bw12_446_Fq("2");
    bw12_446_Fq::root_of_unity = bw12_446_Fq("1");
    bw12_446_Fq::nqr = bw12_446_Fq("3"); 
    bw12_446_Fq::nqr_to_t = bw12_446_Fq("218297830370226601612193514776502382704221475011035725792624279635059316315225569375126760584265176077347911285544462735949991646330886");

    /* parameters for twist field Fq2 */
    bw12_446_Fq2::euler = bigint<2*bw12_446_q_limbs>("");
    bw12_446_Fq2::s = ;
    bw12_446_Fq2::t = bigint<2*bw12_446_q_limbs>("");
    bw12_446_Fq2::t_minus_1_over_2 = bigint<2*bw12_446_q_limbs>("");
    bw12_446_Fq2::non_residue = bw12_446_Fq(""); 
    bw12_446_Fq2::nqr = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq("")); 
    bw12_446_Fq2::nqr_to_t = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq2::Frobenius_coeffs_c1[0] = bw12_446_Fq("");
    bw12_446_Fq2::Frobenius_coeffs_c1[1] = bw12_446_Fq("");

    /* parameters for Fq6 */

    bw12_446_Fq6::non_residue = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq6::Frobenius_coeffs_c1[0] = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq6::Frobenius_coeffs_c1[1] = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq6::Frobenius_coeffs_c1[2] = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq6::Frobenius_coeffs_c1[3] = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq6::Frobenius_coeffs_c1[4] = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq6::Frobenius_coeffs_c1[5] = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq6::Frobenius_coeffs_c2[0] = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq6::Frobenius_coeffs_c2[1] = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq6::Frobenius_coeffs_c2[2] = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq6::Frobenius_coeffs_c2[3] = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq6::Frobenius_coeffs_c2[4] = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq6::Frobenius_coeffs_c2[5] = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));

    /* parameters for Fq12 */
    bw12_446_Fq12::non_residue = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq12::Frobenius_coeffs_c1[0] = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq12::Frobenius_coeffs_c1[1] = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq12::Frobenius_coeffs_c1[2] = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq12::Frobenius_coeffs_c1[3] = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq12::Frobenius_coeffs_c1[4] = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq12::Frobenius_coeffs_c1[5] = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq12::Frobenius_coeffs_c1[6] = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq12::Frobenius_coeffs_c1[7] = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq12::Frobenius_coeffs_c1[8] = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq12::Frobenius_coeffs_c1[9] = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq12::Frobenius_coeffs_c1[10] = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));
    bw12_446_Fq12::Frobenius_coeffs_c1[11] = bw12_446_Fq2(bw12_446_Fq(""),bw12_446_Fq(""));

    /* choice of short Weierstrass curve and its twist */
    bw12_446_coeff_b = bw12_446_Fq("");
    bw12_446_twist = bw12_446_Fq2(bw12_446_Fq(""), bw12_446_Fq(""));
    bw12_446_twist_coeff_b = bw12_446_coeff_b * bw12_446_twist.inverse();
    // bw12_446_twist_coeff_b = bw12_446_Fq2(bw12_446_Fq(""), bw12_446_Fq(""));
    bw12_446_twist_mul_by_b_c0 = bw12_446_coeff_b * bw12_446_Fq2::non_residue;
    bw12_446_twist_mul_by_b_c1 = bw12_446_coeff_b * bw12_446_Fq2::non_residue;
    bw12_446_twist_mul_by_q_X = bw12_446_Fq2(bw12_446_Fq(""),
                                             bw12_446_Fq(""));
    bw12_446_twist_mul_by_q_Y = bw12_446_Fq2(bw12_446_Fq(""),
                                             bw12_446_Fq(""));

    /* choice of group G1 */
    bw12_446_G1::G1_zero = bw12_446_G1(bw12_446_Fq::zero(),
                                       bw12_446_Fq::one(),
                                       bw12_446_Fq::zero());
    bw12_446_G1::G1_one = bw12_446_G1(bw12_446_Fq(""),
                                      bw12_446_Fq(""),
                                      bw12_446_Fq::one());

    // TODO: wNAF window table
    bw12_446_G1::wnaf_window_table.resize(0);
    bw12_446_G1::wnaf_window_table.push_back(11);
    bw12_446_G1::wnaf_window_table.push_back(24);
    bw12_446_G1::wnaf_window_table.push_back(60);
    bw12_446_G1::wnaf_window_table.push_back(127);

    // TODO: fixed-base exponentiation table
    bw12_446_G1::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.99]
    bw12_446_G1::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.99, 10.99]
    bw12_446_G1::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.99, 32.29]
    bw12_446_G1::fixed_base_exp_window_table.push_back(11);
    // window 4 is unbeaten in [32.29, 55.23]
    bw12_446_G1::fixed_base_exp_window_table.push_back(32);
    // window 5 is unbeaten in [55.23, 162.03]
    bw12_446_G1::fixed_base_exp_window_table.push_back(55);
    // window 6 is unbeaten in [162.03, 360.15]
    bw12_446_G1::fixed_base_exp_window_table.push_back(162);
    // window 7 is unbeaten in [360.15, 815.44]
    bw12_446_G1::fixed_base_exp_window_table.push_back(360);
    // window 8 is unbeaten in [815.44, 2373.07]
    bw12_446_G1::fixed_base_exp_window_table.push_back(815);
    // window 9 is unbeaten in [2373.07, 6977.75]
    bw12_446_G1::fixed_base_exp_window_table.push_back(2373);
    // window 10 is unbeaten in [6977.75, 7122.23]
    bw12_446_G1::fixed_base_exp_window_table.push_back(6978);
    // window 11 is unbeaten in [7122.23, 57818.46]
    bw12_446_G1::fixed_base_exp_window_table.push_back(7122);
    // window 12 is never the best
    bw12_446_G1::fixed_base_exp_window_table.push_back(0);
    // window 13 is unbeaten in [57818.46, 169679.14]
    bw12_446_G1::fixed_base_exp_window_table.push_back(57818);
    // window 14 is never the best
    bw12_446_G1::fixed_base_exp_window_table.push_back(0);
    // window 15 is unbeaten in [169679.14, 439758.91]
    bw12_446_G1::fixed_base_exp_window_table.push_back(169679);
    // window 16 is unbeaten in [439758.91, 936073.41]
    bw12_446_G1::fixed_base_exp_window_table.push_back(439759);
    // window 17 is unbeaten in [936073.41, 4666554.74]
    bw12_446_G1::fixed_base_exp_window_table.push_back(936073);
    // window 18 is never the best
    bw12_446_G1::fixed_base_exp_window_table.push_back(0);
    // window 19 is unbeaten in [4666554.74, 7580404.42]
    bw12_446_G1::fixed_base_exp_window_table.push_back(4666555);
    // window 20 is unbeaten in [7580404.42, 34552892.20]
    bw12_446_G1::fixed_base_exp_window_table.push_back(7580404);
    // window 21 is never the best
    bw12_446_G1::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [34552892.20, inf]
    bw12_446_G1::fixed_base_exp_window_table.push_back(34552892);


    /* choice of group G2 */
    bw12_446_G2::G2_zero = bw12_446_G2(bw12_446_Fq2::zero(),
                                         bw12_446_Fq2::one(),
                                         bw12_446_Fq2::zero());
    // G2 generator
    bw12_446_G2::G2_one = bw12_446_G2(bw12_446_Fq2(bw12_446_Fq(""),
                                                   bw12_446_Fq("")),
                                      bw12_446_Fq2(bw12_446_Fq(""),
                                                   bw12_446_Fq("")),
                                      bw12_446_Fq2::one());


    // TODO: wNAF window table
    bw12_446_G2::wnaf_window_table.resize(0);
    bw12_446_G2::wnaf_window_table.push_back(5);
    bw12_446_G2::wnaf_window_table.push_back(15);
    bw12_446_G2::wnaf_window_table.push_back(39);
    bw12_446_G2::wnaf_window_table.push_back(109);

    // TODO: fixed-base exponentiation table
    bw12_446_G2::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 5.10]
    bw12_446_G2::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [5.10, 10.43]
    bw12_446_G2::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.43, 25.28]
    bw12_446_G2::fixed_base_exp_window_table.push_back(10);
    // window 4 is unbeaten in [25.28, 59.00]
    bw12_446_G2::fixed_base_exp_window_table.push_back(25);
    // window 5 is unbeaten in [59.00, 154.03]
    bw12_446_G2::fixed_base_exp_window_table.push_back(59);
    // window 6 is unbeaten in [154.03, 334.25]
    bw12_446_G2::fixed_base_exp_window_table.push_back(154);
    // window 7 is unbeaten in [334.25, 742.58]
    bw12_446_G2::fixed_base_exp_window_table.push_back(334);
    // window 8 is unbeaten in [742.58, 2034.40]
    bw12_446_G2::fixed_base_exp_window_table.push_back(743);
    // window 9 is unbeaten in [2034.40, 4987.56]
    bw12_446_G2::fixed_base_exp_window_table.push_back(2034);
    // window 10 is unbeaten in [4987.56, 8888.27]
    bw12_446_G2::fixed_base_exp_window_table.push_back(4988);
    // window 11 is unbeaten in [8888.27, 26271.13]
    bw12_446_G2::fixed_base_exp_window_table.push_back(8888);
    // window 12 is unbeaten in [26271.13, 39768.20]
    bw12_446_G2::fixed_base_exp_window_table.push_back(26271);
    // window 13 is unbeaten in [39768.20, 106275.75]
    bw12_446_G2::fixed_base_exp_window_table.push_back(39768);
    // window 14 is unbeaten in [106275.75, 141703.40]
    bw12_446_G2::fixed_base_exp_window_table.push_back(106276);
    // window 15 is unbeaten in [141703.40, 462422.97]
    bw12_446_G2::fixed_base_exp_window_table.push_back(141703);
    // window 16 is unbeaten in [462422.97, 926871.84]
    bw12_446_G2::fixed_base_exp_window_table.push_back(462423);
    // window 17 is unbeaten in [926871.84, 4873049.17]
    bw12_446_G2::fixed_base_exp_window_table.push_back(926872);
    // window 18 is never the best
    bw12_446_G2::fixed_base_exp_window_table.push_back(0);
    // window 19 is unbeaten in [4873049.17, 5706707.88]
    bw12_446_G2::fixed_base_exp_window_table.push_back(4873049);
    // window 20 is unbeaten in [5706707.88, 31673814.95]
    bw12_446_G2::fixed_base_exp_window_table.push_back(5706708);
    // window 21 is never the best
    bw12_446_G2::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [31673814.95, inf]
    bw12_446_G2::fixed_base_exp_window_table.push_back(31673815);



    /* pairing parameters */

    bw12_446_ate_loop_count = bigint_q("");
    bw12_446_ate_is_loop_count_neg = ;
    bw12_446_final_exponent = bigint<12*bw12_446_q_limbs>("");
    bw12_446_final_exponent_z = bigint_q("");
    bw12_446_final_exponent_is_z_neg = ;
}
} // libff
