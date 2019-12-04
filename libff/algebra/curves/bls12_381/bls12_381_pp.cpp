/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/bls12_381/bls12_381_pp.hpp>
#include <libff/algebra/curves/bls12.tcc>


namespace libff {


bigint<bls12_381_r_limbs> bls12_381_modulus_r;
bigint<bls12_381_q_limbs> bls12_381_modulus_q;


bls12_381_G1 bls12_381_G1::_one;
bls12_381_G1 bls12_381_G1::_zero;
bls12_381_Fq bls12_381_G1::coeff_b;
std::vector<size_t> bls12_381_G1::wnaf_window_table;
std::vector<size_t> bls12_381_G1::fixed_base_exp_window_table;


bls12_381_G2 bls12_381_G2::_one;
bls12_381_G2 bls12_381_G2::_zero;
bls12_381_Fq2 bls12_381_G2::coeff_b;
std::vector<size_t> bls12_381_G2::wnaf_window_table;
std::vector<size_t> bls12_381_G2::fixed_base_exp_window_table;


void bls12_381_pp::init_public_params()
{
    typedef bigint<bls12_381_r_limbs> bigint_r;
    typedef bigint<bls12_381_q_limbs> bigint_q;

    assert(sizeof(mp_limb_t) == 8 || sizeof(mp_limb_t) == 4); // Montgomery assumes this

    /* parameters for scalar field Fr */

    bls12_381_modulus_r = bigint_r("52435875175126190479447740508185965837690552500527637822603658699938581184513");
    assert(bls12_381_Fr::modulus_is_valid());
    if constexpr(sizeof(mp_limb_t) == 8)
    {
        bls12_381_Fr::Rsquared = bigint_r("3294906474794265442129797520630710739278575682199800681788903916070560242797"); // Rsquared = (W**k)**2 % r where k=4
        bls12_381_Fr::Rcubed = bigint_r("49829253988540319354550742249276084460127446355315915089527227471280320770991");
        bls12_381_Fr::inv = 0xfffffffeffffffff; // (-1/modulus) mod W
    }
    else if constexpr(sizeof(mp_limb_t) == 4)
    {
        bls12_381_Fr::Rsquared = bigint_r("3294906474794265442129797520630710739278575682199800681788903916070560242797");
        bls12_381_Fr::Rcubed = bigint_r("49829253988540319354550742249276084460127446355315915089527227471280320770991");
        bls12_381_Fr::inv = 0xffffffff;
    }
    bls12_381_Fr::num_bits = 255;
    bls12_381_Fr::euler = bigint_r("26217937587563095239723870254092982918845276250263818911301829349969290592256");
    bls12_381_Fr::s = 32; // 2-adic order of modulus-1
    bls12_381_Fr::t = bigint_r("12208678567578594777604504606729831043093128246378069236549469339647"); //(modulus-1)/2^s
    bls12_381_Fr::t_minus_1_over_2 = bigint_r("6104339283789297388802252303364915521546564123189034618274734669823");
    bls12_381_Fr::multiplicative_generator = bls12_381_Fr("7");
    bls12_381_Fr::root_of_unity = bls12_381_Fr("10238227357739495823651030575849232062558860180284477541189508159991286009131");
    bls12_381_Fr::nqr = bls12_381_Fr("5");
    bls12_381_Fr::nqr_to_t = bls12_381_Fr("937917089079007706106976984802249742464848817460758522850752807661925904159");
   
    /* parameters for base field Fq */
    bls12_381_modulus_q = bigint_q("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787");
    assert(bls12_381_Fq::modulus_is_valid());
    if constexpr(sizeof(mp_limb_t) == 8)
    {
        bls12_381_Fq::Rsquared = bigint_q("2708263910654730174793787626328176511836455197166317677006154293982164122222515399004018013397331347120527951271750"); // k=6
        bls12_381_Fq::Rcubed = bigint_q("1639067542774625894236716575548084905938753837211594095883637014582201460755008380976950835174037649440777609978336");

        bls12_381_Fq::inv = 0x89f3fffcfffcfffd;
    }
    else if constexpr(sizeof(mp_limb_t) == 4)
    {
        bls12_381_Fq::Rsquared = bigint_q("2708263910654730174793787626328176511836455197166317677006154293982164122222515399004018013397331347120527951271750");
        bls12_381_Fq::Rcubed = bigint_q("1639067542774625894236716575548084905938753837211594095883637014582201460755008380976950835174037649440777609978336");
        bls12_381_Fq::inv = 0xfffcfffd;
    }
    bls12_381_Fq::num_bits = 381;
    bls12_381_Fq::euler = bigint_q("2001204777610833696708894912867952078278441409969503942666029068062015825245418932221343814564507832018947136279893");
    bls12_381_Fq::s = 1;
    bls12_381_Fq::t = bigint_q("2001204777610833696708894912867952078278441409969503942666029068062015825245418932221343814564507832018947136279893");
    bls12_381_Fq::t_minus_1_over_2 = bigint_q("1000602388805416848354447456433976039139220704984751971333014534031007912622709466110671907282253916009473568139946");
    bls12_381_Fq::multiplicative_generator = bls12_381_Fq("2");
    bls12_381_Fq::root_of_unity = bls12_381_Fq("2");
    bls12_381_Fq::nqr = bls12_381_Fq("2");
    bls12_381_Fq::nqr_to_t = bls12_381_Fq("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559786");

    /* parameters for twist field Fq2 */
    bls12_381_Fq2::euler = bigint<2*bls12_381_q_limbs>("8009641123864852705971874322159486308847560049665276329931192268492988374245678571700328039651096714987477192770085365265551942269853452968100101210518217905546506517135906379008203984748165830709270511838887449985712996744742684");
    bls12_381_Fq2::s = 3;
    bls12_381_Fq2::t = bigint<2*bls12_381_q_limbs>("2002410280966213176492968580539871577211890012416319082482798067123247093561419642925082009912774178746869298192521341316387985567463363242025025302629554476386626629283976594752050996187041457677317627959721862496428249186185671");
    bls12_381_Fq2::t_minus_1_over_2 = bigint<2*bls12_381_q_limbs>("1001205140483106588246484290269935788605945006208159541241399033561623546780709821462541004956387089373434649096260670658193992783731681621012512651314777238193313314641988297376025498093520728838658813979860931248214124593092835");
    bls12_381_Fq2::non_residue = bls12_381_Fq("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559786");
    bls12_381_Fq2::nqr = bls12_381_Fq2(bls12_381_Fq("1"),bls12_381_Fq("1")); // u+1
    bls12_381_Fq2::nqr_to_t = bls12_381_Fq2(bls12_381_Fq("1028732146235106349975324479215795277384839936929757896155643118032610843298655225875571310552543014690878354869257"),bls12_381_Fq("2973677408986561043442465346520108879172042883009249989176415018091420807192182638567116318576472649347015917690530"));
    bls12_381_Fq2::Frobenius_coeffs_c1[0] = bls12_381_Fq("1");
    bls12_381_Fq2::Frobenius_coeffs_c1[1] = bls12_381_Fq("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559786");

    /* parameters for Fq6 */
    bls12_381_Fq6::non_residue = bls12_381_Fq2::nqr; // u + 1
    bls12_381_Fq6::Frobenius_coeffs_c1[0] = bls12_381_Fq2(bls12_381_Fq("1"), bls12_381_Fq("0"));
    bls12_381_Fq6::Frobenius_coeffs_c1[1] = bls12_381_Fq2(bls12_381_Fq("0"), bls12_381_Fq("4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939436"));
    bls12_381_Fq6::Frobenius_coeffs_c1[2] = bls12_381_Fq2(bls12_381_Fq("793479390729215512621379701633421447060886740281060493010456487427281649075476305620758731620350"), bls12_381_Fq("0"));
    bls12_381_Fq6::Frobenius_coeffs_c1[3] = bls12_381_Fq2(bls12_381_Fq("0"), bls12_381_Fq("1"));
    bls12_381_Fq6::Frobenius_coeffs_c1[4] = bls12_381_Fq2(bls12_381_Fq("4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939436"), bls12_381_Fq("0"));
    bls12_381_Fq6::Frobenius_coeffs_c1[5] = bls12_381_Fq2(bls12_381_Fq("0"), bls12_381_Fq("793479390729215512621379701633421447060886740281060493010456487427281649075476305620758731620350"));
    bls12_381_Fq6::Frobenius_coeffs_c2[0] = bls12_381_Fq2(bls12_381_Fq("1"), bls12_381_Fq("0"));
    bls12_381_Fq6::Frobenius_coeffs_c2[1] = bls12_381_Fq2(bls12_381_Fq("4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939437"), bls12_381_Fq("0"));
    bls12_381_Fq6::Frobenius_coeffs_c2[2] = bls12_381_Fq2(bls12_381_Fq("4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939436"), bls12_381_Fq("0"));
    bls12_381_Fq6::Frobenius_coeffs_c2[3] = bls12_381_Fq2(bls12_381_Fq("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559786"), bls12_381_Fq("0"));
    bls12_381_Fq6::Frobenius_coeffs_c2[4] = bls12_381_Fq2(bls12_381_Fq("793479390729215512621379701633421447060886740281060493010456487427281649075476305620758731620350"), bls12_381_Fq("0"));
    bls12_381_Fq6::Frobenius_coeffs_c2[5] = bls12_381_Fq2(bls12_381_Fq("793479390729215512621379701633421447060886740281060493010456487427281649075476305620758731620351"), bls12_381_Fq("0"));


    /* parameters for Fq12 */
    bls12_381_Fq12::non_residue = bls12_381_Fq2::nqr;    // u + 1
    bls12_381_Fq12::Frobenius_coeffs_c1[0] = bls12_381_Fq2(bls12_381_Fq("1"), bls12_381_Fq("0"));
    bls12_381_Fq12::Frobenius_coeffs_c1[1] = bls12_381_Fq2(bls12_381_Fq("3850754370037169011952147076051364057158807420970682438676050522613628423219637725072182697113062777891589506424760"), bls12_381_Fq("151655185184498381465642749684540099398075398968325446656007613510403227271200139370504932015952886146304766135027"));
    bls12_381_Fq12::Frobenius_coeffs_c1[2] = bls12_381_Fq2(bls12_381_Fq("793479390729215512621379701633421447060886740281060493010456487427281649075476305620758731620351"), bls12_381_Fq("0"));
    bls12_381_Fq12::Frobenius_coeffs_c1[3] = bls12_381_Fq2(bls12_381_Fq("2973677408986561043442465346520108879172042883009249989176415018091420807192182638567116318576472649347015917690530"), bls12_381_Fq("1028732146235106349975324479215795277384839936929757896155643118032610843298655225875571310552543014690878354869257"));
    bls12_381_Fq12::Frobenius_coeffs_c1[4] = bls12_381_Fq2(bls12_381_Fq("793479390729215512621379701633421447060886740281060493010456487427281649075476305620758731620350"), bls12_381_Fq("0"));
    bls12_381_Fq12::Frobenius_coeffs_c1[5] = bls12_381_Fq2(bls12_381_Fq("3125332594171059424908108096204648978570118281977575435832422631601824034463382777937621250592425535493320683825557"), bls12_381_Fq("877076961050607968509681729531255177986764537961432449499635504522207616027455086505066378536590128544573588734230"));
    bls12_381_Fq12::Frobenius_coeffs_c1[6] = bls12_381_Fq2(bls12_381_Fq("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559786"), bls12_381_Fq("0"));
    bls12_381_Fq12::Frobenius_coeffs_c1[7] = bls12_381_Fq2(bls12_381_Fq("151655185184498381465642749684540099398075398968325446656007613510403227271200139370504932015952886146304766135027"), bls12_381_Fq("3850754370037169011952147076051364057158807420970682438676050522613628423219637725072182697113062777891589506424760"));
    bls12_381_Fq12::Frobenius_coeffs_c1[8] = bls12_381_Fq2(bls12_381_Fq("4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939436"), bls12_381_Fq("0"));
    bls12_381_Fq12::Frobenius_coeffs_c1[9] = bls12_381_Fq2(bls12_381_Fq("1028732146235106349975324479215795277384839936929757896155643118032610843298655225875571310552543014690878354869257"), bls12_381_Fq("2973677408986561043442465346520108879172042883009249989176415018091420807192182638567116318576472649347015917690530"));
    bls12_381_Fq12::Frobenius_coeffs_c1[10] = bls12_381_Fq2(bls12_381_Fq("4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939437"), bls12_381_Fq("0"));
    bls12_381_Fq12::Frobenius_coeffs_c1[11] = bls12_381_Fq2(bls12_381_Fq("877076961050607968509681729531255177986764537961432449499635504522207616027455086505066378536590128544573588734230"), bls12_381_Fq("3125332594171059424908108096204648978570118281977575435832422631601824034463382777937621250592425535493320683825557"));
    
    /* choice of short Weierstrass curve and its twist */
    bls12_381_G1::coeff_b = bls12_381_Fq("4");
    //bls12_381_twist = bls12_381_Fq2::nqr;
    bls12_381_G2::coeff_b = bls12_381_G1::coeff_b * bls12_381_Fq2::nqr;  // M type twist, Fq2::b should be [4,4]
    //bls12_381_twist_mul_by_b_c0 = bls12_381_coeff_b * bls12_381_Fq2::non_residue;
    //bls12_381_twist_mul_by_b_c1 = bls12_381_coeff_b * bls12_381_Fq2::non_residue;

    // xiToPMinus1Over3 is ξ^((p-1)/3) where ξ = i+9.
    //bls12_381_twist_mul_by_q_X = bls12_381_Fq6::Frobenius_coeffs_c1[1]; /*bls12_381_Fq2(bls12_381_Fq("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559786"),
    //                                                                     bls12_381_Fq("0")); */

    // xiToPMinus1Over2 is ξ^((p-1)/2) where ξ = i+9.
    //bls12_381_twist_mul_by_q_Y = bls12_381_Fq2(bls12_381_Fq("1028732146235106349975324479215795277384839936929757896155643118032610843298655225875571310552543014690878354869257"),
    //                                       bls12_381_Fq("1028732146235106349975324479215795277384839936929757896155643118032610843298655225875571310552543014690878354869257"));
     

    /* choice of group G1 */
    bls12_381_G1::_zero = bls12_381_G1(bls12_381_Fq::zero(), bls12_381_Fq::one(), bls12_381_Fq::zero());
    bls12_381_G1::_one = bls12_381_G1(bls12_381_Fq("3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507"),
                                    bls12_381_Fq("1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569"),
                                    bls12_381_Fq::one());


    // TODO: wNAF window table
    bls12_381_G1::wnaf_window_table.resize(0);
    // bls12_381_G1::wnaf_window_table.push_back(11);
    // bls12_381_G1::wnaf_window_table.push_back(24);
    // bls12_381_G1::wnaf_window_table.push_back(60);
    // bls12_381_G1::wnaf_window_table.push_back(127);

    // // TODO: fixed-base exponentiation table
    bls12_381_G1::fixed_base_exp_window_table.resize(0);
    // // window 1 is unbeaten in [-inf, 4.99]
    // bls12_381_G1::fixed_base_exp_window_table.push_back(1);
    // // window 2 is unbeaten in [4.99, 10.99]
    // bls12_381_G1::fixed_base_exp_window_table.push_back(5);
    // // window 3 is unbeaten in [10.99, 32.29]
    // bls12_381_G1::fixed_base_exp_window_table.push_back(11);
    // // window 4 is unbeaten in [32.29, 55.23]
    // bls12_381_G1::fixed_base_exp_window_table.push_back(32);
    // // window 5 is unbeaten in [55.23, 162.03]
    // bls12_381_G1::fixed_base_exp_window_table.push_back(55);
    // // window 6 is unbeaten in [162.03, 360.15]
    // bls12_381_G1::fixed_base_exp_window_table.push_back(162);
    // // window 7 is unbeaten in [360.15, 815.44]
    // bls12_381_G1::fixed_base_exp_window_table.push_back(360);
    // // window 8 is unbeaten in [815.44, 2373.07]
    // bls12_381_G1::fixed_base_exp_window_table.push_back(815);
    // // window 9 is unbeaten in [2373.07, 6977.75]
    // bls12_381_G1::fixed_base_exp_window_table.push_back(2373);
    // // window 10 is unbeaten in [6977.75, 7122.23]
    // bls12_381_G1::fixed_base_exp_window_table.push_back(6978);
    // // window 11 is unbeaten in [7122.23, 57818.46]
    // bls12_381_G1::fixed_base_exp_window_table.push_back(7122);
    // // window 12 is never the best
    // bls12_381_G1::fixed_base_exp_window_table.push_back(0);
    // // window 13 is unbeaten in [57818.46, 169679.14]
    // bls12_381_G1::fixed_base_exp_window_table.push_back(57818);
    // // window 14 is never the best
    // bls12_381_G1::fixed_base_exp_window_table.push_back(0);
    // // window 15 is unbeaten in [169679.14, 439758.91]
    // bls12_381_G1::fixed_base_exp_window_table.push_back(169679);
    // // window 16 is unbeaten in [439758.91, 936073.41]
    // bls12_381_G1::fixed_base_exp_window_table.push_back(439759);
    // // window 17 is unbeaten in [936073.41, 4666554.74]
    // bls12_381_G1::fixed_base_exp_window_table.push_back(936073);
    // // window 18 is never the best
    // bls12_381_G1::fixed_base_exp_window_table.push_back(0);
    // // window 19 is unbeaten in [4666554.74, 7580404.42]
    // bls12_381_G1::fixed_base_exp_window_table.push_back(4666555);
    // // window 20 is unbeaten in [7580404.42, 34552892.20]
    // bls12_381_G1::fixed_base_exp_window_table.push_back(7580404);
    // // window 21 is never the best
    // bls12_381_G1::fixed_base_exp_window_table.push_back(0);
    // // window 22 is unbeaten in [34552892.20, inf]
    // bls12_381_G1::fixed_base_exp_window_table.push_back(34552892);


    /* choice of group G2 */
    bls12_381_G2::_zero = bls12_381_G2(
                            bls12_381_Fq2::zero(),
                            bls12_381_Fq2::one(),
                            bls12_381_Fq2::zero());

    // simple G2 generator
    bls12_381_G2::_one = bls12_381_G2(
                            bls12_381_Fq2(
                                bls12_381_Fq("352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160"),
                                bls12_381_Fq("3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758")),
                            bls12_381_Fq2(
                                bls12_381_Fq("1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905"),
                                bls12_381_Fq("927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582")),
                            bls12_381_Fq2::one());


    // // TODO: wNAF window table
    bls12_381_G2::wnaf_window_table.resize(0);
    // bls12_381_G2::wnaf_window_table.push_back(5);
    // bls12_381_G2::wnaf_window_table.push_back(15);
    // bls12_381_G2::wnaf_window_table.push_back(39);
    // bls12_381_G2::wnaf_window_table.push_back(109);

    // // TODO: fixed-base exponentiation table 
    bls12_381_G2::fixed_base_exp_window_table.resize(0);
    // // window 1 is unbeaten in [-inf, 5.10]
    // bls12_381_G2::fixed_base_exp_window_table.push_back(1);
    // // window 2 is unbeaten in [5.10, 10.43]
    // bls12_381_G2::fixed_base_exp_window_table.push_back(5);
    // // window 3 is unbeaten in [10.43, 25.28]
    // bls12_381_G2::fixed_base_exp_window_table.push_back(10);
    // // window 4 is unbeaten in [25.28, 59.00]
    // bls12_381_G2::fixed_base_exp_window_table.push_back(25);
    // // window 5 is unbeaten in [59.00, 154.03]
    // bls12_381_G2::fixed_base_exp_window_table.push_back(59);
    // // window 6 is unbeaten in [154.03, 334.25]
    // bls12_381_G2::fixed_base_exp_window_table.push_back(154);
    // // window 7 is unbeaten in [334.25, 742.58]
    // bls12_381_G2::fixed_base_exp_window_table.push_back(334);
    // // window 8 is unbeaten in [742.58, 2034.40]
    // bls12_381_G2::fixed_base_exp_window_table.push_back(743);
    // // window 9 is unbeaten in [2034.40, 4987.56]
    // bls12_381_G2::fixed_base_exp_window_table.push_back(2034);
    // // window 10 is unbeaten in [4987.56, 8888.27]
    // bls12_381_G2::fixed_base_exp_window_table.push_back(4988);
    // // window 11 is unbeaten in [8888.27, 26271.13]
    // bls12_381_G2::fixed_base_exp_window_table.push_back(8888);
    // // window 12 is unbeaten in [26271.13, 39768.20]
    // bls12_381_G2::fixed_base_exp_window_table.push_back(26271);
    // // window 13 is unbeaten in [39768.20, 106275.75]
    // bls12_381_G2::fixed_base_exp_window_table.push_back(39768);
    // // window 14 is unbeaten in [106275.75, 141703.40]
    // bls12_381_G2::fixed_base_exp_window_table.push_back(106276);
    // // window 15 is unbeaten in [141703.40, 462422.97]
    // bls12_381_G2::fixed_base_exp_window_table.push_back(141703);
    // // window 16 is unbeaten in [462422.97, 926871.84]
    // bls12_381_G2::fixed_base_exp_window_table.push_back(462423);
    // // window 17 is unbeaten in [926871.84, 4873049.17]
    // bls12_381_G2::fixed_base_exp_window_table.push_back(926872);
    // // window 18 is never the best
    // bls12_381_G2::fixed_base_exp_window_table.push_back(0);
    // // window 19 is unbeaten in [4873049.17, 5706707.88]
    // bls12_381_G2::fixed_base_exp_window_table.push_back(4873049);
    // // window 20 is unbeaten in [5706707.88, 31673814.95]
    // bls12_381_G2::fixed_base_exp_window_table.push_back(5706708);
    // // window 21 is never the best
    // bls12_381_G2::fixed_base_exp_window_table.push_back(0);
    // // window 22 is unbeaten in [31673814.95, inf]
    // bls12_381_G2::fixed_base_exp_window_table.push_back(31673815);

    /* 
     * https://eprint.iacr.org/2017/1174.pdf
     * ate_loop_count=t where t in the value chosen in q(t) and r(t) parameterization of BLS curve.
     * for BLS12_377 t=3·2^46·(7·13·499)+1 s.t. t=1 (mod 3·2^46) to have a high 2-adicity for q and r. 
    */
}


bls12_381_pp::G1_precomp_type bls12_381_pp::precompute_G1(const G1_type &P)
{
    return P;
}


bls12_381_pp::G2_precomp_type bls12_381_pp::precompute_G2(const G2_type &Q)
{
    return bls12::G2Prepared<bls12_381_pp>(Q);
}


bls12_381_Fq12 bls12_381_pp::miller_loop(const G1_precomp_type &prec_P,
                                         const G2_precomp_type &prec_Q)
{
    return bls12::miller_loop<bls12_381_pp>({
        bls12::PreparedPair<bls12_381_pp>(prec_P, prec_Q)
    });
}


bls12_381_Fq12 bls12_381_pp::double_miller_loop(const G1_precomp_type &prec_P1,
                                                const G2_precomp_type &prec_Q1,
                                                const G1_precomp_type &prec_P2,
                                                const G2_precomp_type &prec_Q2)
{
    return bls12::miller_loop<bls12_381_pp>({
        bls12::PreparedPair<bls12_381_pp>(prec_P1, prec_Q1),
        bls12::PreparedPair<bls12_381_pp>(prec_P2, prec_Q2)
    });
}


bls12_381_Fq12 bls12_381_pp::final_exponentiation(const bls12_381_Fq12 &elt)
{
    return bls12::final_exponentiation<bls12_381_pp>(elt);
}


bls12_381_Fq12 bls12_381_pp::pairing(const G1_type &P,
                                     const G2_type &Q)
{
    return bls12::miller_loop<bls12_381_pp>(P, Q);
}


bls12_381_Fq12 bls12_381_pp::reduced_pairing(const G1_type &P,
                                             const G2_type &Q)
{
    return final_exponentiation(pairing(P, Q));
}


} // libff
