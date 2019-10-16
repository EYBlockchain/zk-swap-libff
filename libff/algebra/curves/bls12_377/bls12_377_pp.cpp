/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/bls12_377/bls12_377_pp.hpp>
#include <libff/algebra/curves/bls12.tcc>


namespace libff {


bigint<bls12_377_r_limbs> bls12_377_modulus_r;
bigint<bls12_377_q_limbs> bls12_377_modulus_q;


bls12_377_G1 bls12_377_G1::_one;
bls12_377_G1 bls12_377_G1::_zero;
bls12_377_Fq bls12_377_G1::coeff_b;
std::vector<size_t> bls12_377_G1::wnaf_window_table;
std::vector<size_t> bls12_377_G1::fixed_base_exp_window_table;


bls12_377_G2 bls12_377_G2::_one;
bls12_377_G2 bls12_377_G2::_zero;
bls12_377_Fq2 bls12_377_G2::coeff_b;
std::vector<size_t> bls12_377_G2::wnaf_window_table;
std::vector<size_t> bls12_377_G2::fixed_base_exp_window_table;


void bls12_377_pp::init_public_params()
{
    typedef bigint<bls12_377_r_limbs> bigint_r;
    typedef bigint<bls12_377_q_limbs> bigint_q;

    assert(sizeof(mp_limb_t) == 8 || sizeof(mp_limb_t) == 4); // Montgomery assumes this

    /* parameters for scalar field Fr */

    bls12_377_modulus_r = bigint_r("8444461749428370424248824938781546531375899335154063827935233455917409239041");
    assert(bls12_377_Fr::modulus_is_valid());
    if constexpr(sizeof(mp_limb_t) == 8)
    {
        bls12_377_Fr::Rsquared = bigint_r("508595941311779472113692600146818027278633330499214071737745792929336755579"); // Rsquared = (W**k)**2 % r where k=4
        bls12_377_Fr::Rcubed = bigint_r("2717187485423313556320207871216538426353201097398909639086937135091399607628");
        bls12_377_Fr::inv = 0xa117fffffffffff; // (-1/modulus) mod W
    }
    else if constexpr(sizeof(mp_limb_t) == 4)
    {
        bls12_377_Fr::Rsquared = bigint_r("508595941311779472113692600146818027278633330499214071737745792929336755579");
        bls12_377_Fr::Rcubed = bigint_r("2717187485423313556320207871216538426353201097398909639086937135091399607628");
        bls12_377_Fr::inv = 0xffffffff;
    }
    bls12_377_Fr::num_bits = 253;
    bls12_377_Fr::euler = bigint_r("4222230874714185212124412469390773265687949667577031913967616727958704619520");
    bls12_377_Fr::s = 47; // 2-adic order of modulus-1
    bls12_377_Fr::t = bigint_r("60001509534603559531609739528203892656505753216962260608619555"); //(modulus-1)/2^s
    bls12_377_Fr::t_minus_1_over_2 = bigint_r("30000754767301779765804869764101946328252876608481130304309777");
    bls12_377_Fr::multiplicative_generator = bls12_377_Fr("22");
    bls12_377_Fr::root_of_unity = bls12_377_Fr("8065159656716812877374967518403273466521432693661810619979959746626482506078");
    bls12_377_Fr::nqr = bls12_377_Fr("11");
    bls12_377_Fr::nqr_to_t = bls12_377_Fr("6924886788847882060123066508223519077232160750698452411071850219367055984476");
   
    /* parameters for base field Fq */
    bls12_377_modulus_q = bigint_q("258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177");
    assert(bls12_377_Fq::modulus_is_valid());
    if constexpr(sizeof(mp_limb_t) == 8)
    {
        bls12_377_Fq::Rsquared = bigint_q("66127428376872697816332570116866232405230528984664918319606315420233909940404532140033099444330447428417853902114"); // k=6
        bls12_377_Fq::Rcubed = bigint_q("157734475176213061358192738313701451942220138363611391489992831740412033225490229541667992423878570205050777755168");

        bls12_377_Fq::inv = 0x8508bfffffffffff;
    }
    else if constexpr(sizeof(mp_limb_t) == 4)
    {
        bls12_377_Fq::Rsquared = bigint_q("66127428376872697816332570116866232405230528984664918319606315420233909940404532140033099444330447428417853902114");
        bls12_377_Fq::Rcubed = bigint_q("157734475176213061358192738313701451942220138363611391489992831740412033225490229541667992423878570205050777755168");
        bls12_377_Fq::inv = 0xffffffff;
    }
    bls12_377_Fq::num_bits = 377;
    bls12_377_Fq::euler = bigint_q("129332213006484547005326366847446766768196756377457330269942131333360234174170411387484444069786680062220160729088");
    bls12_377_Fq::s = 46;
    bls12_377_Fq::t = bigint_q("3675842578061421676390135839012792950148785745837396071634149488243117337281387659330802195819009059");
    bls12_377_Fq::t_minus_1_over_2 = bigint_q("1837921289030710838195067919506396475074392872918698035817074744121558668640693829665401097909504529");
    bls12_377_Fq::multiplicative_generator = bls12_377_Fq("15");
    bls12_377_Fq::root_of_unity = bls12_377_Fq("146552004846884389553264564610149105174701957497228680529098805315416492923550540437026734404078567406251254115855");
    bls12_377_Fq::nqr = bls12_377_Fq("258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458172"); // -5
    bls12_377_Fq::nqr_to_t = bls12_377_Fq("32863578547254505029601261939868325669770508939375122462904745766352256812585773382134936404344547323199885654433");

    // parameters for twist field Fq2
    bls12_377_Fq2::euler = bigint<2*bls12_377_q_limbs>("33453642642309381258089625946249069288005760010886479253070957453297957116339370141113413635838485065209570299254148838549585056123015878375022724998041828785227090063466658233059433323033772513321990316560167027213559780081664");
    bls12_377_Fq2::s = 47;
    bls12_377_Fq2::t = bigint<2*bls12_377_q_limbs>("475404855284145089315325463221726483993816145966867441829193658311651761271425728823393990805904040047516478740222806302278755994777496288961383541476974255391881599499962735436887347234371823579436839914935817251");
    bls12_377_Fq2::t_minus_1_over_2 = bigint<2*bls12_377_q_limbs>("237702427642072544657662731610863241996908072983433720914596829155825880635712864411696995402952020023758239370111403151139377997388748144480691770738487127695940799749981367718443673617185911789718419957467908625");
    bls12_377_Fq2::non_residue = bls12_377_Fq("258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458172"); // 5 as Fq2[u]/u^2-5
    bls12_377_Fq2::nqr = bls12_377_Fq2(bls12_377_Fq("0"),bls12_377_Fq("1")); // u
    bls12_377_Fq2::nqr_to_t = bls12_377_Fq2(bls12_377_Fq("0"),bls12_377_Fq("257286236321774568987262729980034669694531728092793737444525294935421142460394028155736019924956637466133519652786"));
    bls12_377_Fq2::Frobenius_coeffs_c1[0] = bls12_377_Fq("1");
    bls12_377_Fq2::Frobenius_coeffs_c1[1] = bls12_377_Fq("258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458176");

    // parameters for Fq6
    bls12_377_Fq6::non_residue = bls12_377_Fq2::nqr; //bls12_377_Fq2(bls12_377_Fq("0"),bls12_377_Fq("1"));
    bls12_377_Fq6::Frobenius_coeffs_c1[0] = bls12_377_Fq2(bls12_377_Fq("1"),bls12_377_Fq("0"));
    bls12_377_Fq6::Frobenius_coeffs_c1[1] = bls12_377_Fq2(bls12_377_Fq("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410946"),bls12_377_Fq("0"));
    bls12_377_Fq6::Frobenius_coeffs_c1[2] = bls12_377_Fq2(bls12_377_Fq("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945"),bls12_377_Fq("0"));
    bls12_377_Fq6::Frobenius_coeffs_c1[3] = bls12_377_Fq2(bls12_377_Fq("258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458176"),bls12_377_Fq("0"));
    bls12_377_Fq6::Frobenius_coeffs_c1[4] = bls12_377_Fq2(bls12_377_Fq("258664426012969093929703085429980814127835149614277183275038967946009968870203535512256352201271898244626862047231"),bls12_377_Fq("0"));
    bls12_377_Fq6::Frobenius_coeffs_c1[5] = bls12_377_Fq2(bls12_377_Fq("258664426012969093929703085429980814127835149614277183275038967946009968870203535512256352201271898244626862047232"),bls12_377_Fq("0"));
    bls12_377_Fq6::Frobenius_coeffs_c2[0] = bls12_377_Fq2(bls12_377_Fq("1"),bls12_377_Fq("0"));
    bls12_377_Fq6::Frobenius_coeffs_c2[1] = bls12_377_Fq2(bls12_377_Fq("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945"),bls12_377_Fq("0"));
    bls12_377_Fq6::Frobenius_coeffs_c2[2] = bls12_377_Fq2(bls12_377_Fq("258664426012969093929703085429980814127835149614277183275038967946009968870203535512256352201271898244626862047231"),bls12_377_Fq("0"));
    bls12_377_Fq6::Frobenius_coeffs_c2[3] = bls12_377_Fq2(bls12_377_Fq("1"),bls12_377_Fq("0"));
    bls12_377_Fq6::Frobenius_coeffs_c2[4] = bls12_377_Fq2(bls12_377_Fq("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945"),bls12_377_Fq("0"));
    bls12_377_Fq6::Frobenius_coeffs_c2[5] = bls12_377_Fq2(bls12_377_Fq("258664426012969093929703085429980814127835149614277183275038967946009968870203535512256352201271898244626862047231"),bls12_377_Fq("0"));

    // parameters for Fq12
    bls12_377_Fq12::non_residue = bls12_377_Fq2::nqr; //bls12_377_Fq2(bls12_377_Fq("0"),bls12_377_Fq("1"));
    bls12_377_Fq12::Frobenius_coeffs_c1[0] = bls12_377_Fq2(bls12_377_Fq("1"),bls12_377_Fq("0"));
    bls12_377_Fq12::Frobenius_coeffs_c1[1] = bls12_377_Fq2(bls12_377_Fq("92949345220277864758624960506473182677953048909283248980960104381795901929519566951595905490535835115111760994353"),bls12_377_Fq("0"));
    bls12_377_Fq12::Frobenius_coeffs_c1[2] = bls12_377_Fq2(bls12_377_Fq("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410946"),bls12_377_Fq("0"));
    bls12_377_Fq12::Frobenius_coeffs_c1[3] = bls12_377_Fq2(bls12_377_Fq("216465761340224619389371505802605247630151569547285782856803747159100223055385581585702401816380679166954762214499"),bls12_377_Fq("0"));
    bls12_377_Fq12::Frobenius_coeffs_c1[4] = bls12_377_Fq2(bls12_377_Fq("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945"),bls12_377_Fq("0"));
    bls12_377_Fq12::Frobenius_coeffs_c1[5] = bls12_377_Fq2(bls12_377_Fq("123516416119946754630746545296132064952198520638002533875843642777304321125866014634106496325844844051843001220146"),bls12_377_Fq("0"));
    bls12_377_Fq12::Frobenius_coeffs_c1[6] = bls12_377_Fq2(bls12_377_Fq("258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458176"),bls12_377_Fq("0"));
    bls12_377_Fq12::Frobenius_coeffs_c1[7] = bls12_377_Fq2(bls12_377_Fq("165715080792691229252027773188420350858440463845631411558924158284924566418821255823372982649037525009328560463824"),bls12_377_Fq("0"));
    bls12_377_Fq12::Frobenius_coeffs_c1[8] = bls12_377_Fq2(bls12_377_Fq("258664426012969093929703085429980814127835149614277183275038967946009968870203535512256352201271898244626862047231"),bls12_377_Fq("0"));
    bls12_377_Fq12::Frobenius_coeffs_c1[9] = bls12_377_Fq2(bls12_377_Fq("42198664672744474621281227892288285906241943207628877683080515507620245292955241189266486323192680957485559243678"),bls12_377_Fq("0"));
    bls12_377_Fq12::Frobenius_coeffs_c1[10] = bls12_377_Fq2(bls12_377_Fq("258664426012969093929703085429980814127835149614277183275038967946009968870203535512256352201271898244626862047232"),bls12_377_Fq("0"));
    bls12_377_Fq12::Frobenius_coeffs_c1[11] = bls12_377_Fq2(bls12_377_Fq("135148009893022339379906188398761468584194992116912126664040619889416147222474808140862391813728516072597320238031"),bls12_377_Fq("0"));
    
    
    /* choice of short Weierstrass curve and its twist */
    bls12_377_G1::coeff_b = bls12_377_Fq("1");
    //bls12_377_twist = bls12_377_Fq2(bls12_377_Fq("0"), bls12_377_Fq("1"));
    bls12_377_G2::coeff_b = bls12_377_Fq2(bls12_377_Fq("0"), bls12_377_Fq("155198655607781456406391640216936120121836107652948796323930557600032281009004493664981332883744016074664192874906"));
    //bls12_377_twist_mul_by_b_c0 = bls12_377_coeff_b * bls12_377_Fq2::non_residue;
    //bls12_377_twist_mul_by_b_c1 = bls12_377_coeff_b * bls12_377_Fq2::non_residue;
    //bls12_377_twist_mul_by_q_X = bls12_377_Fq2(bls12_377_Fq("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410946"),
    //                                           bls12_377_Fq("0"));
    //bls12_377_twist_mul_by_q_Y = bls12_377_Fq2(bls12_377_Fq("42198664672744474621281227892288285906241943207628877683080515507620245292955241189266486323192680957485559243678"),
    //                                       bls12_377_Fq("0"));
    // for debug
    // bls12_377_twist_mul_by_q_X = bls12_377_Fq2(bls12_377_Fq("1"),bls12_377_Fq("0"));
    // bls12_377_twist_mul_by_q_Y = bls12_377_Fq2(bls12_377_Fq("1"),bls12_377_Fq("0"));

    /* choice of group G1 */
    bls12_377_G1::_zero = bls12_377_G1(bls12_377_Fq::zero(),
                                     bls12_377_Fq::one(),
                                     bls12_377_Fq::zero());
    bls12_377_G1::_one = bls12_377_G1(bls12_377_Fq("81937999373150964239938255573465948239988671502647976594219695644855304257327692006745978603320413799295628339695"),
                                    bls12_377_Fq("241266749859715473739788878240585681733927191168601896383759122102112907357779751001206799952863815012735208165030"),
                                    bls12_377_Fq::one());

    // simple generator for debug
    // bls12_377_G1::G1_one = bls12_377_G1(bls12_377_Fq::zero(), bls12_377_Fq::one(), bls12_377_Fq::one());

    // TODO: wNAF window table
    // bls12_377_G1::wnaf_window_table.resize(0);
    // bls12_377_G1::wnaf_window_table.push_back(11);
    // bls12_377_G1::wnaf_window_table.push_back(24);
    // bls12_377_G1::wnaf_window_table.push_back(60);
    // bls12_377_G1::wnaf_window_table.push_back(127);

    // // TODO: fixed-base exponentiation table
    // bls12_377_G1::fixed_base_exp_window_table.resize(0);
    // // window 1 is unbeaten in [-inf, 4.99]
    // bls12_377_G1::fixed_base_exp_window_table.push_back(1);
    // // window 2 is unbeaten in [4.99, 10.99]
    // bls12_377_G1::fixed_base_exp_window_table.push_back(5);
    // // window 3 is unbeaten in [10.99, 32.29]
    // bls12_377_G1::fixed_base_exp_window_table.push_back(11);
    // // window 4 is unbeaten in [32.29, 55.23]
    // bls12_377_G1::fixed_base_exp_window_table.push_back(32);
    // // window 5 is unbeaten in [55.23, 162.03]
    // bls12_377_G1::fixed_base_exp_window_table.push_back(55);
    // // window 6 is unbeaten in [162.03, 360.15]
    // bls12_377_G1::fixed_base_exp_window_table.push_back(162);
    // // window 7 is unbeaten in [360.15, 815.44]
    // bls12_377_G1::fixed_base_exp_window_table.push_back(360);
    // // window 8 is unbeaten in [815.44, 2373.07]
    // bls12_377_G1::fixed_base_exp_window_table.push_back(815);
    // // window 9 is unbeaten in [2373.07, 6977.75]
    // bls12_377_G1::fixed_base_exp_window_table.push_back(2373);
    // // window 10 is unbeaten in [6977.75, 7122.23]
    // bls12_377_G1::fixed_base_exp_window_table.push_back(6978);
    // // window 11 is unbeaten in [7122.23, 57818.46]
    // bls12_377_G1::fixed_base_exp_window_table.push_back(7122);
    // // window 12 is never the best
    // bls12_377_G1::fixed_base_exp_window_table.push_back(0);
    // // window 13 is unbeaten in [57818.46, 169679.14]
    // bls12_377_G1::fixed_base_exp_window_table.push_back(57818);
    // // window 14 is never the best
    // bls12_377_G1::fixed_base_exp_window_table.push_back(0);
    // // window 15 is unbeaten in [169679.14, 439758.91]
    // bls12_377_G1::fixed_base_exp_window_table.push_back(169679);
    // // window 16 is unbeaten in [439758.91, 936073.41]
    // bls12_377_G1::fixed_base_exp_window_table.push_back(439759);
    // // window 17 is unbeaten in [936073.41, 4666554.74]
    // bls12_377_G1::fixed_base_exp_window_table.push_back(936073);
    // // window 18 is never the best
    // bls12_377_G1::fixed_base_exp_window_table.push_back(0);
    // // window 19 is unbeaten in [4666554.74, 7580404.42]
    // bls12_377_G1::fixed_base_exp_window_table.push_back(4666555);
    // // window 20 is unbeaten in [7580404.42, 34552892.20]
    // bls12_377_G1::fixed_base_exp_window_table.push_back(7580404);
    // // window 21 is never the best
    // bls12_377_G1::fixed_base_exp_window_table.push_back(0);
    // // window 22 is unbeaten in [34552892.20, inf]
    // bls12_377_G1::fixed_base_exp_window_table.push_back(34552892);


    /* choice of group G2 */
    bls12_377_G2::_zero = bls12_377_G2(bls12_377_Fq2::zero(),
                                         bls12_377_Fq2::one(),
                                         bls12_377_Fq2::zero());

    // simple G2 generator
    bls12_377_G2::_one = bls12_377_G2(bls12_377_Fq2(bls12_377_Fq("233578398248691099356572568220835526895379068987715365179118596935057653620464273615301663571204657964920925606294"),
                                                      bls12_377_Fq("140913150380207355837477652521042157274541796891053068589147167627541651775299824604154852141315666357241556069118")),
                                        bls12_377_Fq2(bls12_377_Fq("63160294768292073209381361943935198908131692476676907196754037919244929611450776219210369229519898517858833747423"),
                                                      bls12_377_Fq("149157405641012693445398062341192467754805999074082136895788947234480009303640899064710353187729182149407503257491")),
                                        bls12_377_Fq2::one());


    // // TODO: wNAF window table
    // bls12_377_G2::wnaf_window_table.resize(0);
    // bls12_377_G2::wnaf_window_table.push_back(5);
    // bls12_377_G2::wnaf_window_table.push_back(15);
    // bls12_377_G2::wnaf_window_table.push_back(39);
    // bls12_377_G2::wnaf_window_table.push_back(109);

    // // TODO: fixed-base exponentiation table 
    // bls12_377_G2::fixed_base_exp_window_table.resize(0);
    // // window 1 is unbeaten in [-inf, 5.10]
    // bls12_377_G2::fixed_base_exp_window_table.push_back(1);
    // // window 2 is unbeaten in [5.10, 10.43]
    // bls12_377_G2::fixed_base_exp_window_table.push_back(5);
    // // window 3 is unbeaten in [10.43, 25.28]
    // bls12_377_G2::fixed_base_exp_window_table.push_back(10);
    // // window 4 is unbeaten in [25.28, 59.00]
    // bls12_377_G2::fixed_base_exp_window_table.push_back(25);
    // // window 5 is unbeaten in [59.00, 154.03]
    // bls12_377_G2::fixed_base_exp_window_table.push_back(59);
    // // window 6 is unbeaten in [154.03, 334.25]
    // bls12_377_G2::fixed_base_exp_window_table.push_back(154);
    // // window 7 is unbeaten in [334.25, 742.58]
    // bls12_377_G2::fixed_base_exp_window_table.push_back(334);
    // // window 8 is unbeaten in [742.58, 2034.40]
    // bls12_377_G2::fixed_base_exp_window_table.push_back(743);
    // // window 9 is unbeaten in [2034.40, 4987.56]
    // bls12_377_G2::fixed_base_exp_window_table.push_back(2034);
    // // window 10 is unbeaten in [4987.56, 8888.27]
    // bls12_377_G2::fixed_base_exp_window_table.push_back(4988);
    // // window 11 is unbeaten in [8888.27, 26271.13]
    // bls12_377_G2::fixed_base_exp_window_table.push_back(8888);
    // // window 12 is unbeaten in [26271.13, 39768.20]
    // bls12_377_G2::fixed_base_exp_window_table.push_back(26271);
    // // window 13 is unbeaten in [39768.20, 106275.75]
    // bls12_377_G2::fixed_base_exp_window_table.push_back(39768);
    // // window 14 is unbeaten in [106275.75, 141703.40]
    // bls12_377_G2::fixed_base_exp_window_table.push_back(106276);
    // // window 15 is unbeaten in [141703.40, 462422.97]
    // bls12_377_G2::fixed_base_exp_window_table.push_back(141703);
    // // window 16 is unbeaten in [462422.97, 926871.84]
    // bls12_377_G2::fixed_base_exp_window_table.push_back(462423);
    // // window 17 is unbeaten in [926871.84, 4873049.17]
    // bls12_377_G2::fixed_base_exp_window_table.push_back(926872);
    // // window 18 is never the best
    // bls12_377_G2::fixed_base_exp_window_table.push_back(0);
    // // window 19 is unbeaten in [4873049.17, 5706707.88]
    // bls12_377_G2::fixed_base_exp_window_table.push_back(4873049);
    // // window 20 is unbeaten in [5706707.88, 31673814.95]
    // bls12_377_G2::fixed_base_exp_window_table.push_back(5706708);
    // // window 21 is never the best
    // bls12_377_G2::fixed_base_exp_window_table.push_back(0);
    // // window 22 is unbeaten in [31673814.95, inf]
    // bls12_377_G2::fixed_base_exp_window_table.push_back(31673815);
}


bls12_377_pp::G1_precomp_type bls12_377_pp::precompute_G1(const bls12_377_G1 &P)
{
    return P;
}


bls12_377_pp::G2_precomp_type bls12_377_pp::precompute_G2(const bls12_377_G2 &Q)
{
    return bls12::G2Prepared<bls12_377_pp>(Q);
}


bls12_377_Fq12 bls12_377_pp::miller_loop(const bls12_377_pp::G1_precomp_type &prec_P,
                                         const bls12_377_pp::G2_precomp_type &prec_Q)
{
    return bls12::miller_loop<bls12_377_pp>({
        bls12::PreparedPair<bls12_377_pp>(prec_P, prec_Q)
    });
}


bls12_377_Fq12 bls12_377_pp::double_miller_loop(const G1_precomp_type &prec_P1,
                                                const G2_precomp_type &prec_Q1,
                                                const G1_precomp_type &prec_P2,
                                                const G2_precomp_type &prec_Q2)
{
    return bls12::miller_loop<bls12_377_pp>({
        bls12::PreparedPair<bls12_377_pp>(prec_P1, prec_Q1),
        bls12::PreparedPair<bls12_377_pp>(prec_P2, prec_Q2)
    });
}


bls12_377_Fq12 bls12_377_pp::final_exponentiation(const bls12_377_Fq12 &elt)
{
    return bls12::final_exponentiation<bls12_377_pp>(elt);
}


bls12_377_Fq12 bls12_377_pp::pairing(const bls12_377_G1 &P,
                                     const bls12_377_G2 &Q)
{
    return bls12::miller_loop<bls12_377_pp>(P, Q);
}


bls12_377_Fq12 bls12_377_pp::reduced_pairing(const bls12_377_G1 &P,
                                             const bls12_377_G2 &Q)
{
    return final_exponentiation(pairing(P, Q));
}


} // libff
