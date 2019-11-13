#include <libff/algebra/curves/test_curve/test_curve_g1.hpp>
#include <libff/algebra/curves/test_curve/test_curve_g2.hpp>
#include <libff/algebra/curves/test_curve/test_curve_init.hpp>

namespace libff {

bigint<test_curve_r_limbs> test_curve_modulus_r;
bigint<test_curve_q_limbs> test_curve_modulus_q;

test_curve_Fq test_curve_coeff_b;
test_curve_Fq2 test_curve_twist;
test_curve_Fq2 test_curve_twist_coeff_b;
test_curve_Fq test_curve_twist_mul_by_b_c0;
test_curve_Fq test_curve_twist_mul_by_b_c1;
test_curve_Fq2 test_curve_twist_mul_by_q_X;
test_curve_Fq2 test_curve_twist_mul_by_q_Y;

bigint<test_curve_q_limbs> test_curve_ate_loop_count;
bool test_curve_ate_is_loop_count_neg;
bigint<12*test_curve_q_limbs> test_curve_final_exponent;
bigint<test_curve_q_limbs> test_curve_final_exponent_z;
bool test_curve_final_exponent_is_z_neg;

void init_test_curve_params()
{
	typedef bigint<test_curve_r_limbs> bigint_r;
	typedef bigint<test_curve_q_limbs> bigint_q;

	assert(sizeof(mp_limb_t) == 8 || sizeof(mp_limb_t) == 4);
	
    /* test_curve Fr parameters */

	test_curve_modulus_r = bigint_r("8444461749428370424248824938781546531375899335154063827935233455917409239041");
	assert(test_curve_Fr::modulus_is_valid());
	if (sizeof(mp_limb_t) == 8)
	{
		test_curve_Fr::Rsquared = bigint_r("508595941311779472113692600146818027278633330499214071737745792929336755579");
		test_curve_Fr::Rcubed = bigint_r("2717187485423313556320207871216538426353201097398909639086937135091399607628");
		test_curve_Fr::inv = 0xa117fffffffffff;
	}
	if (sizeof(mp_limb_t) == 4)
	{
	    test_curve_Fr::Rsquared = bigint_r("508595941311779472113692600146818027278633330499214071737745792929336755579");
	    test_curve_Fr::Rcubed = bigint_r("2717187485423313556320207871216538426353201097398909639086937135091399607628");
	    test_curve_Fr::inv = 0xffffffff;
	}
	test_curve_Fr::num_bits = 253;
	test_curve_Fr::euler = bigint_r("4222230874714185212124412469390773265687949667577031913967616727958704619520");
	test_curve_Fr::s = 47; 
	test_curve_Fr::t = bigint_r("60001509534603559531609739528203892656505753216962260608619555");
	test_curve_Fr::t_minus_1_over_2 = bigint_r("30000754767301779765804869764101946328252876608481130304309777");
	test_curve_Fr::multiplicative_generator = test_curve_Fr("22");
	test_curve_Fr::root_of_unity = test_curve_Fr("8065159656716812877374967518403273466521432693661810619979959746626482506078");
	test_curve_Fr::nqr = test_curve_Fr("11"); 
	test_curve_Fr::nqr_to_t = test_curve_Fr("6924886788847882060123066508223519077232160750698452411071850219367055984476");
	
    /* test_curve Fq parameters */

	test_curve_modulus_q = bigint_q("258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177");
	assert(test_curve_Fq::modulus_is_valid());
	if (sizeof(mp_limb_t) == 8)
	{
		test_curve_Fq::Rsquared = bigint_q("66127428376872697816332570116866232405230528984664918319606315420233909940404532140033099444330447428417853902114");
		test_curve_Fq::Rcubed = bigint_q("157734475176213061358192738313701451942220138363611391489992831740412033225490229541667992423878570205050777755168");
		test_curve_Fq::inv = 0x8508bfffffffffff;
	}
	if (sizeof(mp_limb_t) == 4)
	{
        test_curve_Fq::Rsquared = bigint_q("66127428376872697816332570116866232405230528984664918319606315420233909940404532140033099444330447428417853902114");
        test_curve_Fq::Rcubed = bigint_q("157734475176213061358192738313701451942220138363611391489992831740412033225490229541667992423878570205050777755168");
        test_curve_Fq::inv = 0xffffffff;
	}
	test_curve_Fq::num_bits = 377;
	test_curve_Fq::euler = bigint_q("129332213006484547005326366847446766768196756377457330269942131333360234174170411387484444069786680062220160729088");
	test_curve_Fq::s = 46; 
	test_curve_Fq::t = bigint_q("3675842578061421676390135839012792950148785745837396071634149488243117337281387659330802195819009059");
	test_curve_Fq::t_minus_1_over_2 = bigint_q("1837921289030710838195067919506396475074392872918698035817074744121558668640693829665401097909504529");
	test_curve_Fq::multiplicative_generator = test_curve_Fq("15");
	test_curve_Fq::root_of_unity = test_curve_Fq("32863578547254505029601261939868325669770508939375122462904745766352256812585773382134936404344547323199885654433");
	test_curve_Fq::nqr = test_curve_Fq("5"); 
	test_curve_Fq::nqr_to_t = test_curve_Fq("33774956008227656219775876656288133547078610493828613777258829345740556592044969439504850374928261397247202212840");
	
     /* test_curve Fq2 parameters */

	test_curve_Fq2::euler = bigint<2*test_curve_q_limbs>("33453642642309381258089625946249069288005760010886479253070957453297957116339370141113413635838485065209570299254148838549585056123015878375022724998041828785227090063466658233059433323033772513321990316560167027213559780081664");
	test_curve_Fq2::s = 47;
	test_curve_Fq2::t = bigint<2*test_curve_q_limbs>("475404855284145089315325463221726483993816145966867441829193658311651761271425728823393990805904040047516478740222806302278755994777496288961383541476974255391881599499962735436887347234371823579436839914935817251");
	test_curve_Fq2::t_minus_1_over_2 = bigint<2*test_curve_q_limbs>("237702427642072544657662731610863241996908072983433720914596829155825880635712864411696995402952020023758239370111403151139377997388748144480691770738487127695940799749981367718443673617185911789718419957467908625");
	test_curve_Fq2::non_residue = test_curve_Fq("5");
	test_curve_Fq2::nqr = test_curve_Fq2(test_curve_Fq("0"),test_curve_Fq("1"));
	test_curve_Fq2::nqr_to_t = test_curve_Fq2(test_curve_Fq("0"),test_curve_Fq("1378189691194525023390003714858863841861784662120923095358967731299325887946794619232868214616722658306801805391"));
	test_curve_Fq2::Frobenius_coeffs_c1[0] = test_curve_Fq("1");
	test_curve_Fq2::Frobenius_coeffs_c1[1] = test_curve_Fq("258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458176");
	
     /* test_curve Fq6 parameters */

	test_curve_Fq6::non_residue = test_curve_Fq2(test_curve_Fq("0"),test_curve_Fq("1"));
	test_curve_Fq6::Frobenius_coeffs_c1[0] = test_curve_Fq2(test_curve_Fq("1"),test_curve_Fq("0"));
	test_curve_Fq6::Frobenius_coeffs_c1[1] = test_curve_Fq2(test_curve_Fq("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410946"),test_curve_Fq("0"));
	test_curve_Fq6::Frobenius_coeffs_c1[2] = test_curve_Fq2(test_curve_Fq("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945"),test_curve_Fq("0"));
	test_curve_Fq6::Frobenius_coeffs_c1[3] = test_curve_Fq2(test_curve_Fq("258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458176"),test_curve_Fq("0"));
	test_curve_Fq6::Frobenius_coeffs_c1[4] = test_curve_Fq2(test_curve_Fq("258664426012969093929703085429980814127835149614277183275038967946009968870203535512256352201271898244626862047231"),test_curve_Fq("0"));
	test_curve_Fq6::Frobenius_coeffs_c1[5] = test_curve_Fq2(test_curve_Fq("258664426012969093929703085429980814127835149614277183275038967946009968870203535512256352201271898244626862047232"),test_curve_Fq("0"));
	test_curve_Fq6::Frobenius_coeffs_c2[0] = test_curve_Fq2(test_curve_Fq("1"),test_curve_Fq("0"));
	test_curve_Fq6::Frobenius_coeffs_c2[1] = test_curve_Fq2(test_curve_Fq("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945"),test_curve_Fq("0"));
	test_curve_Fq6::Frobenius_coeffs_c2[2] = test_curve_Fq2(test_curve_Fq("258664426012969093929703085429980814127835149614277183275038967946009968870203535512256352201271898244626862047231"),test_curve_Fq("0"));
	test_curve_Fq6::Frobenius_coeffs_c2[3] = test_curve_Fq2(test_curve_Fq("1"),test_curve_Fq("0"));
	test_curve_Fq6::Frobenius_coeffs_c2[4] = test_curve_Fq2(test_curve_Fq("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945"),test_curve_Fq("0"));
	test_curve_Fq6::Frobenius_coeffs_c2[5] = test_curve_Fq2(test_curve_Fq("258664426012969093929703085429980814127835149614277183275038967946009968870203535512256352201271898244626862047231"),test_curve_Fq("0"));
	
     /* test_curve Fq12 parameters */

	test_curve_Fq12::non_residue = test_curve_Fq2(test_curve_Fq("0"),test_curve_Fq("1"));
	test_curve_Fq12::Frobenius_coeffs_c1[0] = test_curve_Fq2(test_curve_Fq("1"),test_curve_Fq("0"));
	test_curve_Fq12::Frobenius_coeffs_c1[1] = test_curve_Fq2(test_curve_Fq("92949345220277864758624960506473182677953048909283248980960104381795901929519566951595905490535835115111760994353"),test_curve_Fq("0"));
	test_curve_Fq12::Frobenius_coeffs_c1[2] = test_curve_Fq2(test_curve_Fq("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410946"),test_curve_Fq("0"));
	test_curve_Fq12::Frobenius_coeffs_c1[3] = test_curve_Fq2(test_curve_Fq("216465761340224619389371505802605247630151569547285782856803747159100223055385581585702401816380679166954762214499"),test_curve_Fq("0"));
	test_curve_Fq12::Frobenius_coeffs_c1[4] = test_curve_Fq2(test_curve_Fq("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945"),test_curve_Fq("0"));
	test_curve_Fq12::Frobenius_coeffs_c1[5] = test_curve_Fq2(test_curve_Fq("123516416119946754630746545296132064952198520638002533875843642777304321125866014634106496325844844051843001220146"),test_curve_Fq("0"));
	test_curve_Fq12::Frobenius_coeffs_c1[6] = test_curve_Fq2(test_curve_Fq("258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458176"),test_curve_Fq("0"));
	test_curve_Fq12::Frobenius_coeffs_c1[7] = test_curve_Fq2(test_curve_Fq("165715080792691229252027773188420350858440463845631411558924158284924566418821255823372982649037525009328560463824"),test_curve_Fq("0"));
	test_curve_Fq12::Frobenius_coeffs_c1[8] = test_curve_Fq2(test_curve_Fq("258664426012969093929703085429980814127835149614277183275038967946009968870203535512256352201271898244626862047231"),test_curve_Fq("0"));
	test_curve_Fq12::Frobenius_coeffs_c1[9] = test_curve_Fq2(test_curve_Fq("42198664672744474621281227892288285906241943207628877683080515507620245292955241189266486323192680957485559243678"),test_curve_Fq("0"));
	test_curve_Fq12::Frobenius_coeffs_c1[10] = test_curve_Fq2(test_curve_Fq("258664426012969093929703085429980814127835149614277183275038967946009968870203535512256352201271898244626862047232"),test_curve_Fq("0"));
	test_curve_Fq12::Frobenius_coeffs_c1[11] = test_curve_Fq2(test_curve_Fq("135148009893022339379906188398761468584194992116912126664040619889416147222474808140862391813728516072597320238031"),test_curve_Fq("0"));
	
    /* choice of short Weierstrass curve and its twist */

	test_curve_coeff_b = test_curve_Fq("1");
	test_curve_twist = test_curve_Fq2(test_curve_Fq("0"),test_curve_Fq("1"));
	test_curve_twist_coeff_b = test_curve_Fq2(test_curve_Fq("0"),test_curve_Fq("103465770405187637604261093477957413414557405101965864215953705066688187339336329109987555255829344049776128583271"));
	test_curve_twist_mul_by_b_c0 = test_curve_coeff_b * test_curve_Fq2::non_residue;
	test_curve_twist_mul_by_b_c1 = test_curve_coeff_b * test_curve_Fq2::non_residue;
	test_curve_twist_mul_by_q_X = test_curve_Fq2(test_curve_Fq("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410946"),test_curve_Fq("0"));
	test_curve_twist_mul_by_q_Y = test_curve_Fq2(test_curve_Fq("216465761340224619389371505802605247630151569547285782856803747159100223055385581585702401816380679166954762214499"),test_curve_Fq("0"));
	
    /* choice of group G1 */

	test_curve_G1::G1_zero = test_curve_G1(test_curve_Fq::zero(),   
                                                   test_curve_Fq::one(),
                                                   test_curve_Fq::zero());
	test_curve_G1::G1_one = test_curve_G1(test_curve_Fq("81937999373150964239938255573465948239988671502647976594219695644855304257327692006745978603320413799295628339695"),
                                                  test_curve_Fq("241266749859715473739788878240585681733927191168601896383759122102112907357779751001206799952863815012735208165030"),
                                                  test_curve_Fq::one());

    // TODO
	test_curve_G1::wnaf_window_table.resize(0);
    test_curve_G1::wnaf_window_table.push_back(11);
    test_curve_G1::wnaf_window_table.push_back(24);
    test_curve_G1::wnaf_window_table.push_back(60);
    test_curve_G1::wnaf_window_table.push_back(127);

    test_curve_G1::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.99]
    test_curve_G1::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.99, 10.99]
    test_curve_G1::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.99, 32.29]
    test_curve_G1::fixed_base_exp_window_table.push_back(11);
    // window 4 is unbeaten in [32.29, 55.23]
    test_curve_G1::fixed_base_exp_window_table.push_back(32);
    // window 5 is unbeaten in [55.23, 162.03]
    test_curve_G1::fixed_base_exp_window_table.push_back(55);
    // window 6 is unbeaten in [162.03, 360.15]
    test_curve_G1::fixed_base_exp_window_table.push_back(162);
    // window 7 is unbeaten in [360.15, 815.44]
    test_curve_G1::fixed_base_exp_window_table.push_back(360);
    // window 8 is unbeaten in [815.44, 2373.07]
    test_curve_G1::fixed_base_exp_window_table.push_back(815);
    // window 9 is unbeaten in [2373.07, 6977.75]
    test_curve_G1::fixed_base_exp_window_table.push_back(2373);
    // window 10 is unbeaten in [6977.75, 7122.23]
    test_curve_G1::fixed_base_exp_window_table.push_back(6978);
    // window 11 is unbeaten in [7122.23, 57818.46]
    test_curve_G1::fixed_base_exp_window_table.push_back(7122);
    // window 12 is never the best
    test_curve_G1::fixed_base_exp_window_table.push_back(0);
    // window 13 is unbeaten in [57818.46, 169679.14]
    test_curve_G1::fixed_base_exp_window_table.push_back(57818);
    // window 14 is never the best
    test_curve_G1::fixed_base_exp_window_table.push_back(0);
    // window 15 is unbeaten in [169679.14, 439758.91]
    test_curve_G1::fixed_base_exp_window_table.push_back(169679);
    // window 16 is unbeaten in [439758.91, 936073.41]
    test_curve_G1::fixed_base_exp_window_table.push_back(439759);
    // window 17 is unbeaten in [936073.41, 4666554.74]
    test_curve_G1::fixed_base_exp_window_table.push_back(936073);
    // window 18 is never the best
    test_curve_G1::fixed_base_exp_window_table.push_back(0);
    // window 19 is unbeaten in [4666554.74, 7580404.42]
    test_curve_G1::fixed_base_exp_window_table.push_back(4666555);
    // window 20 is unbeaten in [7580404.42, 34552892.20]
    test_curve_G1::fixed_base_exp_window_table.push_back(7580404);
    // window 21 is never the best
    test_curve_G1::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [34552892.20, inf]
    test_curve_G1::fixed_base_exp_window_table.push_back(34552892);

    /* choice of group G2 */

	test_curve_G2::G2_zero = test_curve_G2(test_curve_Fq2::zero(),
                                 test_curve_Fq2::one(),
                                 test_curve_Fq2::zero());
	test_curve_G2::G2_one = test_curve_G2(test_curve_Fq2(test_curve_Fq("229333750276549292773392407107592235223117097523216664105847024235830162889413156552545986841119901481759443034658"),
                                                                     test_curve_Fq("138488675805149526409586699643852536836155818985450873344679280398096052441132846915951790515164142044356717177747")),
                                                  test_curve_Fq2(test_curve_Fq("46928676784594546038825359633204275549652239625881023595005107328853903064771615718508047529230056308887605558782"),
                                                                     test_curve_Fq("199044010096646683085962364522264780781592913855437531311176701976142857176410696121640365883643668237137795880430")),
                                                  test_curve_Fq2::one());
    // TODO
	test_curve_G2::wnaf_window_table.resize(0);
    test_curve_G2::wnaf_window_table.push_back(5);
    test_curve_G2::wnaf_window_table.push_back(15);
    test_curve_G2::wnaf_window_table.push_back(39);
    test_curve_G2::wnaf_window_table.push_back(109);

    test_curve_G2::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 5.10]
    test_curve_G2::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [5.10, 10.43]
    test_curve_G2::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.43, 25.28]
    test_curve_G2::fixed_base_exp_window_table.push_back(10);
    // window 4 is unbeaten in [25.28, 59.00]
    test_curve_G2::fixed_base_exp_window_table.push_back(25);
    // window 5 is unbeaten in [59.00, 154.03]
    test_curve_G2::fixed_base_exp_window_table.push_back(59);
    // window 6 is unbeaten in [154.03, 334.25]
    test_curve_G2::fixed_base_exp_window_table.push_back(154);
    // window 7 is unbeaten in [334.25, 742.58]
    test_curve_G2::fixed_base_exp_window_table.push_back(334);
    // window 8 is unbeaten in [742.58, 2034.40]
    test_curve_G2::fixed_base_exp_window_table.push_back(743);
    // window 9 is unbeaten in [2034.40, 4987.56]
    test_curve_G2::fixed_base_exp_window_table.push_back(2034);
    // window 10 is unbeaten in [4987.56, 8888.27]
    test_curve_G2::fixed_base_exp_window_table.push_back(4988);
    // window 11 is unbeaten in [8888.27, 26271.13]
    test_curve_G2::fixed_base_exp_window_table.push_back(8888);
    // window 12 is unbeaten in [26271.13, 39768.20]
    test_curve_G2::fixed_base_exp_window_table.push_back(26271);
    // window 13 is unbeaten in [39768.20, 106275.75]
    test_curve_G2::fixed_base_exp_window_table.push_back(39768);
    // window 14 is unbeaten in [106275.75, 141703.40]
    test_curve_G2::fixed_base_exp_window_table.push_back(106276);
    // window 15 is unbeaten in [141703.40, 462422.97]
    test_curve_G2::fixed_base_exp_window_table.push_back(141703);
    // window 16 is unbeaten in [462422.97, 926871.84]
    test_curve_G2::fixed_base_exp_window_table.push_back(462423);
    // window 17 is unbeaten in [926871.84, 4873049.17]
    test_curve_G2::fixed_base_exp_window_table.push_back(926872);
    // window 18 is never the best
    test_curve_G2::fixed_base_exp_window_table.push_back(0);
    // window 19 is unbeaten in [4873049.17, 5706707.88]
    test_curve_G2::fixed_base_exp_window_table.push_back(4873049);
    // window 20 is unbeaten in [5706707.88, 31673814.95]
    test_curve_G2::fixed_base_exp_window_table.push_back(5706708);
    // window 21 is never the best
    test_curve_G2::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [31673814.95, inf]
    test_curve_G2::fixed_base_exp_window_table.push_back(31673815);

    /* choice of pairing */

	test_curve_ate_loop_count = bigint_q("9586122913090633729");
	test_curve_ate_is_loop_count_neg = false;
	test_curve_final_exponent = bigint<12*test_curve_q_limbs>("10623521018019860488254031663707568428798032905123811199571213965079129114663661236359849629341526275899063345613340067081670062620727617884137487754739150147491204559514205186492385590272208934467461444944652711005169371168250068790820776124772095630237102189827733019989835063334551453893534663070786533932633573962932272563471643288531959637300817070265537429506484880990981069041269405383502889677357082012807298529931118124428569059822346289745077401570134157444973271520981774047146918354408632568723153146248333028827919406785654402107153546667815607201488590832478225403444136409349877481268154817904541340614173261949772403060924324366861723245182619859389254985008236007465814273361497134138868945580557938161335670207544906643574043606819537336472235809927599628123275314288006170804044560238676463931639339711913111080974582593228138704154320599775683095604041309000197025419968125718018311805959315220036948621879242495199408833915486421612374480018459896018440926235261824654956932384859260479372776022979736734221629097297890154692194441528462770218811795624471108972377573690833913231260547835550851256817740247389770320334698430697237343583761719223414894063451411431859122738488311580005412765070251810159991897110936324943232526870280724876946523218213525646968094720");
	test_curve_final_exponent_z = bigint_q("9586122913090633729");
	test_curve_final_exponent_is_z_neg = false;

}
} // libff
