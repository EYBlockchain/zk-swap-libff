#include <libff/algebra/curves/pendulum/pendulum_g1.hpp>
#include <libff/algebra/curves/pendulum/pendulum_g2.hpp>
#include <libff/algebra/curves/pendulum/pendulum_init.hpp>

namespace libff {

bigint<pendulum_r_limbs> pendulum_modulus_r;
bigint<pendulum_q_limbs> pendulum_modulus_q;

pendulum_Fq pendulum_coeff_a;
pendulum_Fq pendulum_coeff_b;
pendulum_Fq3 pendulum_twist;
pendulum_Fq3 pendulum_twist_coeff_a;
pendulum_Fq3 pendulum_twist_coeff_b;
pendulum_Fq pendulum_twist_mul_by_a_c0;
pendulum_Fq pendulum_twist_mul_by_a_c1;
pendulum_Fq pendulum_twist_mul_by_a_c2;
pendulum_Fq pendulum_twist_mul_by_b_c0;
pendulum_Fq pendulum_twist_mul_by_b_c1;
pendulum_Fq pendulum_twist_mul_by_b_c2;
pendulum_Fq pendulum_twist_mul_by_q_X;
pendulum_Fq pendulum_twist_mul_by_q_Y;

bigint<pendulum_q_limbs> pendulum_ate_loop_count;
bool pendulum_ate_is_loop_count_neg;
bigint<6*pendulum_q_limbs> pendulum_final_exponent;
bigint<pendulum_q_limbs> pendulum_final_exponent_last_chunk_abs_of_w0;
bool pendulum_final_exponent_last_chunk_is_w0_neg;
bigint<pendulum_q_limbs> pendulum_final_exponent_last_chunk_w1;


void init_pendulum_params()
{
    typedef bigint<pendulum_r_limbs> bigint_r;
    typedef bigint<pendulum_q_limbs> bigint_q;

    assert(sizeof(mp_limb_t) == 8 || sizeof(mp_limb_t) == 4); // Montgomery assumes this

    /* parameters for scalar field Fr */
    pendulum_modulus_r = bigint_r("475922286169261325753349249653048451545124878552823515553267735739164647307408490559963137");
    assert(pendulum_Fr::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        pendulum_Fr::Rsquared = bigint_r("163983144722506446826715124368972380525894397127205577781234305496325861831001705438796139");
        pendulum_Fr::Rcubed = bigint_r("207236281459091063710247635236340312578688659363066707916716212805695955118593239854980171");
        pendulum_Fr::inv = 0xbb4334a3ffffffff;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        pendulum_Fr::Rsquared = bigint_r("163983144722506446826715124368972380525894397127205577781234305496325861831001705438796139");
        pendulum_Fr::Rcubed = bigint_r("207236281459091063710247635236340312578688659363066707916716212805695955118593239854980171");
        pendulum_Fr::inv = 0xffffffff;
    }
    pendulum_Fr::num_bits = 298;
    pendulum_Fr::euler = bigint_r("237961143084630662876674624826524225772562439276411757776633867869582323653704245279981568");
    pendulum_Fr::s = 34;
    pendulum_Fr::t = bigint_r("27702323054502562488973446286577291993024111641153199339359284829066871159442729");
    pendulum_Fr::t_minus_1_over_2 = bigint_r("13851161527251281244486723143288645996512055820576599669679642414533435579721364");
    pendulum_Fr::multiplicative_generator = pendulum_Fr("10");
    pendulum_Fr::root_of_unity = pendulum_Fr("120638817826913173458768829485690099845377008030891618010109772937363554409782252579816313");
    pendulum_Fr::nqr = pendulum_Fr("5");
    pendulum_Fr::nqr_to_t = pendulum_Fr("406220604243090401056429458730298145937262552508985450684842547562990900634752279902740880");
  
    /* parameters for base field Fq */
    pendulum_modulus_q = bigint_q("19050022797317891600939264904924934656417895081121634056186244048763811669585984032184028629480644260294123843823582617865870693473572190965725707704312821545976965077621486794922414287");
    assert(pendulum_Fq::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        pendulum_Fq::Rsquared = bigint_q("8578827046306886914027642114663788257945766401669528758478475985831297403454046497330331843669617003368544670659754777554745652937396507224870436411421385726201339868370359556399726891"); // k=10 
        pendulum_Fq::Rcubed = bigint_q("3711622197431008092504984114265758823383736404180765249918451588406242080653182878515083753572244107488387081156702012017861504011357458301066605577664563595709701483509240497734732206");

        pendulum_Fq::inv = 0x3db7b35378c349d1;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        pendulum_Fq::Rsquared = bigint_q("8578827046306886914027642114663788257945766401669528758478475985831297403454046497330331843669617003368544670659754777554745652937396507224870436411421385726201339868370359556399726891");
        pendulum_Fq::Rcubed = bigint_q("3711622197431008092504984114265758823383736404180765249918451588406242080653182878515083753572244107488387081156702012017861504011357458301066605577664563595709701483509240497734732206");
        pendulum_Fq::inv = 0x78c349d1;
    }
    pendulum_Fq::num_bits = 613;
    pendulum_Fq::euler = bigint_q("9525011398658945800469632452462467328208947540560817028093122024381905834792992016092014314740322130147061921911791308932935346736786095482862853852156410772988482538810743397461207143");
    pendulum_Fq::s = 1;
    pendulum_Fq::t = bigint_q("9525011398658945800469632452462467328208947540560817028093122024381905834792992016092014314740322130147061921911791308932935346736786095482862853852156410772988482538810743397461207143");
    pendulum_Fq::t_minus_1_over_2 = bigint_q("4762505699329472900234816226231233664104473770280408514046561012190952917396496008046007157370161065073530960955895654466467673368393047741431426926078205386494241269405371698730603571");
    pendulum_Fq::multiplicative_generator = pendulum_Fq("3");
    pendulum_Fq::root_of_unity = pendulum_Fq("19050022797317891600939264904924934656417895081121634056186244048763811669585984032184028629480644260294123843823582617865870693473572190965725707704312821545976965077621486794922414286");
    pendulum_Fq::nqr = pendulum_Fq("3");
    pendulum_Fq::nqr_to_t = pendulum_Fq("19050022797317891600939264904924934656417895081121634056186244048763811669585984032184028629480644260294123843823582617865870693473572190965725707704312821545976965077621486794922414286");

    /* parameters for twist field Fq3 */
    pendulum_Fq3::euler = bigint<3*pendulum_q_limbs>("3456658722320335174353551207932785579644752512910401406680016736181084292870281340317233717232774313005965367789218718893970022322131118477117545529247397233536018152366440805049179283934783621509201378722205089278922068612054841699043112497084824541424223630998363455549086211925192885477366859918144938687339934699921603675415165657156051498271739873235479474457281360122254082709882805427591951436821089707675571700256304510341212328830846433132430016893957416054679062299747314098566512239537008235914634050155823111780235496289082751225792467968951");
    pendulum_Fq3::s = 1;
    pendulum_Fq3::t = bigint<3*pendulum_q_limbs>("3456658722320335174353551207932785579644752512910401406680016736181084292870281340317233717232774313005965367789218718893970022322131118477117545529247397233536018152366440805049179283934783621509201378722205089278922068612054841699043112497084824541424223630998363455549086211925192885477366859918144938687339934699921603675415165657156051498271739873235479474457281360122254082709882805427591951436821089707675571700256304510341212328830846433132430016893957416054679062299747314098566512239537008235914634050155823111780235496289082751225792467968951");
    pendulum_Fq3::t_minus_1_over_2 = bigint<3*pendulum_q_limbs>("1728329361160167587176775603966392789822376256455200703340008368090542146435140670158616858616387156502982683894609359446985011161065559238558772764623698616768009076183220402524589641967391810754600689361102544639461034306027420849521556248542412270712111815499181727774543105962596442738683429959072469343669967349960801837707582828578025749135869936617739737228640680061127041354941402713795975718410544853837785850128152255170606164415423216566215008446978708027339531149873657049283256119768504117957317025077911555890117748144541375612896233984475");
    pendulum_Fq3::non_residue = pendulum_Fq("3");
    pendulum_Fq3::nqr = pendulum_Fq3(pendulum_Fq("0"),pendulum_Fq("1"),pendulum_Fq("0"));
    pendulum_Fq3::nqr_to_t = pendulum_Fq3(pendulum_Fq("19050022797317891600939264904924934656417895081121634056186244048763811669585984032184028629480644260294123843823582617865870693473572190965725707704312821545976965077621486794922414286"),pendulum_Fq("0"),pendulum_Fq("0"));
    pendulum_Fq3::Frobenius_coeffs_c1[0] = pendulum_Fq("1");
    pendulum_Fq3::Frobenius_coeffs_c1[1] = pendulum_Fq("17900163747628920354847454606245864377330997902231447778281388313919659811991993446336892445839751089711383680525708491913537773742066921127626000186860542518496145894403248164444394781");
    pendulum_Fq3::Frobenius_coeffs_c1[2] = pendulum_Fq("1149859049688971246091810298679070279086897178890186277904855734844151857593990585847136183640893170582740163297874125952332919731505269838099707517452279027480819183218238630478019505");
    pendulum_Fq3::Frobenius_coeffs_c2[0] = pendulum_Fq("1");
    pendulum_Fq3::Frobenius_coeffs_c2[1] = pendulum_Fq("1149859049688971246091810298679070279086897178890186277904855734844151857593990585847136183640893170582740163297874125952332919731505269838099707517452279027480819183218238630478019505");
    pendulum_Fq3::Frobenius_coeffs_c2[2] = pendulum_Fq("17900163747628920354847454606245864377330997902231447778281388313919659811991993446336892445839751089711383680525708491913537773742066921127626000186860542518496145894403248164444394781");


    /* parameters for Fq6 */
    pendulum_Fq6::non_residue = pendulum_Fq("3");
    pendulum_Fq6::Frobenius_coeffs_c1[0] = pendulum_Fq("1");
    pendulum_Fq6::Frobenius_coeffs_c1[1] = pendulum_Fq("17900163747628920354847454606245864377330997902231447778281388313919659811991993446336892445839751089711383680525708491913537773742066921127626000186860542518496145894403248164444394782");
    pendulum_Fq6::Frobenius_coeffs_c1[2] = pendulum_Fq("17900163747628920354847454606245864377330997902231447778281388313919659811991993446336892445839751089711383680525708491913537773742066921127626000186860542518496145894403248164444394781");
    pendulum_Fq6::Frobenius_coeffs_c1[3] = pendulum_Fq("19050022797317891600939264904924934656417895081121634056186244048763811669585984032184028629480644260294123843823582617865870693473572190965725707704312821545976965077621486794922414286");
    pendulum_Fq6::Frobenius_coeffs_c1[4] = pendulum_Fq("1149859049688971246091810298679070279086897178890186277904855734844151857593990585847136183640893170582740163297874125952332919731505269838099707517452279027480819183218238630478019505");
    pendulum_Fq6::Frobenius_coeffs_c1[5] = pendulum_Fq("1149859049688971246091810298679070279086897178890186277904855734844151857593990585847136183640893170582740163297874125952332919731505269838099707517452279027480819183218238630478019506");
    pendulum_Fq6::my_Fp2::non_residue = pendulum_Fq3::non_residue;

    /* choice of short Weierstrass curve and its twist */
    pendulum_coeff_a = pendulum_Fq("0");
    pendulum_coeff_b = pendulum_Fq("3779136");
    pendulum_twist = pendulum_Fq3(pendulum_Fq::zero(), pendulum_Fq::one(), pendulum_Fq::zero()); // from zexe
    pendulum_twist_coeff_a = pendulum_Fq3(pendulum_Fq::zero(), pendulum_Fq::zero(), pendulum_G1::coeff_a);
    pendulum_twist_coeff_b = pendulum_Fq3(pendulum_G1::coeff_b * pendulum_Fq3::non_residue, pendulum_Fq::zero(), pendulum_Fq::zero());
    pendulum_G2::twist = pendulum_twist;
    pendulum_G2::coeff_a = pendulum_twist_coeff_a;
    pendulum_G2::coeff_b = pendulum_twist_coeff_b;
    pendulum_twist_mul_by_a_c0 = pendulum_G1::coeff_a * pendulum_Fq3::non_residue;
    pendulum_twist_mul_by_a_c1 = pendulum_G1::coeff_a * pendulum_Fq3::non_residue;
    pendulum_twist_mul_by_a_c2 = pendulum_G1::coeff_a;
    pendulum_twist_mul_by_b_c0 = pendulum_G1::coeff_b * pendulum_Fq3::non_residue;
    pendulum_twist_mul_by_b_c1 = pendulum_G1::coeff_b * pendulum_Fq3::non_residue;
    pendulum_twist_mul_by_b_c2 = pendulum_G1::coeff_b * pendulum_Fq3::non_residue;
    pendulum_twist_mul_by_q_X = pendulum_Fq("1149859049688971246091810298679070279086897178890186277904855734844151857593990585847136183640893170582740163297874125952332919731505269838099707517452279027480819183218238630478019505");
    pendulum_twist_mul_by_q_Y = pendulum_Fq("19050022797317891600939264904924934656417895081121634056186244048763811669585984032184028629480644260294123843823582617865870693473572190965725707704312821545976965077621486794922414286");


    /* choice of group G1 */
    pendulum_G1::G1_zero = pendulum_G1(pendulum_Fq::zero(),
                                       pendulum_Fq::one(),
                                       pendulum_Fq::zero());
    pendulum_G1::G1_one = pendulum_G1(pendulum_Fq("10429529130963884009088672788347332381112680208340155261267801535021515284160681139607251964825852124795332292830211013740612118630760054643157598813437729109347416973526597847768430196"),
                                      pendulum_Fq("16905957489427084020812201420874004139445860435619110387337972510865682281055685071559898890725067816998514457554471177636644195808052025100386087604283564319339703394765085347212298098"),
                                      pendulum_Fq::one());


    // TODO: wNAF window table
    // pendulum_G1::wnaf_window_table.resize(0);
    // pendulum_G1::wnaf_window_table.push_back(11);
    // pendulum_G1::wnaf_window_table.push_back(24);
    // pendulum_G1::wnaf_window_table.push_back(60);
    // pendulum_G1::wnaf_window_table.push_back(127);

    // // TODO: fixed-base exponentiation table
    // pendulum_G1::fixed_base_exp_window_table.resize(0);
    // // window 1 is unbeaten in [-inf, 4.99]
    // pendulum_G1::fixed_base_exp_window_table.push_back(1);
    // // window 2 is unbeaten in [4.99, 10.99]
    // pendulum_G1::fixed_base_exp_window_table.push_back(5);
    // // window 3 is unbeaten in [10.99, 32.29]
    // pendulum_G1::fixed_base_exp_window_table.push_back(11);
    // // window 4 is unbeaten in [32.29, 55.23]
    // pendulum_G1::fixed_base_exp_window_table.push_back(32);
    // // window 5 is unbeaten in [55.23, 162.03]
    // pendulum_G1::fixed_base_exp_window_table.push_back(55);
    // // window 6 is unbeaten in [162.03, 360.15]
    // pendulum_G1::fixed_base_exp_window_table.push_back(162);
    // // window 7 is unbeaten in [360.15, 815.44]
    // pendulum_G1::fixed_base_exp_window_table.push_back(360);
    // // window 8 is unbeaten in [815.44, 2373.07]
    // pendulum_G1::fixed_base_exp_window_table.push_back(815);
    // // window 9 is unbeaten in [2373.07, 6977.75]
    // pendulum_G1::fixed_base_exp_window_table.push_back(2373);
    // // window 10 is unbeaten in [6977.75, 7122.23]
    // pendulum_G1::fixed_base_exp_window_table.push_back(6978);
    // // window 11 is unbeaten in [7122.23, 57818.46]
    // pendulum_G1::fixed_base_exp_window_table.push_back(7122);
    // // window 12 is never the best
    // pendulum_G1::fixed_base_exp_window_table.push_back(0);
    // // window 13 is unbeaten in [57818.46, 169679.14]
    // pendulum_G1::fixed_base_exp_window_table.push_back(57818);
    // // window 14 is never the best
    // pendulum_G1::fixed_base_exp_window_table.push_back(0);
    // // window 15 is unbeaten in [169679.14, 439758.91]
    // pendulum_G1::fixed_base_exp_window_table.push_back(169679);
    // // window 16 is unbeaten in [439758.91, 936073.41]
    // pendulum_G1::fixed_base_exp_window_table.push_back(439759);
    // // window 17 is unbeaten in [936073.41, 4666554.74]
    // pendulum_G1::fixed_base_exp_window_table.push_back(936073);
    // // window 18 is never the best
    // pendulum_G1::fixed_base_exp_window_table.push_back(0);
    // // window 19 is unbeaten in [4666554.74, 7580404.42]
    // pendulum_G1::fixed_base_exp_window_table.push_back(4666555);
    // // window 20 is unbeaten in [7580404.42, 34552892.20]
    // pendulum_G1::fixed_base_exp_window_table.push_back(7580404);
    // // window 21 is never the best
    // pendulum_G1::fixed_base_exp_window_table.push_back(0);
    // // window 22 is unbeaten in [34552892.20, inf]
    // pendulum_G1::fixed_base_exp_window_table.push_back(34552892);


    /* choice of group G2 */
    pendulum_G2::G2_zero = pendulum_G2(pendulum_Fq3::zero(),
                             pendulum_Fq3::one(),
                             pendulum_Fq3::zero());

    // simple G2 generator
    pendulum_G2::G2_zero = pendulum_G2(pendulum_Fq3::zero(),
                               pendulum_Fq3::one(),
                               pendulum_Fq3::zero());
    pendulum_G2::G2_one = pendulum_G2(pendulum_Fq3(pendulum_Fq("13912037883745548354885384080735939137135972624544456324778703816799088874865621783205940029962037178108763091296730317064882342261056217993081826235556952843127200569255758772006603222"),
                                    pendulum_Fq("0"),
                                    pendulum_Fq("0")),
                            pendulum_Fq3(pendulum_Fq("17848446077687999456842554535611770905438459561283772120727453120807955998955357528566232686711238095565974536067108118163952231407922320805847265871073335951065323478695430783920346311"),
                                    pendulum_Fq("0"),
                                    pendulum_Fq("0")),
                            pendulum_Fq3::one());


    // // TODO: wNAF window table
    // pendulum_G2::wnaf_window_table.resize(0);
    // pendulum_G2::wnaf_window_table.push_back(5);
    // pendulum_G2::wnaf_window_table.push_back(15);
    // pendulum_G2::wnaf_window_table.push_back(39);
    // pendulum_G2::wnaf_window_table.push_back(109);

    // // TODO: fixed-base exponentiation table 
    // pendulum_G2::fixed_base_exp_window_table.resize(0);
    // // window 1 is unbeaten in [-inf, 5.10]
    // pendulum_G2::fixed_base_exp_window_table.push_back(1);
    // // window 2 is unbeaten in [5.10, 10.43]
    // pendulum_G2::fixed_base_exp_window_table.push_back(5);
    // // window 3 is unbeaten in [10.43, 25.28]
    // pendulum_G2::fixed_base_exp_window_table.push_back(10);
    // // window 4 is unbeaten in [25.28, 59.00]
    // pendulum_G2::fixed_base_exp_window_table.push_back(25);
    // // window 5 is unbeaten in [59.00, 154.03]
    // pendulum_G2::fixed_base_exp_window_table.push_back(59);
    // // window 6 is unbeaten in [154.03, 334.25]
    // pendulum_G2::fixed_base_exp_window_table.push_back(154);
    // // window 7 is unbeaten in [334.25, 742.58]
    // pendulum_G2::fixed_base_exp_window_table.push_back(334);
    // // window 8 is unbeaten in [742.58, 2034.40]
    // pendulum_G2::fixed_base_exp_window_table.push_back(743);
    // // window 9 is unbeaten in [2034.40, 4987.56]
    // pendulum_G2::fixed_base_exp_window_table.push_back(2034);
    // // window 10 is unbeaten in [4987.56, 8888.27]
    // pendulum_G2::fixed_base_exp_window_table.push_back(4988);
    // // window 11 is unbeaten in [8888.27, 26271.13]
    // pendulum_G2::fixed_base_exp_window_table.push_back(8888);
    // // window 12 is unbeaten in [26271.13, 39768.20]
    // pendulum_G2::fixed_base_exp_window_table.push_back(26271);
    // // window 13 is unbeaten in [39768.20, 106275.75]
    // pendulum_G2::fixed_base_exp_window_table.push_back(39768);
    // // window 14 is unbeaten in [106275.75, 141703.40]
    // pendulum_G2::fixed_base_exp_window_table.push_back(106276);
    // // window 15 is unbeaten in [141703.40, 462422.97]
    // pendulum_G2::fixed_base_exp_window_table.push_back(141703);
    // // window 16 is unbeaten in [462422.97, 926871.84]
    // pendulum_G2::fixed_base_exp_window_table.push_back(462423);
    // // window 17 is unbeaten in [926871.84, 4873049.17]
    // pendulum_G2::fixed_base_exp_window_table.push_back(926872);
    // // window 18 is never the best
    // pendulum_G2::fixed_base_exp_window_table.push_back(0);
    // // window 19 is unbeaten in [4873049.17, 5706707.88]
    // pendulum_G2::fixed_base_exp_window_table.push_back(4873049);
    // // window 20 is unbeaten in [5706707.88, 31673814.95]
    // pendulum_G2::fixed_base_exp_window_table.push_back(5706708);
    // // window 21 is never the best
    // pendulum_G2::fixed_base_exp_window_table.push_back(0);
    // // window 22 is unbeaten in [31673814.95, inf]
    // pendulum_G2::fixed_base_exp_window_table.push_back(31673815);



    /* pairing parameters */
    pendulum_ate_loop_count = bigint_q("506464946133393486072777102926336625944849939610982267859828541006717966526573193706126370441346337661774335955699621");
    pendulum_ate_is_loop_count_neg = true;
    pendulum_final_exponent = bigint<6*pendulum_q_limbs>("100423870617765377462288921991705058064066317441544284495241123579182211596199015851611264231095298339998088763771082746874787805299346743855491714431237761533371754609061885918820013688992199567176740528257075196642520203571395720349075208345116351733410271370221200176246986072404830843366450512893766369356307021152474431528599741498202226004015201673575698774526383309793439778791410813617073875322701411074298623735325647821678998811323097221570247014577699527521316660579585310467502764262886643629279956136589108065749949275502904571054732187489208977384998875130114367813643869701439127571469962838383471719004134141664250354239990346715730264937727465514553571373542401086719861657200622270483688837039467865500870083321548750233771612931970133032178527755565626224514504506456192754852318102059860644159431186965577474145619581368938985932531405729469296392643918817935332526021874413186293408331487941093239399997803417973633666435650866228341242449813102571910623252432186420977607250534445625335582474784");
    pendulum_final_exponent_last_chunk_abs_of_w0 = bigint_q("7000705447348627246181409558336018323010329260726930841638672011287206690002601216854775649561085256265269640040570922609783227469279331691880282815325569032149343779036142830666859805506518426649197067288711084398033");
    pendulum_final_exponent_last_chunk_is_w0_neg = true;
    pendulum_final_exponent_last_chunk_w1 = bigint_q("86482221941698704497288378992285180119495364068003923046442785886272123124361700722982503222189455144364945735564951562986");
}
} // libff
