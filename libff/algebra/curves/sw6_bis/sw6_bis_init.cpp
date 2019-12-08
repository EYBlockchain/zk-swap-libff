#include <libff/algebra/curves/sw6_bis/sw6_bis_g1.hpp>
#include <libff/algebra/curves/sw6_bis/sw6_bis_g2.hpp>
#include <libff/algebra/curves/sw6_bis/sw6_bis_init.hpp>

namespace libff {

bigint<sw6_bis_r_limbs> sw6_bis_modulus_r;
bigint<sw6_bis_q_limbs> sw6_bis_modulus_q;

sw6_bis_Fq sw6_bis_coeff_a;
sw6_bis_Fq sw6_bis_coeff_b;
sw6_bis_Fq3 sw6_bis_twist;
sw6_bis_Fq3 sw6_bis_twist_coeff_a;
sw6_bis_Fq3 sw6_bis_twist_coeff_b;
sw6_bis_Fq sw6_bis_twist_mul_by_a_c0;
sw6_bis_Fq sw6_bis_twist_mul_by_a_c1;
sw6_bis_Fq sw6_bis_twist_mul_by_a_c2;
sw6_bis_Fq sw6_bis_twist_mul_by_b_c0;
sw6_bis_Fq sw6_bis_twist_mul_by_b_c1;
sw6_bis_Fq sw6_bis_twist_mul_by_b_c2;
sw6_bis_Fq sw6_bis_twist_mul_by_q_X;
sw6_bis_Fq sw6_bis_twist_mul_by_q_Y;

bigint<sw6_bis_q_limbs> sw6_bis_ate_loop_count;
bool sw6_bis_ate_is_loop_count_neg;
bigint<6*sw6_bis_q_limbs> sw6_bis_final_exponent;
bigint<sw6_bis_q_limbs> sw6_bis_final_exponent_z;
bool sw6_bis_final_exponent_is_z_neg;


void init_sw6_bis_params()
{
    typedef bigint<sw6_bis_r_limbs> bigint_r;
    typedef bigint<sw6_bis_q_limbs> bigint_q;

    assert(sizeof(mp_limb_t) == 8 || sizeof(mp_limb_t) == 4); // Montgomery assumes this

    /* parameters for scalar field Fr */

    sw6_bis_modulus_r = bigint_r("258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177");
    assert(sw6_bis_Fr::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        sw6_bis_Fr::Rsquared = bigint_r("66127428376872697816332570116866232405230528984664918319606315420233909940404532140033099444330447428417853902114"); // Rsquared = (W**k)**2 % r where k=6
        sw6_bis_Fr::Rcubed = bigint_r("157734475176213061358192738313701451942220138363611391489992831740412033225490229541667992423878570205050777755168");
        sw6_bis_Fr::inv = 0x8508bfffffffffff; // (-1/modulus) mod W
    }
    if (sizeof(mp_limb_t) == 4)
    {
        sw6_bis_Fr::Rsquared = bigint_r("66127428376872697816332570116866232405230528984664918319606315420233909940404532140033099444330447428417853902114");
        sw6_bis_Fr::Rcubed = bigint_r("157734475176213061358192738313701451942220138363611391489992831740412033225490229541667992423878570205050777755168");
        sw6_bis_Fr::inv = 0xffffffff;
    }
    sw6_bis_Fr::num_bits = 377;
    sw6_bis_Fr::euler = bigint_r("129332213006484547005326366847446766768196756377457330269942131333360234174170411387484444069786680062220160729088");
    sw6_bis_Fr::s = 46;
    sw6_bis_Fr::t = bigint_r("3675842578061421676390135839012792950148785745837396071634149488243117337281387659330802195819009059");
    sw6_bis_Fr::t_minus_1_over_2 = bigint_r("1837921289030710838195067919506396475074392872918698035817074744121558668640693829665401097909504529");
    sw6_bis_Fr::multiplicative_generator = sw6_bis_Fr("15");
    sw6_bis_Fr::root_of_unity = sw6_bis_Fr("32863578547254505029601261939868325669770508939375122462904745766352256812585773382134936404344547323199885654433");
    sw6_bis_Fr::nqr = sw6_bis_Fr("5");
    sw6_bis_Fr::nqr_to_t = sw6_bis_Fr("33774956008227656219775876656288133547078610493828613777258829345740556592044969439504850374928261397247202212840");
    // Done with Fr

    /* parameters for base field Fq */
    sw6_bis_modulus_q = bigint_q("6891450384315732539396789682275657542479668912536150109513790160209623422243491736087683183289411687640864567753786613451161759120554247759349511699125301598951605099378508850372543631423596795951899700429969112842764913119068299");
    assert(sw6_bis_Fq::modulus_is_valid());
    if (sizeof(mp_limb_t) == 8)
    {
        sw6_bis_Fq::Rsquared = bigint_q("4101737105507298352442561313393192324180371814155294089883586780083371310025435312104187656671185260872966272843049570295923422980866771377818994384387830909209154498924545983803406507410808360495749428678951279422657716620863065");
        sw6_bis_Fq::Rcubed = bigint_q("297415514632118086526834853439131925910825545338380841392089322640624198677430665521626961318971291885148318712859327336480040754936765047304839676586453060891603903912680899124818868997690823101153644786484398736894438670450156");

        sw6_bis_Fq::inv = 0xa5593568fa798dd;
    }
    if (sizeof(mp_limb_t) == 4)
    {
        sw6_bis_Fq::Rsquared = bigint_q("4101737105507298352442561313393192324180371814155294089883586780083371310025435312104187656671185260872966272843049570295923422980866771377818994384387830909209154498924545983803406507410808360495749428678951279422657716620863065");
        sw6_bis_Fq::Rcubed = bigint_q("297415514632118086526834853439131925910825545338380841392089322640624198677430665521626961318971291885148318712859327336480040754936765047304839676586453060891603903912680899124818868997690823101153644786484398736894438670450156");
        sw6_bis_Fq::inv = 0x8fa798dd;
    }
    sw6_bis_Fq::num_bits = 761;
    sw6_bis_Fq::euler = bigint_q("3445725192157866269698394841137828771239834456268075054756895080104811711121745868043841591644705843820432283876893306725580879560277123879674755849562650799475802549689254425186271815711798397975949850214984556421382456559534149");
    sw6_bis_Fq::s = 1;
    sw6_bis_Fq::t = bigint_q("3445725192157866269698394841137828771239834456268075054756895080104811711121745868043841591644705843820432283876893306725580879560277123879674755849562650799475802549689254425186271815711798397975949850214984556421382456559534149");
    sw6_bis_Fq::t_minus_1_over_2 = bigint_q("1722862596078933134849197420568914385619917228134037527378447540052405855560872934021920795822352921910216141938446653362790439780138561939837377924781325399737901274844627212593135907855899198987974925107492278210691228279767074");
    sw6_bis_Fq::multiplicative_generator = sw6_bis_Fq("2");
    sw6_bis_Fq::root_of_unity = sw6_bis_Fq("6891450384315732539396789682275657542479668912536150109513790160209623422243491736087683183289411687640864567753786613451161759120554247759349511699125301598951605099378508850372543631423596795951899700429969112842764913119068298");
    sw6_bis_Fq::nqr = sw6_bis_Fq("2");
    sw6_bis_Fq::nqr_to_t = sw6_bis_Fq("6891450384315732539396789682275657542479668912536150109513790160209623422243491736087683183289411687640864567753786613451161759120554247759349511699125301598951605099378508850372543631423596795951899700429969112842764913119068298");

    /* parameters for twist field Fq3 */
    sw6_bis_Fq3::euler = bigint<3*sw6_bis_q_limbs>("163644685426295400324914169940954888955920619637363322051209614326898453970148805686419011719286049425180559084625546731806777278679106679984663099143697116161447047621868304705520956489172590015220902479424161563352136860292139566900393097133389363042061971461347807883724603987331244933262704461518541983713459799376468603514421848466717682565239322945704382553665794724881767934477627356011124130848888786345848916593956748157692324278157203986274300190783511909765368289882962935850660388936997934686011687910278075743164270849285352627899889251291640002032048904566729137025370148302740362549750397891505412715221417910833483274288268810463328915955262815220574933205354950574767449");
    sw6_bis_Fq3::s = 1;
    sw6_bis_Fq3::t = bigint<3*sw6_bis_q_limbs>("163644685426295400324914169940954888955920619637363322051209614326898453970148805686419011719286049425180559084625546731806777278679106679984663099143697116161447047621868304705520956489172590015220902479424161563352136860292139566900393097133389363042061971461347807883724603987331244933262704461518541983713459799376468603514421848466717682565239322945704382553665794724881767934477627356011124130848888786345848916593956748157692324278157203986274300190783511909765368289882962935850660388936997934686011687910278075743164270849285352627899889251291640002032048904566729137025370148302740362549750397891505412715221417910833483274288268810463328915955262815220574933205354950574767449");
    sw6_bis_Fq3::t_minus_1_over_2 = bigint<3*sw6_bis_q_limbs>("81822342713147700162457084970477444477960309818681661025604807163449226985074402843209505859643024712590279542312773365903388639339553339992331549571848558080723523810934152352760478244586295007610451239712080781676068430146069783450196548566694681521030985730673903941862301993665622466631352230759270991856729899688234301757210924233358841282619661472852191276832897362440883967238813678005562065424444393172924458296978374078846162139078601993137150095391755954882684144941481467925330194468498967343005843955139037871582135424642676313949944625645820001016024452283364568512685074151370181274875198945752706357610708955416741637144134405231664457977631407610287466602677475287383724");
    sw6_bis_Fq3::non_residue = sw6_bis_Fq("2");
    sw6_bis_Fq3::nqr = sw6_bis_Fq3(sw6_bis_Fq("0"),sw6_bis_Fq("1"),sw6_bis_Fq("0"));
    sw6_bis_Fq3::nqr_to_t = sw6_bis_Fq3(sw6_bis_Fq("6891450384315732539396789682275657542479668912536150109513790160209623422243491736087683183289411687640864567753786613451161759120554247759349511699125301598951605099378508850372543631423596795951899700429969112842764913119068298"),sw6_bis_Fq("0"),sw6_bis_Fq("0"));
    sw6_bis_Fq3::Frobenius_coeffs_c1[0] = sw6_bis_Fq("1");
    sw6_bis_Fq3::Frobenius_coeffs_c1[1] = sw6_bis_Fq("1968985824090209297278610739700577151397666382303825728450741611566800370218827257750865013421937292370006175842381275743914023380727582819905021229583192207421122272650305267822868639090213645505120388400344940985710520836292650");
    sw6_bis_Fq3::Frobenius_coeffs_c1[2] = sw6_bis_Fq("4922464560225523242118178942575080391082002530232324381063048548642823052024664478336818169867474395270858391911405337707247735739826664939444490469542109391530482826728203582549674992333383150446779312029624171857054392282775648");
    sw6_bis_Fq3::Frobenius_coeffs_c2[0] = sw6_bis_Fq("1");
    sw6_bis_Fq3::Frobenius_coeffs_c2[1] = sw6_bis_Fq("4922464560225523242118178942575080391082002530232324381063048548642823052024664478336818169867474395270858391911405337707247735739826664939444490469542109391530482826728203582549674992333383150446779312029624171857054392282775648");
    sw6_bis_Fq3::Frobenius_coeffs_c2[2] = sw6_bis_Fq("1968985824090209297278610739700577151397666382303825728450741611566800370218827257750865013421937292370006175842381275743914023380727582819905021229583192207421122272650305267822868639090213645505120388400344940985710520836292650");


    /* parameters for Fq6 */
    sw6_bis_Fq6::non_residue = sw6_bis_Fq("2");
    sw6_bis_Fq6::Frobenius_coeffs_c1[0] = sw6_bis_Fq("1");
    sw6_bis_Fq6::Frobenius_coeffs_c1[1] = sw6_bis_Fq("1968985824090209297278610739700577151397666382303825728450741611566800370218827257750865013421937292370006175842381275743914023380727582819905021229583192207421122272650305267822868639090213645505120388400344940985710520836292651");
    sw6_bis_Fq6::Frobenius_coeffs_c1[2] = sw6_bis_Fq("1968985824090209297278610739700577151397666382303825728450741611566800370218827257750865013421937292370006175842381275743914023380727582819905021229583192207421122272650305267822868639090213645505120388400344940985710520836292650");
    sw6_bis_Fq6::Frobenius_coeffs_c1[3] = sw6_bis_Fq("6891450384315732539396789682275657542479668912536150109513790160209623422243491736087683183289411687640864567753786613451161759120554247759349511699125301598951605099378508850372543631423596795951899700429969112842764913119068298");
    sw6_bis_Fq6::Frobenius_coeffs_c1[4] = sw6_bis_Fq("4922464560225523242118178942575080391082002530232324381063048548642823052024664478336818169867474395270858391911405337707247735739826664939444490469542109391530482826728203582549674992333383150446779312029624171857054392282775648");
    sw6_bis_Fq6::Frobenius_coeffs_c1[5] = sw6_bis_Fq("4922464560225523242118178942575080391082002530232324381063048548642823052024664478336818169867474395270858391911405337707247735739826664939444490469542109391530482826728203582549674992333383150446779312029624171857054392282775649");
    sw6_bis_Fq6::my_Fp2::non_residue = sw6_bis_Fq3::non_residue;

    /* choice of short Weierstrass curve and its twist */
    sw6_bis_coeff_a = sw6_bis_Fq("0");
    sw6_bis_coeff_b = sw6_bis_Fq("5428247903343207843304490009542442997117969973913823318164330064320104021081180430153151788629347606889122435541645149581251622306618937027915190165889600280602116819284463614893841480114928247406042232431052188770336204942290254");
    sw6_bis_G1::coeff_a = sw6_bis_coeff_a;
    sw6_bis_G1::coeff_b = sw6_bis_coeff_b;
    sw6_bis_twist = sw6_bis_Fq3(sw6_bis_Fq::zero(), sw6_bis_Fq::one(), sw6_bis_Fq::zero());
    sw6_bis_twist_coeff_a = sw6_bis_Fq3(sw6_bis_Fq::zero(), sw6_bis_Fq::zero(), sw6_bis_G1::coeff_a);
    sw6_bis_twist_coeff_b = sw6_bis_Fq3(sw6_bis_G1::coeff_b * sw6_bis_Fq3::non_residue, sw6_bis_Fq::zero(), sw6_bis_Fq::zero());
    sw6_bis_G2::twist = sw6_bis_twist;
    sw6_bis_G2::coeff_a = sw6_bis_twist_coeff_a;
    sw6_bis_G2::coeff_b = sw6_bis_twist_coeff_b;
    sw6_bis_twist_mul_by_a_c0 = sw6_bis_G1::coeff_a * sw6_bis_Fq3::non_residue;
    sw6_bis_twist_mul_by_a_c1 = sw6_bis_G1::coeff_a * sw6_bis_Fq3::non_residue;
    sw6_bis_twist_mul_by_a_c2 = sw6_bis_G1::coeff_a;
    sw6_bis_twist_mul_by_b_c0 = sw6_bis_G1::coeff_b * sw6_bis_Fq3::non_residue;
    sw6_bis_twist_mul_by_b_c1 = sw6_bis_G1::coeff_b * sw6_bis_Fq3::non_residue;
    sw6_bis_twist_mul_by_b_c2 = sw6_bis_G1::coeff_b * sw6_bis_Fq3::non_residue;
    sw6_bis_twist_mul_by_q_X = sw6_bis_Fq("4922464560225523242118178942575080391082002530232324381063048548642823052024664478336818169867474395270858391911405337707247735739826664939444490469542109391530482826728203582549674992333383150446779312029624171857054392282775648");
    sw6_bis_twist_mul_by_q_Y = sw6_bis_Fq("6891450384315732539396789682275657542479668912536150109513790160209623422243491736087683183289411687640864567753786613451161759120554247759349511699125301598951605099378508850372543631423596795951899700429969112842764913119068298");


    /* choice of group G1 */
    sw6_bis_G1::G1_zero = sw6_bis_G1(sw6_bis_Fq::zero(),
                             sw6_bis_Fq::one(),
                             sw6_bis_Fq::zero());
    sw6_bis_G1::G1_one = sw6_bis_G1(sw6_bis_Fq("5579967068336365913530314810395365133780834399862286099844214670112509916823249335294962213691537977152362219416893517101994883989493084838175435762233462039369475294964563970855229509232628421518123427796600636650128165679151508"),
                            sw6_bis_Fq("3286106293257681699728332709973922945771547429394479993180335107651349218966829706081783660531045380298579922041512096833422398940161774037026391058614654750408470238881887954061918258846716359094162226636558919540815399302189033"),
                            sw6_bis_Fq::one());


    // TODO: wNAF window table
    sw6_bis_G1::wnaf_window_table.resize(0);
    sw6_bis_G1::wnaf_window_table.push_back(11);
    sw6_bis_G1::wnaf_window_table.push_back(24);
    sw6_bis_G1::wnaf_window_table.push_back(60);
    sw6_bis_G1::wnaf_window_table.push_back(127);

    // TODO: fixed-base exponentiation table
    sw6_bis_G1::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.99]
    sw6_bis_G1::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.99, 10.99]
    sw6_bis_G1::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.99, 32.29]
    sw6_bis_G1::fixed_base_exp_window_table.push_back(11);
    // window 4 is unbeaten in [32.29, 55.23]
    sw6_bis_G1::fixed_base_exp_window_table.push_back(32);
    // window 5 is unbeaten in [55.23, 162.03]
    sw6_bis_G1::fixed_base_exp_window_table.push_back(55);
    // window 6 is unbeaten in [162.03, 360.15]
    sw6_bis_G1::fixed_base_exp_window_table.push_back(162);
    // window 7 is unbeaten in [360.15, 815.44]
    sw6_bis_G1::fixed_base_exp_window_table.push_back(360);
    // window 8 is unbeaten in [815.44, 2373.07]
    sw6_bis_G1::fixed_base_exp_window_table.push_back(815);
    // window 9 is unbeaten in [2373.07, 6977.75]
    sw6_bis_G1::fixed_base_exp_window_table.push_back(2373);
    // window 10 is unbeaten in [6977.75, 7122.23]
    sw6_bis_G1::fixed_base_exp_window_table.push_back(6978);
    // window 11 is unbeaten in [7122.23, 57818.46]
    sw6_bis_G1::fixed_base_exp_window_table.push_back(7122);
    // window 12 is never the best
    sw6_bis_G1::fixed_base_exp_window_table.push_back(0);
    // window 13 is unbeaten in [57818.46, 169679.14]
    sw6_bis_G1::fixed_base_exp_window_table.push_back(57818);
    // window 14 is never the best
    sw6_bis_G1::fixed_base_exp_window_table.push_back(0);
    // window 15 is unbeaten in [169679.14, 439758.91]
    sw6_bis_G1::fixed_base_exp_window_table.push_back(169679);
    // window 16 is unbeaten in [439758.91, 936073.41]
    sw6_bis_G1::fixed_base_exp_window_table.push_back(439759);
    // window 17 is unbeaten in [936073.41, 4666554.74]
    sw6_bis_G1::fixed_base_exp_window_table.push_back(936073);
    // window 18 is never the best
    sw6_bis_G1::fixed_base_exp_window_table.push_back(0);
    // window 19 is unbeaten in [4666554.74, 7580404.42]
    sw6_bis_G1::fixed_base_exp_window_table.push_back(4666555);
    // window 20 is unbeaten in [7580404.42, 34552892.20]
    sw6_bis_G1::fixed_base_exp_window_table.push_back(7580404);
    // window 21 is never the best
    sw6_bis_G1::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [34552892.20, inf]
    sw6_bis_G1::fixed_base_exp_window_table.push_back(34552892);


    /* choice of group G2 */
    sw6_bis_G2::G2_zero = sw6_bis_G2(sw6_bis_Fq3::zero(),
                             sw6_bis_Fq3::one(),
                             sw6_bis_Fq3::zero());

    sw6_bis_G2::G2_one = sw6_bis_G2(sw6_bis_Fq3(sw6_bis_Fq("0"),
                                    sw6_bis_Fq("0"),
                                    sw6_bis_Fq("2762333537634180879188497831180460021045024768964216860077326831768282319295307213190468027456390916014688653832658199933857161399195589808528083537985014526923490693736135387257722421714776391535494333278083320743942589072320465")),
                            sw6_bis_Fq3(sw6_bis_Fq("6381154894827903503049202516448541497527668437176080964419566768101131101581568586920051373882212504439521771700091985568521381031252030757451435249323352844301649858691839209117065545550058565983491476754520157763374643574911714"),
                                    sw6_bis_Fq("0"),
                                    sw6_bis_Fq("0")),
                            sw6_bis_Fq3::one());


    // TODO: wNAF window table
    sw6_bis_G2::wnaf_window_table.resize(0);
    sw6_bis_G2::wnaf_window_table.push_back(5);
    sw6_bis_G2::wnaf_window_table.push_back(15);
    sw6_bis_G2::wnaf_window_table.push_back(39);
    sw6_bis_G2::wnaf_window_table.push_back(109);

    // TODO: fixed-base exponentiation table
    sw6_bis_G2::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 5.10]
    sw6_bis_G2::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [5.10, 10.43]
    sw6_bis_G2::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.43, 25.28]
    sw6_bis_G2::fixed_base_exp_window_table.push_back(10);
    // window 4 is unbeaten in [25.28, 59.00]
    sw6_bis_G2::fixed_base_exp_window_table.push_back(25);
    // window 5 is unbeaten in [59.00, 154.03]
    sw6_bis_G2::fixed_base_exp_window_table.push_back(59);
    // window 6 is unbeaten in [154.03, 334.25]
    sw6_bis_G2::fixed_base_exp_window_table.push_back(154);
    // window 7 is unbeaten in [334.25, 742.58]
    sw6_bis_G2::fixed_base_exp_window_table.push_back(334);
    // window 8 is unbeaten in [742.58, 2034.40]
    sw6_bis_G2::fixed_base_exp_window_table.push_back(743);
    // window 9 is unbeaten in [2034.40, 4987.56]
    sw6_bis_G2::fixed_base_exp_window_table.push_back(2034);
    // window 10 is unbeaten in [4987.56, 8888.27]
    sw6_bis_G2::fixed_base_exp_window_table.push_back(4988);
    // window 11 is unbeaten in [8888.27, 26271.13]
    sw6_bis_G2::fixed_base_exp_window_table.push_back(8888);
    // window 12 is unbeaten in [26271.13, 39768.20]
    sw6_bis_G2::fixed_base_exp_window_table.push_back(26271);
    // window 13 is unbeaten in [39768.20, 106275.75]
    sw6_bis_G2::fixed_base_exp_window_table.push_back(39768);
    // window 14 is unbeaten in [106275.75, 141703.40]
    sw6_bis_G2::fixed_base_exp_window_table.push_back(106276);
    // window 15 is unbeaten in [141703.40, 462422.97]
    sw6_bis_G2::fixed_base_exp_window_table.push_back(141703);
    // window 16 is unbeaten in [462422.97, 926871.84]
    sw6_bis_G2::fixed_base_exp_window_table.push_back(462423);
    // window 17 is unbeaten in [926871.84, 4873049.17]
    sw6_bis_G2::fixed_base_exp_window_table.push_back(926872);
    // window 18 is never the best
    sw6_bis_G2::fixed_base_exp_window_table.push_back(0);
    // window 19 is unbeaten in [4873049.17, 5706707.88]
    sw6_bis_G2::fixed_base_exp_window_table.push_back(4873049);
    // window 20 is unbeaten in [5706707.88, 31673814.95]
    sw6_bis_G2::fixed_base_exp_window_table.push_back(5706708);
    // window 21 is never the best
    sw6_bis_G2::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [31673814.95, inf]
    sw6_bis_G2::fixed_base_exp_window_table.push_back(31673815);



    /* pairing parameters */
    sw6_bis_ate_loop_count = bigint_q("6891450384315732539396789682275657542479668912536150109513790160209623422243491736087683183289411687640864567753786354786735746151460237106615816805591765205438850184717968966109876910955248455129124731541829539482640472797610122");
    sw6_bis_ate_is_loop_count_neg = false;
    sw6_bis_final_exponent = bigint<6*sw6_bis_q_limbs>("414120851190081857259378940145350857957559993168243292101436486535347107467806118274557843557325028910702591282977055677717179784934937158667585739736274135957119947179274672229952137675351215358786414154809447021601389903058849249956180740529241718879212085712503766864947946026605663395381450048308422396237123218798227193109519343275683314075860379350896921476723748271268458882462692568332631752289700097002045301296223372458293764434701950420210421715000131013980234881438585729277162712365054040927527479267564012274462728109348776756595766165627114213448725413038547637569495361079248349231297636997111755969684221288359818271693104882768483818996611430102893654553413082460832591443507134248207058064236465653945427795555589686925914139903467405281128720290325448429710297756756712477520844506970433659101005435238626369875843069961517840661170857661265550121820426065701082282309656905184495939454532818839387869011877715974614992299600773124822953681729343211729158584094701247610275434353440208893937178948454546809228173671015895702635689980844340167854010716240140417687462852252840200256248355635867968532779336190285724662667718441288196836598576766668991221618805280949647788766988939944453753742556685039123899111760770924773432347181567842600");
    sw6_bis_final_exponent_z = bigint_q("9586122913090633729");
    sw6_bis_final_exponent_is_z_neg = false;

}
} // libff
