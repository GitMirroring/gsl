#include <config.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>

#include "../test_sf_new.h"

int
test_lngamma(void)
{
  int s = 0;

  /* Near x = 1 rational region: 0.8 <= x <= 1.23395 */
  TEST_SF(s, gsl_new_sf_lngamma_e, (0.8, &r),
          0.152059678399837588778292602291, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (0.9, &r),
          0.0663762397347429711887167398671, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (1.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (1.1, &r),
          -0.0498724412598397241482898083159, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (1.23395, &r),
          -0.09446606702381254116508568555, TEST_TOL1, GSL_SUCCESS);

  /* Middle rational region: 1.23395 < x <= 1.67895 */
  TEST_SF(s, gsl_new_sf_lngamma_e, (1.25, &r),
          -0.0982718364218131614638538026966, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (1.5, &r),
          -0.120782237635245222345518445782, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (1.67895, &r),
          -0.1000209420168271613054818956, TEST_TOL1, GSL_SUCCESS);

  /* Near x = 2 rational region: 1.67895 < x < 2.3 */
  TEST_SF(s, gsl_new_sf_lngamma_e, (1.75, &r),
          -0.0844011210204855559577860341314, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (2.0, &r),
          0.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (2.25, &r),
          0.124871714892396594302441287613, TEST_TOL1, GSL_SUCCESS);

  /* Recurrence region: 2.3 <= x < 8 */
  TEST_SF(s, gsl_new_sf_lngamma_e, (2.3, &r),
          0.154189454959630581089917911489, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (2.5, &r),
          0.284682870472919159632494669683, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (2.67895, &r),
          0.418147656004422021021414196095, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (3.0, &r),
          0.693147180559945309417232121458, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (3.3, &r),
          0.987098577894734587878679288615, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (3.67895, &r),
          1.4035725827122472653951624762, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (4.5, &r),
          2.45373657084244222050414250344, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (5.5, &r),
          3.95781396761871629387740085582, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (6.5, &r),
          5.66256205985714152852211231233, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (7.5, &r),
          7.53436423675873295515836763244, TEST_TOL2, GSL_SUCCESS);

  /* Stirling U4 region: x >= 8 */
  TEST_SF(s, gsl_new_sf_lngamma_e, (8.0, &r),
          8.52516136106541430016553103635, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (10.0, &r),
          12.8018274800814696112077178746, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (12.0, &r),
          17.5023078458738858392876529072, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (32.0, &r),
          78.0922235533153106314168080587, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (64.0, &r),
          201.009316399281526679282039157, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (100.0, &r),
          359.13420536957539877604401046, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (1000.0, &r),
          5905.22042320918121182607691236, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (1000000.0, &r),
          12815504.569147611659976971785, TEST_TOL2, GSL_SUCCESS);

  /* Negative zero-relative region: -28 < x < -3 */
  TEST_SF(s, gsl_new_sf_lngamma_e, (-3.0499999999999998, &r),
          1.1449295208528985707438027720705, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-3.1435808883499798, &r),
          1.6978655906121084302227284796428e-15, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-3.25, &r),
          -0.62315345132099864662544730510872, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-3.75, &r),
          -1.3172679424463636738500787165944, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-3.9552942848585979, &r),
          -4.1438275075770499507244016834057e-16, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-4.0393618397405371, &r),
          -5.6645780740603349450285078420119e-15, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-4.25, &r),
          -2.0700724342573241080405176800655, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-4.75, &r),
          -2.8754125604929135150246419055659, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-4.9915446405600479, &r),
          1.7683619350849613451468704448523e-14, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-5.0082181683225935, &r),
          5.4188509265538102571959551633162e-15, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-5.25, &r),
          -3.7283005108608564857066514175148, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-5.75, &r),
          -4.6246124153021725869969304944598, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-5.9986074800808753, &r),
          -2.3721063667118472752237831005407e-13, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-6.0013852944531552, &r),
          -4.6433836942838006596799446326157e-14, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-6.25, &r),
          -5.5608819746091666160737058410509, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-6.75, &r),
          -6.534154920186611042348201962311, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-6.9998015078906377, &r),
          5.3130110657359020089733710643346e-14, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-7.000198333407325, &r),
          -1.2632037434939769640001792161819e-12, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-7.25, &r),
          -7.5418834434757500244225136304964, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-7.75, &r),
          -8.581847763551866669442902043937, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-8.0000248002706815, &r),
          1.7559556198603900467075421783668e-11, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-8.25, &r),
          -9.6520966438223396410452382024677, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-8.75, &r),
          -10.75090146392138973031454987769, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-9.0000027557148226, &r),
          3.444263328391509080169646312889e-11, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-9.25, &r),
          -11.876720195346673466578869630583, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-9.75, &r),
          -13.028168748931145538928818313261, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-10.000000275573013, &r),
          0.0000000012668051387565236605006265979856, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-10.25, &r),
          -14.203997900931090651611168760704, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-10.75, &r),
          -15.403074503504817343567196583691, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-11.000000025052106, &r),
          0.00000002734640389948086188539909967139, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-11.25, &r),
          -16.624366029581519790167954324859, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-11.75, &r),
          -17.866927744094985311553683010547, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-12.000000002087676, &r),
          -0.00000005854619992113373154438814384638, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-12.25, &r),
          -19.129891966572255781544195568829, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-12.75, &r),
          -20.412459015699420464363998622426, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-13.00000000016059, &r),
          0.0000042007069212685969408352038562843, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-13.25, &r),
          -21.713889519004486996854200464941, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-13.75, &r),
          -23.033497839812000764192237290701, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-14.00000000001147, &r),
          0.000070572041497083638297434379593237, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-14.25, &r),
          -24.370646425719146529424008890835, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-14.75, &r),
          -25.724740922597829595973823421504, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-15.000000000000764, &r),
          0.0011552549067268102218571173970878, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-15.25, &r),
          -27.095225928772567159340933751345, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-15.75, &r),
          -28.481581287869471665035202395876, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-16.000000000000046, &r),
          0.034255203335035831173333412019739, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-16.25, &r),
          -29.88331883754831365116071628322, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-16.75, &r),
          -31.299979546140547105870810149324, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-17.000000000000004, &r),
          -0.23400878325952420104406836235182, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-17.25, &r),
          -32.731130981025682414528250109036, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-17.75, &r),
          -34.17636506206197190836580043892, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-18.25, &r),
          -35.635296061054182924785894723558, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-18.75, &r),
          -37.107558814478391730128100099379, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-19.25, &r),
          -38.59280712178797615511872680205, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-19.75, &r),
          -40.090712305825522605466581397944, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-20.25, &r),
          -41.600961915340524301865243506824, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-20.75, &r),
          -43.123258552502229910107589378319, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-21.25, &r),
          -44.657318810710950137881073215007, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-21.75, &r),
          -46.202872310036923009851642404687, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-22.25, &r),
          -47.759660819323199357364424512761, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-22.75, &r),
          -49.327437455433882432176018346779, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-23.25, &r),
          -50.905965951356564675854369831309, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-23.75, &r),
          -52.495019985914532647951340868977, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-24.25, &r),
          -54.094382568740056879136627210097, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-24.75, &r),
          -55.703845474929231955969310677871, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-25.25, &r),
          -57.323208724461425711186361234093, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-25.75, &r),
          -58.952280102038977107903448750008, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-26.25, &r),
          -60.590874713499058463453254304769, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-26.75, &r),
          -62.238814575380992662373383324115, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-27.25, &r),
          -63.895928234608311544996106514062, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-27.75, &r),
          -65.562050415573436179302259989152, TEST_TOL6, GSL_SUCCESS);

  /* Large-negative reflection region: x <= -28 */
  TEST_SF(s, gsl_new_sf_lngamma_e, (-28.000000000000004, &r),
          -34.61867847030417203229627399158509, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-28.25, &r),
          -67.23702169220076149474377376476146, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-28.5, &r),
          -68.42435103497423400849630744022931, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-28.75, &r),
          -68.92068818281679562587530791127238, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-29.000000000000004, &r),
          -37.98597430029064618198691425844703, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-29.25, &r),
          -70.61290126587862699475328743725548, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-29.5, &r),
          -71.80874129832000814969512569249083, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-29.75, &r),
          -72.31351731480843439239573102967233, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-31.25, &r),
          -77.46441682653788796901052410703159, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-31.5, &r),
          -78.67645552776496146750789377028832, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-31.75, &r),
          -79.19730003439966192303625196151585, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-32.25, &r),
          -80.93793486977966946504414761438361, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-32.5, &r),
          -82.15769561710065326874490842362165, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-32.75, &r),
          -82.68620299648092284835121172976902, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-40.25, &r),
          -109.7554380653398315321672354251256, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-40.5, &r),
          -111.029647157328477335218671620959, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-40.75, &r),
          -111.6122333928907274943029661919722, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-64.25, &r),
          -204.7190486779166394283223660722472, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-64.5, &r),
          -206.1087401737804492419808172070866, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-64.75, &r),
          -206.8062460109525578353494408059117, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-100.25, &r),
          -363.4009232278215407065774773916365, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-100.5, &r),
          -364.9009683094273518227565704629958, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-100.75, &r),
          -365.7084850179312873651435587991706, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-256.25, &r),
          -1167.152879203871903426365914845991, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-256.5, &r),
          -1168.886600338462295513478808825679, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-256.75, &r),
          -1169.927417482885745012975903429177, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-1024.25, &r),
          -6078.45360177596590585229854225924, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-1024.5, &r),
          -6080.533256858772292409362888896278, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-1024.75, &r),
          -6081.919825736624257857551787262157, TEST_TOL6, GSL_SUCCESS);

  /* Special glibc polynomial region: -3 < x < -2 */
  TEST_SF(s, gsl_new_sf_lngamma_e, (-2.9999999999999996, &r),
          33.5587467393291563372946636885678129, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-2.9500000000000002, &r),
          1.27053795251190077044012793771420479, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-2.875, &r),
          0.46828841650663106816726946888859951, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-2.75, &r),
          0.00448789753595577331146161063794197027, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-2.7476826467274127, &r),
          1.73350924402450086109664880252406536e-16, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-2.625, &r),
          -0.117511100308904859994666567899607071, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-2.5000000000000004, &r),
          -0.0562437164976745405725554150922302611, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-2.5, &r),
          -0.0562437164976740506725945300976542841, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-2.4570247382208006, &r),
          5.61919235895009645091256851694733931e-17, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-2.375, &r),
          0.158207579895527270473830532170860059, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-2.25, &r),
          0.555501545020647470593575893540250059, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-2.125, &r),
          1.29374541087240664521742958517489887, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-2.0499999999999998, &r),
          2.26007111147221876781876045016499848, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-2.0000000000000004, &r),
          34.6573590279972650610630398370154388, TEST_TOL6, GSL_SUCCESS);

  /* Negative ordinary-reflection region: -2 < x < -1 */
  TEST_SF(s, gsl_new_sf_lngamma_e, (-1.9999999999999998, &r),
          35.3505062085572109851781213123139139, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-1.95, &r),
          2.35234312286362508400073915460134233, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-1.875, &r),
          1.52434109075594483072232593632426592, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-1.75, &r),
          1.01608880921443569853894094568671813, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-1.625, &r),
          0.847569795734682208254235048091568659, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-1.5, &r),
          0.860047015376481014510932681670356787, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-1.375, &r),
          1.02320501738213180223116159968418389, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-1.25, &r),
          1.36643176123697623454960212446894833, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-1.125, &r),
          2.0475172132487867972152678386734957, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-1.05, &r),
          2.97791090462253098636448674523657155, TEST_TOL6, GSL_SUCCESS);

  TEST_SF(s, gsl_new_sf_lngamma_e, (-1.0000000000000002, &r),
          36.0436533891171559958190896603940632, TEST_TOL6, GSL_SUCCESS);

  return s;
}
