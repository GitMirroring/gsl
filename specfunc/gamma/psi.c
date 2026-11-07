/* specfunc/gamma/psi.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2005, 2006 Gerard Jungman
 * Copyright (C) 2026 Patrick Alken
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*
 * References: the reflection, recurrence and asymptotic formulas are
 * standard (NIST DLMF 5.5, 5.11, 5.15).  The coefficients of the series
 * about the zeros of psi(x) and Re psi(1+iy) and the Chebyshev fits of
 * Re psi(1+iy) are mathematical constants computed for GSL.  The tables
 * for integer arguments and the complex and polygamma routines are
 * carried over from the previous GSL implementation (G. Jungman).
 */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_dd.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_sincos_pi.h>
#include <gsl/gsl_complex_math.h>

#include "constants.h"
#include "error.h"
#include "eval.h"
#include "chebyshev.h"
#include "cheb_eval.c"

/* zero of psi(x) on x > 0, as a double-double pair */
#define PSI_X0_HI  1.46163214496836225e+00
#define PSI_X0_LO  9.54999542996569743e-17

/* zero of Re psi(1+iy) on y > 0, as a double-double pair */
#define PSI_1PIY_Y0_HI  8.83820637894184746e-01
#define PSI_1PIY_Y0_LO  -9.31171737642724284e-19

/*
 * psi(x) for x < 0 is recomputed in double-double when its two terms
 * cancel by more than this factor (near the zeros of psi)
 */
#define PSI_REFINE_AMP  2.0

/* smallest argument of the double precision asymptotic series */
#define PSI_ASYMP_MIN  10.0

/* smallest argument of the double-double asymptotic series */
#define PSI_DD_ASYMP_MIN  20.0

/* absolute error of the double-double evaluation relative to its terms */
#define PSI_DD_REL_ERR  4.0e-32

static double psi_x0_series(const double th, const double tl);
static int psi_pos(const double x, gsl_sf_result * result);
static int psi_neg(const double x, gsl_sf_result * result);
static double psi1_pos(const double x);
static gsl_dd psi_dd(const gsl_dd x, double * mag);
static gsl_dd psi_dd_pi_cot_pi(const double x);

#define PSI_TABLE_NMAX 100
static double psi_table[PSI_TABLE_NMAX+1] = {
  0.0,  /* Infinity */              /* psi(0) */
 -M_EULER,                          /* psi(1) */
  0.42278433509846713939348790992,  /* ...    */
  0.92278433509846713939348790992,
  1.25611766843180047272682124325,
  1.50611766843180047272682124325,
  1.70611766843180047272682124325,
  1.87278433509846713939348790992,
  2.01564147795560999653634505277,
  2.14064147795560999653634505277,
  2.25175258906672110764745616389,
  2.35175258906672110764745616389,
  2.44266167997581201673836525479,
  2.52599501330914535007169858813,
  2.60291809023222227314862166505,
  2.67434666166079370172005023648,
  2.74101332832746036838671690315,
  2.80351332832746036838671690315,
  2.86233685773922507426906984432,
  2.91789241329478062982462539988,
  2.97052399224214905087725697883,
  3.02052399224214905087725697883,
  3.06814303986119666992487602645,
  3.11359758531574212447033057190,
  3.15707584618530734186163491973,
  3.1987425128519740085283015864,
  3.2387425128519740085283015864,
  3.2772040513135124700667631249,
  3.3142410883505495071038001619,
  3.3499553740648352213895144476,
  3.3844381326855248765619282407,
  3.4177714660188582098952615740,
  3.4500295305349872421533260902,
  3.4812795305349872421533260902,
  3.5115825608380175451836291205,
  3.5409943255438998981248055911,
  3.5695657541153284695533770196,
  3.5973435318931062473311547974,
  3.6243705589201332743581818244,
  3.6506863483938174848844976139,
  3.6763273740348431259101386396,
  3.7013273740348431259101386396,
  3.7257176179372821503003825420,
  3.7495271417468059598241920658,
  3.7727829557002943319172153216,
  3.7955102284275670591899425943,
  3.8177324506497892814121648166,
  3.8394715810845718901078169905,
  3.8607481768292527411716467777,
  3.8815815101625860745049801110,
  3.9019896734278921969539597029,
  3.9219896734278921969539597029,
  3.9415975165651470989147440166,
  3.9608282857959163296839747858,
  3.9796962103242182164764276160,
  3.9982147288427367349949461345,
  4.0163965470245549168131279527,
  4.0342536898816977739559850956,
  4.0517975495308205809735289552,
  4.0690389288411654085597358518,
  4.0859880813835382899156680552,
  4.1026547480502049565823347218,
  4.1190481906731557762544658694,
  4.1351772229312202923834981274,
  4.1510502388042361653993711433,
  4.1666752388042361653993711433,
  4.1820598541888515500147557587,
  4.1972113693403667015299072739,
  4.2121367424746950597388624977,
  4.2268426248276362362094507330,
  4.2413353784508246420065521823,
  4.2556210927365389277208378966,
  4.2697055997787924488475984600,
  4.2835944886676813377364873489,
  4.2972931188046676391063503626,
  4.3108066323181811526198638761,
  4.3241399656515144859531972094,
  4.3372978603883565912163551041,
  4.3502848733753695782293421171,
  4.3631053861958823987421626300,
  4.3757636140439836645649474401,
  4.3882636140439836645649474401,
  4.4006092930563293435772931191,
  4.4128044150075488557724150703,
  4.4248526077786331931218126607,
  4.4367573696833950978837174226,
  4.4485220755657480390601880108,
  4.4601499825424922251066996387,
  4.4716442354160554434975042364,
  4.4830078717796918071338678728,
  4.4942438268358715824147667492,
  4.5053549379469826935258778603,
  4.5163439489359936825368668713,
  4.5272135141533849868846929582,
  4.5379662023254279976373811303,
  4.5486045001977684231692960239,
  4.5591308159872421073798223397,
  4.5695474826539087740464890064,
  4.5798567610044242379640147796,
  4.5900608426370772991885045755,
  4.6001618527380874001986055856
};


#define PSI_1_TABLE_NMAX 100
static double psi_1_table[PSI_1_TABLE_NMAX+1] = {
  0.0,  /* Infinity */              /* psi(1,0) */
  M_PI*M_PI/6.0,                    /* psi(1,1) */
  0.644934066848226436472415,       /* ...      */
  0.394934066848226436472415,
  0.2838229557371153253613041,
  0.2213229557371153253613041,
  0.1813229557371153253613041,
  0.1535451779593375475835263,
  0.1331370146940314251345467,
  0.1175120146940314251345467,
  0.1051663356816857461222010,
  0.0951663356816857461222010,
  0.0869018728717683907503002,
  0.0799574284273239463058557,
  0.0740402686640103368384001,
  0.0689382278476838062261552,
  0.0644937834032393617817108,
  0.0605875334032393617817108,
  0.0571273257907826143768665,
  0.0540409060376961946237801,
  0.0512708229352031198315363,
  0.0487708229352031198315363,
  0.0465032492390579951149830,
  0.0444371335365786562720078,
  0.0425467743683366902984728,
  0.0408106632572255791873617,
  0.0392106632572255791873617,
  0.0377313733163971768204978,
  0.0363596312039143235969038,
  0.0350841209998326909438426,
  0.0338950603577399442137594,
  0.0327839492466288331026483,
  0.0317433665203020901265817,
  0.03076680402030209012658168,
  0.02984853037475571730748159,
  0.02898347847164153045627052,
  0.02816715194102928555831133,
  0.02739554700275768062003973,
  0.02666508681283803124093089,
  0.02597256603721476254286995,
  0.02531510384129102815759710,
  0.02469010384129102815759710,
  0.02409521984367056414807896,
  0.02352832641963428296894063,
  0.02298749353699501850166102,
  0.02247096461137518379091722,
  0.02197713745088135663042339,
  0.02150454765882086513703965,
  0.02105185413233829383780923,
  0.02061782635456051606003145,
  0.02020133322669712580597065,
  0.01980133322669712580597065,
  0.01941686571420193164987683,
  0.01904704322899483105816086,
  0.01869104465298913508094477,
  0.01834810912486842177504628,
  0.01801753061247172756017024,
  0.01769865306145131939690494,
  0.01739086605006319997554452,
  0.01709360088954001329302371,
  0.01680632711763538818529605,
  0.01652854933985761040751827,
  0.01625980437882562975715546,
  0.01599965869724394401313881,
  0.01574770606433893015574400,
  0.01550356543933893015574400,
  0.01526687904880638577704578,
  0.01503731063741979257227076,
  0.01481454387422086185273411,
  0.01459828089844231513993134,
  0.01438824099085987447620523,
  0.01418415935820681325171544,
  0.01398578601958352422176106,
  0.01379288478501562298719316,
  0.01360523231738567365335942,
  0.01342261726990576130858221,
  0.01324483949212798353080444,
  0.01307170929822216635628920,
  0.01290304679189732236910755,
  0.01273868124291638877278934,
  0.01257845051066194236996928,
  0.01242220051066194236996928,
  0.01226978472038606978956995,
  0.01212106372098095378719041,
  0.01197590477193174490346273,
  0.01183418141592267460867815,
  0.01169577311142440471248438,
  0.01156056489076458859566448,
  0.01142844704164317229232189,
  0.01129931481023821361463594,
  0.01117306812421372175754719,
  0.01104961133409026496742374,
  0.01092885297157366069257770,
  0.01081070552355853781923177,
  0.01069508522063334415522437,
  0.01058191183901270133041676,
  0.01047110851491297833872701,
  0.01036260157046853389428257,
  0.01025632035036012704977199,  /* ...        */
  0.01015219706839427948625679,  /* psi(1,99)  */
  0.01005016666333357139524567   /* psi(1,100) */
};

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_psi_int_e(const int n, gsl_sf_result * result)
{
  if (n <= 0)
    {
      DOMAIN_ERROR(result);
    }
  else if (n <= PSI_TABLE_NMAX)
    {
      result->val = psi_table[n];
      result->err = GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
  else
    {
      /* Abramowitz+Stegun 6.3.18 */
      const double c2 = -1.0 / 12.0;
      const double c3 =  1.0 / 120.0;
      const double c4 = -1.0 / 252.0;
      const double c5 =  1.0 / 240.0;
      const double ni2 = (1.0 / n) * (1.0 / n);
      const double ser = ni2 * (c2 + ni2 * (c3 + ni2 * (c4 + ni2 * c5)));
      result->val  = log(n) - 0.5 / n + ser;
      result->err  = GSL_DBL_EPSILON * (fabs(log(n)) + fabs(0.5 / n) + fabs(ser));
      result->err += GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
}

/*
 * Digamma function psi(x) = Gamma'(x)/Gamma(x).
 *
 * x > 0: psi_pos(), with relative accuracy up to the zero x0 = 1.4616...
 * x < 0: the reflection formula (psi_neg()), with a double-double
 * refinement near the zeros of psi on the negative axis.
 * Zero and the negative integers are poles (GSL_EDOM).
 */
int
gsl_sf_psi_e(const double x, gsl_sf_result * result)
{
  if (x <= 0.0 && x == floor(x))
    {
      DOMAIN_ERROR(result);
    }
  else if (x > 0.0)
    {
      return psi_pos(x, result);
    }
  else
    {
      return psi_neg(x, result);
    }
}

/*
 * Re psi(1+iy), an even function of y with a single zero y0 = 0.8838...
 * on y > 0:
 *
 *   |y| <= 1/2:     -Euler + sum_{k>=1} (-1)^(k+1) zeta(2k+1) y^(2k)
 *   1/2 < |y| <= 1.25: Taylor series about y0 (relative accuracy at the zero)
 *   1.25 < |y| <= 10:  Chebyshev fits on [1.25,4] and [4,10]
 *   |y| > 10:       log y + sum_{k>=1} (-1)^(k+1) B_2k/(2k) y^(-2k)
 */
int
gsl_sf_psi_1piy_e(const double y, gsl_sf_result * result)
{
  /* zeta(2k+1), k = 1..30 */
  static const double z[30] = {
       1.20205690315959424e+00,    1.03692775514336999e+00,
       1.00834927738192293e+00,    1.00200839282608212e+00,
       1.00049418860411943e+00,    1.00012271334757852e+00,
       1.00003058823630697e+00,    1.00000763719763786e+00,
       1.00000190821271651e+00,    1.00000047693298688e+00,
       1.00000011921992593e+00,    1.00000002980350344e+00,
       1.00000000745071183e+00,    1.00000000186265980e+00,
       1.00000000046566284e+00,    1.00000000011641554e+00,
       1.00000000002910383e+00,    1.00000000000727596e+00,
       1.00000000000181899e+00,    1.00000000000045475e+00,
       1.00000000000011369e+00,    1.00000000000002842e+00,
       1.00000000000000711e+00,    1.00000000000000178e+00,
       1.00000000000000044e+00,    1.00000000000000022e+00,
       1.00000000000000000e+00,    1.00000000000000000e+00,
       1.00000000000000000e+00,    1.00000000000000000e+00
  };
  /* Taylor coefficients of Re psi(1+iy) about y0 */
  static const double d[42] = {
       8.33935793911211531e-01,   -1.52724357576886210e-01,
      -1.34504164716683944e-01,    2.18280065091297104e-01,
      -1.58186897945267246e-01,    4.17875934153908862e-02,
       4.72611482587214707e-02,   -7.17194637991090261e-02,
       4.51955196435194426e-02,   -4.52814042033327334e-03,
      -2.10174163571579843e-02,    2.34403486271278945e-02,
      -1.14492991178855134e-02,   -1.81118212711013659e-03,
       8.22773730817463918e-03,   -7.14653257078361973e-03,
       2.47180620344362870e-03,    1.55929358805749507e-03,
      -2.93500381679082187e-03,    2.03721931089553875e-03,
      -3.73971942945517473e-04,   -7.72609928028310071e-04,
       9.76710984280452707e-04,   -5.35539413885981392e-04,
      -1.68796959685303266e-05,    3.17423968531961608e-04,
      -3.05541996599619447e-04,    1.25012733665633862e-04,
       4.74776804891770444e-05,   -1.17304900915367834e-04,
       8.97601696460443855e-05,   -2.32204642477431650e-05,
      -2.73503242603964208e-05,    4.01799248260262628e-05,
      -2.45199111979405567e-05,    1.77553894564031722e-06,
       1.20043389481891067e-05,   -1.29102247753547317e-05,
       6.07269150157943055e-06,    1.22162543865081971e-06,
      -4.62181076524322502e-06,    3.90091865195996891e-06
  };
  /* B_2k/(2k), k = 1..10 */
  static const double b[10] = {
       8.33333333333333287e-02,   -8.33333333333333322e-03,
       3.96825396825396803e-03,   -4.16666666666666661e-03,
       7.57575757575757597e-03,   -2.10927960927960942e-02,
       8.33333333333333287e-02,   -4.43259803921568607e-01,
       3.05395433027011975e+00,   -2.64562121212121220e+01
  };
  /* Chebyshev fits of Re psi(1+iy) on [1.25,4] and [4,10] */
  static const double cs1[34] = {
       1.81738914446025457e+00,    5.43558920211707708e-01,
      -7.04763135026604676e-02,    1.13529360012503099e-02,
      -1.81200864383459950e-03,    2.27237378961392031e-04,
       7.70441499780008230e-07,   -1.53833069028488833e-05,
       7.65790091226951562e-06,   -2.75614219229542468e-06,
       8.31004710964271945e-07,   -2.15620548913666935e-07,
       4.71164610674761877e-08,   -7.82724755666526243e-09,
       5.15497025994041669e-10,    2.87592058369800939e-10,
      -1.71088767559150771e-10,    6.14746906478973073e-11,
      -1.75457520710836177e-11,    4.17540358000594500e-12,
      -8.04374864860851738e-13,    1.03572965857006065e-13,
       3.95169514813429030e-15,   -8.57857643820814240e-15,
       3.75400257634763967e-15,   -1.19852741980750832e-15,
       3.16802304068277538e-16,   -7.03189738346474269e-17,
       1.22913997969262573e-17,   -1.16979560794429331e-18,
      -2.55539085475484025e-19,    1.96714460060975807e-19,
      -7.56728721304322104e-20,    2.26555302111407989e-20
  };
  static const double cs2[29] = {
       3.79755293667786109e+00,    4.48307006828077059e-01,
      -5.00288017769490931e-02,    7.41056420137154500e-03,
      -1.22902481308985283e-03,    2.16298416384768050e-04,
      -3.94273419880513290e-05,    7.34494405295170155e-06,
      -1.38653697282011811e-06,    2.63602948427913426e-07,
      -5.02146394868363793e-08,    9.53798502677051861e-09,
      -1.79674752443235192e-09,    3.33373175836459239e-10,
      -6.03142838208757194e-11,    1.04630514460546013e-11,
      -1.68393172257496783e-12,    2.31492614745617143e-13,
      -1.90506417159770772e-14,   -3.25328340772303036e-15,
       2.48094832028597363e-15,   -9.69234214195195628e-16,
       3.15482282936879289e-16,   -9.41203374949085338e-17,
       2.66460222399949866e-17,   -7.27118827278595164e-18,
       1.92699940502824338e-18,   -4.97717861125068927e-19,
       1.25442535589647278e-19
  };
  const double ay = fabs(y);

  if (ay <= 0.5)
    {
      const double y2 = ay * ay;
      double s = z[29];
      int k;

      for (k = 28; k >= 0; --k)
        s = z[k] - y2 * s;

      result->val = y2 * s - M_EULER;
      result->err = 2.0 * GSL_DBL_EPSILON * (M_EULER + fabs(result->val));
      return GSL_SUCCESS;
    }
  else if (ay <= 1.25)
    {
      const double t = (ay - PSI_1PIY_Y0_HI) - PSI_1PIY_Y0_LO;   /* ay - Y0_HI is exact */
      double p = d[41];
      int k;

      for (k = 40; k >= 0; --k)
        p = p * t + d[k];

      result->val = t * p;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
  else if (ay <= 10.0)
    {
      static const cheb_series cs_lo = { (double *) cs1, 33, 1.25, 4.0, 33 };
      static const cheb_series cs_hi = { (double *) cs2, 28, 4.0, 10.0, 28 };
      const int status = cheb_eval_e(ay <= 4.0 ? &cs_lo : &cs_hi, ay, result);

      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return status;
    }
  else
    {
      const double w = 1.0 / (ay * ay);
      const double ly = log(ay);
      double s = b[9];
      int k;

      for (k = 8; k >= 0; --k)
        s = b[k] - w * s;

      s *= w;
      result->val = ly + s;
      result->err = 2.0 * GSL_DBL_EPSILON * (fabs(ly) + fabs(s));
      return GSL_SUCCESS;
    }
}

int
gsl_sf_psi_1_int_e(const int n, gsl_sf_result * result)
{
  if (n <= 0)
    {
      DOMAIN_ERROR(result);
    }
  else if (n <= PSI_1_TABLE_NMAX)
    {
      result->val = psi_1_table[n];
      result->err = GSL_DBL_EPSILON * result->val;
      return GSL_SUCCESS;
    }
  else
    {
      /* the positive-argument trigamma kernel (a sum of positive terms) */
      result->val = psi1_pos((double) n);
      result->err = 3.0 * GSL_DBL_EPSILON * result->val;
      return GSL_SUCCESS;
    }
}

/*
 * Trigamma function psi_1(x) = psi'(x).
 *
 * x > 0: psi1_pos(), a sum of positive terms.
 * x < 0: psi_1(x) = pi^2/sin^2(pi x) - psi_1(1-x), with sin(pi x) from
 * the exact argument; the result is at least 5 pi^2/6, so the
 * subtraction loses at most a factor ~1.2.
 */
int
gsl_sf_psi_1_e(const double x, gsl_sf_result * result)
{
  if (x <= 0.0 && x == floor(x))
    {
      DOMAIN_ERROR(result);
    }
  else if (x > 0.0)
    {
      const double v = psi1_pos(x);

      if (!gsl_finite(v))
        OVERFLOW_ERROR(result);

      result->val = v;
      result->err = 3.0 * GSL_DBL_EPSILON * v;
      return GSL_SUCCESS;
    }
  else
    {
      const double s = gsl_sf_sin_pi(x);
      const double q = M_PI / s;
      const double d = q * q;
      const double r = psi1_pos(1.0 - x);

      if (!gsl_finite(d))
        OVERFLOW_ERROR(result);

      result->val = d - r;
      result->err = 3.0 * GSL_DBL_EPSILON * (d + r);
      return GSL_SUCCESS;
    }
}

int
gsl_sf_psi_n_e(const int n, const double x, gsl_sf_result * result)
{
  if (n == 0)
    {
      return gsl_sf_psi_e(x, result);
    }
  else if (n == 1)
    {
      return gsl_sf_psi_1_e(x, result);
    }
  else if (n < 0 || x <= 0.0)
    {
      DOMAIN_ERROR(result);
    }
  else
    {
      /* Abramowitz + Stegun 6.4.10 */
      gsl_sf_result ln_nf;
      gsl_sf_result hzeta;
      int stat_hz = gsl_sf_hzeta_e(n + 1.0, x, &hzeta);
      int stat_nf = gsl_sf_lnfact_e((unsigned int) n, &ln_nf);
      int stat_e  = gsl_sf_exp_mult_err_e(ln_nf.val, ln_nf.err,
                                          hzeta.val, hzeta.err,
                                          result);
      if (GSL_IS_EVEN(n))
        result->val = -result->val;
      return GSL_ERROR_SELECT_3(stat_e, stat_nf, stat_hz);
    }
}

/* psi(z) for large |z| in the right half-plane; [Abramowitz + Stegun, 6.3.18] */
static
gsl_complex
psi_complex_asymp(gsl_complex z)
{
  /* coefficients in the asymptotic expansion for large z;
   * let w = z^(-2) and write the expression in the form
   *
   *   ln(z) - 1/(2z) - 1/12 w (1 + c1 w + c2 w + c3 w + ... )
   */
  static const double c1 = -0.1;
  static const double c2 =  1.0/21.0;
  static const double c3 = -0.05;

  gsl_complex zi = gsl_complex_inverse(z);
  gsl_complex w  = gsl_complex_mul(zi, zi);
  gsl_complex cs;

  /* Horner method evaluation of term in parentheses */
  gsl_complex sum;
  sum = gsl_complex_mul_real(w, c3/c2);
  sum = gsl_complex_add_real(sum, 1.0);
  sum = gsl_complex_mul_real(sum, c2/c1);
  sum = gsl_complex_mul(sum, w);
  sum = gsl_complex_add_real(sum, 1.0);
  sum = gsl_complex_mul_real(sum, c1);
  sum = gsl_complex_mul(sum, w);
  sum = gsl_complex_add_real(sum, 1.0);

  /* correction added to log(z) */
  cs = gsl_complex_mul(sum, w);
  cs = gsl_complex_mul_real(cs, -1.0/12.0);
  cs = gsl_complex_add(cs, gsl_complex_mul_real(zi, -0.5));

  return gsl_complex_add(gsl_complex_log(z), cs);
}



/* psi(z) for complex z in the right half-plane */
static int
psi_complex_rhp(
  gsl_complex z,
  gsl_sf_result * result_re,
  gsl_sf_result * result_im
  )
{
  int n_recurse = 0;
  int i;
  gsl_complex a;

  if(GSL_REAL(z) == 0.0 && GSL_IMAG(z) == 0.0)
  {
    result_re->val = 0.0;
    result_im->val = 0.0;
    result_re->err = 0.0;
    result_im->err = 0.0;
    return GSL_EDOM;
  }

  /* compute the number of recurrences to apply */
  if(GSL_REAL(z) < 20.0 && fabs(GSL_IMAG(z)) < 20.0)
  {
    const double sp = sqrt(20.0 + GSL_IMAG(z));
    const double sn = sqrt(20.0 - GSL_IMAG(z));
    const double rhs = sp*sn - GSL_REAL(z);
    if(rhs > 0.0) n_recurse = ceil(rhs);
  }

  /* compute asymptotic at the large value z + n_recurse */
  a = psi_complex_asymp(gsl_complex_add_real(z, n_recurse));

  result_re->err = 2.0 * GSL_DBL_EPSILON * fabs(GSL_REAL(a));
  result_im->err = 2.0 * GSL_DBL_EPSILON * fabs(GSL_IMAG(a));

  /* descend recursively, if necessary */
  for(i = n_recurse; i >= 1; --i)
  {
    gsl_complex zn = gsl_complex_add_real(z, i - 1.0);
    gsl_complex zn_inverse = gsl_complex_inverse(zn);
    a = gsl_complex_sub(a, zn_inverse);

    /* accumulate the error, to catch cancellations */
    result_re->err += 2.0 * GSL_DBL_EPSILON * fabs(GSL_REAL(zn_inverse));
    result_im->err += 2.0 * GSL_DBL_EPSILON * fabs(GSL_IMAG(zn_inverse));
  }

  result_re->val = GSL_REAL(a);
  result_im->val = GSL_IMAG(a);

  result_re->err += 2.0 * GSL_DBL_EPSILON * fabs(result_re->val);
  result_im->err += 2.0 * GSL_DBL_EPSILON * fabs(result_im->val);

  return GSL_SUCCESS;
}

int
gsl_sf_complex_psi_e(const double x, const double y,
                     gsl_sf_result * result_re, gsl_sf_result * result_im)
{
  if (x >= 0.0)
    {
      gsl_complex z = gsl_complex_rect(x, y);
      return psi_complex_rhp(z, result_re, result_im);
    }
  else
    {
      /* reflection formula [Abramowitz+Stegun, 6.3.7] */
      gsl_complex z = gsl_complex_rect(x, y);
      gsl_complex omz = gsl_complex_rect(1.0 - x, -y);
      gsl_complex zpi = gsl_complex_mul_real(z, M_PI);
      gsl_complex cotzpi = gsl_complex_cot(zpi);
      int ret_val = psi_complex_rhp(omz, result_re, result_im);

      if (GSL_IS_REAL(GSL_REAL(cotzpi)) && GSL_IS_REAL(GSL_IMAG(cotzpi)))
        {
          result_re->val -= M_PI * GSL_REAL(cotzpi);
          result_im->val -= M_PI * GSL_IMAG(cotzpi);
          return ret_val;
        }
      else
        {
          GSL_ERROR("singularity", GSL_EDOM);
        }
    }
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_psi_int(const int n) { EVAL_RESULT(gsl_sf_psi_int_e(n, &result)); }
double gsl_sf_psi(const double x) { EVAL_RESULT(gsl_sf_psi_e(x, &result)); }
double gsl_sf_psi_1piy(const double x) { EVAL_RESULT(gsl_sf_psi_1piy_e(x, &result)); }
double gsl_sf_psi_1_int(const int n) { EVAL_RESULT(gsl_sf_psi_1_int_e(n, &result)); }
double gsl_sf_psi_1(const double x) { EVAL_RESULT(gsl_sf_psi_1_e(x, &result)); }
double gsl_sf_psi_n(const int n, const double x) { EVAL_RESULT(gsl_sf_psi_n_e(n, x, &result)); }

/*********************************************************
 * INTERNAL ROUTINES: DIGAMMA                             *
 *********************************************************/

/*
 * psi(x0 + t) for t = th + tl, |t| <= 0.55, from the Taylor series about
 * the zero x0 of psi,
 *
 *   psi(x0 + t) = sum_{k=1}^{40} c_k t^k,   c_k = psi^(k)(x0)/k!
 *                                              = (-1)^(k+1) zeta(k+1, x0),
 *
 * which converges for |t| < x0 (pole at 0); the truncation error is
 * below 1e-18 relative for |t| <= 0.55.  The low part tl enters to first
 * order through the derivative of the series.  Relative accuracy is
 * retained up to the zero.
 */
static double
psi_x0_series(const double th, const double tl)
{
  static const double c[40] = {
       9.67672245447621204e-01,   -4.42763168983592081e-01,
       2.58499760955651026e-01,   -1.63942705442406522e-01,
       1.07824050691262371e-01,   -7.21995612564547140e-02,
       4.88042881641431101e-02,   -3.31611264748473619e-02,
       2.25976482322181038e-02,   -1.54247659049489595e-02,
       1.05387916166121750e-02,   -7.20453438635686866e-03,
       4.92678139572985327e-03,   -3.36980165543932821e-03,
       2.30512632673492797e-03,   -1.57693677143019720e-03,
       1.07882520191629667e-03,   -7.38070938996005150e-04,
       5.04953265834601987e-04,   -3.45468025106307692e-04,
       2.36356015640270530e-04,   -1.61706220919748030e-04,
       1.10633727687474103e-04,   -7.56917958219506607e-05,
       5.17857579522208093e-05,   -3.54300709476596036e-05,
       2.42400661186013184e-05,   -1.65842422718541346e-05,
       1.13463845846638497e-05,   -7.76281766846209433e-06,
       5.31106092088986356e-06,   -3.63365078980104556e-06,
       2.48602273312953779e-06,   -1.70085388543326074e-06,
       1.16366753635488430e-06,   -7.96142543124196994e-07,
       5.44694193066944566e-07,   -3.72661612834382273e-07,
       2.54962655202155422e-07,   -1.74436951177277454e-07
  };
  double p = c[39], dp = 40.0 * c[39];
  int k;

  for (k = 38; k >= 0; --k)
    {
      p = p * th + c[k];
      dp = dp * th + (k + 1) * c[k];
    }

  return th * p + tl * dp;
}

/*
 * psi(x) for x > 0:
 *
 *   x < 1:        psi(x) = psi(1+x) - 1/x, with psi(1+x) from the series
 *                 about x0 at the exact offset t = x + 1 - x0
 *   1 <= x <= 2:  the series about x0
 *   2 < x < 10:   psi(x) = psi(y) + sum_{k=0}^{n-1} 1/(y+k), y = x - n in [1,2)
 *   x >= 10:      log x - 1/(2x) - sum_{k=1}^{10} B_2k/(2k) x^(-2k)
 */
static int
psi_pos(const double x, gsl_sf_result * result)
{
  if (x < 1.0)
    {
      /* t = x + 1 - x0 = x - (x0_hi - 1) - x0_lo, with x0_hi - 1 exact */
      const gsl_dd t = gsl_dd_two_sum(x, 1.0 - PSI_X0_HI);
      const double p = psi_x0_series(t.hi, t.lo - PSI_X0_LO);
      const double r = 1.0 / x;

      if (!gsl_finite(r))
        OVERFLOW_ERROR(result);

      result->val = p - r;
      result->err = 2.0 * GSL_DBL_EPSILON * (fabs(p) + r + fabs(result->val));
      return GSL_SUCCESS;
    }
  else if (x <= 2.0)
    {
      result->val = psi_x0_series(x - PSI_X0_HI, -PSI_X0_LO);    /* x - x0_hi is exact */
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
  else if (x < PSI_ASYMP_MIN)
    {
      const int n = (int) x - 1;
      const double y = x - n;               /* exact, in [1,2) */
      const double p = psi_x0_series(y - PSI_X0_HI, -PSI_X0_LO);
      double s = 0.0;
      int k;

      for (k = n - 1; k >= 0; --k)
        s += 1.0 / (y + k);

      result->val = p + s;
      result->err = 2.0 * GSL_DBL_EPSILON * (fabs(p) + s + fabs(result->val));
      return GSL_SUCCESS;
    }
  else
    {
      /* B_2k/(2k), k = 1..10 */
      static const double b[10] = {
           8.33333333333333287e-02,   -8.33333333333333322e-03,
           3.96825396825396803e-03,   -4.16666666666666661e-03,
           7.57575757575757597e-03,   -2.10927960927960942e-02,
           8.33333333333333287e-02,   -4.43259803921568607e-01,
           3.05395433027011975e+00,   -2.64562121212121220e+01
      };
      const double w = 1.0 / (x * x);
      const double lx = log(x);
      double s = b[9];
      int k;

      for (k = 8; k >= 0; --k)
        s = s * w + b[k];

      s = s * w + 0.5 / x;
      result->val = lx - s;
      result->err = 2.0 * GSL_DBL_EPSILON * (fabs(lx) + fabs(s));
      return GSL_SUCCESS;
    }
}

/*
 * psi(x) for x < 0, x not an integer, by reflection,
 *
 *   psi(x) = psi(1-x) - pi cot(pi x),
 *
 * with cot(pi x) from sin_pi/cos_pi of the exact argument and 1-x formed
 * exactly (its low part enters through psi_1).  Near the zeros of psi
 * the two terms cancel; when they do by more than PSI_REFINE_AMP the
 * result is recomputed in double-double.
 */
static int
psi_neg(const double x, gsl_sf_result * result)
{
  const gsl_dd u = gsl_dd_two_sum(1.0, -x);
  const double c = M_PI * gsl_sf_cos_pi(x) / gsl_sf_sin_pi(x);
  gsl_sf_result p;
  double t;
  int status;

  status = psi_pos(u.hi, &p);
  if (status)
    {
      result->val = GSL_NAN;
      result->err = GSL_NAN;
      return status;
    }

  if (u.lo != 0.0)
    p.val += psi1_pos(u.hi) * u.lo;

  if (!gsl_finite(c))
    OVERFLOW_ERROR(result);

  result->val = p.val - c;
  t = fabs(p.val) + fabs(c);

  if (t > PSI_REFINE_AMP * fabs(result->val))
    {
      double mag;
      const gsl_dd v = gsl_dd_sub(psi_dd(u, &mag), psi_dd_pi_cot_pi(x));

      result->val = v.hi + v.lo;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val) + PSI_DD_REL_ERR * (mag + fabs(c));
    }
  else
    {
      result->err = p.err + 3.0 * GSL_DBL_EPSILON * (fabs(c) + fabs(result->val));
    }

  return GSL_SUCCESS;
}

/*
 * psi(x) in double-double for x = x.hi + x.lo > 1: shift to
 * z = x + n >= PSI_DD_ASYMP_MIN, psi(x) = psi(z) - sum_{k<n} 1/(x+k), with
 *
 *   psi(z) = log z - 1/(2z) - sum_{k=1}^{15} B_2k/(2k) z^(-2k)
 *
 * (first omitted term below 4e-34 relative for z >= 20).  The shift
 * sum is formed in pairs, 1/a + 1/(a+1) = (2a+1)/(a(a+1)), halving the
 * number of divisions.  *mag receives the magnitude of the terms.
 */
static gsl_dd
psi_dd(const gsl_dd x, double * mag)
{
  /* B_2k/(2k), k = 1..15, as double-double pairs */
  static const double bdd[15][2] = {
    { 8.33333333333333287e-02, 4.62592926927148533e-18 },
    { -8.33333333333333322e-03, -1.15648231731787138e-19 },
    { 3.96825396825396803e-03, 2.20282346155785016e-19 },
    { -4.16666666666666661e-03, -5.78241158658935690e-20 },
    { 7.57575757575757597e-03, -2.10269512239612987e-19 },
    { -2.10927960927960942e-02, 1.39116773995307320e-18 },
    { 8.33333333333333287e-02, 4.62592926927148533e-18 },
    { -4.43259803921568607e-01, -2.04629341793656316e-17 },
    { 3.05395433027011975e+00, -1.08827208206086071e-17 },
    { -2.64562121212121220e+01, 7.44993292645438342e-16 },
    { 2.81460144927536248e+02, -1.64763532929878310e-14 },
    { -3.60751054639804624e+03, -1.53470290335798157e-13 },
    { 5.48275833333333358e+04, -2.42531920472780850e-12 },
    { -9.74936823850574670e+05, -4.22841858589786249e-11 },
    { 2.00526957966880798e+07, -8.87840912727505469e-10 }
  };
  const gsl_dd one = { 1.0, 0.0 };
  const int n = (x.hi < PSI_DD_ASYMP_MIN) ? (int) ceil(PSI_DD_ASYMP_MIN - x.hi) : 0;
  gsl_dd s = { 0.0, 0.0 }, z, w, w2, a, bk;
  int k;

  for (k = 0; k + 1 < n; k += 2)
    {
      const gsl_dd a = gsl_dd_add_d(x, (double) k);
      s = gsl_dd_add(s, gsl_dd_div(gsl_dd_add_d(gsl_dd_mul_d(a, 2.0), 1.0),
                                   gsl_dd_mul(a, gsl_dd_add_d(a, 1.0))));
    }
  if (k < n)
    s = gsl_dd_add(s, gsl_dd_div(one, gsl_dd_add_d(x, (double) k)));

  z = gsl_dd_add_d(x, (double) n);
  w = gsl_dd_div(one, z);
  w2 = gsl_dd_mul(w, w);

  a.hi = bdd[14][0];
  a.lo = bdd[14][1];
  for (k = 13; k >= 0; --k)
    {
      bk.hi = bdd[k][0];
      bk.lo = bdd[k][1];
      a = gsl_dd_add(gsl_dd_mul(a, w2), bk);
    }

  a = gsl_dd_add(gsl_dd_mul(a, w2), gsl_dd_mul_d(w, 0.5));     /* 1/(2z) + series */

  {
    const gsl_dd lz = gsl_dd_log(z);

    *mag = fabs(lz.hi) + fabs(s.hi);
    return gsl_dd_sub(gsl_dd_sub(lz, a), s);
  }
}

/*
 * pi cot(pi x) in double-double for x not an integer, from sin(pi x)
 * and cos(pi x).
 */
static gsl_dd
psi_dd_pi_cot_pi(const double x)
{
  const gsl_dd pi = { M_PI, GSL_SF_PI_LO };
  gsl_dd xd, sn, cs;

  xd.hi = x;
  xd.lo = 0.0;
  gsl_dd_sincos_pi(xd, &sn, &cs);

  return gsl_dd_div(gsl_dd_mul(pi, cs), sn);
}

/*********************************************************
 * INTERNAL ROUTINES: TRIGAMMA                            *
 *********************************************************/

/*
 * psi_1(x) for x > 0: for x < 10, psi_1(x) = psi_1(x+n) + sum_{k<n} 1/(x+k)^2
 * with x + n >= 10, and
 *
 *   psi_1(z) = 1/z + 1/(2z^2) + sum_{k=1}^{10} B_2k z^(-2k-1)
 *
 * for z >= 10.  All terms are positive; the smallest are added first.
 * Returns inf on overflow (x below ~1e-154).
 */
static double
psi1_pos(const double x)
{
  /* B_2k, k = 1..10 */
  static const double b1[10] = {
       1.66666666666666657e-01,   -3.33333333333333329e-02,
       2.38095238095238082e-02,   -3.33333333333333329e-02,
       7.57575757575757597e-02,   -2.53113553113553102e-01,
       1.16666666666666674e+00,   -7.09215686274509771e+00,
       5.49711779448621556e+01,   -5.29124242424242425e+02
  };
  double z = x, s = 0.0, w, a;
  int n = 0, k;

  if (x < PSI_ASYMP_MIN)
    {
      n = (int) ceil(PSI_ASYMP_MIN - x);
      z = x + n;
    }

  w = 1.0 / (z * z);
  a = b1[9];
  for (k = 8; k >= 0; --k)
    a = a * w + b1[k];

  s = (a * w + 0.5 / z + 1.0) / z;

  for (k = n - 1; k >= 0; --k)
    {
      const double y = x + k;
      s += 1.0 / (y * y);
    }

  return s;
}
