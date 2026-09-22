/* specfunc/gamma/gamma.c
 *
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

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>

#include "check.h"
#include "error.h"
#include "eval.h"

#include "chebyshev.h"
#include "cheb_eval.c"

#include "gamma_kernels.c"
#include "gamma_fallback.c"
#include "gammastar_kernels.c"

#ifdef HAVE_TGAMMA
static int gamma_system_e(const double x, gsl_sf_result * result);
#endif

/*
 * Scaled factorial table
 *
 * Represents
 *
 *   n! = fact_mant[n] * 2^fact_exp[n]
 *
 * with 0.5 <= fact_mant[n] < 1, for 0 <= n <= GSL_SF_FACT_TABLE_NMAX.
 * The scaled representation permits factorial information beyond the
 * double-precision overflow limit and is shared by the factorial,
 * double-factorial, logarithmic-factorial, and Taylor coefficient routines.
 */

#define GSL_SF_FACT_TABLE_NMAX 300
static const double fact_mant[GSL_SF_FACT_TABLE_NMAX + 1] = {
  5.00000000000000000000e-01, 5.00000000000000000000e-01, 5.00000000000000000000e-01,
  7.50000000000000000000e-01, 7.50000000000000000000e-01, 9.37500000000000000000e-01,
  7.03125000000000000000e-01, 6.15234375000000000000e-01, 6.15234375000000000000e-01,
  6.92138671875000000000e-01, 8.65173339843750000000e-01, 5.94806671142578125000e-01,
  8.92210006713867187500e-01, 7.24920630455017089844e-01, 6.34305551648139953613e-01,
  5.94661454670131206512e-01, 5.94661454670131206512e-01, 6.31827795587014406919e-01,
  7.10806270035391207784e-01, 8.44082445667027059244e-01, 5.27551528541891912027e-01,
  6.92411381211233134536e-01, 9.52065649165445559987e-01, 6.84297185337664037874e-01,
  5.13222889003247972894e-01, 8.01910764067574999281e-01, 6.51552495804904707732e-01,
  5.49747418335388315924e-01, 9.62057982086929608379e-01, 8.71865046266279919429e-01,
  8.17373480874637459159e-01, 7.91830559597304972641e-01, 7.91830559597304972641e-01,
  8.16575264584720739158e-01, 8.67611218621265889439e-01, 9.48949770367009559635e-01,
  5.33784245831442905050e-01, 6.17188034242605820801e-01, 7.32910790663094391384e-01,
  8.93235026120646313785e-01, 5.58271891325403890605e-01, 7.15285860760673797287e-01,
  9.38812692248384372817e-01, 6.30764777604383208853e-01, 8.67301569206026967684e-01,
  6.09821415847987657877e-01, 8.76618285281482334526e-01, 6.43766553253588558192e-01,
  9.65649829880382837288e-01, 7.39325651002168071635e-01, 5.77598164845443795556e-01,
  9.20547075222426069985e-01, 7.47944498618221209618e-01, 6.19391537918214440950e-01,
  5.22611610118493441490e-01, 8.98238704891160644195e-01, 7.85958866779765563670e-01,
  6.99994615725728630551e-01, 6.34370120501441592253e-01, 5.84809954837266499084e-01,
  5.48259332659937315135e-01, 5.22559676441502807265e-01, 5.06229686552705837599e-01,
  9.96639695400639569201e-01, 9.96639695400639569201e-01, 5.06106095320637328072e-01,
  5.21921910799407196002e-01, 5.46387000368129371886e-01, 5.80536187891137478445e-01,
  6.25890577570132555785e-01, 6.84567819217332562687e-01, 7.59442424444228270097e-01,
  8.54372727499756789982e-01, 9.74518892304410133676e-01, 5.63393734613487073837e-01,
  6.60227032750180242715e-01, 7.84019601390839038224e-01, 9.43273582923353171026e-01,
  5.74807339593918276144e-01, 7.09527809811242970994e-01, 8.86909762264053713743e-01,
  5.61247583932721449962e-01, 7.19098466913799416744e-01, 9.32580824278833597774e-01,
  6.12006165932984513844e-01, 8.12820689129745033163e-01, 5.46113900509047511811e-01,
  7.42373583504486456164e-01, 5.10381838659334396979e-01, 7.09749744385636871513e-01,
  9.98085578042301913015e-01, 7.09576465639448983325e-01, 5.10008084678354012276e-01,
  7.41105498048233113373e-01, 5.44249350129171194368e-01, 8.07870129097988498579e-01,
  6.05902596823491346179e-01, 9.18321123310604114920e-01, 7.03089610034681289363e-01,
  5.43795870261198799334e-01, 8.49681047283123103142e-01, 6.70451451371839368676e-01,
  5.34266000311934496914e-01, 8.59834344252019544541e-01, 6.98615404704765907695e-01,
  5.73082949171878297534e-01, 9.49168634565923374780e-01, 7.93445655457451604065e-01,
  6.69469771792224777052e-01, 5.70095352541816424718e-01, 9.79851387181246979985e-01,
  8.49714874821237531322e-01, 7.43500515468582867662e-01, 6.56371548812108351889e-01,
  5.84580910660783992228e-01, 5.25209411921798152711e-01, 9.51942059108259019951e-01,
  8.70134538403643076343e-01, 8.02155277590858473147e-01, 7.45753734635251186447e-01,
  6.99144126220548001172e-01, 6.60909681817861804909e-01, 6.29929540482649485966e-01,
  6.05322917807545990421e-01, 5.86406576626060216384e-01, 5.72662672486386914450e-01,
  5.63714818228787106769e-01, 5.59310796211374694487e-01, 5.59310796211374694487e-01,
  5.63680411806776082706e-01, 5.72487918241256998897e-01, 5.85905603825036425114e-01,
  6.04215153944568816868e-01, 6.27817308395528517195e-01, 6.57246244726568895622e-01,
  6.93189398735053141642e-01, 7.36513736155994025445e-01, 7.88299858229462357961e-01,
  8.49885784653639064778e-01, 9.22922844272311193592e-01, 5.04723430461420208282e-01,
  5.55984403867658127929e-01, 6.16795198040683279039e-01, 6.89075885311075886364e-01,
  7.75210370974960261137e-01, 8.78167998370072200309e-01, 5.00830186570431767912e-01,
  5.75172167389480293309e-01, 6.65042818544086533628e-01, 7.74151405961475802542e-01,
  9.07208678861104456104e-01, 5.35111369171979545101e-01, 6.35444750891725695929e-01,
  7.59555053800265889130e-01, 9.13839674103444910003e-01, 5.53301365179820203721e-01,
  6.74336038812905824713e-01, 8.27115297606454857871e-01, 5.10485222741483890552e-01,
  6.34118362624186970855e-01, 7.92647953280233741324e-01, 9.97002503735293976739e-01,
  6.30915646894990689297e-01, 8.03431644092839736260e-01, 5.14698396996975393591e-01,
  6.63478402378913645343e-01, 8.60448553085153666764e-01, 5.61308235801643240492e-01,
  7.36717059489656711513e-01, 9.72696742607437325745e-01, 6.45931430637751380353e-01,
  8.62924020617620968210e-01, 5.79777076352464115772e-01, 7.83604954757627258921e-01,
  5.32606492686824806171e-01, 7.28172939220268267135e-01, 5.00618895713934364267e-01,
  6.92262066729424918599e-01, 9.62676936545606598550e-01, 6.73121764225248386815e-01,
  9.46577480941755422528e-01, 6.69259859572100546110e-01, 9.51603862829080493491e-01,
  6.80248073819225473891e-01, 9.77856606115136584023e-01, 7.06654188012891704140e-01,
  5.13428433478116619249e-01, 7.50086852034436035730e-01, 5.50845031962789000168e-01,
  8.13357117507555571478e-01, 6.03663485650138942340e-01, 9.00779107493566599807e-01,
  6.75584330620174977611e-01, 5.09327249256616321205e-01, 7.71949112154559080601e-01,
  5.88008112773980573529e-01, 9.00387422685157701174e-01, 6.92876258863187710624e-01,
  5.35896481464496754477e-01, 8.33151561026834852619e-01, 6.50899657052214752895e-01,
  5.11057933857402946209e-01, 8.06513301868714060916e-01, 6.39539844841206850390e-01,
  5.09633313857836722782e-01, 8.16209604225441642811e-01, 6.56793665900160084092e-01,
  5.31079253286457553251e-01, 8.63003786590493482400e-01, 7.04561685146145033620e-01,
  5.77960757346447140392e-01, 9.52732185938283904214e-01, 7.88981341480141429301e-01,
  6.56457131778398927935e-01, 5.48757133596005308718e-01, 9.21740497837040151374e-01,
  7.77718545050002707519e-01, 6.59237985452541308540e-01, 5.61382346986929681165e-01,
  9.60490109297950089839e-01, 8.25421187677925827231e-01, 7.12570634675084413345e-01,
  6.17932347257299730003e-01, 5.38277005618663428876e-01, 9.41984759832661056045e-01,
  8.27916292821674781877e-01, 7.30894852256634752763e-01, 6.48098169774437793400e-01,
  5.77212432455358692707e-01, 5.16334558719832581453e-01, 9.27788660199699188880e-01,
  8.37184298852072350527e-01, 7.58698270834690546849e-01, 6.90533973064386352725e-01,
  6.31191209754165649670e-01, 5.79413805829019246119e-01, 5.34147102248627048127e-01,
  9.89006744007223548287e-01, 9.19467207319215718009e-01, 8.58408838083174030587e-01,
  8.04758285702975584286e-01, 7.57604479900066851705e-01, 7.16172984905532050703e-01,
  6.79804825515797905844e-01, 6.47938974319744898089e-01, 6.20097846516943329576e-01,
  5.95875274387375308827e-01, 5.74926534272194089681e-01, 5.56960080076188024378e-01,
  5.41730702886604742474e-01, 5.29033889537700030559e-01, 5.18701196382666762652e-01,
  5.10596490189187623976e-01, 5.04612937569783093217e-01, 5.00670648995019118566e-01,
  9.97429808544764706646e-01, 9.97429808544764706646e-01, 5.00663009367196321087e-01,
  5.04574439127877583111e-01, 5.10487420836407412850e-01, 5.18463786786976243981e-01,
  5.28590032622659355965e-01, 5.40978861512252939825e-01, 5.55771252256728631735e-01,
  5.73139103889751422294e-01, 5.93288525510875408564e-01, 6.16463858538644049467e-01,
  6.42952539960226387983e-01, 6.73090940270862025940e-01, 7.07271339581491687021e-01,
  7.45950240964854516257e-01, 7.89658262896388962737e-01, 8.39011904327413349236e-01,
  8.94727538599155547416e-01, 9.57638068656908703069e-01, 5.14356384532519306241e-01,
  5.54540477074122351020e-01, 6.00030125584109064185e-01, 6.51595214501493380510e-01,
  7.10136972054361881312e-01, 7.76712313184458325033e-01, 8.52563125018878098516e-01,
  9.39151567403607945295e-01, 5.19101354639103562683e-01, 5.75878065302755537402e-01,
  6.41114252387833372104e-01, 7.16244828839532554632e-01, 8.02977601081819702955e-01,
  9.03349801217047110313e-01, 5.09898618265091085888e-01, 5.77619528503423462773e-01,
  6.56590948416000852461e-01, 7.48924050537001062544e-01, 8.57166979716176946624e-01,
  9.84402703267796996833e-01, 5.67185151296875189608e-01, 6.55807831187011913698e-01,
  7.60839554150556818257e-01, 8.85664793503382563600e-01, 5.17214400893576953244e-01,
  6.06110626047160483409e-01
};

static const short fact_exp[GSL_SF_FACT_TABLE_NMAX + 1] = {
  1, 1, 2, 3, 5, 7, 10, 13, 16, 19, 22, 26,
  29, 33, 37, 41, 45, 49, 53, 57, 62, 66, 70, 75,
  80, 84, 89, 94, 98, 103, 108, 113, 118, 123, 128, 133,
  139, 144, 149, 154, 160, 165, 170, 176, 181, 187, 192, 198,
  203, 209, 215, 220, 226, 232, 238, 243, 249, 255, 261, 267,
  273, 279, 285, 290, 296, 303, 309, 315, 321, 327, 333, 339,
  345, 351, 358, 364, 370, 376, 383, 389, 395, 402, 408, 414,
  421, 427, 434, 440, 447, 453, 459, 466, 473, 479, 486, 492,
  499, 505, 512, 519, 525, 532, 539, 545, 552, 559, 565, 572,
  579, 586, 592, 599, 606, 613, 620, 627, 633, 640, 647, 654,
  661, 668, 675, 682, 689, 696, 703, 710, 717, 724, 731, 738,
  745, 752, 759, 766, 773, 780, 787, 794, 802, 809, 816, 823,
  830, 837, 845, 852, 859, 866, 873, 881, 888, 895, 902, 910,
  917, 924, 932, 939, 946, 953, 961, 968, 976, 983, 990, 998,
  1005, 1012, 1020, 1027, 1035, 1042, 1050, 1057, 1065, 1072, 1079, 1087,
  1094, 1102, 1109, 1117, 1124, 1132, 1140, 1147, 1155, 1162, 1170, 1177,
  1185, 1193, 1200, 1208, 1215, 1223, 1231, 1238, 1246, 1254, 1261, 1269,
  1277, 1284, 1292, 1300, 1307, 1315, 1323, 1330, 1338, 1346, 1354, 1361,
  1369, 1377, 1385, 1392, 1400, 1408, 1416, 1424, 1431, 1439, 1447, 1455,
  1463, 1471, 1478, 1486, 1494, 1502, 1510, 1518, 1526, 1533, 1541, 1549,
  1557, 1565, 1573, 1581, 1589, 1597, 1605, 1613, 1621, 1629, 1637, 1645,
  1653, 1661, 1669, 1676, 1684, 1693, 1701, 1709, 1717, 1725, 1733, 1741,
  1749, 1757, 1765, 1773, 1781, 1789, 1797, 1805, 1813, 1821, 1829, 1838,
  1846, 1854, 1862, 1870, 1878, 1886, 1894, 1903, 1911, 1919, 1927, 1935,
  1943, 1952, 1960, 1968, 1976, 1984, 1992, 2001, 2009, 2017, 2025, 2034,
  2042
};

/*
 * Gamma function
 *
 * Computes
 *
 *   Gamma(x)
 *
 * over the real line.  Negative integers and zero are poles and return
 * GSL_EDOM.  Results which overflow or underflow the double-precision
 * range return GSL_EOVRFLW or GSL_EUNDRFLW, respectively.
 *
 * If the system tgamma() function is available, it is used directly.
 * Otherwise gamma_fallback_e() provides a portable implementation using
 * the Gamma kernels in gamma_kernels.c.
 */

int
gsl_new_sf_gamma_e(const double x, gsl_sf_result * result)
{
#ifdef HAVE_TGAMMA

  /* use system tgamma() */
  return gamma_system_e(x, result);

#else

  return gamma_fallback_e(x, result);

#endif
}

/*
 * Reciprocal gamma function
 *
 * Computes
 *
 *   1 / Gamma(x)
 *
 * over the real line.  Several representations are used according to
 * the argument:
 *
 *   x < -10
 *     The reflection formula
 *
 *       1/Gamma(x) = (-x) sin(pi*x) Gamma(-x) / pi
 *
 *     is evaluated using sin(pi*x) and a scaled positive Gamma kernel.
 *     Negative integer poles of Gamma are returned as exact zeros.
 *
 *   -10 <= x < -1/2
 *     The argument is shifted into the central interval and the
 *     recurrence relation is used.
 *
 *   -1/2 <= x <= 2
 *     A rational approximation to the auxiliary function H(x) is used
 *     to evaluate reciprocal Gamma near unity without loss of accuracy.
 *
 *   x > 2
 *     A scaled positive Gamma representation is inverted directly,
 *     avoiding intermediate overflow of Gamma(x).
 */

int
gsl_new_sf_gammainv_e(const double x, gsl_sf_result * result)
{
  int status;
  gsl_sf_result tmp;
  gsl_sf_result gm;
  double r, base, p;
  double d, ed, term;
  double ebase, eprod;
  double qmant, eqmant;
  int n, k;
  int ge, qexp;

  if (x < -10.0)
    {
      /* Reflection: 1/Gamma(x) = (-x) sin(pi*x) Gamma(-x) / pi */
      status = gsl_sf_sin_pi_e(x, &tmp);
      if (status)
        return status;

      /* Gamma poles: 1/Gamma(-n) = 0 exactly */
      if (tmp.val == 0.0)
        {
          result->val = 0.0;
          result->err = 0.0;
          return GSL_SUCCESS;
        }

      /* The remaining nonintegral negative tail eventually overflows */
      if (-x > GAMMA_SCALED_X_MAX)
        {
          result->val = tmp.val < 0.0 ? GSL_NEGINF : GSL_POSINF;
          result->err = GSL_POSINF;
          return GSL_EOVRFLW;
        }

      status = gamma_2_inf_scaled(-x, &gm, &ge);
      if (status)
        return status;

      /*
       * qmant = (-x) sin(pi*x) Gamma(-x) / pi, with the Gamma
       * factor represented in scaled form.
       */
      qmant = ((-x) * tmp.val / M_PI) * gm.val;
      qmant = frexp(qmant, &qexp);
      qexp += ge;

      result->val = ldexp(qmant, qexp);

      if (gsl_isinf(result->val))
        {
          result->err = GSL_POSINF;
          return GSL_EOVRFLW;
        }

      /* Propagate the sin(pi*x) and scaled Gamma errors */
      eqmant = fabs((-x) / M_PI) * tmp.err * fabs(gm.val) +
               fabs(((-x) * tmp.val / M_PI)) * gm.err +
               3.0 * GSL_DBL_EPSILON * fabs(qmant);

      result->err = ldexp(eqmant, qexp);

      return GSL_SUCCESS;
    }
  else if (x < -0.5)
    {
      /* Moderate negative x: 1/Gamma(x) = x...(x+n-1) / Gamma(x+n) */
      n = (int) ceil(-x - 0.5);
      if (n < 1)
        n = 1;

      r = x + (double) n;

      status = gamma1p_H(r, &tmp);
      if (status)
        return status;

      d = 1.0 + r * (r - 1.0) * tmp.val;
      ed = fabs(r * (r - 1.0)) * tmp.err +
           2.0 * GSL_DBL_EPSILON *
           fabs(r * (r - 1.0) * tmp.val) +
           GSL_DBL_EPSILON * fabs(d);

      base = r / d;
      ebase = fabs(base) * (GSL_DBL_EPSILON + ed / fabs(d));

      p = 1.0;
      for (k = 0; k < n; ++k)
        p *= x + (double) k;

      eprod = (double) n * GSL_DBL_EPSILON * fabs(p);

      result->val = p * base;
      result->err = fabs(base) * eprod +
                    fabs(p) * ebase +
                    GSL_DBL_EPSILON * fabs(result->val);

      return GSL_SUCCESS;
    }
  else if (x <= 1.0)
    {
      /* Central H representation: 1/Gamma(x) = x / [1 + x(x-1) H(x)] */
      status = gamma1p_H(x, &tmp);
      if (status)
        return status;

      d = 1.0 + x * (x - 1.0) * tmp.val;
      ed = fabs(x * (x - 1.0)) * tmp.err +
           2.0 * GSL_DBL_EPSILON *
           fabs(x * (x - 1.0) * tmp.val) +
           GSL_DBL_EPSILON * fabs(d);

      result->val = x / d;
      result->err = fabs(result->val) * (GSL_DBL_EPSILON + ed / fabs(d));

      return GSL_SUCCESS;
    }
  else if (x <= 2.0)
    {
      /* 1/Gamma(x) = 1 / [1 + r(r-1)H(r)], r = x-1 */
      r = x - 1.0;

      status = gamma1p_H(r, &tmp);
      if (status)
        return status;

      term = r * (r - 1.0) * tmp.val;
      d = 1.0 + term;
      ed = fabs(r * (r - 1.0)) * tmp.err +
           2.0 * GSL_DBL_EPSILON * fabs(term) +
           GSL_DBL_EPSILON * fabs(d);

      result->val = 1.0 / d;
      result->err = fabs(result->val) * (GSL_DBL_EPSILON + ed / fabs(d));

      return GSL_SUCCESS;
    }
  else if (x <= GAMMA_SCALED_X_MAX)
    {
      status = gamma_2_inf_scaled(x, &gm, &ge);
      if (status)
        return status;

      result->val = ldexp(1.0 / gm.val, -ge);

      if (fabs(result->val) < GSL_DBL_MIN)
        result->err = GSL_SF_TRUE_MIN;
      else
        result->err = fabs(result->val) * (GSL_DBL_EPSILON + gm.err / fabs(gm.val));

      return GSL_SUCCESS;
    }
  else
    {
      /* 1/Gamma(x) has underflowed to zero. */
      result->val = 0.0;
      result->err = GSL_SF_TRUE_MIN;
      return GSL_SUCCESS;
    }
}

/*
 * Scaled Gamma function
 *
 * Computes
 *
 *   Gamma*(x) = Gamma(x) / [sqrt(2*pi/x) (x/e)^x]
 *
 * for x > 0.  The implementation uses several representations according
 * to the argument:
 *
 *   0 < x < 0.4
 *     Remove the square-root singularity at the origin and evaluate
 *
 *       Gamma*(x) = exp[log(Gamma(1+x)) + x(1-log(x))] / [sqrt(2*pi) sqrt(x)].
 *
 *     log(Gamma(1+x)) is evaluated by the small-argument rational kernel,
 *     and expm1() is used to reconstruct the exponential accurately.
 *
 *   0.4 <= x < 2
 *     Use a Chebyshev approximation to Gamma*(x).
 *
 *   2 <= x < 16
 *     Use a second Chebyshev approximation to Gamma*(x).
 *
 *   x >= 16
 *     Use the direct Stirling asymptotic expansion for Gamma*(x),
 *     truncated at its least term and evaluated with compensated
 *     summation.
 *
 * Gamma*(x) is defined here only for x > 0; nonpositive arguments
 * return GSL_EDOM.
 */

int
gsl_new_sf_gammastar_e(const double x, gsl_sf_result * result)
{
  if (x <= 0.0)
    DOMAIN_ERROR(result);

  if (x < GAMMASTAR_X_SMALL)
    return gammastar_small_e(x, result);
  else if (x < GAMMASTAR_X_MID)
    return gammastar_mid1_e(x, result);
  else if (x < GAMMASTAR_X_ASYM)
    return gammastar_mid2_e(x, result);
  else
    return gammastar_asym_e(x, result);
}

/*
 * Taylor coefficient
 *
 * Computes
 *
 *   x^n / n!
 *
 * for integer n >= 0 and x >= 0.
 *
 * The cases n = 0, n = 1, and x = 0 are returned exactly.  For
 * 2 <= n <= TAYLORCOEFF_POW_NMAX, a scaled power representation is
 * used to reduce accumulated rounding error and avoid intermediate
 * overflow in x^n.  Larger orders use the general product formulation,
 * with logarithmic tests for underflow and overflow.
 */

int
gsl_new_sf_taylorcoeff_e(const int n, const double x, gsl_sf_result * result)
{
  if (x < 0.0 || n < 0)
    DOMAIN_ERROR(result);
  else if (n == 0)
    {
      result->val = 1.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
  else if (n == 1)
    {
      result->val = x;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
  else if (x == 0.0)
    {
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
  else if (n <= GSL_SF_FACT_TABLE_NMAX)
    {
      int e, scale;
      double m, p, y;

      m = frexp(x, &e);
      p = pow(m, (double) n);
      p /= fact_mant[n];
      scale = n * e - (int) fact_exp[n];
      y = ldexp(p, scale);

      if (gsl_isinf(y))
        {
          result->val = y;
          result->err = GSL_POSINF;
          return GSL_EOVRFLW;
        }
      else if (y == 0.0)
        {
          result->val = 0.0;
          result->err = GSL_SF_TRUE_MIN;
          return GSL_EUNDRFLW;
        }
      else
        {
          result->val = y;

          if (fabs(y) < GSL_DBL_MIN)
            result->err = 2.0 * GSL_SF_TRUE_MIN;
          else
            result->err = 8.0 * GSL_DBL_EPSILON * fabs(y);

          return GSL_SUCCESS;
        }
    }
  else
    {
      const double log2pi = M_LNPI + M_LN2;
      const double ln_test = n*(log(x) + 1.0) +
                             1.0 -
                             (n + 0.5)*log(n + 1.0) +
                             0.5*log2pi;

      if (ln_test < GSL_LOG_DBL_MIN + 1.0)
        UNDERFLOW_ERROR(result);
      else if (ln_test > GSL_LOG_DBL_MAX - 1.0)
        OVERFLOW_ERROR(result);
      else
        {
          double product = 1.0;
          int k;

          for (k = 1; k <= n; ++k)
            product *= (x / (double) k);

          result->val = product;
          result->err = n * GSL_DBL_EPSILON * product;
          CHECK_UNDERFLOW(result);
          return GSL_SUCCESS;
        }
    }
}

/*
 * Factorial
 *
 * Computes
 *
 *   n!
 *
 * for integer n >= 0.  Factorials through GSL_SF_FACT_NMAX are
 * reconstructed from the scaled factorial table.  Larger values
 * overflow the double-precision range and return GSL_EOVRFLW.
 */

int
gsl_new_sf_fact_e(const unsigned int n, gsl_sf_result * result)
{
  if (n > GSL_SF_FACT_NMAX)
    OVERFLOW_ERROR(result);
  else if (n > GSL_SF_FACT_TABLE_NMAX)
    {
      return GSL_EOVRFLW;
    }
  else
    {
      result->val = ldexp(fact_mant[n], (int) fact_exp[n]);

      if (n <= 22U)
        result->err = 0.0;
      else
        result->err = GSL_DBL_EPSILON * fabs(result->val);

      return GSL_SUCCESS;
    }
}

/*
 * Double factorial
 *
 * Computes
 *
 *   n!!
 *
 * for integer n >= 0.  Even and odd double factorials are reconstructed
 * from the scaled factorial table using
 *
 *   (2k)!!   = 2^k k!
 *   (2k-1)!! = (2k)! / (2^k k!).
 *
 * Values above GSL_SF_DOUBLEFACT_NMAX overflow the double-precision
 * range and return GSL_EOVRFLW.
 */

int
gsl_new_sf_doublefact_e(const unsigned int n, gsl_sf_result * result)
{
  unsigned int k;
  double mant;
  int e, de;

  if (n > GSL_SF_DOUBLEFACT_NMAX)
    OVERFLOW_ERROR(result);
  else if (n == 0U || n == 1U)
    {
      result->val = 1.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
  else if ((n & 1U) == 0U)
    {
      k = n / 2U;
      mant = fact_mant[k];
      e = (int) fact_exp[k] + (int) k;
    }
  else
    {
      k = (n + 1U) / 2U;
      mant = fact_mant[2U*k] / fact_mant[k];
      e = (int) fact_exp[2U*k] - (int) fact_exp[k] - (int) k;
      mant = frexp(mant, &de);
      e += de;
    }

  result->val = ldexp(mant, e);
  result->err = 4.0 * GSL_DBL_EPSILON * fabs(result->val);

  return GSL_SUCCESS;
}

/*
 * Logarithm of factorial
 *
 * Computes
 *
 *   log(n!)
 *
 * for integer n >= 0.  For n <= GSL_SF_FACT_TABLE_NMAX the logarithm
 * is reconstructed directly from the scaled factorial table.  Larger
 * arguments are evaluated using the log-Gamma function.
 */

int
gsl_new_sf_lnfact_e(const unsigned int n, gsl_sf_result * result)
{
  if (n <= GSL_SF_FACT_TABLE_NMAX)
    {
      const double e2 = (double) fact_exp[n] * M_LN2;
      result->val = log(fact_mant[n]) + e2;
      result->err = 3.0 * GSL_DBL_EPSILON * fabs(result->val) +
                    2.0 * GSL_DBL_EPSILON * fabs(e2);
      return GSL_SUCCESS;
    }
  else
    return gsl_sf_lngamma_e(n + 1.0, result);
}

/*
 * Logarithm of double factorial
 *
 * Computes
 *
 *   log(n!!)
 *
 * for integer n >= 0.  For n <= GSL_SF_FACT_TABLE_NMAX the result is
 * reconstructed from the scaled factorial table using the even and odd
 * double-factorial identities.  Larger arguments are evaluated using
 * the corresponding log-Gamma representations.
 */

int
gsl_new_sf_lndoublefact_e(const unsigned int n, gsl_sf_result * result)
{
  unsigned int k;
  double l2 = M_LN2;
  double v;

  if (n <= GSL_SF_FACT_TABLE_NMAX)
    {
      if (n == 0U || n == 1U)
        {
          result->val = 0.0;
          result->err = 0.0;
          return GSL_SUCCESS;
        }
      else if ((n & 1U) == 0U)
        {
          k = n / 2U;
          v = (double) k * l2 + log(fact_mant[k]) + (double) fact_exp[k] * l2;
        }
      else
        {
          k = (n + 1U) / 2U;
          v = log(fact_mant[2U*k]) - log(fact_mant[k])
              + (double) ((int) fact_exp[2U*k] - (int) fact_exp[k] - (int) k) * l2;
        }

      result->val = v;
      result->err = 5.0 * GSL_DBL_EPSILON * fabs(v)
                   + 2.0 * GSL_DBL_EPSILON * (double) n * l2;

      return GSL_SUCCESS;
    }
  else if (GSL_IS_ODD(n))
    {
      gsl_sf_result lg;
      gsl_sf_lngamma_e(0.5 * (n + 2.0), &lg);
      result->val = 0.5 * (n + 1.0) * l2 - 0.5 * M_LNPI + lg.val;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val) + lg.err;
      return GSL_SUCCESS;
    }
  else
    {
      gsl_sf_result lg;
      gsl_sf_lngamma_e(0.5 * n + 1.0, &lg);
      result->val = 0.5 * n * l2 + lg.val;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val) + lg.err;
      return GSL_SUCCESS;
    }
}

double gsl_new_sf_gamma(const double x) { EVAL_RESULT(gsl_new_sf_gamma_e(x, &result)); }
double gsl_new_sf_gammainv(const double x) { EVAL_RESULT(gsl_new_sf_gammainv_e(x, &result)); }
double gsl_new_sf_gammastar(const double x) { EVAL_RESULT(gsl_new_sf_gammastar_e(x, &result)); }
double gsl_new_sf_taylorcoeff(const int n, const double x) { EVAL_RESULT(gsl_new_sf_taylorcoeff_e(n, x, &result)); }
double gsl_new_sf_fact(const unsigned int n) { EVAL_RESULT(gsl_new_sf_fact_e(n, &result)); }
double gsl_new_sf_doublefact(const unsigned int n) { EVAL_RESULT(gsl_new_sf_doublefact_e(n, &result)); }
double gsl_new_sf_lnfact(const unsigned int n) { EVAL_RESULT(gsl_new_sf_lnfact_e(n, &result)); }
double gsl_new_sf_lndoublefact(const unsigned int n) { EVAL_RESULT(gsl_new_sf_lndoublefact_e(n, &result)); }

/********************************************
 * INTERNAL ROUTINES FOR GAMMA              *
 ********************************************/

#ifdef HAVE_TGAMMA

/*
 * Evaluate Gamma(x) using the system tgamma() function.
 *
 * The system function supplies only the value, so an absolute error
 * estimate is constructed here.  For normal finite results the estimate
 * is
 *
 *   12 eps |Gamma(x)|,
 *
 * while nonzero subnormal results use 2*GSL_SF_TRUE_MIN.  Exact zero
 * results are reported as underflow and infinite results as overflow.
 * Gamma poles at zero and the negative integers are detected before
 * calling tgamma().
 */

static int
gamma_system_e(const double x, gsl_sf_result * result)
{
  if (x == 0.0 || (x < 0.0 && x == floor(x)))
    DOMAIN_ERROR(result);

  result->val = tgamma(x);

  if (gsl_isinf(result->val))
    {
      result->err = GSL_POSINF;
      return GSL_EOVRFLW;
    }
  else if (result->val == 0.0)
    {
      result->err = GSL_SF_TRUE_MIN;
      return GSL_EUNDRFLW;
    }
  else if (fabs(result->val) < GSL_DBL_MIN)
    {
      result->err = 2.0 * GSL_SF_TRUE_MIN;
      return GSL_SUCCESS;
    }
  else
    {
      result->err = 12.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
}

#endif /* HAVE_TGAMMA */
