/* roots/test_aps.c
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

/* The test problems of

   G. E. Alefeld, F. A. Potra and Y. Shi, "Algorithm 748: Enclosing
   zeros of continuous functions", ACM Transactions on Mathematical
   Software 21, 327-344 (1995), Table I:

   15 functions with the brackets and parameter values of the paper,
   154 cases in all. The reference roots were computed in high
   precision and rounded to double precision.

   The paper gives no starting points for methods without a bracket.
   For the fdfsolvers we start from the root perturbed by 20%, x0 = 1.2
   root (x0 = 0.2 for a root at 0), or from the midpoint between the
   root and the upper end of the bracket if 1.2 root lies outside it.
   Problem 13 is not used for the fdfsolvers: its root is flat to all
   orders, so that Newton-type steps x -> 2x/(x^2 + 2) converge
   arbitrarily slowly.

   This file is included by test.c after test_funcs.c. */

#define APS_N(p) (((const root_test_params *) (p))->n)
#define APS_A(p) (((const root_test_params *) (p))->a)
#define APS_B(p) (((const root_test_params *) (p))->b)

static double aps1_f  (double x, void * p) { U(p); return sin (x) - x / 2.0; }
static double aps1_df (double x, void * p) { U(p); return cos (x) - 0.5; }

static double aps2_f (double x, void * p)
{
  double sum = 0.0;
  int i;

  U(p);

  for (i = 1; i <= 20; i++)
    {
      const double c = 2.0 * i - 5.0;
      sum += c * c / gsl_pow_3 (x - (double) (i * i));
    }

  return -2.0 * sum;
}

static double aps2_df (double x, void * p)
{
  double sum = 0.0;
  int i;

  U(p);

  for (i = 1; i <= 20; i++)
    {
      const double c = 2.0 * i - 5.0;
      sum += c * c / gsl_pow_4 (x - (double) (i * i));
    }

  return 6.0 * sum;
}

static double aps3_f  (double x, void * p) { return APS_A(p) * x * exp (APS_B(p) * x); }
static double aps3_df (double x, void * p) { return APS_A(p) * (1.0 + APS_B(p) * x) * exp (APS_B(p) * x); }

static double aps4_f  (double x, void * p) { return pow (x, APS_N(p)) - APS_A(p); }
static double aps4_df (double x, void * p) { return APS_N(p) * pow (x, APS_N(p) - 1.0); }

static double aps5_f  (double x, void * p) { U(p); return sin (x) - 0.5; }
static double aps5_df (double x, void * p) { U(p); return cos (x); }

static double aps6_f (double x, void * p)
{
  const double n = APS_N(p);
  return 2.0 * x * exp (-n) - 2.0 * exp (-n * x) + 1.0;
}

static double aps6_df (double x, void * p)
{
  const double n = APS_N(p);
  return 2.0 * exp (-n) + 2.0 * n * exp (-n * x);
}

static double aps7_f (double x, void * p)
{
  const double n = APS_N(p);
  return (1.0 + gsl_pow_2 (1.0 - n)) * x - gsl_pow_2 (1.0 - n * x);
}

static double aps7_df (double x, void * p)
{
  const double n = APS_N(p);
  return (1.0 + gsl_pow_2 (1.0 - n)) + 2.0 * n * (1.0 - n * x);
}

static double aps8_f  (double x, void * p) { return x * x - pow (1.0 - x, APS_N(p)); }
static double aps8_df (double x, void * p) { return 2.0 * x + APS_N(p) * pow (1.0 - x, APS_N(p) - 1.0); }

static double aps9_f (double x, void * p)
{
  const double n = APS_N(p);
  return (1.0 + gsl_pow_4 (1.0 - n)) * x - gsl_pow_4 (1.0 - n * x);
}

static double aps9_df (double x, void * p)
{
  const double n = APS_N(p);
  return (1.0 + gsl_pow_4 (1.0 - n)) + 4.0 * n * gsl_pow_3 (1.0 - n * x);
}

static double aps10_f (double x, void * p)
{
  const double n = APS_N(p);
  return exp (-n * x) * (x - 1.0) + pow (x, n);
}

static double aps10_df (double x, void * p)
{
  const double n = APS_N(p);
  return exp (-n * x) * (1.0 - n * (x - 1.0)) + n * pow (x, n - 1.0);
}

static double aps11_f (double x, void * p)
{
  const double n = APS_N(p);
  return (n * x - 1.0) / ((n - 1.0) * x);
}

static double aps11_df (double x, void * p)
{
  const double n = APS_N(p);
  return 1.0 / ((n - 1.0) * x * x);
}

static double aps12_f (double x, void * p)
{
  const double n = APS_N(p);
  return pow (x, 1.0 / n) - pow (n, 1.0 / n);
}

static double aps12_df (double x, void * p)
{
  const double n = APS_N(p);
  return pow (x, 1.0 / n - 1.0) / n;
}

/* x exp(-1/x^2) is exactly zero in double precision for |x| < 0.037,
   so any point in that interval is a root of the computed function */
static double aps13_f (double x, void * p)
{
  U(p);
  return (x == 0.0) ? 0.0 : x * exp (-1.0 / (x * x));
}

static double aps13_df (double x, void * p)
{
  U(p);
  return (x == 0.0) ? 0.0 : (1.0 + 2.0 / (x * x)) * exp (-1.0 / (x * x));
}

static double aps14_f (double x, void * p)
{
  const double n = APS_N(p);
  return (x >= 0.0) ? n / 20.0 * (x / 1.5 + sin (x) - 1.0) : -n / 20.0;
}

static double aps14_df (double x, void * p)
{
  const double n = APS_N(p);
  return (x >= 0.0) ? n / 20.0 * (1.0 / 1.5 + cos (x)) : 0.0;
}

static double aps15_f (double x, void * p)
{
  const double n = APS_N(p);

  if (x > 2.0e-3 / (1.0 + n))
    return M_E - 1.859;
  else if (x >= 0.0)
    return exp ((n + 1.0) * x / 2.0 * 1000.0) - 1.859;
  else
    return -0.859;
}

static double aps15_df (double x, void * p)
{
  const double n = APS_N(p);

  if (x >= 0.0 && x <= 2.0e-3 / (1.0 + n))
    return (n + 1.0) / 2.0 * 1000.0 * exp ((n + 1.0) * x / 2.0 * 1000.0);
  else
    return 0.0;
}

static const root_test_t aps1_test  = { aps1_f,  aps1_df };
static const root_test_t aps2_test  = { aps2_f,  aps2_df };
static const root_test_t aps3_test  = { aps3_f,  aps3_df };
static const root_test_t aps4_test  = { aps4_f,  aps4_df };
static const root_test_t aps5_test  = { aps5_f,  aps5_df };
static const root_test_t aps6_test  = { aps6_f,  aps6_df };
static const root_test_t aps7_test  = { aps7_f,  aps7_df };
static const root_test_t aps8_test  = { aps8_f,  aps8_df };
static const root_test_t aps9_test  = { aps9_f,  aps9_df };
static const root_test_t aps10_test = { aps10_f, aps10_df };
static const root_test_t aps11_test = { aps11_f, aps11_df };
static const root_test_t aps12_test = { aps12_f, aps12_df };
static const root_test_t aps13_test = { aps13_f, aps13_df };
static const root_test_t aps14_test = { aps14_f, aps14_df };
static const root_test_t aps15_test = { aps15_f, aps15_df };

static const root_test_function test_aps_functions[] =
{
  /* function   params {n, a, b}    xmin  xmax  x0  root  flags */
  { &aps1_test, { 0.0, 0.0, 0.0 }, M_PI_2, M_PI, 2.274593120440777, 1.895494267033981, 0, "sin(x) - x/2" },
  { &aps2_test, { 1.0, 0.0, 0.0 }, 1.0 + 1.0e-9, 4.0 - 1.0e-9, 3.6274984167276685, 3.0229153472730568, 0, "-2 sum (2i-5)^2/(x-i^2)^3, n=1" },
  { &aps2_test, { 2.0, 0.0, 0.0 }, 4.0 + 1.0e-9, 9.0 - 1.0e-9, 8.020504272969694, 6.6837535608080785, 0, "-2 sum (2i-5)^2/(x-i^2)^3, n=2" },
  { &aps2_test, { 3.0, 0.0, 0.0 }, 9.0 + 1.0e-9, 16.0 - 1.0e-9, 13.486441986002655, 11.238701655002211, 0, "-2 sum (2i-5)^2/(x-i^2)^3, n=3" },
  { &aps2_test, { 4.0, 0.0, 0.0 }, 16.0 + 1.0e-9, 25.0 - 1.0e-9, 23.611200096748092, 19.67600008062341, 0, "-2 sum (2i-5)^2/(x-i^2)^3, n=4" },
  { &aps2_test, { 5.0, 0.0, 0.0 }, 25.0 + 1.0e-9, 36.0 - 1.0e-9, 35.793872791805704, 29.828227326504756, 0, "-2 sum (2i-5)^2/(x-i^2)^3, n=5" },
  { &aps2_test, { 6.0, 0.0, 0.0 }, 36.0 + 1.0e-9, 49.0 - 1.0e-9, 45.453058097144705, 41.906116195289414, 0, "-2 sum (2i-5)^2/(x-i^2)^3, n=6" },
  { &aps2_test, { 7.0, 0.0, 0.0 }, 49.0 + 1.0e-9, 64.0 - 1.0e-9, 59.97679789957155, 55.95359580014309, 0, "-2 sum (2i-5)^2/(x-i^2)^3, n=7" },
  { &aps2_test, { 8.0, 0.0, 0.0 }, 64.0 + 1.0e-9, 81.0 - 1.0e-9, 76.49283279279389, 71.9856655865878, 0, "-2 sum (2i-5)^2/(x-i^2)^3, n=8" },
  { &aps2_test, { 9.0, 0.0, 0.0 }, 81.0 + 1.0e-9, 100.0 - 1.0e-9, 95.00443426908333, 90.00886853916667, 0, "-2 sum (2i-5)^2/(x-i^2)^3, n=9" },
  { &aps2_test, { 10.0, 0.0, 0.0 }, 100.0 + 1.0e-9, 121.0 - 1.0e-9, 115.51326637366509, 110.0265327483302, 0, "-2 sum (2i-5)^2/(x-i^2)^3, n=10" },
  { &aps3_test, { 0.0, -40.0, -1.0 }, -9.0, 31.0, 0.2, 0.0, 0, "a x exp(b x), a=-40 b=-1" },
  { &aps3_test, { 0.0, -100.0, -2.0 }, -9.0, 31.0, 0.2, 0.0, 0, "a x exp(b x), a=-100 b=-2" },
  { &aps3_test, { 0.0, -200.0, -3.0 }, -9.0, 31.0, 0.2, 0.0, 0, "a x exp(b x), a=-200 b=-3" },
  { &aps4_test, { 4.0, 0.2, 0.0 }, 0.0, 5.0, 0.8024883659717065, 0.668740304976422, 0, "x^n - a, n=4 a=0.2" },
  { &aps4_test, { 6.0, 0.2, 0.0 }, 0.0, 5.0, 0.917669389598076, 0.76472449133173, 0, "x^n - a, n=6 a=0.2" },
  { &aps4_test, { 8.0, 0.2, 0.0 }, 0.0, 5.0, 0.981318520749531, 0.8177654339579425, 0, "x^n - a, n=8 a=0.2" },
  { &aps4_test, { 10.0, 0.2, 0.0 }, 0.0, 5.0, 1.0216079070249415, 0.8513399225207846, 0, "x^n - a, n=10 a=0.2" },
  { &aps4_test, { 12.0, 0.2, 0.0 }, 0.0, 5.0, 1.0493823266654014, 0.8744852722211678, 0, "x^n - a, n=12 a=0.2" },
  { &aps4_test, { 4.0, 1.0, 0.0 }, 0.0, 5.0, 1.2, 1.0, 0, "x^n - a, n=4 a=1" },
  { &aps4_test, { 6.0, 1.0, 0.0 }, 0.0, 5.0, 1.2, 1.0, 0, "x^n - a, n=6 a=1" },
  { &aps4_test, { 8.0, 1.0, 0.0 }, 0.0, 5.0, 1.2, 1.0, 0, "x^n - a, n=8 a=1" },
  { &aps4_test, { 10.0, 1.0, 0.0 }, 0.0, 5.0, 1.2, 1.0, 0, "x^n - a, n=10 a=1" },
  { &aps4_test, { 12.0, 1.0, 0.0 }, 0.0, 5.0, 1.2, 1.0, 0, "x^n - a, n=12 a=1" },
  { &aps4_test, { 8.0, 1.0, 0.0 }, -0.95, 4.05, 1.2, 1.0, 0, "x^n - a, n=8 a=1" },
  { &aps4_test, { 10.0, 1.0, 0.0 }, -0.95, 4.05, 1.2, 1.0, 0, "x^n - a, n=10 a=1" },
  { &aps4_test, { 12.0, 1.0, 0.0 }, -0.95, 4.05, 1.2, 1.0, 0, "x^n - a, n=12 a=1" },
  { &aps4_test, { 14.0, 1.0, 0.0 }, -0.95, 4.05, 1.2, 1.0, 0, "x^n - a, n=14 a=1" },
  { &aps5_test, { 0.0, 0.0, 0.0 }, 0.0, 1.5, 0.6283185307179586, 0.5235987755982989, 0, "sin(x) - 0.5" },
  { &aps6_test, { 1.0, 0.0, 0.0 }, 0.0, 1.0, 0.506973251569484, 0.42247770964123665, 0, "2x exp(-n) - 2exp(-nx) + 1, n=1" },
  { &aps6_test, { 2.0, 0.0, 0.0 }, 0.0, 1.0, 0.36803929257984447, 0.3066994104832037, 0, "2x exp(-n) - 2exp(-nx) + 1, n=2" },
  { &aps6_test, { 3.0, 0.0, 0.0 }, 0.0, 1.0, 0.2684465491855956, 0.22370545765466296, 0, "2x exp(-n) - 2exp(-nx) + 1, n=3" },
  { &aps6_test, { 4.0, 0.0, 0.0 }, 0.0, 1.0, 0.20606297702341006, 0.1717191475195084, 0, "2x exp(-n) - 2exp(-nx) + 1, n=4" },
  { &aps6_test, { 5.0, 0.0, 0.0 }, 0.0, 1.0, 0.1659085860681889, 0.13825715505682407, 0, "2x exp(-n) - 2exp(-nx) + 1, n=5" },
  { &aps6_test, { 20.0, 0.0, 0.0 }, 0.0, 1.0, 0.04158883082502462, 0.03465735902085385, 0, "2x exp(-n) - 2exp(-nx) + 1, n=20" },
  { &aps6_test, { 40.0, 0.0, 0.0 }, 0.0, 1.0, 0.02079441541679836, 0.01732867951399863, 0, "2x exp(-n) - 2exp(-nx) + 1, n=40" },
  { &aps6_test, { 60.0, 0.0, 0.0 }, 0.0, 1.0, 0.013862943611198907, 0.011552453009332421, 0, "2x exp(-n) - 2exp(-nx) + 1, n=60" },
  { &aps6_test, { 80.0, 0.0, 0.0 }, 0.0, 1.0, 0.01039720770839918, 0.008664339756999316, 0, "2x exp(-n) - 2exp(-nx) + 1, n=80" },
  { &aps6_test, { 100.0, 0.0, 0.0 }, 0.0, 1.0, 0.008317766166719344, 0.006931471805599453, 0, "2x exp(-n) - 2exp(-nx) + 1, n=100" },
  { &aps7_test, { 5.0, 0.0, 0.0 }, 0.0, 1.0, 0.04608306220874628, 0.0384025518406219, 0, "(1 + (1-n)^2)x - (1-nx)^2, n=5" },
  { &aps7_test, { 10.0, 0.0, 0.0 }, 0.0, 1.0, 0.0118800119976006, 0.0099000099980005, 0, "(1 + (1-n)^2)x - (1-nx)^2, n=10" },
  { &aps7_test, { 20.0, 0.0, 0.0 }, 0.0, 1.0, 0.002992500046874414, 0.0024937500390620117, 0, "(1 + (1-n)^2)x - (1-nx)^2, n=20" },
  { &aps8_test, { 2.0, 0.0, 0.0 }, 0.0, 1.0, 0.6, 0.5, 0, "x^2 - (1-x)^n, n=2" },
  { &aps8_test, { 5.0, 0.0, 0.0 }, 0.0, 1.0, 0.4151457790178904, 0.345954815848242, 0, "x^2 - (1-x)^n, n=5" },
  { &aps8_test, { 10.0, 0.0, 0.0 }, 0.0, 1.0, 0.2941468005039687, 0.24512233375330725, 0, "x^2 - (1-x)^n, n=10" },
  { &aps8_test, { 15.0, 0.0, 0.0 }, 0.0, 1.0, 0.23465714824387873, 0.1955476235365656, 0, "x^2 - (1-x)^n, n=15" },
  { &aps8_test, { 20.0, 0.0, 0.0 }, 0.0, 1.0, 0.19790514873172915, 0.16492095727644096, 0, "x^2 - (1-x)^n, n=20" },
  { &aps9_test, { 1.0, 0.0, 0.0 }, 0.0, 1.0, 0.33060964919938124, 0.2755080409994844, 0, "(1 + (1-n)^4)x - (1-nx)^4, n=1" },
  { &aps9_test, { 2.0, 0.0, 0.0 }, 0.0, 1.0, 0.16530482459969062, 0.1377540204997422, 0, "(1 + (1-n)^4)x - (1-nx)^4, n=2" },
  { &aps9_test, { 4.0, 0.0, 0.0 }, 0.0, 1.0, 0.012366340533787733, 0.010305283778156444, 0, "(1 + (1-n)^4)x - (1-nx)^4, n=4" },
  { &aps9_test, { 5.0, 0.0, 0.0 }, 0.0, 1.0, 0.004340529814684876, 0.0036171081789040634, 0, "(1 + (1-n)^4)x - (1-nx)^4, n=5" },
  { &aps9_test, { 8.0, 0.0, 0.0 }, 0.0, 1.0, 0.0004930475021956745, 0.00041087291849639543, 0, "(1 + (1-n)^4)x - (1-nx)^4, n=8" },
  { &aps9_test, { 15.0, 0.0, 0.0 }, 0.0, 1.0, 3.1187491071489155e-05, 2.5989575892907626e-05, 0, "(1 + (1-n)^4)x - (1-nx)^4, n=15" },
  { &aps9_test, { 20.0, 0.0, 0.0 }, 0.0, 1.0, 9.202314146622404e-06, 7.668595122185337e-06, 0, "(1 + (1-n)^4)x - (1-nx)^4, n=20" },
  { &aps10_test, { 1.0, 0.0, 0.0 }, 0.0, 1.0, 0.48126976504985647, 0.401058137541547, 0, "exp(-nx)(x-1) + x^n, n=1" },
  { &aps10_test, { 5.0, 0.0, 0.0 }, 0.0, 1.0, 0.6193842225095203, 0.5161535187579336, 0, "exp(-nx)(x-1) + x^n, n=5" },
  { &aps10_test, { 10.0, 0.0, 0.0 }, 0.0, 1.0, 0.647426672290099, 0.5395222269084159, 0, "exp(-nx)(x-1) + x^n, n=10" },
  { &aps10_test, { 15.0, 0.0, 0.0 }, 0.0, 1.0, 0.6578187532087864, 0.5481822943406552, 0, "exp(-nx)(x-1) + x^n, n=15" },
  { &aps10_test, { 20.0, 0.0, 0.0 }, 0.0, 1.0, 0.6632456000141853, 0.5527046666784878, 0, "exp(-nx)(x-1) + x^n, n=20" },
  { &aps11_test, { 2.0, 0.0, 0.0 }, 0.01, 1.0, 0.6, 0.5, 0, "(nx-1)/((n-1)x), n=2" },
  { &aps11_test, { 5.0, 0.0, 0.0 }, 0.01, 1.0, 0.24, 0.2, 0, "(nx-1)/((n-1)x), n=5" },
  { &aps11_test, { 15.0, 0.0, 0.0 }, 0.01, 1.0, 0.08, 0.06666666666666667, 0, "(nx-1)/((n-1)x), n=15" },
  { &aps11_test, { 20.0, 0.0, 0.0 }, 0.01, 1.0, 0.06, 0.05, 0, "(nx-1)/((n-1)x), n=20" },
  { &aps12_test, { 2.0, 0.0, 0.0 }, 1.0, 100.0, 2.4, 2.0, 0, "x^(1/n) - n^(1/n), n=2" },
  { &aps12_test, { 3.0, 0.0, 0.0 }, 1.0, 100.0, 3.6, 3.0, 0, "x^(1/n) - n^(1/n), n=3" },
  { &aps12_test, { 4.0, 0.0, 0.0 }, 1.0, 100.0, 4.8, 4.0, 0, "x^(1/n) - n^(1/n), n=4" },
  { &aps12_test, { 5.0, 0.0, 0.0 }, 1.0, 100.0, 6.0, 5.0, 0, "x^(1/n) - n^(1/n), n=5" },
  { &aps12_test, { 6.0, 0.0, 0.0 }, 1.0, 100.0, 7.2, 6.0, 0, "x^(1/n) - n^(1/n), n=6" },
  { &aps12_test, { 7.0, 0.0, 0.0 }, 1.0, 100.0, 8.4, 7.0, 0, "x^(1/n) - n^(1/n), n=7" },
  { &aps12_test, { 9.0, 0.0, 0.0 }, 1.0, 100.0, 10.8, 9.0, 0, "x^(1/n) - n^(1/n), n=9" },
  { &aps12_test, { 11.0, 0.0, 0.0 }, 1.0, 100.0, 13.2, 11.0, 0, "x^(1/n) - n^(1/n), n=11" },
  { &aps12_test, { 13.0, 0.0, 0.0 }, 1.0, 100.0, 15.6, 13.0, 0, "x^(1/n) - n^(1/n), n=13" },
  { &aps12_test, { 15.0, 0.0, 0.0 }, 1.0, 100.0, 18.0, 15.0, 0, "x^(1/n) - n^(1/n), n=15" },
  { &aps12_test, { 17.0, 0.0, 0.0 }, 1.0, 100.0, 20.4, 17.0, 0, "x^(1/n) - n^(1/n), n=17" },
  { &aps12_test, { 19.0, 0.0, 0.0 }, 1.0, 100.0, 22.8, 19.0, 0, "x^(1/n) - n^(1/n), n=19" },
  { &aps12_test, { 21.0, 0.0, 0.0 }, 1.0, 100.0, 25.2, 21.0, 0, "x^(1/n) - n^(1/n), n=21" },
  { &aps12_test, { 23.0, 0.0, 0.0 }, 1.0, 100.0, 27.6, 23.0, 0, "x^(1/n) - n^(1/n), n=23" },
  { &aps12_test, { 25.0, 0.0, 0.0 }, 1.0, 100.0, 30.0, 25.0, 0, "x^(1/n) - n^(1/n), n=25" },
  { &aps12_test, { 27.0, 0.0, 0.0 }, 1.0, 100.0, 32.4, 27.0, 0, "x^(1/n) - n^(1/n), n=27" },
  { &aps12_test, { 29.0, 0.0, 0.0 }, 1.0, 100.0, 34.8, 29.0, 0, "x^(1/n) - n^(1/n), n=29" },
  { &aps12_test, { 31.0, 0.0, 0.0 }, 1.0, 100.0, 37.2, 31.0, 0, "x^(1/n) - n^(1/n), n=31" },
  { &aps12_test, { 33.0, 0.0, 0.0 }, 1.0, 100.0, 39.6, 33.0, 0, "x^(1/n) - n^(1/n), n=33" },
  { &aps13_test, { 0.0, 0.0, 0.0 }, -1.0, 4.0, 0.0, 0.0, ROOT_TEST_NO_X0 | ROOT_TEST_MULTIPLE, "x exp(-1/x^2)" },
  { &aps14_test, { 1.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=1" },
  { &aps14_test, { 2.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=2" },
  { &aps14_test, { 3.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=3" },
  { &aps14_test, { 4.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=4" },
  { &aps14_test, { 5.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=5" },
  { &aps14_test, { 6.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=6" },
  { &aps14_test, { 7.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=7" },
  { &aps14_test, { 8.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=8" },
  { &aps14_test, { 9.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=9" },
  { &aps14_test, { 10.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=10" },
  { &aps14_test, { 11.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=11" },
  { &aps14_test, { 12.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=12" },
  { &aps14_test, { 13.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=13" },
  { &aps14_test, { 14.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=14" },
  { &aps14_test, { 15.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=15" },
  { &aps14_test, { 16.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=16" },
  { &aps14_test, { 17.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=17" },
  { &aps14_test, { 18.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=18" },
  { &aps14_test, { 19.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=19" },
  { &aps14_test, { 20.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=20" },
  { &aps14_test, { 21.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=21" },
  { &aps14_test, { 22.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=22" },
  { &aps14_test, { 23.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=23" },
  { &aps14_test, { 24.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=24" },
  { &aps14_test, { 25.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=25" },
  { &aps14_test, { 26.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=26" },
  { &aps14_test, { 27.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=27" },
  { &aps14_test, { 28.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=28" },
  { &aps14_test, { 29.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=29" },
  { &aps14_test, { 30.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=30" },
  { &aps14_test, { 31.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=31" },
  { &aps14_test, { 32.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=32" },
  { &aps14_test, { 33.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=33" },
  { &aps14_test, { 34.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=34" },
  { &aps14_test, { 35.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=35" },
  { &aps14_test, { 36.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=36" },
  { &aps14_test, { 37.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=37" },
  { &aps14_test, { 38.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=38" },
  { &aps14_test, { 39.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=39" },
  { &aps14_test, { 40.0, 0.0, 0.0 }, -1.0e4, M_PI_2, 0.7485678227539347, 0.6238065189616123, 0, "piecewise (n/20)(x/1.5 + sin(x) - 1), n=40" },
  { &aps15_test, { 20.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 7.086156671306365e-05, 5.905130559421971e-05, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=20" },
  { &aps15_test, { 21.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 6.76405864079244e-05, 5.6367155339937e-05, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=21" },
  { &aps15_test, { 22.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 6.469969134671029e-05, 5.391640945559191e-05, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=22" },
  { &aps15_test, { 23.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 6.200387087393069e-05, 5.166989239494225e-05, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=23" },
  { &aps15_test, { 24.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 5.952371603897347e-05, 4.9603096699144554e-05, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=24" },
  { &aps15_test, { 25.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 5.72343423451668e-05, 4.7695285287638995e-05, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=25" },
  { &aps15_test, { 26.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 5.51145518879384e-05, 4.5928793239948666e-05, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=26" },
  { &aps15_test, { 27.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 5.314617503479774e-05, 4.4288479195664784e-05, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=27" },
  { &aps15_test, { 28.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 5.131354830945989e-05, 4.276129025788324e-05, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=28" },
  { &aps15_test, { 29.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 4.9603096699144554e-05, 4.1335913915953796e-05, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=29" },
  { &aps15_test, { 30.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 4.8002996805623766e-05, 4.00024973380198e-05, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=30" },
  { &aps15_test, { 31.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 4.650290315544802e-05, 3.875241929620669e-05, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=31" },
  { &aps15_test, { 32.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 4.5093724271949595e-05, 3.7578103559958e-05, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=32" },
  { &aps15_test, { 33.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 4.376743826395108e-05, 3.6472865219959236e-05, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=33" },
  { &aps15_test, { 34.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 4.2516940027838194e-05, 3.543078335653183e-05, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=34" },
  { &aps15_test, { 35.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 4.1335913915953796e-05, 3.44465949299615e-05, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=35" },
  { &aps15_test, { 36.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 4.021872705336045e-05, 3.351560587780038e-05, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=36" },
  { &aps15_test, { 37.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 3.9160339499324654e-05, 3.2633616249437206e-05, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=37" },
  { &aps15_test, { 38.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 3.81562282301112e-05, 3.1796856858426e-05, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=38" },
  { &aps15_test, { 39.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 3.720232252435842e-05, 3.1001935436965345e-05, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=39" },
  { &aps15_test, { 40.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 3.6294948804252113e-05, 3.0245790670210097e-05, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=40" },
  { &aps15_test, { 100.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 1.4733593078953828e-05, 1.2277994232461525e-05, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=100" },
  { &aps15_test, { 200.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 7.403447268529038e-06, 6.169539390440865e-06, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=200" },
  { &aps15_test, { 300.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 4.943830235795138e-06, 4.1198585298292825e-06, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=300" },
  { &aps15_test, { 400.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 3.71095486527266e-06, 3.092462387727217e-06, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=400" },
  { &aps15_test, { 500.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 2.9702453113260214e-06, 2.475204426105018e-06, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=500" },
  { &aps15_test, { 600.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 2.4760281214215255e-06, 2.063356767851271e-06, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=600" },
  { &aps15_test, { 700.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 2.122814409378512e-06, 1.7690120078154265e-06, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=700" },
  { &aps15_test, { 800.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 1.857793883863092e-06, 1.54816156988591e-06, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=800" },
  { &aps15_test, { 900.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 1.6516014439226823e-06, 1.3763345366022351e-06, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=900" },
  { &aps15_test, { 1000.0, 0.0, 0.0 }, -1.0e4, 1.0e-4, 1.486606294679657e-06, 1.2388385788997142e-06, 0, "piecewise exp((n+1)x/2 1000) - 1.859, n=1000" },
  { NULL, { 0.0, 0.0, 0.0 }, 0.0, 0.0, 0.0, 0.0, 0, NULL }
};
