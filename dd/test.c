/* dd/test.c
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
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_dd.h>
#include <gsl/gsl_ieee_utils.h>

#include "internal_common.c"

/* Always compile the portable static implementation into this test TU. */
#include "internal_p.c"

/* When implemented, make the FMA static implementation available too. */
#ifdef HAVE_DECL_FMA
#include "internal_fma.c"
#endif

typedef struct { gsl_dd a, b, ref; } dd_binary_case;
typedef struct { gsl_dd a; double b; gsl_dd ref; } dd_mixed_case;
typedef struct { gsl_dd a, ref; } dd_unary_case;

static const dd_binary_case add_cases[] = {
  { {1, 5.5511151231257827e-17}, {3.1415926535897931, 1.2246467991473532e-16}, {4.1415926535897931, 1.7797583114599315e-16} },
  { {3.1415926535897931, 1.2246467991473532e-16}, {-7, 2e-16}, {-3.8584073464102064, -1.216245299353273e-16} },
  { {1e+20, 4096}, {9.9999999999999995e-21, 9.9999999999999994e-37}, {1e+20, 4096} },
  { {9.9999999999999995e-21, 9.9999999999999994e-37}, {0.5, -2.7755575615628914e-17}, {0.5, -2.7745575615628914e-17} },
  { {-7, 2e-16}, {1, 5.5511151231257827e-17}, {-6, 2.5551115123125782e-16} },
};

static const dd_mixed_case add_d_cases[] = {
  { {1, 5.5511151231257827e-17}, 0.75, {1.75, 5.5511151231257827e-17} },
  { {3.1415926535897931, 1.2246467991473532e-16}, -2.5, {0.64159265358979323, 1.1442377452219667e-17} },
  { {1e+20, 4096}, 1e-10, {1e+20, 4096.0000000001} },
  { {9.9999999999999995e-21, 9.9999999999999994e-37}, 10000000000, {10000000000, 1.0000000000000001e-20} },
  { {-7, 2e-16}, 3, {-4, 2e-16} },
};

static const dd_binary_case sub_cases[] = {
  { {1, 5.5511151231257827e-17}, {3.1415926535897931, 1.2246467991473532e-16}, {-2.1415926535897931, -6.6953528683477494e-17} },
  { {3.1415926535897931, 1.2246467991473532e-16}, {-7, 2e-16}, {10.141592653589793, -7.7535320085264675e-17} },
  { {1e+20, 4096}, {9.9999999999999995e-21, 9.9999999999999994e-37}, {1e+20, 4096} },
  { {9.9999999999999995e-21, 9.9999999999999994e-37}, {0.5, -2.7755575615628914e-17}, {-0.49999999999999994, -2.7745575615628914e-17} },
  { {-7, 2e-16}, {1, 5.5511151231257827e-17}, {-8, 1.4448884876874217e-16} },
};

static const dd_mixed_case sub_d_cases[] = {
  { {1, 5.5511151231257827e-17}, 0.75, {0.25000000000000006, 0} },
  { {3.1415926535897931, 1.2246467991473532e-16}, -2.5, {5.6415926535897931, 1.2246467991473532e-16} },
  { {1e+20, 4096}, 1e-10, {1e+20, 4095.9999999999} },
  { {9.9999999999999995e-21, 9.9999999999999994e-37}, 10000000000, {-10000000000, 1.0000000000000001e-20} },
  { {-7, 2e-16}, 3, {-10, 2e-16} },
};


static const dd_binary_case mul_cases[] = {
  { {1, 5.5511151231257827e-17}, {3.1415926535897931, 1.2246467991473532e-16}, {3.1415926535897936, -1.4723110503489571e-16} },
  { {3.1415926535897931, 1.2246467991473532e-16}, {-7, 2e-16}, {-21.991148575128552, -2.289342286851886e-16} },
  { {1e+20, 4096}, {9.9999999999999995e-21, 9.9999999999999994e-37}, {1, 8.6113271454209573e-17} },
  { {9.9999999999999995e-21, 9.9999999999999994e-37}, {0.5, -2.7755575615628914e-17}, {4.9999999999999997e-21, 2.2244424384371081e-37} },
  { {-7, 2e-16}, {1, 5.5511151231257827e-17}, {-7, -1.8857805861880479e-16} },
};

static const dd_mixed_case mul_d_cases[] = {
  { {1, 5.5511151231257827e-17}, 0.75, {0.75, 4.163336342344337e-17} },
  { {3.1415926535897931, 1.2246467991473532e-16}, -2.5, {-7.8539816339744828, -3.0616169978683831e-16} },
  { {1e+20, 4096}, 1e-10, {10000000000, 7.7392197315497748e-07} },
  { {9.9999999999999995e-21, 9.9999999999999994e-37}, 10000000000, {1e-10, 8.721074138711825e-28} },
  { {-7, 2e-16}, 3, {-21, 5.9999999999999999e-16} },
};

static const dd_binary_case div_cases[] = {
  { {1, 5.5511151231257827e-17}, {3.1415926535897931, 1.2246467991473532e-16}, {0.31830988618379069, -2.0089284448296168e-18} },
  { {3.1415926535897931, 1.2246467991473532e-16}, {-7, 2e-16}, {-0.44879895051282764, 2.519336980021485e-17} },
  { {1e+20, 4096}, {9.9999999999999995e-21, 9.9999999999999994e-37}, {1e+40, -3.457187429690993e+23} },
  { {9.9999999999999995e-21, 9.9999999999999994e-37}, {0.5, -2.7755575615628914e-17}, {2.0000000000000002e-20, 1.009574865201005e-37} },
  { {-7, 2e-16}, {1, 5.5511151231257827e-17}, {-6.9999999999999991, -2.996003610813205e-16} },
};

static const dd_mixed_case div_d_cases[] = {
  { {1, 5.5511151231257827e-17}, 0.75, {1.3333333333333335, -7.401486830834377e-17} },
  { {3.1415926535897931, 1.2246467991473532e-16}, -2.5, {-1.2566370614359172, -4.8985871965894131e-17} },
  { {1e+20, 4096}, 1e-10, {1e+30, -15356822154153.742} },
  { {9.9999999999999995e-21, 9.9999999999999994e-37}, 10000000000, {1e-30, -3.818314915337642e-47} },
  { {-7, 2e-16}, 3, {-2.3333333333333335, 2.146964032833542e-16} },
};


static const dd_unary_case exp_cases[] = {
  { { 0, 0 }, { 1, 0 } },                                                                                      /* zero */
  { { 1, 0 }, { 2.7182818284590451, 1.4456468917292502e-16 } },                                                /* one */
  { { -1, 0 }, { 0.36787944117144233, -1.2428753672788363e-17 } },                                             /* minus_one */
  { { 0.5, 0 }, { 1.6487212707001282, -4.7315684794358332e-17 } },                                             /* half */
  { { -0.5, 0 }, { 0.60653065971263342, -6.5931784154914137e-19 } },                                           /* minus_half */
  { { 2, 0 }, { 7.3890560989306504, -1.7971139497839148e-16 } },                                               /* two */
  { { -2, 0 }, { 0.1353352832366127, -1.042381423288669e-17 } },                                               /* minus_two */
  { { 10, 0 }, { 22026.465794806718, -1.3780134700517372e-12 } },                                              /* ten */
  { { -10, 0 }, { 4.5399929762484854e-05, -2.6375540553275309e-21 } },                                         /* minus_ten */
  { { 100, 0 }, { 2.6881171418161356e+43, -1.6101271449201627e+27 } },                                         /* hundred */
  { { -100, 0 }, { 3.7200759760208361e-44, -1.5705024907732008e-60 } },                                        /* minus_hundred */
  { { 600, 0 }, { 3.7730203009299397e+260, 1.6116934109232247e+244 } },                                        /* six_hundred */
  { { -600, 0 }, { 2.6503965530043108e-261, 6.3773428174913948e-278 } },                                       /* minus_six_hundred */
  { { 700, 0 }, { 1.0142320547350045e+304, 1.6666571920734673e+287 } },                                        /* seven_hundred */
  { { -650, 0 }, { 5.111951948651156e-283, 2.849121073604364e-299 } },                                         /* minus_650 */
  { { 3.1415926535897931, 1.2246467991473532e-16 }, { 23.14069263277927, -1.3488747091995788e-15 } },          /* pi */
  { { -3.1415926535897931, -1.2246467991473532e-16 }, { 0.043213918263772251, -1.2638616556157167e-18 } },     /* minus_pi */
  { { 0.69314718055994529, 2.3190468138462996e-17 }, { 2, -1.1415416876832424e-33 } },                         /* ln2 */
  { { -0.69314718055994529, -2.3190468138462996e-17 }, { 0.5, 2.853854219208106e-34 } },                       /* minus_ln2 */
  { { 1.4142135623730951, -9.6672933134529135e-17 }, { 4.1132503787829275, -2.3836812819073011e-17 } },        /* sqrt2 */
  { { 7.8886090522101181e-31, 4.3790577010150533e-47 }, { 1, 7.8886090522101181e-31 } },                       /* tiny_pos */
  { { -7.8886090522101181e-31, 4.3790577010150533e-47 }, { 1, -7.8886090522101172e-31 } },                     /* tiny_neg */
  { { 1, 5.5511151231257827e-17 }, { 2.7182818284590455, -1.486295670083675e-16 } },                           /* one_low_plus */
  { { 1, -5.5511151231257827e-17 }, { 2.7182818284590451, -6.3302644958450886e-18 } },                         /* one_low_minus */
  { { 0.69314718055994529, 2.3190468138462996e-17 }, { 2, -1.1415416876832424e-33 } },                         /* center_1 */
  { { 11.78350206951907, 2.8321565589135527e-16 }, { 131072, -8.6790850927499448e-28 } },                      /* center_17 */
  { { 69.314718055994533, -1.6777560748042639e-15 }, { 1.2676506002282294e+30, 0.068271199721139597 } },       /* center_100 */
  { { 346.57359027997268, -2.2599635089223323e-14 }, { 3.2733906078961419e+150, 1.5270289553055257e+120 } },   /* center_500 */
  { { 693.14718055994535, -4.5199270178446646e-14 }, { 1.0715086071862673e+301, 9.9971244805651305e+270 } },   /* center_1000 */
  { { -0.69314718055994529, -2.3190468138462996e-17 }, { 0.5, 2.853854219208106e-34 } },                       /* center_-1 */
  { { -11.78350206951907, -2.8321565589135527e-16 }, { 7.62939453125e-06, 5.0518924211791861e-38 } },          /* center_-17 */
  { { -69.314718055994533, 1.6777560748042639e-15 }, { 7.8886090522101181e-31, -4.2485271890256103e-62 } },    /* center_-100 */
  { { -346.57359027997268, 2.2599635089223323e-14 }, { 3.0549363634996047e-151, -1.425120568387625e-181 } },   /* center_-500 */
  { { -623.83246250395075, -2.7532759472367637e-14 }, { 1.1830521861667747e-271, -1.7400147651507838e-301 } }, /* center_-900 */
  { { 0.34657359027997264, 1.1595234069231498e-17 }, { 1.4142135623730951, -9.6672933134529135e-17 } },        /* halfbd_0 */
  { { 12.130075659799044, -6.4887868097079624e-16 }, { 185363.80004736633, -1.2671114691808995e-11 } },        /* halfbd_17 */
  { { 69.661291646274506, -2.6098504116664153e-15 }, { 1.7927286711931566e+30, -122547501713809.06 } },        /* halfbd_100 */
  { { 346.92016387025262, 4.8899800043185329e-15 }, { 4.6292733926314343e+150, -3.1644827136033766e+134 } },   /* halfbd_500 */
  { { 693.4937541502253, -1.770965508490479e-14 }, { 1.5153420044823246e+301, -1.035858799355892e+285 } },     /* halfbd_1000 */
  { { -0.34657359027997264, -1.1595234069231498e-17 }, { 0.70710678118654757, -4.8336466567264567e-17 } },     /* halfbd_-1 */
  { { -11.436928479239098, 5.6104684664674358e-16 }, { 1.0789593218788873e-05, -7.3755594737647388e-22 } },    /* halfbd_-17 */
  { { -68.968144465714559, 7.4566173794211241e-16 }, { 1.1156177909894717e-30, -7.6261497542874941e-47 } },    /* halfbd_-100 */
  { { -346.22701668969268, -6.7541686780428354e-15 }, { 4.3203324374478846e-151, -2.9532965879883789e-167 } }, /* halfbd_-500 */
  { { -623.48588891367081, -4.3144378825779911e-17 }, { 1.6730884466721926e-271, -1.1436912488795913e-287 } }, /* halfbd_-900 */
};

static const dd_unary_case log_cases[] = {
  { { 1.00000000000000000e+00,  0.00000000000000000e+00 }, { 0.00000000000000000e+00,  0.00000000000000000e+00 } },  /* one */
  { { 2.00000000000000000e+00,  0.00000000000000000e+00 }, { 6.93147180559945286e-01,  2.31904681384629956e-17 } },  /* two */
  { { 5.00000000000000000e-01,  0.00000000000000000e+00 }, {-6.93147180559945286e-01, -2.31904681384629956e-17 } },  /* half */
  { { 3.14159265358979312e+00,  1.22464679914735321e-16 }, { 1.14472988584940016e+00,  1.02659511627078279e-17 } },  /* pi */
  { { 1.41421356237309515e+00, -9.66729331345291345e-17 }, { 3.46573590279972643e-01,  1.15952340692314947e-17 } },  /* sqrt2 */
  { { 2.71828182845904509e+00,  1.44564689172925016e-16 }, { 1.00000000000000000e+00,  7.82743380676001721e-34 } },  /* e */
  { { 1.00000000000000000e+00,  5.55111512312578270e-17 }, { 5.55111512312578270e-17, -1.54074395550978868e-33 } },  /* near1_plus */
  { { 1.00000000000000000e+00, -5.55111512312578270e-17 }, {-5.55111512312578270e-17, -1.54074395550978868e-33 } },  /* near1_minus */
  { { 8.75000000000000000e-01,  0.00000000000000000e+00 }, {-1.33531392624522627e-01,  3.66445766366008474e-18 } },  /* seven_eighths */
  { { 1.12500000000000000e+00,  0.00000000000000000e+00 }, { 1.17783035656383456e-01, -1.19716857475936773e-18 } },  /* nine_eighths */
  { { 2.50000000000000000e-01,  0.00000000000000000e+00 }, {-1.38629436111989057e+00, -4.63809362769259912e-17 } },  /* quarter */
  { { 4.00000000000000000e+00,  0.00000000000000000e+00 }, { 1.38629436111989057e+00,  4.63809362769259912e-17 } },  /* four */
  { { 1.00000000000000000e+01,  0.00000000000000000e+00 }, { 2.30258509299404590e+00, -2.17075622338224935e-16 } },  /* ten */
  { { 3.05493636349960468e-151,  0.00000000000000000e+00 }, {-3.46573590279972677e+02,  2.25996350892233229e-14 } }, /* 2^-500 */
  { { 1.26790687472559659e+271,  0.00000000000000000e+00 }, { 6.24237927612058911e+02,  3.17487258279172694e-14 } }, /* 1.5*2^900 */

  /* Exact reduction centers. */
  { { 1.00390625000000000e+00,  0.00000000000000000e+00 }, { 3.89864041565732289e-03,  1.25416590383049728e-19 } },  /* center_0 */
  { { 1.01171875000000000e+00,  0.00000000000000000e+00 }, { 1.16506172199752744e-02, -2.36187885155090345e-19 } },  /* center_1 */
  { { 1.24609375000000000e+00,  0.00000000000000000e+00 }, { 2.20013658305282106e-01, -1.00795744224419992e-17 } },  /* center_31 */
  { { 1.49609375000000000e+00,  0.00000000000000000e+00 }, { 4.02857544701083536e-01, -2.17651298256827917e-17 } },  /* center_63 */
  { { 1.50390625000000000e+00,  0.00000000000000000e+00 }, { 4.08065889808221727e-01,  2.16879664336615209e-17 } },  /* center_64 */
  { { 1.69921875000000000e+00,  0.00000000000000000e+00 }, { 5.30168586609121584e-01,  3.40999644338310668e-17 } },  /* center_89 */
  { { 1.84765625000000000e+00,  0.00000000000000000e+00 }, { 6.13917944012370542e-01, -5.00424804679914907e-17 } },  /* center_108 */
  { { 1.99609375000000000e+00,  0.00000000000000000e+00 }, { 6.91192145724142004e-01, -4.52888001650340699e-17 } },  /* center_127 */

  /* Exact center-selection boundaries. */
  { { 1.00781250000000000e+00,  0.00000000000000000e+00 }, { 7.78214044205494896e-03, -1.28191791233438450e-20 } },  /* boundary_1 */
  { { 1.25000000000000000e+00,  0.00000000000000000e+00 }, { 2.23143551314209765e-01, -9.09127059732479905e-18 } },  /* boundary_32 */
  { { 1.50000000000000000e+00,  0.00000000000000000e+00 }, { 4.05465108108164385e-01, -2.88113802596264264e-18 } },  /* boundary_64 */
  { { 1.75000000000000000e+00,  0.00000000000000000e+00 }, { 5.59615787935422659e-01,  2.68549258021230796e-17 } },  /* boundary_96 */
  { { 1.99218750000000000e+00,  0.00000000000000000e+00 }, { 6.89233281238809004e-01, -2.32077938372292054e-17 } }   /* boundary_127 */
};

static const dd_unary_case log1p_cases[] = {
  { { 0.00000000000000000e+00,  0.00000000000000000e+00 }, { 0.00000000000000000e+00,  0.00000000000000000e+00 } },  /* zero */

  /* Tiny arguments: direct Taylor-36 path. */
  { { 7.88860905221011805e-31,  2.18952885050752667e-47 }, { 7.88860905221011805e-31,  2.18952885050749556e-47 } },  /* tiny_pos */
  { {-7.88860905221011805e-31,  2.18952885050752667e-47 }, {-7.88860905221011805e-31,  2.18952885050749556e-47 } },  /* tiny_neg */

  /* Direct-kernel boundaries |x| = 1/8. */
  { { 1.25000000000000000e-01,  0.00000000000000000e+00 }, { 1.17783035656383456e-01, -1.19716857475936773e-18 } },  /* one_eighth */
  { {-1.25000000000000000e-01,  0.00000000000000000e+00 }, {-1.33531392624522627e-01,  3.66445766366008474e-18 } },  /* minus_one_eighth */

  /* Just outside the direct-kernel interval. */
  { { 1.25000000000000000e-01,  6.93889390390722838e-18 }, { 1.17783035656383456e-01,  4.97073711760261309e-18 } },  /* above_eighth */
  { {-1.25000000000000000e-01, -6.93889390390722838e-18 }, {-1.33531392624522627e-01, -4.26570679794817659e-18 } },  /* below_minus_eighth */

  /* Representative Candidate-C core cases. */
  { { 5.00000000000000000e-01,  0.00000000000000000e+00 }, { 4.05465108108164385e-01, -2.88113802596264264e-18 } },  /* half */
  { {-5.00000000000000000e-01,  0.00000000000000000e+00 }, {-6.93147180559945286e-01, -2.31904681384629956e-17 } },  /* minus_half */
  { { 7.50000000000000000e-01,  0.00000000000000000e+00 }, { 5.59615787935422659e-01,  2.68549258021230796e-17 } },  /* three_quarters */
  { {-7.50000000000000000e-01,  0.00000000000000000e+00 }, {-1.38629436111989057e+00, -4.63809362769259912e-17 } },  /* minus_three_quarters */

  /* DD low-limb sensitivity in the core path. */
  { { 5.00000000000000000e-01,  2.77555756156289135e-17 }, { 4.05465108108164385e-01,  1.56225790511233002e-17 } },  /* half_low_plus */
  { {-5.00000000000000000e-01,  2.77555756156289135e-17 }, {-6.93147180559945286e-01,  3.23206830927948284e-17 } },  /* minus_half_low_plus */

  /*
   * Near the singular point x = -1.  These specifically test DD inputs
   * whose distance from -1 resides entirely in the low component.
   */
  { {-1.00000000000000000e+00,  1.11022302462515654e-16 }, {-3.67368005696771007e+01, -6.73983299025960560e-16 } },  /* -1 + 2^-53 */
  { {-1.00000000000000000e+00,  5.55111512312578270e-17 }, {-3.74299477502370479e+01,  1.19020537469834266e-15 } },  /* -1 + 2^-54 */
  { { 1.00000000000000000e+00,  0.00000000000000000e+00 }, { 6.93147180559945286e-01,  2.31904681384629956e-17 } },  /* one */

  /*
   * Large-x transition.  x=7 uses the Candidate-C core; x=8 enters
   * log1p_large().
   */
  { { 7.00000000000000000e+00,  0.00000000000000000e+00 }, { 2.07944154167983575e+00,  1.80593706877904650e-16 } },  /* seven */
  { { 8.00000000000000000e+00,  0.00000000000000000e+00 }, { 2.19722457733621956e+00, -1.81425944700030599e-16 } },  /* eight */
  { { 1.00000000000000000e+01,  0.00000000000000000e+00 }, { 2.39789527279837067e+00, -1.25358421142316092e-16 } },  /* ten */
  { { 1.00000000000000000e+02,  0.00000000000000000e+00 }, { 4.61512051684125968e+00, -2.28213021947403639e-16 } },  /* hundred */

  /* Very large arguments exercise log(x) + K(1/x). */
  { { 1.00000000000000002e+100,  0.00000000000000000e+00 }, { 2.30258509299404579e+02, -1.10335183063112318e-14 } }, /* 1e100 */
  { { 1.79769313486231571e+308,  0.00000000000000000e+00 }, { 7.09782712893383973e+02,  2.36360170713235918e-14 } }  /* DBL_MAX */
};

#define NELEM(a) (sizeof(a) / sizeof((a)[0]))

/*
 * DD relative error test.
 *
 * The reference is the canonical (hi,lo) pair generated at high precision.
 * Form the error as a double-double subtraction so that the low component is
 * not discarded.  Scaling by |ref.hi| is appropriate for normalized nonzero
 * DD references.  Near zero, use an absolute DD error test.
 *
 * Relative tolerance for the initial DD unit tests.
 * The underlying binary64 unit roundoff is u = 2^-53.  Accurate
 * double-word arithmetic has errors of order u^2 = 2^-106.
 */
#define DD_U2      1.232595164407830945955825883254353483864e-32
#define DD_TOL     (16.0 * DD_U2)
#define DD_LOG_TOL (1.0 * DD_U2)

static double
dd_abs_value(const gsl_dd x)
{
  return fabs(x.hi) + fabs(x.lo);
}

static gsl_dd
dd_test_diff(const gsl_dd a, const gsl_dd b)
{
  gsl_dd s, t, u;
  double e;

  /*
   * This is test machinery only.  Form
   *
   *   (a.hi + a.lo) - (b.hi + b.lo)
   *
   * without collapsing either DD number to binary64.
   */
  s = two_sum(a.hi, -b.hi);
  t = two_sum(a.lo, -b.lo);

  e = s.lo + t.hi;
  u = two_sum(s.hi, e);
  e = u.lo + t.lo;

  return two_sum(u.hi, e);
}

static void
test_dd_rel(const gsl_dd got, const gsl_dd ref, const double tol,
            const char *name, const size_t i)
{
  gsl_dd d;
  double err;
  double scale;
  double rel;

  d = dd_test_diff(got, ref);
  err = dd_abs_value(d);
  scale = dd_abs_value(ref);

  if (scale == 0.0)
    {
      gsl_test(err > tol,
               "%s[%u] absolute DD error %.18g > %.18g",
               name, (unsigned int) i, err, tol);
    }
  else
    {
      rel = err / scale;

      gsl_test(rel > tol,
               "%s[%u] relative DD error %.18g > %.18g",
               name, (unsigned int) i, rel, tol);
    }
}

/*
 * Check the usual nonoverlap condition for a normalized DD pair.  For a
 * normal finite high word, |lo| should not exceed approximately half an ulp
 * of hi.  The tolerance factor accommodates exact halfway cases.
 */
static void
test_dd_normalized(const gsl_dd x, const char *name, const size_t i)
{
  double ah;
  int e;
  double half_ulp;
  int status;

  if (x.hi == 0.0)
    {
      status = (x.lo != 0.0);
    }
  else if (!gsl_finite(x.hi))
    {
      status = 0;
    }
  else
    {
      ah = fabs(x.hi);
      (void) frexp(ah, &e);
      half_ulp = ldexp(1.0, e - 54);
      status = fabs(x.lo) > half_ulp;
    }

  gsl_test(status, "%s[%u] normalized DD result",
           name, (unsigned int) i);
}

/*
 * Construct a function name from prefix + operation + suffix.
 *
 * Examples:
 *
 *   DD_FN(,     dd_add, _p)   -> dd_add_p
 *   DD_FN(,     dd_add, _fma) -> dd_add_fma
 *   DD_FN(gsl_, dd_add, )     -> gsl_dd_add
 */
#define DD_CAT3_I(a,b,c) a##b##c
#define DD_CAT3(a,b,c) DD_CAT3_I(a,b,c)
#define DD_FN(prefix,name,suffix) DD_CAT3(prefix,name,suffix)

#define TEST_BINARY(prefix,suffix,opname,cases)                         \
  do {                                                                  \
    size_t i_;                                                          \
    for (i_ = 0; i_ < NELEM(cases); ++i_)                              \
      {                                                                 \
        gsl_dd r_ = DD_FN(prefix,opname,suffix)                         \
                      ((cases)[i_].a, (cases)[i_].b);                   \
        test_dd_rel(r_, (cases)[i_].ref, DD_TOL,                       \
                    #prefix #opname #suffix, i_);                       \
        test_dd_normalized(r_, #prefix #opname #suffix, i_);           \
      }                                                                 \
  } while (0)

#define TEST_MIXED(prefix,suffix,opname,cases)                          \
  do {                                                                  \
    size_t i_;                                                          \
    for (i_ = 0; i_ < NELEM(cases); ++i_)                              \
      {                                                                 \
        gsl_dd r_ = DD_FN(prefix,opname,suffix)                         \
                      ((cases)[i_].a, (cases)[i_].b);                   \
        test_dd_rel(r_, (cases)[i_].ref, DD_TOL,                       \
                    #prefix #opname #suffix, i_);                       \
        test_dd_normalized(r_, #prefix #opname #suffix, i_);           \
      }                                                                 \
  } while (0)


#define TEST_UNARY(prefix,suffix,opname,cases,tol)                      \
  do {                                                                  \
    size_t i_;                                                          \
    for (i_ = 0; i_ < NELEM(cases); ++i_)                               \
      {                                                                 \
        gsl_dd r_ = DD_FN(prefix,opname,suffix)((cases)[i_].a);         \
        test_dd_rel(r_, (cases)[i_].ref, tol,                           \
                    #prefix #opname #suffix, i_);                       \
        test_dd_normalized(r_, #prefix #opname #suffix, i_);            \
      }                                                                 \
  } while (0)

#define TEST_EXP_SPECIAL(prefix,suffix)                                 \
  do {                                                                  \
    gsl_dd a_, r_;                                                      \
    a_.hi = -710.0;                                                     \
    a_.lo = 0.0;                                                        \
    r_ = DD_FN(prefix,dd_exp,suffix)(a_);                              \
    gsl_test(!gsl_finite(r_.hi),                                       \
             "%s nonfinite at -710", #prefix "dd_exp" #suffix);       \
    a_.hi = -746.0;                                                     \
    a_.lo = 0.0;                                                        \
    r_ = DD_FN(prefix,dd_exp,suffix)(a_);                              \
    gsl_test(r_.hi != 0.0 || r_.lo != 0.0,                             \
             "%s nonzero at -746", #prefix "dd_exp" #suffix);         \
    a_.hi = 709.0;                                                      \
    a_.lo = 0.0;                                                        \
    r_ = DD_FN(prefix,dd_exp,suffix)(a_);                              \
    gsl_test(!gsl_finite(r_.hi),                                       \
             "%s nonfinite at +709", #prefix "dd_exp" #suffix);       \
    a_.hi = 710.0;                                                      \
    a_.lo = 0.0;                                                        \
    r_ = DD_FN(prefix,dd_exp,suffix)(a_);                              \
    gsl_test(gsl_finite(r_.hi),                                        \
             "%s finite at +710", #prefix "dd_exp" #suffix);          \
  } while (0)

#define TEST_COMMON(prefix,suffix)                                      \
  do {                                                                  \
    TEST_BINARY(prefix,suffix,dd_add,add_cases);                        \
    TEST_MIXED(prefix,suffix,dd_add_d,add_d_cases);                     \
    TEST_BINARY(prefix,suffix,dd_sub,sub_cases);                        \
    TEST_MIXED(prefix,suffix,dd_sub_d,sub_d_cases);                     \
  } while (0)

#define TEST_BACKEND(prefix,suffix)                                     \
  do {                                                                  \
    TEST_BINARY(prefix,suffix,dd_mul,mul_cases);                        \
    TEST_MIXED(prefix,suffix,dd_mul_d,mul_d_cases);                     \
    TEST_BINARY(prefix,suffix,dd_div,div_cases);                        \
    TEST_MIXED(prefix,suffix,dd_div_d,div_d_cases);                     \
    TEST_UNARY(prefix,suffix,dd_exp,exp_cases,DD_TOL);                  \
    TEST_EXP_SPECIAL(prefix,suffix);                                    \
    TEST_UNARY(prefix,suffix,dd_log,log_cases,DD_LOG_TOL);              \
    TEST_UNARY(prefix,suffix,dd_log1p,log1p_cases,DD_LOG_TOL);          \
  } while (0)

#define TEST_API(prefix,suffix)                                         \
  do {                                                                  \
    TEST_COMMON(prefix,suffix);                                         \
    TEST_BACKEND(prefix,suffix);                                        \
  } while (0)

int
main(void)
{
  gsl_ieee_env_setup();

  /* Operations common to all implementations. */
  TEST_COMMON(, );

  /* Portable C89 backend is always tested. */
  TEST_BACKEND(, _p);

#ifdef HAVE_DECL_FMA
  /* Test the FMA backend independently when fma() is available. */
  TEST_BACKEND(, _fma);
#endif

  /* Test the public API against the same high-precision references. */
  TEST_API(gsl_, );

  exit(gsl_test_summary());
}
