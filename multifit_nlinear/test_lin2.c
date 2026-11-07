#define lin2_N         20  /* can be anything >= p */
#define lin2_P         5

static double lin2_x0[lin2_P] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
static double lin2_epsrel = 1.0e-9;

static void
lin2_checksol(const double x[], const double sumsq,
              const double epsrel, const char *sname,
              const char *pname)
{
  size_t i;
  const double n = (double) lin2_N;
  const double sumsq_exact = 0.5 * (n*(n - 1.0)) / (2.0*n + 1.0);
  const double sum_exact = 3.0 / (2.0*n + 1.0);
  const double c2 = n * (n + 1.0) * (2.0*n + 1.0) / 6.0;  /* sum_i c_i^2 */
  double sum = 0.0, sum_tol;

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < lin2_P; ++i)
    sum += (i + 1.0) * x[i];

  /*
   * sumsq = sumsq* + 1/2 sum_i c_i^2 (s - s*)^2, with c_i the multiplier
   * of the coefficient sum s in residual i (f_i = c_i s - 1), depends on
   * s only to second order: a solver stopping with sumsq within 64 eps of
   * the minimum (by the step size test, on this rank-deficient problem)
   * determines s only to sqrt(2 * 64 eps sumsq* / sum_i c_i^2), which
   * depends on rounding (e.g. FMA contraction) and can exceed epsrel
   */
  sum_tol = GSL_MAX (epsrel, sqrt (128.0 * GSL_DBL_EPSILON * sumsq_exact / c2) / sum_exact);

  gsl_test_rel(sum, sum_exact, sum_tol, "%s/%s coeff sum",
               sname, pname);
}

static int
lin2_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  size_t i, j;

  for (i = 0; i < lin2_N; ++i)
    {
      double fi = 0.0;

      for (j = 0; j < lin2_P; ++j)
        {
          double xj = gsl_vector_get(x, j);
          fi += (j + 1) * xj;
        }

      fi = (i + 1) * fi - 1.0;
      gsl_vector_set(f, i, fi);
    }

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
lin2_df (const gsl_vector * x, void *params, gsl_matrix * J)
{
  size_t i, j;

  for (i = 0; i < lin2_N; ++i)
    {
      for (j = 0; j < lin2_P; ++j)
        {
          gsl_matrix_set(J, i, j, (i + 1.0) * (j + 1.0));
        }
    }

  (void)x;      /* avoid unused parameter warning */
  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
lin2_fvv (const gsl_vector * x, const gsl_vector * v,
          void *params, gsl_vector * fvv)
{
  (void)x; /* avoid unused parameter warnings */
  (void)v;
  (void)params;

  gsl_vector_set_zero(fvv);

  return GSL_SUCCESS;
}

static gsl_multifit_nlinear_fdf lin2_func =
{
  lin2_f,
  lin2_df,
  lin2_fvv,
  lin2_N,
  lin2_P,
  NULL,
  0,
  0,
  0
};

static test_fdf_problem lin2_problem =
{
  "linear_rank1",
  lin2_x0,
  NULL,
  NULL,
  &lin2_epsrel,
  &lin2_checksol,
  &lin2_func
};
