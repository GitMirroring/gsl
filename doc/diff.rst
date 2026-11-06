.. index::
   single: differentiation of functions, numeric
   single: functions, numerical differentiation
   single: derivatives, calculating numerically
   single: numerical derivatives
   single: slope, see numerical derivative

*************************
Numerical Differentiation
*************************

The functions described in this chapter compute numerical derivatives by
finite differencing.  Both fixed-step and adaptive algorithms are
provided: the adaptive algorithm selects a finite-difference step size
and estimates the error in the derivative.
These functions are declared in the header file :file:`gsl_deriv.h`.

The basic idea is to approximate the :math:`m`-th derivative of
a function at a point :math:`x_0` as

.. math:: D_h^{(m)} f(x_0) \approx \frac{1}{h^m} \sum_{k=1}^N w_k f(x_0 + h s_k)

where :math:`h` is the step size, :math:`w_k` are the finite difference weights,
:math:`s_k` are dimensionless stencil offsets, and :math:`N` is the stencil size.
The stencil size, weights, and offsets can be chosen to produce a derivative estimate
which has a leading truncation error of :math:`O(h^p)`. Here, :math:`p` is known
as the *method order*, so that

.. math:: D_h^{(m)} f(x_0) = f^{(m)}(x_0) + C f^{(m+p)}(x_0) h^p + O(h^{p+1})

where the constant :math:`C` is

.. math:: C = \frac{1}{(m+p)!} \sum_{k=1}^N w_k s_k^{m+p}

For a given derivative order :math:`m`, method order :math:`p`,
and stencil type, GSL selects an appropriate stencil size and offsets.
The corresponding finite-difference weights :math:`w_k`
are then computed using the Fornberg algorithm.

By far, the most difficult problem in numerical differentiation is to determine
the optimal step size :math:`h`. For a sufficiently smooth function and sufficiently
small :math:`h`, the truncation error decreases as :math:`h` decreases.
However, roundoff error, caused by catastrophic cancellation
of the terms in the weighted sum, grows as :math:`h \rightarrow 0`. Therefore,
there is typically an intermediate range of useful step sizes, often visible
as a *plateau* in the derivative estimates or a minimum in the error curve.
If :math:`h` is too large, truncation error dominates, but if :math:`h` is too small,
roundoff error dominates. Ignoring problem-dependent and stencil-dependent proportionality
constants, the errors scale as

.. math::

   E_{\mathrm{trunc}}(h) \sim h^p \\
   E_{\mathrm{round}}(h) \sim \frac{\epsilon}{h^m}

where :math:`\epsilon` is the machine epsilon (:macro:`GSL_DBL_EPSILON`).
A useful theoretical estimate for the optimal step size, designed to balance
truncation and roundoff error is

.. math:: h_{\mathrm{opt}} \sim \epsilon^{\frac{1}{m+p}}

The omitted proportionality constant depends on the finite-difference stencil
and on derivatives and characteristic values of the function near :math:`x_0`.
For a typical black box function, these are not known in advance and so
:math:`\epsilon^{1/(m+p)}` is often used as an initial step size.

The adaptive algorithm used by GSL, which is based on the DERIVEST
algorithm, computes derivative estimates :math:`D_h^{(m)} f(x_0)`
for multiple step sizes :math:`h` spanning several decades around the value
:math:`h_{\mathrm{opt}}`. The finite difference derivative estimates are grouped into
successive overlapping windows and extrapolated using the expected truncation
error expansion. For each window, the derivative estimates are fitted by least
squares to a model consisting of the derivative value at :math:`h=0` together with the
leading truncation error terms. For central stencils, symmetry causes alternating
powers to vanish, so these terms have powers :math:`p,p+2,p+4,\ldots`, while for forward
and backward stencils they generally have powers :math:`p,p+1,p+2,\ldots`. The fitted
intercept at :math:`h=0` provides the Romberg-extrapolated derivative estimate for the
window. Since the least-squares solution is linear in the data, this estimate is a
fixed linear combination :math:`\sum_k g_k D_{h_k}^{(m)} f(x_0)` of the derivative
estimates in the window, where the extrapolation weights :math:`g_k` depend only on the
ratio between successive step sizes and on the truncation error powers.

Each Romberg window produces a candidate derivative estimate together with an error
estimate, which is the largest of the following three quantities:

1. The 95% uncertainty of the fitted intercept. The residual norm of the least-squares
   fit is used to estimate the standard error of the intercept, which is multiplied by the
   appropriate Student-t factor. This measures how closely the finite difference
   estimates conform to the expected truncation error behavior.

2. The roundoff error propagated through the extrapolation,

   .. math:: \sum_k |g_k| E_{\mathrm{round}}(h_k), \qquad E_{\mathrm{round}}(h) = \frac{\epsilon}{h^m} \sum_{j=1}^N \left| w_j f(x_0 + h s_j) \right|

   where :math:`E_{\mathrm{round}}(h)` is an estimate of the roundoff error of a single
   finite difference. An extrapolated value cannot be more accurate than the roundoff
   errors of the estimates from which it is computed allow.

3. The smaller of the differences between the candidate and the candidates of the two
   neighboring windows. Neighboring windows share all but one finite difference estimate,
   so two candidates on the plateau agree closely, while a candidate computed from
   roundoff-dominated estimates generally differs from both of its neighbors.

The candidate with the smallest error estimate is selected. The second and third
quantities are needed because each least-squares fit has only one residual degree of
freedom (see :data:`romberg_terms` below): its residual, and hence the first quantity,
can be very small by chance. When the minimum is taken over many windows, a window in
the roundoff-dominated region whose residual happens to be small would otherwise be
selected, with an error estimate much smaller than its actual error.

If the sampled range of step sizes encompasses the plateau region, selecting the
candidate with the smallest estimated error will typically provide an accurate
estimate of the derivative, with an error estimate that bounds the actual error.
No error estimate can compensate for a range of step sizes which does not contain the
plateau, however. The default range is centered on :math:`h_{\mathrm{opt}}`, which assumes
that the function varies on a length scale of order one. If the function has a
different characteristic length scale :math:`L` (for example :math:`L = 1/a` for
:math:`\sin(ax)` or :math:`e^{ax}`), the plateau lies near :math:`L h_{\mathrm{opt}}`,
and the step sizes should be scaled accordingly with :macro:`GSL_DERIV_STEP_USER_SCALE`.
For other functions, the user may need to tune the search range of :math:`h` values
to locate the plateau; the diagnostic output (see below) can be used for this purpose.

Parameters for Differentiation
==============================

This section details the parameters that users can tune for their
specific differentiation problem.

.. type:: gsl_deriv_stencil_t

   This specifies the type of stencil

   .. macro:: GSL_DERIV_STENCIL_CENTRAL

      This specifies a central stencil, e.g. :math:`\dots, -2, -1, 0, 1, 2, \dots`.

   .. macro:: GSL_DERIV_STENCIL_FORWARD

      This specifies a forward stencil, e.g. :math:`0, 1, 2, \dots`.

   .. macro:: GSL_DERIV_STENCIL_BACKWARD

      This specifies a backward stencil, e.g. :math:`\dots, -2, -1, 0`.

.. type:: gsl_deriv_step_policy_t

   The :code:`min_step` and :code:`max_step` of :type:`gsl_deriv_parameters` specify
   the range of base step sizes :math:`h`. The step policy determines how these
   base steps are converted to the actual step :math:`h_i` used for each input
   :math:`x_i`.

   .. macro:: GSL_DERIV_STEP_ABSOLUTE

      For a given point :math:`x_i`, the step size is :math:`h_i = h`.

   .. macro:: GSL_DERIV_STEP_RELATIVE

      For a given point :math:`x_i`, the step size is :math:`h_i = \max(|x_i|,1) h`.

   .. macro:: GSL_DERIV_STEP_USER_SCALE

      For a given point :math:`x_i`, the step size is :math:`h_i = u_i h`,
      where :math:`u` is a vector of user-supplied characteristic scales. The user-supplied
      value :math:`u_i` represents a characteristic magnitude or scale length
      for :math:`x_i` and should have the same physical units as :math:`x_i`.

.. type:: gsl_deriv_parameters

   This structure contains the parameters for adaptive and fixed-step differentiation::

     typedef struct
     {
       unsigned int deriv_order;            /* derivative order */
       unsigned int method_order;           /* method order */
       gsl_deriv_stencil_t stencil_type;    /* stencil type */
       unsigned int romberg_terms;          /* number of truncation terms in Romberg fit */
       gsl_deriv_step_policy_t step_policy; /* step policy */
       double min_step;                     /* minimum step size */
       double max_step;                     /* maximum step size */
       size_t nstep;                        /* number of step sizes */
       const gsl_vector * scales;           /* user-supplied scales, same length as x */
     } gsl_deriv_parameters;

   The parameter :data:`romberg_terms` specifies the number of truncation-error
   terms included in each Romberg extrapolation. If :data:`romberg_terms` is
   :math:`r`, each window of finite-difference estimates is fitted to a model of the form

   .. math:: D(h) = D_0 + \sum_{j=1}^{r} C_j h^{q_j},

   where :math:`D_0` is the extrapolated derivative at :math:`h=0`. For central
   stencils, the exponents are :math:`q_j=p+2(j-1)`, while for forward and backward
   stencils they are :math:`q_j=p+(j-1)`. Thus, for example, :code:`romberg_terms = 2`
   fits the leading two truncation-error terms in addition to the extrapolated derivative.
   Each fit contains :math:`r+1` fitted parameters and uses :math:`r+2` finite-difference
   estimates, leaving one residual degree of
   freedom that is used to estimate the uncertainty in :math:`D_0` (this uncertainty is
   combined with the propagated roundoff error and the consistency with neighboring
   windows, as described above). Increasing
   :data:`romberg_terms` allows additional truncation-error terms to be modeled and
   extrapolated away, but requires a
   larger window of derivative estimates and may increase sensitivity to roundoff error
   or to departures from the assumed asymptotic truncation-error expansion. For example,
   :code:`romberg_terms = 3` often gives more accurate results for analytic functions,
   while the default of 2 is better suited to functions of limited smoothness, such as
   splines, whose higher derivatives are discontinuous.

   The number of step sizes :data:`nstep` must be large enough to form at least
   one Romberg window, and increasing :data:`romberg_terms` reduces the number of
   overlapping Romberg candidates available for a fixed :data:`nstep`. The relation is

   .. math:: N_{\mathrm{cand}} = N_{\mathrm{step}} - r - 1

.. function:: gsl_deriv_parameters gsl_deriv_default_parameters(const unsigned int deriv_order, const unsigned int method_order)

   This function returns default parameters given the
   derivative order :math:`m`, specified by :data:`deriv_order` and
   the method order :math:`p`, specified by :data:`method_order`.
   The default values are shown in the table below.

   ============= ==================================
   Parameter     Default Value
   ============= ==================================
   stencil_type  :macro:`GSL_DERIV_STENCIL_CENTRAL`
   romberg_terms :math:`2`
   step_policy   :macro:`GSL_DERIV_STEP_ABSOLUTE`
   min_step      :math:`10^{-2} h_{\mathrm{opt}}`
   max_step      :math:`10^3 h_{\mathrm{opt}}`
   nstep         :math:`25`
   scales        :code:`NULL`
   ============= ==================================

Providing Functions for Differentiation
=======================================

.. type:: gsl_deriv_function

   This data type defines a general function with vector input and scalar output, with parameters,

   .. math:: f: \mathbb{R}^n \rightarrow \mathbb{R}

   :code:`double (* f) (const gsl_vector * x, void * params)`

      this function should return the value
      :math:`f(x,params)` for argument :data:`x` and parameters :data:`params`

   :code:`size_t n`

      the size of the input vector :data:`x`.

   :code:`void * params`

      a pointer to the parameters of the function

.. type:: gsl_deriv_vec_function

   This data type defines a general function with vector input and vector output, with parameters,

   .. math:: F: \mathbb{R}^n \rightarrow \mathbb{R}^m

   :code:`int (* f) (const gsl_vector * x, void * params, gsl_vector * f)`

      This function should store the value :math:`F(x)` into the vector
      :data:`f` for argument :data:`x` and parameters :data:`params`.

   :code:`size_t n`

      the size of the input vector :data:`x`.

   :code:`size_t m`

      the size of the output vector :data:`f`.

   :code:`void * params`

      a pointer to the parameters of the function

Fixed Step Differentiation
==========================

These functions take a fixed step size :math:`h` provided by the user and
estimate the function derivative.

.. function:: int gsl_deriv_fixed (const gsl_deriv_parameters * params, const gsl_function * f, const double x, const double h, double * result, double * esterr, unsigned int * neval)
              int gsl_deriv_fixed_v (const gsl_deriv_parameters * params, const gsl_function * f, const gsl_vector * x, const double h, gsl_vector * result, gsl_vector * esterr, unsigned int * neval)

   These functions estimate :math:`D_h^{(m)} f(x)` given the step size
   :math:`h` and point :math:`x`. The :data:`params` input specifies
   the derivative order, method order, and stencil type. The step
   size :math:`h` supplied to these functions is used directly and is
   not modified according to :code:`step_policy`.
   The result is returned in :data:`result`.
   The function additionally computes :math:`D_{h/2}^{(m)} f(x)` in order to
   calculate an absolute error estimate (which includes both truncation and roundoff),
   which is returned in :data:`esterr`. The explicit formula is

   .. math:: \mathrm{esterr} = \left| D_h^{(m)} f(x) - D_{h/2}^{(m)} f(x) \right| + E_{\mathrm{round}}(h)

   with

   .. math:: E_{\mathrm{round}}(h) = \frac{\epsilon}{h^m} \sum_{k=1}^N \left| w_k f(x + h s_k) \right|

   The total number of function evaluations is returned in :data:`neval`.
   The function :func:`gsl_deriv_fixed` operates on a single point,
   while :func:`gsl_deriv_fixed_v` applies the same differentiation
   operation independently to each element of the input vector.

.. function:: int gsl_deriv_directional_fixed(const gsl_deriv_parameters * params, const gsl_deriv_function * f, const gsl_vector * x, const gsl_vector * v, const double h, double * result, double * esterr, unsigned int * neval, gsl_vector * work)

   This function estimates the :math:`m`-th directional derivative of the function
   :math:`f(x)` in the direction of the vector :math:`v`,
   :math:`D_v^{(m)} f(x) = \frac{d^m}{dt^m} f(x + t v) \left. \right|_{t=0}`.
   The fixed absolute step size :math:`h` applies to the scalar parameter
   :math:`t`, so finite difference evaluations are performed at
   :math:`x + h s_k v`. The vector :math:`v` is used as supplied and
   is not normalized. The estimated
   directional derivative is returned in :data:`result` and the
   estimated absolute error is returned in :data:`esterr`. The
   number of function evaluations is returned in :data:`neval`.
   Additional workspace of length :math:`n` is required in :data:`work`,
   where :math:`n` is the length of the vectors :math:`x` and :math:`v`.

.. function:: int gsl_deriv_gradient_fixed(const gsl_deriv_parameters * params, const gsl_deriv_function * f, const gsl_vector * x, const gsl_vector * h, gsl_vector * result, gsl_vector * esterr, unsigned int * neval, gsl_vector * work)

   This function computes the vector of partial derivatives
   :math:`\partial^m f(x) / \partial x_j^m` using fixed absolute
   step sizes provided in the vector :math:`h`. The element
   :math:`h_j` specifies the fixed absolute step used to
   compute the derivative with respect to :math:`x_j` (the
   remaining components of :math:`x` are held fixed).
   The output vector is stored in :data:`result` and the absolute
   estimated errors for each partial derivative are stored
   in the output :data:`esterr`. The number of function evaluations
   is returned in :data:`neval`.
   Additional workspace of length :math:`n` is required in :data:`work`,
   where :math:`n` is the length of the vector :math:`x`.

Adaptive Differentiation
========================

These functions use the adaptive algorithm described above
to select a finite-difference step size and compute the
derivative estimate and associated absolute error.

.. function:: int gsl_deriv_adaptive(const gsl_deriv_parameters * params, const gsl_function * f, const double x, double * result, double * esterr, double * final_step, unsigned int * neval, gsl_deriv_diag * diag)
              int gsl_deriv_adaptive_v(const gsl_deriv_parameters * params, const gsl_function * f, const gsl_vector * x, gsl_vector * result, gsl_vector * esterr, gsl_vector * final_step, unsigned int * neval, gsl_deriv_diag * diag)

   These functions estimate :math:`f^{(m)}(x)` using the
   adaptive algorithm described above to select a finite-difference
   step size. The :data:`params` input specifies
   the parameters controlling the adaptive algorithm.
   The estimated derivative and absolute error are returned in :data:`result`
   and :data:`esterr`, respectively, and the selected step size (after
   applying the step policy) is returned in :data:`final_step`.
   The total number of function evaluations is returned in :data:`neval`.
   The function :func:`gsl_deriv_adaptive` operates on a single point,
   while :func:`gsl_deriv_adaptive_v` applies the same differentiation
   operation independently to each element of the input vector; for the
   latter, :data:`result`, :data:`esterr`, and :data:`final_step` contain one
   value for each input point.

   If :data:`diag` is non-null, diagnostic information is returned for
   the derivative estimate. For :func:`gsl_deriv_adaptive_v`, the diagnostic
   information corresponds only to the last element processed.
   This information can be used to assess the reliability of the estimate
   and whether the sampled step range adequately encompasses the plateau region.

.. function:: int gsl_deriv_directional(const gsl_deriv_parameters * params, const gsl_deriv_function * f, const gsl_vector * x, const gsl_vector * v, double * result, double * esterr, double * final_step, unsigned int * neval, gsl_deriv_diag * diag, gsl_vector * work)

   This function estimates the :math:`m`-th directional derivative of the function
   :math:`f(x)` in the direction of the vector :math:`v`,
   :math:`D_v^{(m)} f(x) = \frac{d^m}{dt^m} f(x + t v) \left. \right|_{t=0}`,
   using the adaptive algorithm described above. The adaptive step sizes
   apply to the scalar parameter :math:`t`, so finite-difference evaluations
   are performed along the line :math:`x + t v`. Since the derivative is
   evaluated at :math:`t=0`, :macro:`GSL_DERIV_STEP_RELATIVE` is equivalent
   to :macro:`GSL_DERIV_STEP_ABSOLUTE` for this function. :macro:`GSL_DERIV_STEP_USER_SCALE`
   is not supported for directional derivatives.
   The vector :math:`v` is used as supplied and is not normalized. The estimated
   directional derivative is returned in :data:`result` and the
   estimated absolute error is returned in :data:`esterr`. The
   selected step size is returned in :data:`final_step`.
   The number of function evaluations is returned in :data:`neval`.
   Additional workspace of length :math:`n` is required in :data:`work`,
   where :math:`n` is the length of the vectors :math:`x` and :math:`v`.

   If :data:`diag` is non-null, diagnostic information is returned for
   the derivative estimate. This information can be used to assess the
   reliability of the estimate and whether the sampled step range adequately
   encompasses the plateau region.

.. function:: int gsl_deriv_gradient(const gsl_deriv_parameters * params, const gsl_deriv_function * f, const gsl_vector * x, gsl_vector * result, gsl_vector * esterr, gsl_vector * final_step, unsigned int * neval, gsl_vector * work)

   This function estimates the vector of partial derivatives
   :math:`\partial^m f(x) / \partial x_j^m` using the adaptive
   algorithm described above to select finite-difference
   step sizes independently for each input :math:`x_j`, with
   the remaining components of :math:`x` held fixed.
   The derivative estimates and corresponding estimated absolute
   errors are returned in :data:`result` and :data:`esterr`,
   respectively.  The selected step sizes (after
   applying the step policy) are returned in :data:`final_step`.
   The number of function evaluations is returned in :data:`neval`.
   Additional workspace of length :math:`n` is required in :data:`work`,
   where :math:`n` is the length of the vector :math:`x`.

Diagnostic Output
-----------------

Some of the adaptive functions provide a diagnostic output
structure,

.. type:: gsl_deriv_diag

which can be analyzed to determine if the adaptive algorithm
is correctly locating the plateau region between where truncation
and roundoff error dominate. An example program which analyzes
this structure is provided in the :ref:`sec_examples-diag` section below.

.. function:: int gsl_deriv_diag_print(FILE * fp, const gsl_deriv_diag * diag)

   This function outputs the diagnostic information contained in
   :data:`diag` to the file pointer :data:`fp`. The printed output
   is in ASCII format suitable for use with the :code:`gnuplot` program.

Jacobian Calculation
====================

For the map

.. math:: F: \mathbb{R}^n \rightarrow \mathbb{R}^l

with

.. math:: F(x) = \left\{ f_1(x), f_2(x), \dots, f_l(x) \right\}

GSL computes the :math:`l`-by-:math:`n` matrix

.. math:: J_{ij}^{(m)} = \frac{\partial^m f_i(x)}{\partial x_j^m}

The classical Jacobian corresponds to :math:`m = 1`, but GSL supports
the calculation of higher order derivatives by specifying the
:code:`deriv_order` parameter.

.. function:: gsl_deriv_jacobian_workspace * gsl_deriv_jacobian_alloc(const gsl_deriv_parameters * params, const size_t nrows, const size_t ncols)

   This function allocates a workspace suitable for calculations of
   :math:`J^{(m)}` for both fixed-step and adaptive routines. The input :data:`nrows`
   specifies :math:`l`, the dimension of the output of :math:`F(x)`,
   and the number of rows of the derivative matrix. The input :data:`ncols`
   specifies :math:`n`, the dimension of the input :math:`x` and 
   the number of columns of the derivative matrix.

.. function:: void gsl_deriv_jacobian_free(gsl_deriv_jacobian_workspace * w)

   This function frees the memory associated with :data:`w`.

.. function:: int gsl_deriv_jacobian_fixed(const gsl_deriv_vec_function * F, const gsl_vector * x, const gsl_vector * h, gsl_matrix * J, gsl_matrix * esterr, unsigned int * neval, gsl_deriv_jacobian_workspace * w)

   This function estimates the :math:`l`-by-:math:`n` matrix :math:`J^{(m)}`
   associated with the function :math:`F(x)`. The function uses fixed step
   sizes, specified in :math:`h`, to compute the derivatives. The element
   :math:`h_j` specifies the absolute step size used to compute all derivatives
   in column :math:`j` of the derivative matrix, corresponding to differentiation with
   respect to :math:`x_j`. The estimated absolute error of each element of
   the derivative matrix is returned in the :math:`l`-by-:math:`n` matrix :data:`esterr`.
   The total number of function evaluations is provided in :data:`neval` on output.
   The workspace :data:`w` must be allocated by calling
   :func:`gsl_deriv_jacobian_alloc`.

.. function:: int gsl_deriv_jacobian(const gsl_deriv_vec_function * F, const gsl_vector * x, gsl_matrix * J, gsl_matrix * esterr, gsl_matrix * final_step, unsigned int * neval, gsl_deriv_jacobian_workspace * w)

   This function estimates the :math:`l`-by-:math:`n` matrix :math:`J^{(m)}`
   associated with the function :math:`F(x)`, using the adaptive algorithm described
   above to select finite-difference step sizes. For each input :math:`x_j`,
   the function computes finite-difference estimates for the entire :math:`j`-th
   column using a common sequence of step sizes. Romberg extrapolation and
   candidate selection are then performed independently for each matrix
   element :math:`J_{ij}^{(m)}`, so different elements within a column may select different
   final step sizes.  The estimated absolute error of each element of
   the derivative matrix is returned in the :math:`l`-by-:math:`n` matrix :data:`esterr`.
   The selected step sizes (after applying the step policy) are returned
   in the :math:`l`-by-:math:`n` matrix :data:`final_step`.
   The total number of function evaluations is provided in :data:`neval` on output.
   The workspace :data:`w` must be allocated by calling
   :func:`gsl_deriv_jacobian_alloc`.

Obsolete Functions
==================

.. function:: int gsl_deriv_central (const gsl_function * f, double x, double h, double * result, double * abserr)

   This function computes the numerical derivative of the function :data:`f`
   at the point :data:`x` using an adaptive central difference algorithm with
   a step-size of :data:`h`.   The derivative is returned in :data:`result` and an
   estimate of its absolute error is returned in :data:`abserr`.

   The initial value of :data:`h` is used to estimate an optimal step-size,
   based on the scaling of the truncation error and round-off error in the
   derivative calculation.  The derivative is computed using a 5-point rule
   for equally spaced abscissae at :math:`x - h`, :math:`x - h/2`, :math:`x`,
   :math:`x + h/2`, :math:`x+h`, with an error estimate taken from the difference
   between the 5-point rule and the corresponding 3-point rule :math:`x-h`,
   :math:`x`, :math:`x+h`.  Note that the value of the function at :math:`x`
   does not contribute to the derivative calculation, so only 4-points are
   actually used.

.. function:: int gsl_deriv_forward (const gsl_function * f, double x, double h, double * result, double * abserr)

   This function computes the numerical derivative of the function :data:`f`
   at the point :data:`x` using an adaptive forward difference algorithm with
   a step-size of :data:`h`. The function is evaluated only at points greater
   than :data:`x`, and never at :data:`x` itself.  The derivative is returned in
   :data:`result` and an estimate of its absolute error is returned in
   :data:`abserr`.  This function should be used if :math:`f(x)` has a
   discontinuity at :data:`x`, or is undefined for values less than :data:`x`.

   The initial value of :data:`h` is used to estimate an optimal step-size,
   based on the scaling of the truncation error and round-off error in the
   derivative calculation.  The derivative at :math:`x` is computed using an
   "open" 4-point rule for equally spaced abscissae at :math:`x+h/4`,
   :math:`x + h/2`, :math:`x + 3h/4`, :math:`x+h`, with an error estimate taken
   from the difference between the 4-point rule and the corresponding
   2-point rule :math:`x+h/2`, :math:`x+h`. 

.. function:: int gsl_deriv_backward (const gsl_function * f, double x, double h, double * result, double * abserr)

   This function computes the numerical derivative of the function :data:`f`
   at the point :data:`x` using an adaptive backward difference algorithm
   with a step-size of :data:`h`. The function is evaluated only at points
   less than :data:`x`, and never at :data:`x` itself.  The derivative is
   returned in :data:`result` and an estimate of its absolute error is
   returned in :data:`abserr`.  This function should be used if :math:`f(x)`
   has a discontinuity at :data:`x`, or is undefined for values greater than
   :data:`x`.

   This function is equivalent to calling :func:`gsl_deriv_forward` with a
   negative step-size.

Examples
========

Ordinary Derivative
-------------------

The following code estimates the derivative of the function 
:math:`f(x) = \log{(1+x)}`
at :math:`x = 5` and at :math:`x = -0.9`.  Since the function :math:`f(x)` is
undefined for :math:`x \leq -1` a forward stencil is used at :math:`x=-0.9`
to ensure that all function evaluations remain within the domain.

.. include:: examples/diff.c
   :code:

Here is the output of the program,

.. include:: examples/diff.txt
   :code:

Directional Derivative
----------------------

The following code estimates the directional derivative of the function 
:math:`f(x) = \| x \|^2`. The derivative in the direction
:math:`v` is given by :math:`D_v f(x) = 2 x \cdot v`.

.. include:: examples/diff_d.c
   :code:

Here is the output of the program,

.. include:: examples/diff_d.txt
   :code:

Gradient Example
----------------

The following code estimates the gradient of the function 
:math:`f(x) = x^T A x` for a symmetric matrix :math:`A`.
Assuming a Cartesian input vector, the gradient is given by
:math:`\nabla f(x) = 2 A x`.

.. include:: examples/diff_g.c
   :code:

Here is the output of the program,

.. include:: examples/diff_g.txt
   :code:

Jacobian Example
----------------

The following code estimates the Jacobian of the function 
:math:`F(x_1,x_2) = \left( x_1^2, x_1 x_2, x_2^2 \right)`.
The Jacobian matrix is given by

.. math:: J = \begin{pmatrix}
                2 x_1 & 0 \\
                x_2 & x_1 \\
                0 & 2 x_2
              \end{pmatrix}

.. include:: examples/diff_j.c
   :code:

Here is the output of the program,

.. include:: examples/diff_j.txt
   :code:

Higher Order Derivatives
------------------------

This example computes the first and tenth derivatives
(:math:`m=1` and :math:`m=10`) of the function :math:`f(x) = e^x`.
Since :math:`f^{(m)}(x) = e^x` for all :math:`m`, the exact result
is known for any derivative order.
The derivative order can be changed by using the :code:`deriv_order` argument
to :func:`gsl_deriv_default_parameters`. Higher order numerical
derivatives are increasingly sensitive to roundoff error, since the
roundoff contribution scales approximately as :math:`\epsilon/h^m`. The
adaptive algorithm therefore generally selects different step sizes as
:math:`m` changes.

The program source code is below.

.. include:: examples/diff_exp.c
   :code:

Here is the output of the program,

.. include:: examples/diff_exp.txt
   :code:

.. _sec_examples-diag:

Diagnostic Example
------------------

Here we illustrate the diagnostic output of the
adaptive differentiation routines on the function

.. math:: f(x) = x^{3/2}

which has derivatives

.. math::

   f'(x) &= \frac{3}{2} \sqrt{x} \\
   f''(x) &= \frac{3}{4 \sqrt{x}}

The domain of this function is :math:`x \geq 0`, and the one-sided
derivative :math:`f'(0) = 0` exists. However, :math:`f''(x)` diverges
as :math:`x \rightarrow 0^{+}`.  Therefore, at the point :math:`x=0`,
the function does not have the higher-order smoothness required for the
nominal finite-difference truncation-error expansion. In particular,
the finite-difference estimates do not follow the integer-power error
model assumed by the Romberg extrapolation.

We illustrate the situation by first computing the derivative at
:math:`x=1` which does satisfy smoothness requirements, and then at
:math:`x=0` which does not. We use a forward stencil for both examples
for two reasons: (1) since the domain of the function is :math:`x \geq 0`,
we want to avoid taking steps into the :math:`x < 0` region, and (2)
we want to hold the stencil type fixed when comparing the behavior at
the smooth point :math:`x=1` and the endpoint :math:`x=0`.

The source code is below.

.. include:: examples/diff_diag.c
   :code:

Here is the output for applying the adaptive algorithm to :math:`x=1`.

.. code-block:: console

   $ ./diff_diag 1 > diff_diag_x1.txt

.. include:: examples/diff_diag_x1_output.txt
   :code:

The diagnostic output is organized into 3 blocks. The
first block contains the initial derivative estimates
and roundoff error estimates computed by applying the
finite-difference formula at each step size. The second
block is the result of applying the Romberg extrapolation
and error estimation, and is usually the most useful data
to plot. The third block is the result after applying
the candidate filter (which removes any :code:`nan` values).
The figures below plot the data from the second block
(Romberg extrapolated derivative estimates and associated
error estimates).

.. figure:: /images/deriv_1.png

We see that the left figure exhibits a clear plateau in the
derivative estimates, and the right figure shows a clear
minimum in the error curve. At smaller step sizes the error
curve eventually increases and becomes irregular as roundoff
error becomes important, while at larger step sizes it rises
rapidly as truncation error dominates.
The derivative candidate corresponding to the minimum of the
estimated error curve is selected as the final derivative estimate.

Now we compute the result for :math:`x=0`.

.. code-block:: console

   $ ./diff_diag 0 > diff_diag_x0.txt

.. include:: examples/diff_diag_x0_output.txt
   :code:

We see that the adaptive algorithm estimates the derivative as
:math:`f'(0) \approx 1.269 \cdot 10^{-3}`, which is not satisfactory, as the
correct value is :math:`f'(0) = 0`. The diagnostic output below explains
the problem.

.. figure:: /images/deriv_2.png

In contrast to the :math:`x=1` case, the derivative estimates show no
plateau and the estimated error continues to decrease over the sampled
range as :math:`h` decreases. Consequently, the adaptive algorithm is driven toward the
smallest available steps rather than identifying an interior region
in which truncation and roundoff errors are balanced.

At :math:`x=0`, the finite-difference approximation becomes

.. math:: D_h f(0) = \frac{1}{h} \sum_k w_k \left( h s_k \right)^{3/2} = \left( \sum_k w_k s_k^{3/2} \right) \sqrt{h}

and so :math:`D_h f(0) \propto \sqrt{h}`. Thus, the observed convergence order
at :math:`x=0` is :math:`1/2`, regardless of the nominal method order :math:`p` of the
forward stencil. However, for a sufficiently smooth function, the Romberg model
for a forward stencil is

.. math:: D(h) = D_0 + C_1 h^p + C_2 h^{p+1} + \cdots

Therefore, the finite-difference approximations do not follow the asymptotic model
assumed by the Romberg extrapolation at :math:`x=0`.
Extending the search to smaller step sizes may reduce the raw finite-difference
error, but it does not correct the mismatch between the observed :math:`\sqrt{h}`
convergence and the integer-power truncation model assumed by the Romberg extrapolation.
For this particular problem, a sufficiently small fixed step can give a useful
approximation because the raw finite-difference estimate still converges to
:math:`f'(0) = 0`, albeit only as :math:`O(\sqrt{h})`.

References and Further Reading
==============================

The algorithms used by these functions are described in the following sources:

* Abramowitz and Stegun, *Handbook of Mathematical Functions*,
  Section 25.3.4, and Table 25.5 (Coefficients for Differentiation).

* J. D'Errico, *Adaptive robust numerical differentiation (DERIVEST)*,
  https://www.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation, MATLAB Central File Exchange, 2007.

* B. Fornberg, *Calculation of weights in finite difference formulas*,
  SIAM Rev., Vol 40, No. 3, pp. 685-691, 1998.

* S.D. Conte and Carl de Boor, *Elementary Numerical Analysis: An
  Algorithmic Approach*, McGraw-Hill, 1972.
