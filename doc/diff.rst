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
window. The residual norm of the least-squares fit is used to estimate the standard
error of the fitted intercept, from which a 95% uncertainty estimate is calculated
using the appropriate Student-t factor. Thus, each Romberg window produces both a
candidate derivative estimate and associated error estimate, based on how
closely the finite difference estimates conform to the expected truncation error behavior.

If the sampled range of step sizes encompasses the plateau region, selecting the
candidate with the smallest estimated error will typically provide an accurate
estimate of the derivative.  For some functions,
the user may need to tune the search range of :math:`h` values to attempt to
locate the plateau.

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
   freedom that is used to estimate the uncertainty in :math:`D_0`. Increasing
   :data:`romberg_terms` allows additional truncation-error terms to be modeled and
   extrapolated away, but requires a
   larger window of derivative estimates and may increase sensitivity to roundoff error
   or to departures from the assumed asymptotic truncation-error expansion.

   The number of step sizes :data:`nstep` must be large enough to form at least
   one Romberg window, and increasing :data:`romberg_terms` reduces the number of
   overlapping Romberg candidates available for a fixed :data:`nstep`. The relation is

   .. math:: N_{\mathrm{cand}} = N_{\mathrm{step}} - r - 1

.. function:: gsl_deriv_parameters gsl_deriv_default_parameters(const unsigned int deriv_order, const unsigned int method_order)

   This function calculates default parameters given the
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

   This function computes the :math:`m`-th directional derivative of the function
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

   This function allocates a workspace suitable for Jacobian calculations
   for both fixed-step and adaptive routines. The input :data:`nrows`
   specifies :math:`l`, the dimension of the output of :math:`F(x)`,
   and the number of rows of the Jacobian. The input :data:`ncols`
   specifies :math:`n`, the dimension of the input :math:`x` and 
   the number of columns of the Jacobian.

.. function:: void gsl_deriv_jacobian_free(gsl_deriv_jacobian_workspace * w)

   This function frees the memory associated with :data:`w`.

.. function:: int gsl_deriv_jacobian_fixed(const gsl_deriv_vec_function * F, const gsl_vector * x, const gsl_vector * h, gsl_matrix * J, gsl_matrix * esterr, unsigned int * neval, gsl_deriv_jacobian_workspace * w)

   This function computes the :math:`l`-by-:math:`n` Jacobian matrix :math:`J`
   associated with the function :math:`F(x)`. The function uses fixed step
   sizes, specified in :math:`h`, to compute the derivatives. The element
   :math:`h_j` specifies the absolute step size used to compute all derivatives
   in column :math:`j` of the Jacobian, corresponding to differentiation with
   respect to :math:`x_j`. The estimated absolute error of each element of
   the Jacobian is returned in the :math:`l`-by-:math:`n` matrix :data:`esterr`.
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

The following code estimates the derivative of the function 
:math:`f(x) = x^{3/2}`
at :math:`x = 2` and at :math:`x = 0`.  The function :math:`f(x)` is
undefined for :math:`x < 0` so the derivative at :math:`x=0` is computed
using :func:`gsl_deriv_forward`.

.. include:: examples/diff.c
   :code:

Here is the output of the program,

.. include:: examples/diff.txt
   :code:

References and Further Reading
==============================

The algorithms used by these functions are described in the following sources:

* Abramowitz and Stegun, *Handbook of Mathematical Functions*,
  Section 25.3.4, and Table 25.5 (Coefficients for Differentiation).

* J. D'Errico, *Adaptive robust numerical differentiation (DERIVEST)*,
  https://www.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation, MATLAB Central File Exchange, 2026.

* B. Fornberg, *Calculation of weights in finite difference formulas*,
  SIAM Rev., Vol 40, No. 3, pp. 685-691, 1998.

* S.D. Conte and Carl de Boor, *Elementary Numerical Analysis: An
  Algorithmic Approach*, McGraw-Hill, 1972.
