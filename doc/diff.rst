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
finite differencing.  An adaptive algorithm is used to find the best
choice of finite difference and to estimate the error in the derivative.
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

For a given derivative order :math:`m` and method order :math:`p`,
the stencil size and offsets can be determined, and the weights :math:`w_k`
can be computed to high accuracy using the Fornberg algorithm.

By far, the most difficult problem in numerical differentiation is to determine
the optimal step size :math:`h`. As :math:`h \rightarrow 0`, the truncation
error shrinks. However, roundoff error, caused by catastrophic cancellation
of the terms in the weighted sum, grows as :math:`h \rightarrow 0`. Therefore,
there is typically a *plateau* of optimal step sizes in between these two regimes.
If :math:`h` is too large, truncation error dominates, but if :math:`h` is too small,
roundoff error dominates. A classic theoretical step size, designed to balance
truncation and roundoff error is

.. math:: h_{opt} \sim \epsilon^{\frac{1}{m+p}}

where :math:`\epsilon` is the machine epsilon value (:code:`GSL_DBL_EPSILON`).
Note that there are proportionality constants also which depend on
:math:`m`, :math:`p`, and function values at :math:`x_0`. For a typical
black box function, these are not known in advance and so
:math:`\epsilon^{1/(m+p)}` is often used as an initial step size.

The adaptive algorithm used by GSL, which is based on the DERIVEST
algorithm, computes derivative estimates :math:`D_h^{(m)} f(x_0)`
for multiple step sizes :math:`h` spanning several decades around the value
:math:`h_{opt}`. Then, Richardson extrapolation is used to cancel
the leading :math:`O(h^p)` error term to improve the derivative estimate.
The Richardson extrapolated derivative estimates are then grouped into
successive windows. For each window, the estimates are fitted by least
squares to the expected truncation error expansion. The residual norm
of the least-squares fit provides an estimate of the error in the derivative
at each step size.

If the original range of step sizes :math:`h` encompasses the optimal plateau,
then the derivative estimate corresponding to the least calculated error will
typically provide the best estimate of the function derivative. For some functions,
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
   a dimensionless range of step sizes :math:`h`. But for a given input :math:`x`,
   it may be desirable to add units to the final step size, especially if the
   input variables have characteristic scale lengths which are not :math:`O(1)`.
   This parameter specifies the type of steps taken by the adaptive differentiation algorithm.

   .. macro:: GSL_DERIV_STEP_ABSOLUTE

      For a given point :math:`x_i`, the step size is calculated as :math:`h`

   .. macro:: GSL_DERIV_STEP_RELATIVE

      For a given point :math:`x_i`, the step size is calculated as :math:`\max(|x_i|,1) h`

   .. macro:: GSL_DERIV_STEP_USER_SCALE

      For a given point :math:`x_i`, the step size is calculated as :math:`u_i h`,
      where :math:`u` is a vector of user-supplied scale factors.

.. type:: gsl_deriv_parameters

   This structure contains the parameters for adaptive and fixed-step differentiation::

     typedef struct
     {
       unsigned int deriv_order;            /* derivative order */
       unsigned int method_order;           /* method order */
       gsl_deriv_stencil_t stencil_type;    /* stencil type */
       unsigned int romberg_terms;          /* number of Romberg terms */
       gsl_deriv_step_policy_t step_policy; /* step policy */
       double min_step;                     /* minimum step size */
       double max_step;                     /* maximum step size */
       size_t nstep;                        /* number of step sizes */
       gsl_vector * scales;                 /* user-supplied scales, same length as x */
     } gsl_deriv_parameters;

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
   min_step      :math:`10^{-2} h_{opt}`
   max_step      :math:`10^3 h_{opt}`
   nstep         :math:`25`
   scales        :code:`NULL`
   ============= ==================================

Fixed Step Differentiation
==========================

These functions take a fixed step size :math:`h` provided by the user and
and estimate the function derivative.

.. function:: int gsl_deriv_fixed (const gsl_deriv_parameters * params, const gsl_function * f, const double x, const double h, double * result, double * esterr, unsigned int * neval)
              int gsl_deriv_fixed_v (const gsl_deriv_parameters * params, const gsl_function * f, const gsl_vector * x, const double h, gsl_vector * result, gsl_vector * esterr, unsigned int * neval)

   These functions estimate :math:`D_h^{(m)} f(x)` given the step size
   :math:`h` and point :math:`x`. The :data:`params` input specifies
   the derivative order, method order, and stencil type.
   The result is returned in :data:`result`.
   The function additionally computes :math:`D_{h/2}^{(m)} f(x)` in order to
   calculate an error estimate (which includes both truncation and roundoff),
   which is returned in :data:`esterr`. The total number of function evaluations
   is returned in :data:`neval`. Functions which operate on both scalar and vector
   inputs are provided.

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

* S.D. Conte and Carl de Boor, *Elementary Numerical Analysis: An
  Algorithmic Approach*, McGraw-Hill, 1972.
