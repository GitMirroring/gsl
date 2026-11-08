.. index::
   single: Legendre polynomials
   single: Legendre functions
   single: spherical harmonics
   single: conical functions
   single: hyperbolic space

The Legendre Functions and Legendre Polynomials are described in
Abramowitz & Stegun, Chapter 8.  These functions are declared in 
the header file :file:`gsl_sf_legendre.h`.

Legendre Polynomials
--------------------

.. function:: double gsl_sf_legendre_P1 (double x)
              double gsl_sf_legendre_P2 (double x)
              double gsl_sf_legendre_P3 (double x)
              int gsl_sf_legendre_P1_e (double x, gsl_sf_result * result)
              int gsl_sf_legendre_P2_e (double x, gsl_sf_result * result)
              int gsl_sf_legendre_P3_e (double x, gsl_sf_result * result)

   These functions evaluate the Legendre polynomials
   :math:`P_l(x)` using explicit representations for :math:`l = 1, 2, 3`.
.. Exceptional Return Values: none

.. function:: double gsl_sf_legendre_Pl (int l, double x)
              int gsl_sf_legendre_Pl_e (int l, double x, gsl_sf_result * result)

   These functions evaluate the Legendre polynomial :math:`P_l(x)`
   for a specific value of :data:`l`, :data:`x` subject to :math:`l \ge 0` and
   :math:`|x| \le 1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: int gsl_sf_legendre_Pl_array (int lmax, double x, double result_array[])
              int gsl_sf_legendre_Pl_deriv_array (int lmax, double x, double result_array[], double result_deriv_array[])

   These functions compute arrays of Legendre polynomials
   :math:`P_l(x)` and derivatives :math:`dP_l(x)/dx`
   for :math:`l = 0, \dots, lmax` and :math:`|x| \le 1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_legendre_Q0 (double x)
              int gsl_sf_legendre_Q0_e (double x, gsl_sf_result * result)

   These routines compute the Legendre function :math:`Q_0(x)` for
   :math:`x > -1` and :math:`x \ne 1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_legendre_Q1 (double x)
              int gsl_sf_legendre_Q1_e (double x, gsl_sf_result * result)

   These routines compute the Legendre function :math:`Q_1(x)` for
   :math:`x > -1` and :math:`x \ne 1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_legendre_Ql (int l, double x)
              int gsl_sf_legendre_Ql_e (int l, double x, gsl_sf_result * result)

   These routines compute the Legendre function :math:`Q_l(x)` for
   :math:`x > -1`, :math:`x \ne 1` and :math:`l \ge 0`.
.. Exceptional Return Values: GSL_EDOM

Conical Functions
-----------------

The Conical Functions :math:`P^\mu_{-(1/2)+i\lambda}(x)`
and :math:`Q^\mu_{-(1/2)+i\lambda}`
are described in Abramowitz & Stegun, Section 8.12.

.. function:: double gsl_sf_conicalP_half (double lambda, double x)
              int gsl_sf_conicalP_half_e (double lambda, double x, gsl_sf_result * result)

   These routines compute the irregular Spherical Conical Function
   :math:`P^{1/2}_{-1/2 + i \lambda}(x)` for :math:`x > -1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_conicalP_mhalf (double lambda, double x)
              int gsl_sf_conicalP_mhalf_e (double lambda, double x, gsl_sf_result * result)

   These routines compute the regular Spherical Conical Function
   :math:`P^{-1/2}_{-1/2 + i \lambda}(x)` for :math:`x > -1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_conicalP_0 (double lambda, double x)
              int gsl_sf_conicalP_0_e (double lambda, double x, gsl_sf_result * result)

   These routines compute the conical function
   :math:`P^0_{-1/2 + i \lambda}(x)` for :math:`x > -1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_conicalP_1 (double lambda, double x)
              int gsl_sf_conicalP_1_e (double lambda, double x, gsl_sf_result * result)

   These routines compute the conical function 
   :math:`P^1_{-1/2 + i \lambda}(x)` for :math:`x > -1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_conicalP_sph_reg (int l, double lambda, double x)
              int gsl_sf_conicalP_sph_reg_e (int l, double lambda, double x, gsl_sf_result * result)

   These routines compute the Regular Spherical Conical Function
   :math:`P^{-1/2-l}_{-1/2 + i \lambda}(x)`
   for :math:`x > -1` and :math:`l \ge -1`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_conicalP_cyl_reg (int m, double lambda, double x)
              int gsl_sf_conicalP_cyl_reg_e (int m, double lambda, double x, gsl_sf_result * result)

   These routines compute the Regular Cylindrical Conical Function
   :math:`P^{-m}_{-1/2 + i \lambda}(x)`
   for :math:`x > -1` and :math:`m \ge -1`.
.. Exceptional Return Values: GSL_EDOM

Radial Functions for Hyperbolic Space
-------------------------------------

The following spherical functions are specializations of Legendre
functions which give the regular eigenfunctions of the Laplacian on a
3-dimensional hyperbolic space :math:`H^3`.  Of particular interest is
the flat limit, :math:`\lambda \to \infty`, :math:`\eta \to 0`,
:math:`\lambda\eta` fixed.
  
.. function:: double gsl_sf_legendre_H3d_0 (double lambda, double eta)
              int gsl_sf_legendre_H3d_0_e (double lambda, double eta, gsl_sf_result * result)

   These routines compute the zeroth radial eigenfunction of the Laplacian on the
   3-dimensional hyperbolic space,

   .. math:: L^{H3d}_0(\lambda,\eta) := {\sin(\lambda\eta) \over \lambda\sinh(\eta)}

   for :math:`\eta \ge 0`.
   In the flat limit this takes the form
   :math:`L^{H3d}_0(\lambda,\eta) = j_0(\lambda\eta)`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_legendre_H3d_1 (double lambda, double eta)
              int gsl_sf_legendre_H3d_1_e (double lambda, double eta, gsl_sf_result * result)

   These routines compute the first radial eigenfunction of the Laplacian on
   the 3-dimensional hyperbolic space,

   .. math:: L^{H3d}_1(\lambda,\eta) := {1\over\sqrt{\lambda^2 + 1}} {\left(\sin(\lambda \eta)\over \lambda \sinh(\eta)\right)} \left(\coth(\eta) - \lambda \cot(\lambda\eta)\right)

   for :math:`\eta \ge 0`
   In the flat limit this takes the form 
   :math:`L^{H3d}_1(\lambda,\eta) = j_1(\lambda\eta)`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_legendre_H3d (int l, double lambda, double eta)
              int gsl_sf_legendre_H3d_e (int l, double lambda, double eta, gsl_sf_result * result)

   These routines compute the :data:`l`-th radial eigenfunction of the
   Laplacian on the 3-dimensional hyperbolic space :math:`\eta \ge 0` and
   :math:`l \ge 0`.
   In the flat limit this takes the form
   :math:`L^{H3d}_l(\lambda,\eta) = j_l(\lambda\eta)`.
.. Exceptional Return Values: GSL_EDOM

.. function:: int gsl_sf_legendre_H3d_array (int lmax, double lambda, double eta, double result_array[])

   This function computes an array of radial eigenfunctions
   :math:`L^{H3d}_l( \lambda, \eta)`
   for :math:`0 \le l \le lmax`.
.. Exceptional Return Values:
