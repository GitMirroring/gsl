.. index::
   single: error function
   single: erf(x)
   single: erfc(x)
   single: erfcx(x)

The error function is described in Abramowitz & Stegun, Chapter 7.  The
functions in this section are declared in the header file
:file:`gsl_sf_erf.h`.

Error Function
--------------

.. function:: double gsl_sf_erf (const double x)
              int gsl_sf_erf_e (const double x, gsl_sf_result * result)

   These routines compute the error function :math:`\erf(x)`,
   where
   :math:`\erf(x) = (2/\sqrt{\pi}) \int_0^x dt \exp(-t^2)`.
.. Exceptional Return Values: none

Complementary Error Function
----------------------------

.. function:: double gsl_sf_erfc (const double x)
              int gsl_sf_erfc_e (const double x, gsl_sf_result * result)

   These routines compute the complementary error function 
   :math:`\erfc(x) = 1 - \erf(x) = (2/\sqrt{\pi}) \int_x^\infty \exp(-t^2)`
.. Exceptional Return Values: none

Scaled Complementary Error Function
-----------------------------------

.. function:: double gsl_sf_erfcx (const double x)
              int gsl_sf_erfcx_e (const double x, gsl_sf_result * result)

   These routines compute the scaled complementary error function

   .. math:: \textrm{erfcx}(x) = \exp{(x^2)} \erfc(x)

   The scaled complementary error function avoids the rapid
   underflow of :math:`\erfc(x)` for large positive :math:`x`.
   For large positive arguments it has the asymptotic behavior

   .. math:: \textrm{erfcx}(x) \sim \frac{1}{\sqrt{\pi} x}

Log Complementary Error Function
--------------------------------

.. function:: double gsl_sf_log_erfc (const double x)
              int gsl_sf_log_erfc_e (const double x, gsl_sf_result * result)

   These routines compute the logarithm of the complementary error function
   :math:`\log(\erfc(x))`.
.. Exceptional Return Values: none

Probability functions
---------------------

The probability functions for the Normal or Gaussian distribution are
described in Abramowitz & Stegun, Section 26.2.

.. function:: double gsl_sf_erf_Z (const double x)
              int gsl_sf_erf_Z_e (const double x, gsl_sf_result * result)

   These routines compute the Gaussian probability density function 
   :math:`Z(x) = (1/\sqrt{2\pi}) \exp(-x^2/2)`

.. function:: double gsl_sf_erf_Q (const double x)
              int gsl_sf_erf_Q_e (const double x, gsl_sf_result * result)

   These routines compute the upper tail of the Gaussian probability function 
   :math:`Q(x) = (1/\sqrt{2\pi}) \int_x^\infty dt \exp(-t^2/2)`

.. Exceptional Return Values: none

.. index::
   single: hazard function, normal distribution
   single:  Mills' ratio, inverse

The *hazard function* for the normal distribution, 
also known as the inverse Mills' ratio, is defined as,

.. only:: not texinfo

   .. math:: h(x) = {Z(x) \over Q(x)} = \sqrt{2 \over \pi} {\exp(-x^2 / 2) \over \erfc(x/\sqrt 2)}

.. only:: texinfo

   ::

      h(x) = Z(x)/Q(x) = \sqrt{2/\pi} \exp(-x^2 / 2) / \erfc(x/\sqrt 2)

It decreases rapidly as :math:`x` approaches :math:`-\infty` and asymptotes
to :math:`h(x) \sim x` as :math:`x` approaches :math:`+\infty`.

.. function:: double gsl_sf_hazard (const double x)
              int gsl_sf_hazard_e (const double x, gsl_sf_result * result)

   These routines compute the hazard function for the normal distribution.
.. Exceptional Return Values: GSL_EUNDRFLW
