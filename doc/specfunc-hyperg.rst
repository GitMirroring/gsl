.. index::
   single: hypergeometric functions
   single: confluent hypergeometric functions

Hypergeometric functions are described in Abramowitz & Stegun, Chapters
13 and 15.  These functions are declared in the header file
:file:`gsl_sf_hyperg.h`.

.. function:: double gsl_sf_hyperg_0F1 (double c, double x)
              int gsl_sf_hyperg_0F1_e (double c, double x, gsl_sf_result * result)

   These routines compute the hypergeometric function
   
   .. only:: not texinfo
      
      .. math:: {}_0F_1(c,x)

   .. only:: texinfo

      .. math:: 0F1(c,x)

.. It is related to Bessel functions
.. 0F1[c,x] =
..   Gamma[c]    x^(1/2(1-c)) I_(c-1)(2 Sqrt[x])
..   Gamma[c] (-x)^(1/2(1-c)) J_(c-1)(2 Sqrt[-x])
.. exceptions: GSL_EOVRFLW, GSL_EUNDRFLW

.. function:: double gsl_sf_hyperg_1F1_int (const int m, const int n, const double x)
              int gsl_sf_hyperg_1F1_int_e (const int m, const int n, const double x, gsl_sf_result * result)

   These routines compute the confluent hypergeometric function

   .. only:: not texinfo

      .. math:: {}_1F_1(m,n,x) = M(m,n,x)

   .. only:: texinfo

      .. math:: 1F1(m,n,x) = M(m,n,x)

   for integer parameters :data:`m`, :data:`n`, as
   :func:`gsl_sf_hyperg_1F1`.
.. exceptions: GSL_EDOM, GSL_ELOSS, GSL_EOVRFLW, GSL_EUNDRFLW

.. function:: double gsl_sf_hyperg_1F1 (const double a, const double b, const double x)
              int gsl_sf_hyperg_1F1_e (const double a, const double b, const double x, gsl_sf_result * result)

   These routines compute the confluent hypergeometric function

   .. only:: not texinfo

      .. math:: {}_1F_1(a,b,x) = M(a,b,x) = \sum_{k=0}^\infty {(a)_k \over (b)_k} {x^k \over k!}

   .. only:: texinfo

      .. math:: 1F1(a,b,x) = M(a,b,x) = \sum_{k=0}^\infty (a)_k / (b)_k x^k / k!

   for general parameters :data:`a`, :data:`b`.  When :math:`b` is a
   non-positive integer the function is defined only if :math:`a` is an
   integer with :math:`b \le a \le 0`, where the series terminates before
   the zero of :math:`(b)_k`; it is then the polynomial
   :math:`\sum_{k=0}^{-a} (a)_k/(b)_k \, x^k/k!`, and otherwise a domain
   error is returned.  (Releases before 2.9 returned
   :math:`e^x M(b-a,b,-x)` for some of these arguments, which is a
   different polynomial.)  :math:`M(0,b,x) = 1` for all :math:`x`.

   The function is evaluated by the series (with its cancellation
   tracked, and a double-double summation when needed), Kummer's
   transformation :math:`M(a,b,x) = e^x M(b-a,b,-x)` for :math:`x < 0`,
   the asymptotic expansion for large :math:`x`, and the recurrences in
   :math:`a` and :math:`b` otherwise.  The error estimate accounts for
   the cancellation and the growth of rounding errors in the recurrences.
   If no significant digit is left, as can happen for very large
   parameters of opposite signs, the error code :macro:`GSL_ELOSS` is
   returned.
.. exceptions: GSL_EDOM, GSL_ELOSS, GSL_EOVRFLW, GSL_EUNDRFLW

.. function:: double gsl_sf_hyperg_1F1_renorm (const double a, const double b, const double x)
              int gsl_sf_hyperg_1F1_renorm_e (const double a, const double b, const double x, gsl_sf_result * result)

   These routines compute the renormalized confluent hypergeometric
   function

   .. only:: not texinfo

      .. math:: {}_1F_1(a,b,x) / \Gamma(b) = {\bf M}(a,b,x)

   .. only:: texinfo

      .. math:: 1F1(a,b,x) / \Gamma(b)

   which is defined for all :math:`b`.  For :math:`b = -n` a
   non-positive integer the limit

   .. math:: {(a)_{n+1} \over (n+1)!} x^{n+1} M(a+n+1, n+2, x)

   is returned; it is zero when :math:`a` is an integer with
   :math:`b \le a \le 0`.
.. exceptions: GSL_EDOM, GSL_ELOSS, GSL_EOVRFLW, GSL_EUNDRFLW

.. function:: double gsl_sf_hyperg_U_int (int m, int n, double x)
              int gsl_sf_hyperg_U_int_e (int m, int n, double x, gsl_sf_result * result)

   These routines compute the confluent hypergeometric function
   :math:`U(m,n,x)` for integer parameters :data:`m`, :data:`n`.
.. exceptions:

.. function:: int gsl_sf_hyperg_U_int_e10_e (int m, int n, double x, gsl_sf_result_e10 * result)

   This routine computes the confluent hypergeometric function
   :math:`U(m,n,x)` for integer parameters :data:`m`, :data:`n` using the
   :type:`gsl_sf_result_e10` type to return a result with extended range.

.. function:: double gsl_sf_hyperg_U (double a, double b, double x)
              int gsl_sf_hyperg_U_e (double a, double b, double x, gsl_sf_result * result)

   These routines compute the confluent hypergeometric function :math:`U(a,b,x)`.
.. exceptions:

.. function:: int gsl_sf_hyperg_U_e10_e (double a, double b, double x, gsl_sf_result_e10 * result)

   This routine computes the confluent hypergeometric function
   :math:`U(a,b,x)` using the :type:`gsl_sf_result_e10` type to return a
   result with extended range. 
.. exceptions:

.. function:: double gsl_sf_hyperg_2F1 (const double a, const double b, const double c, const double x)
              int gsl_sf_hyperg_2F1_e (const double a, const double b, const double c, const double x, gsl_sf_result * result)

   These routines compute the Gauss hypergeometric function

   .. only:: not texinfo

      .. math:: {}_2F_1(a,b,c,x) = F(a,b,c,x)

   .. only:: texinfo

      .. math:: 2F1(a,b,c,x) = F(a,b,c,x)
         
   for real :math:`x < 1`.  At :math:`x = 1` the value
   :math:`\Gamma(c)\Gamma(c-a-b)/(\Gamma(c-a)\Gamma(c-b))` is returned when
   :math:`c - a - b > 0`.  When :math:`a` or :math:`b` is a non-positive
   integer the series terminates and the function, a polynomial in :math:`x`,
   is computed for all real :math:`x`; similarly when :math:`c-a` or
   :math:`c-b` is a non-positive integer, where
   :math:`F = (1-x)^{c-a-b}` times a polynomial.  For :math:`x > 1` the
   function is otherwise complex and a domain error is returned, as it is
   when :math:`c` is a non-positive integer and the series does not
   terminate before the pole.

   The function is evaluated with the transformation of :math:`x` that gives
   the smallest variable (:math:`x`, :math:`x/(x-1)`, :math:`1-x` or
   :math:`1/(1-x)`); the connection formulas are evaluated uniformly in
   the distance of :math:`c-a-b` (or :math:`b-a`) from an integer, so that
   there is no loss of accuracy near the degenerate cases.  The error
   estimate accounts for cancellation between terms.  For very large
   parameters the cancellation can exceed what is recoverable in double
   precision (with double-double summation as a fallback); if no
   significant digit is left the error code :macro:`GSL_ELOSS` is returned.
.. exceptions: GSL_EDOM, GSL_ELOSS, GSL_EOVRFLW

.. function:: double gsl_sf_hyperg_2F1_conj (const double aR, const double aI, const double c, const double x)
              int gsl_sf_hyperg_2F1_conj_e (const double aR, const double aI, const double c, const double x, gsl_sf_result * result)

   These routines compute the Gauss hypergeometric function

   .. only:: not texinfo

      .. math:: {}_2F_1(a_R + i a_I, aR - i aI, c, x)

   .. only:: texinfo

      .. math:: 2F1(a_R + i a_I, aR - i aI, c, x)

   with complex parameters for real :math:`x < 1`, and for :math:`x = 1`
   when :math:`c - 2 a_R > 0`.
.. exceptions: GSL_EDOM, GSL_ELOSS, GSL_EOVRFLW

.. function:: double gsl_sf_hyperg_2F1_renorm (const double a, const double b, const double c, const double x)
              int gsl_sf_hyperg_2F1_renorm_e (const double a, const double b, const double c, const double x, gsl_sf_result * result)

   These routines compute the renormalized Gauss hypergeometric function

   .. only:: not texinfo

      .. math:: {}_2F_1(a,b,c,x) / \Gamma(c)

   .. only:: texinfo

      .. math:: 2F1(a,b,c,x) / \Gamma(c)

   on the domain of :func:`gsl_sf_hyperg_2F1`.  For :math:`c = -n` a
   non-positive integer the limit

   .. math:: {(a)_{n+1} (b)_{n+1} \over (n+1)!} x^{n+1} F(a+n+1, b+n+1; n+2; x)

   is returned.
.. exceptions: GSL_EDOM, GSL_ELOSS, GSL_EOVRFLW

.. function:: double gsl_sf_hyperg_2F1_conj_renorm (const double aR, const double aI, const double c, const double x)
              int gsl_sf_hyperg_2F1_conj_renorm_e (const double aR, const double aI, const double c, const double x, gsl_sf_result * result)

   These routines compute the renormalized Gauss hypergeometric function

   .. only:: not texinfo

      .. math:: {}_2F_1(a_R + i a_I, a_R - i a_I, c, x) / \Gamma(c)

   .. only:: texinfo

      .. math:: 2F1(a_R + i a_I, a_R - i a_I, c, x) / \Gamma(c)

   on the domain of :func:`gsl_sf_hyperg_2F1_conj`, with the limit for
   non-positive integer :math:`c` as for :func:`gsl_sf_hyperg_2F1_renorm`.
.. exceptions:

.. function:: double gsl_sf_hyperg_2F0 (double a, double b, double x)
              int gsl_sf_hyperg_2F0_e (double a, double b, double x, gsl_sf_result * result)

   These routines compute the hypergeometric function
   
   .. only:: not texinfo
      
      .. math:: {}_2F_0(a,b,x)

   .. only:: texinfo

      .. math:: 2F0(a,b,x)

   The series representation is a divergent hypergeometric series.
   However, for :math:`x < 0` we have 

   .. only:: not texinfo

      .. math:: {}_2F_0(a,b,x) = (-1/x)^a U(a,1+a-b,-1/x)

   .. only:: texinfo

      .. math:: 2F0(a,b,x) = (-1/x)^a U(a,1+a-b,-1/x)

.. exceptions: GSL_EDOM
