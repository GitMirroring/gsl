.. index::
   single: associated Legendre functions
   single: Legendre functions, associated
   single: spherical harmonics

The associated Legendre functions are described in
Abramowitz & Stegun, Chapter 8.  These functions are declared in 
the header file :file:`gsl_sf_alf.h`.

The associated Legendre functions :math:`P_l^m(x)` are solutions
of the differential equation

.. only:: not texinfo

   .. math:: (1 - x^2) {d^2 \over dx^2} P_l^m(x) - 2x {d \over dx} P_l^m(x) +
             \left( l(l+1) - {m^2 \over 1 - x^2} \right) P_l^m(x) = 0

.. only:: texinfo

   ::

      (1 - x^2) d^2/dx^2 P_l^m(x) - 2x d/dx P_l^m(x) +
      ( l(l+1) - m^2 / (1 - x^2) ) P_l^m(x) = 0

where the degree :math:`l` and order :math:`m` satisfy :math:`0 \le l` and
:math:`0 \le m \le l`.
The functions :math:`P_l^m(x)` grow combinatorially with
:math:`l` and can overflow in double precision for :math:`l` larger than about 150.
The routines below therefore compute normalized associated Legendre functions

.. math:: \bar{P}_l^m(x) = K_l^m P_l^m(x)

where :math:`P_l^m(x)` is defined without the Condon-Shortley phase factor
:math:`(-1)^m`; the factor can be included by setting the flag
:macro:`GSL_SF_ALF_FLG_CSPHASE` in :func:`gsl_sf_alf_precompute`.

The normalizations are those used in spherical harmonic expansions, either
with complex spherical harmonics

.. math:: \bar{P}_l^m(\cos{\theta}) e^{i m \phi}

or with real spherical harmonics

.. math:: \bar{P}_l^m(\cos{\theta}) \cos{m \phi}, \quad \bar{P}_l^m(\cos{\theta}) \sin{m \phi}

The two require different factors for :math:`m > 0`: the integral of
:math:`\cos^2{m \phi}` or :math:`\sin^2{m \phi}` over :math:`[0, 2\pi]` is
:math:`\pi` for :math:`m > 0`, but the integral of :math:`|e^{i m \phi}|^2` is
:math:`2 \pi` (and for :math:`m = 0` there is the single real function, with
integral :math:`2 \pi`). Equivalently, :math:`\cos{m\phi}` and
:math:`\sin{m\phi}` combine the complex functions of orders :math:`m` and
:math:`-m`. Normalizations for real spherical harmonics therefore contain the
factor

.. math:: c_m = 2 - \delta_{m0}

under the square root, which is :math:`c_m = 1` for complex spherical harmonics.
In the functions below, a parameter of type :type:`gsl_sf_alf_t` specifies the
normalization:

.. type:: gsl_sf_alf_t

   .. list-table::
      :header-rows: 1

      * - Value
        - :math:`K_l^m`
        - :math:`\int_{-1}^1 \bar{P}_k^m(x) \bar{P}_l^m(x) dx`
      * - :code:`GSL_SF_ALF_SCHMIDT_REAL`
        - :math:`\sqrt{c_m \frac{(l-m)!}{(l+m)!}}`, :math:`c_m = 2 - \delta_{m0}`
        - :math:`\frac{2 c_m}{2l+1} \delta_{kl}`
      * - :code:`GSL_SF_ALF_SCHMIDT_COMPLEX`
        - :math:`\sqrt{\frac{(l-m)!}{(l+m)!}}`
        - :math:`\frac{2}{2l+1} \delta_{kl}`
      * - :code:`GSL_SF_ALF_FOURPI_REAL`
        - :math:`\sqrt{c_m (2l+1) \frac{(l-m)!}{(l+m)!}}`, :math:`c_m = 2 - \delta_{m0}`
        - :math:`2 c_m \delta_{kl}`
      * - :code:`GSL_SF_ALF_FOURPI_COMPLEX`
        - :math:`\sqrt{(2l+1) \frac{(l-m)!}{(l+m)!}}`
        - :math:`2 \delta_{kl}`
      * - :code:`GSL_SF_ALF_ORTHO_REAL`
        - :math:`\sqrt{c_m \frac{2l+1}{4\pi} \frac{(l-m)!}{(l+m)!}}`, :math:`c_m = 2 - \delta_{m0}`
        - :math:`\frac{c_m}{2\pi} \delta_{kl}`
      * - :code:`GSL_SF_ALF_ORTHO_COMPLEX`
        - :math:`\sqrt{\frac{2l+1}{4\pi} \frac{(l-m)!}{(l+m)!}}`
        - :math:`\frac{1}{2\pi} \delta_{kl}`
      * - :code:`GSL_SF_ALF_FULL`
        - :math:`\sqrt{(l + \frac{1}{2}) \frac{(l-m)!}{(l+m)!}}`
        - :math:`\delta_{kl}`
      * - :code:`GSL_SF_ALF_NONE`
        - :math:`1`
        - :math:`\frac{2}{2l+1} \frac{(l+m)!}{(l-m)!} \delta_{kl}`

With the factor :math:`\pi c_m` from the integral over :math:`\phi` for real,
and :math:`2 \pi` for complex spherical harmonics, the spherical harmonics have
the following squared norms over the unit sphere, independent of :math:`m`:

* :code:`GSL_SF_ALF_SCHMIDT_REAL`, :code:`GSL_SF_ALF_SCHMIDT_COMPLEX`:
  :math:`4 \pi / (2l+1)` (Schmidt semi-normalization, common in geomagnetism)

* :code:`GSL_SF_ALF_FOURPI_REAL`, :code:`GSL_SF_ALF_FOURPI_COMPLEX`:
  :math:`4 \pi` (:math:`4 \pi` normalization, common in geodesy; see
  Hofmann-Wellenhof and Moritz, 2006)

* :code:`GSL_SF_ALF_ORTHO_REAL`, :code:`GSL_SF_ALF_ORTHO_COMPLEX`:
  :math:`1` (orthonormal spherical harmonics)

:code:`GSL_SF_ALF_FULL` gives functions which are orthonormal on :math:`[-1,1]`,
and :code:`GSL_SF_ALF_NONE` the unnormalized functions :math:`P_l^m(x)`.

The normalized associated Legendre routines below use the recurrence
relations of Holmes and Featherstone (2002). The sectoral terms
:math:`P_m^m(x)`, which contain the factor :math:`\sin^m{\theta}`, are carried
with a separate power-of-two exponent, so they do not underflow at high
order; values below the smallest normal double precision number are
returned as zero. Near the poles, the relative rounding error of the
recurrence grows roughly like :math:`\epsilon L^2` for maximum degree :math:`L`.
Routines are provided to compute
first and second derivatives
:math:`dP_l^m(x)/dx` and :math:`d^2 P_l^m(x)/dx^2` as well as their alternate
versions :math:`d P_l^m(\cos{\theta})/d\theta` and
:math:`d^2 P_l^m(\cos{\theta})/d\theta^2`. While there is a simple
scaling relationship between the two forms, the derivatives
involving :math:`\theta` are heavily used in spherical harmonic
expansions, and also do not suffer from singularities at the poles,
:math:`x = \pm 1`, and so these routines are also provided.

The routines below which return an array of ALFs organize the outputs into blocks
of the same order :math:`m`.  This scheme uses the following indexing function to
locate an ALF of degree :math:`l` and order :math:`m`,

.. math:: \mathcal{I}_m(l,m,L) = m L - \frac{m(m-1)}{2} + l

where :math:`L \geq 0` is the maximum degree, :code:`lmax`. This corresponds to the following memory layout,

.. math:: l \quad \overbrace{0 \; 1 \; 2 \; \cdots \; L}^{m = 0} \quad \overbrace{1 \; 2 \; \cdots \; L}^{m = 1} \quad \overbrace{2 \; 3 \; \cdots \; L}^{m = 2} \quad \cdots \quad \overbrace{M \; \cdots \; L}^{m = M}

Here, :math:`M` is the maximum order, :code:`mmax`, which satisfies :math:`0 \leq M \leq L`.
This memory layout corresponds to the order that ALFs are computed in their
recurrence relations, and is therefore cache-efficient. The following
code demonstrates how to access the array elements in order,

.. code::

   idx = 0;
   for (m = 0; m <= mmax; ++m) {
     for (l = m; l <= lmax; ++l) {
       double value = Plm[idx]; /* (l,m) element */
       ++idx;
     }
   }

.. function:: int gsl_sf_alf_precompute (const gsl_sf_alf_t norm, const size_t lmax, const size_t mmax, const size_t flags, double result_array[])

   This function precomputes the multiplicative factors needed for the
   associated Legendre recurrence relations. The input :data:`norm`
   specifies the ALF normalization. The input :data:`lmax` specifies
   the maximum ALF degree. The input :data:`mmax` specifies the maximum
   ALF order. The input :data:`flags` is a bitmask which
   specifies how the ALFs are computed and stored. It can contain the
   following values,

   .. macro:: GSL_SF_ALF_FLG_CSPHASE

      This flag will include the Condon-Shortley phase factor in the calculation
      of the ALFs

   The output array :data:`result_array` should have a length
   as returned by the function :func:`gsl_sf_alf_array_size`. The computed
   recurrence factors are stored at the end of :data:`result_array`, leaving
   room at the front for the calculation of the ALFs.

   This routine must be called prior to calling ALF array functions. The
   normalization and phase convention chosen here apply to all values the
   array functions subsequently compute in :data:`result_array`. The
   array functions check that :data:`result_array` was initialized by this
   routine with the same :data:`lmax` and :data:`mmax`, and return
   :macro:`GSL_EINVAL` otherwise. The array functions overwrite only the
   first :code:`gsl_sf_alf_nlm(lmax, mmax)` elements of :data:`result_array`
   (the ALF values), so the precomputed factors can be reused for any
   number of points :math:`x`.

   All array functions below return :macro:`GSL_EDOM` if :math:`x` is outside
   :math:`[-1,1]` or is a NaN.

.. function:: int gsl_sf_alf_array (const size_t lmax, const size_t mmax, const double x, double result_array[])

   This function calculates all associated Legendre polynomials
   for :math:`0 \le l \le lmax` and :math:`0 \le m \le \min{(l,mmax)}` for :math:`|x| \le 1`,
   in the normalization selected in :func:`gsl_sf_alf_precompute`.
   The :math:`P_l^m(x)` values are stored in :data:`result_array`, whose
   minimum size can be obtained from calling :func:`gsl_sf_alf_array_size`.
   The array index of :math:`P_l^m(x)` is obtained from calling
   :code:`gsl_sf_alf_array_index(l, m, lmax)`.

   The function :func:`gsl_sf_alf_precompute` must be called first
   using the same :data:`lmax` and :data:`mmax` inputs
   to initialize :data:`result_array` with the multiplicative factors
   used in the recurrence relations.

.. function:: int gsl_sf_alf_deriv_array (const size_t lmax, const size_t mmax, const double x, double result_array[], double result_deriv_array[])

   This function calculates all associated Legendre
   functions and their first derivatives for :math:`0 \leq l \leq lmax`
   and :math:`0 \leq m \leq \min{(l,mmax)}` for :math:`|x| < 1`.
   The :math:`P_l^m(x)` values and their derivatives
   :math:`dP_l^m(x)/dx` are stored in :data:`result_array` and
   :data:`result_deriv_array` respectively. The array :data:`result_deriv_array`
   must have length at least :code:`gsl_sf_alf_nlm(lmax, mmax)`.

   Note that for some orders :math:`m`, the derivatives :math:`dP_l^m(x)/dx` have
   singularities at the end points :math:`x = \pm 1`, and so this function only
   accepts interior points as input, :math:`x \in (-1,1)`.

   The function :func:`gsl_sf_alf_precompute` must be called first
   using the same :data:`lmax` and :data:`mmax` inputs
   to initialize :data:`result_array` with the multiplicative factors
   used in the recurrence relations.

.. function:: int gsl_sf_alf_deriv2_array (const size_t lmax, const size_t mmax, const double x, double result_array[], double result_deriv_array[], double result_deriv2_array[])

   This function calculates all associated Legendre
   functions, their first derivatives, and their second
   derivatives for :math:`0 \leq l \leq lmax`
   and :math:`0 \leq m \leq \min{(l,mmax)}` for :math:`|x| < 1`.
   The :math:`P_l^m(x)` values, their derivatives
   :math:`dP_l^m(x)/dx`, and second derivatives
   :math:`d^2 P_l^m(x) / dx^2` are stored in :data:`result_array`,
   :data:`result_deriv_array`, and :data:`result_deriv2_array` respectively.
   The arrays :data:`result_deriv_array` and :data:`result_deriv2_array`
   must have length at least :code:`gsl_sf_alf_nlm(lmax, mmax)`.

   Note that for some orders :math:`m`, the derivatives have
   singularities at the end points :math:`x = \pm 1`, and so this function only
   accepts interior points as input, :math:`x \in (-1,1)`.

   The function :func:`gsl_sf_alf_precompute` must be called first
   using the same :data:`lmax` and :data:`mmax` inputs
   to initialize :data:`result_array` with the multiplicative factors
   used in the recurrence relations.

.. function:: int gsl_sf_alf_vsh_array (const size_t lmax, const size_t mmax, const double x, double result_array[], double result_deriv_array[])

   In vector spherical harmonic expansions, it is often necessary to simultaneously
   compute terms of the form

   .. math:: \frac{d}{d\theta} P_l^m(\cos{\theta}) e^{i m \phi}

   and

   .. math:: \frac{i m}{\sin{\theta}} P_l^m(\cos{\theta}) e^{i m \phi}

   Analogous expressions to those above also arise when using real-valued spherical
   harmonics. The factor :math:`\frac{1}{\sin{\theta}}` in the second term could
   present problems at the poles, :math:`\theta = 0,\pi` (e.g. :math:`x = \pm 1`). However,
   for :math:`m \neq 0`, :math:`P_l^m(\cos{\theta}) / \sin{\theta}` is well defined
   and can be computed stably for all :math:`\theta`.

   This function computes and stores the following quantities in the output arrays,

   .. math::

      \textrm{result\_array[index(l,0,lmax)]} &= P_l^0(x) \\
      \textrm{result\_array[index(l,m,lmax)]} &= \frac{P_l^m(x)}{\sin{\theta}}, \qquad m > 0 \\
      \textrm{result\_deriv\_array[index(l,m,lmax)]} &= \frac{d}{d\theta} P_l^m(x), \qquad m \geq 0
   
   The :code:`index` notation above refers to the :func:`gsl_sf_alf_array_index` function.
   All associated Legendre functions and their first :math:`\theta` derivatives are computed
   for :math:`0 \leq l \leq lmax` and :math:`0 \leq m \leq \min{(l,mmax)}` for :math:`|x| \leq 1`.
   The array :data:`result_deriv_array` must have length at least :code:`gsl_sf_alf_nlm(lmax, mmax)`.

   The function :func:`gsl_sf_alf_precompute` must be called first
   using the same :data:`lmax` and :data:`mmax` inputs
   to initialize :data:`result_array` with the multiplicative factors
   used in the recurrence relations.

.. function:: int gsl_sf_alf_vsh2_array (const size_t lmax, const size_t mmax, const double x, double result_array[], double result_deriv_array[], double result_deriv2_array[])

   This function computes the same quantities as :func:`gsl_sf_alf_vsh_array`, stored in
   :data:`result_array` and :data:`result_deriv_array` with the same conventions, together with
   the second derivatives

   .. math:: \textrm{result\_deriv2\_array[index(l,m,lmax)]} = \frac{d^2}{d\theta^2} P_l^m(x), \qquad m \geq 0

   for :math:`0 \leq l \leq lmax` and :math:`0 \leq m \leq \min{(l,mmax)}` for :math:`|x| \leq 1`.
   The arrays :data:`result_deriv_array` and :data:`result_deriv2_array` must have length
   at least :code:`gsl_sf_alf_nlm(lmax, mmax)`.
   The derivatives are obtained by differentiating the recurrence relations with respect to
   :math:`\theta`, which involves only multiplications by :math:`\sin{\theta}` and :math:`\cos{\theta}`,
   so the results are well defined at the poles, and any :math:`mmax \leq lmax` may be used.

   The function :func:`gsl_sf_alf_precompute` must be called first
   using the same :data:`lmax` and :data:`mmax` inputs
   to initialize :data:`result_array` with the multiplicative factors
   used in the recurrence relations.

.. function:: int gsl_sf_alf_theta_derivk_array (const size_t lmax, const double Plm[], const double input_array[], double output_array[])

   This function computes derivatives (with respect to :math:`\theta`) of order :math:`k`,
   given derivatives of order :math:`k-1`, using the algorithm of Bosch (2000). The
   derivatives :math:`d^{k-1} P_l^m/d\theta^{k-1}` are provided in the array :data:`input_array`,
   while the derivatives :math:`d^k P_l^m/d\theta^k` are stored in the array :data:`output_array`;
   both have length :code:`gsl_sf_alf_nlm(lmax, lmax)` and use the indexing of
   :func:`gsl_sf_alf_array_index`.
   The array :data:`Plm` provides the recurrence factors: it must have been initialized by
   :func:`gsl_sf_alf_precompute`, which also selects the normalization and phase convention.
   The maximum degree :math:`L` is given in :data:`lmax`,
   and for this algorithm to work, it is required that :math:`L = M`: :data:`Plm` must be
   initialized with :code:`mmax = lmax`, otherwise :macro:`GSL_EINVAL` is returned.

   For the first derivatives (:math:`k = 1`), :data:`input_array` must hold the values
   :math:`P_l^m(x)` themselves, as computed by :func:`gsl_sf_alf_array` (not
   :math:`P_l^m(x)/\sin{\theta}` as stored by :func:`gsl_sf_alf_vsh_array`); since
   :func:`gsl_sf_alf_array` stores them at the front of :data:`Plm`, the same array may be
   passed as :data:`Plm` and :data:`input_array`. For :math:`k = 2`, :data:`input_array`
   may be the :data:`result_deriv_array` of :func:`gsl_sf_alf_vsh_array`. Repeated calls give
   derivatives of any order.

.. function:: size_t gsl_sf_alf_nlm (const size_t lmax, const size_t mmax)

   This function returns the total number of associated Legendre
   functions :math:`P_l^m(x)` for a given :data:`lmax` and :data:`mmax`.

   An inline version of this function is used if :macro:`HAVE_INLINE` is
   defined.

.. function:: size_t gsl_sf_alf_array_size (const size_t lmax, const size_t mmax)

   This function returns the minimum array size for maximum degree :data:`lmax`
   and maximum order :data:`mmax` needed for the array versions of the associated Legendre functions.
   Size is calculated as the total number of :math:`P_l^m(x)` functions
   (see :func:`gsl_sf_alf_nlm`),
   plus extra space for precomputing multiplicative factors used in the
   recurrence relations.

.. function:: size_t gsl_sf_alf_array_index (const size_t l, const size_t m, const size_t lmax)

   This function returns the index into the ALF arrays corresponding
   to the function of degree :math:`l` and order :math:`m`. The maximum
   degree :math:`L` must also be specified in the parameter :data:`lmax`.

   An inline version of this function is used if :macro:`HAVE_INLINE` is
   defined.

Single Associated Legendre Functions
------------------------------------

The following functions compute a single associated Legendre polynomial
:math:`P_l^m(x)`, unnormalized or with the spherical harmonic
normalization, without precomputation.  They are declared in the header
file :file:`gsl_sf_legendre.h`.  To compute all degrees and orders up to
a maximum, in any normalization, or their derivatives, the array
functions above are more efficient.

.. function:: double gsl_sf_legendre_Plm (int l, int m, double x)
              int gsl_sf_legendre_Plm_e (int l, int m, double x, gsl_sf_result * result)

   These routines compute the associated Legendre polynomial
   :math:`P_l^m(x)` for :math:`m \ge 0`,
   :math:`l \ge m`, and :math:`|x| \le 1`.
.. Exceptional Return Values: GSL_EDOM, GSL_EOVRFLW

.. function:: double gsl_sf_legendre_sphPlm (int l, int m, double x)
              int gsl_sf_legendre_sphPlm_e (int l, int m, double x, gsl_sf_result * result)

   These routines compute the normalized associated Legendre polynomial
   :math:`\sqrt{(2l+1)/(4\pi)} \sqrt{(l-m)!/(l+m)!} P_l^m(x)` suitable
   for use in spherical harmonics.  The parameters must satisfy :math:`m \ge 0`,
   :math:`l \ge m`, and :math:`|x| \le 1`.
   These routines avoid the overflows
   that occur for the standard normalization of :math:`P_l^m(x)`.
.. Exceptional Return Values: GSL_EDOM
