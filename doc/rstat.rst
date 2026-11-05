.. index::
   single: running statistics
   single: online statistics

******************
Running Statistics
******************

This chapter describes routines for computing running statistics,
also known as online statistics, of data. These routines are
suitable for handling large datasets for which it may be
inconvenient or impractical to store in memory all at once.
The data can be processed in a single pass, one point at a time.
Each time a data point is added to the accumulator, internal
parameters are updated in order to compute the current mean, variance,
standard deviation, skewness, and kurtosis. These statistics are
exact, and are updated with numerically stable single-pass algorithms.
The median, arbitrary quantiles, and the median absolute deviation
are also available. These are computed from a summary of the data of
fixed size (see :ref:`sec_rstat-quantile-estimation`), and are exact
until the number of data exceeds a user-selectable accuracy parameter
:math:`k`, and approximations with a bounded rank error beyond.

The functions described in this chapter are declared in the header files
:file:`gsl_rstat.h` and :file:`gsl_wrstat.h` (weighted data).

Initializing the Accumulator
============================

.. type:: gsl_rstat_workspace

   This workspace contains parameters used to calculate various statistics
   and are updated after each data point is added to the accumulator.

.. function:: gsl_rstat_workspace * gsl_rstat_alloc (void)

   This function allocates a workspace for computing running statistics,
   with the default accuracy parameter :math:`k = 200` for the median,
   quantiles and median absolute deviation. It is equivalent to
   :code:`gsl_rstat_alloc2(200)`.

.. function:: gsl_rstat_workspace * gsl_rstat_alloc2 (const size_t k)

   This function allocates a workspace for computing running statistics,
   with accuracy parameter :data:`k` for the median, quantiles and median
   absolute deviation. These are exact as long as at most :data:`k` data have
   been added; beyond, their normalized rank error is :math:`O(1/k)`, about
   :math:`0.4/k` on average and :math:`1.3/k` at the 99th percentile
   (see :ref:`sec_rstat-quantile-estimation`). The size of the workspace
   is :math:`O(k)`, independent of the number of data; about :math:`3k`
   data are retained. The parameter :data:`k` must be at least :math:`8`.

.. function:: void gsl_rstat_free (gsl_rstat_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: int gsl_rstat_reset (gsl_rstat_workspace * w)

   This function resets the workspace :data:`w` to its initial state,
   so it can begin working on a new set of data.

.. function:: int gsl_rstat_merge (const gsl_rstat_workspace * src, gsl_rstat_workspace * dest)

   This function merges the data accumulated in :data:`src` into :data:`dest`,
   so that :data:`dest` holds the statistics of the data added to either
   workspace, as if all of them had been added to :data:`dest`. The workspace
   :data:`src` is unchanged. This allows data to be accumulated in parallel,
   for example in one workspace per thread, and combined at the end. The
   minimum, maximum, mean, variance, skewness and kurtosis are combined
   exactly, using the formulas of Pebay et al. (2016). The median, quantiles,
   median absolute deviation and cumulative distribution function have the
   same error bound after merging as for a single workspace. The workspaces
   must have the same accuracy parameter :math:`k`, and must be different.

Adding Data to the Accumulator
==============================

.. function:: int gsl_rstat_add (const double x, gsl_rstat_workspace * w)

   This function adds the data point :data:`x` to the statistical
   accumulator, updating calculations of the mean, variance,
   standard deviation, skewness, kurtosis, median, quantiles and median
   absolute deviation. A NaN value is counted, but is excluded from the
   median, quantiles and median absolute deviation, since it has no rank.

.. function:: size_t gsl_rstat_n (const gsl_rstat_workspace * w)

   This function returns the number of data so far added to the accumulator.

Current Statistics
==================

.. function:: double gsl_rstat_min (const gsl_rstat_workspace * w)

   This function returns the minimum value added to the accumulator.

.. function:: double gsl_rstat_max (const gsl_rstat_workspace * w)

   This function returns the maximum value added to the accumulator.

.. function:: double gsl_rstat_mean (const gsl_rstat_workspace * w)

   This function returns the mean of all data added to the accumulator,
   defined as

   .. only:: not texinfo

      .. math:: {\Hat\mu} = {1 \over N} \sum x_i

   .. only:: texinfo

      ::

         \Hat\mu = (1/N) \sum x_i

.. function:: double gsl_rstat_variance (const gsl_rstat_workspace * w)

   This function returns the variance of all data added to the accumulator,
   defined as

   .. only:: not texinfo

      .. math:: {\Hat\sigma}^2 = {1 \over (N-1)} \sum (x_i - {\Hat\mu})^2

   .. only:: texinfo

      ::

         \Hat\sigma^2 = (1/(N-1)) \sum (x_i - \Hat\mu)^2

.. function:: double gsl_rstat_sd (const gsl_rstat_workspace * w)

   This function returns the standard deviation of all data added to the
   accumulator, defined as the square root of the variance given above.

.. function:: double gsl_rstat_sd_mean (const gsl_rstat_workspace * w)

   This function returns the standard deviation of the mean, defined as

   .. only:: not texinfo
   
      .. math:: \Hat\sigma_{\Hat\mu} = {\Hat\sigma \over \sqrt{N}}

   .. only:: texinfo

      ::

         sd_mean = \Hat\sigma / \sqrt{N}

.. function:: double gsl_rstat_rms (const gsl_rstat_workspace * w)

   This function returns the root mean square of all data added to the
   accumulator, defined as

   .. math:: rms = \sqrt{{1 \over N} \sum x_i^2}

.. function:: double gsl_rstat_skew (const gsl_rstat_workspace * w)

   This function returns the skewness of all data added to the accumulator,
   defined as

   .. only:: not texinfo

      .. math::

         skew = {1 \over N} \sum 
          {\left( x_i - {\Hat\mu} \over {\Hat\sigma} \right)}^3

   .. only:: texinfo

      ::

         skew = (1/N) \sum ((x_i - \Hat\mu)/\Hat\sigma)^3

.. function:: double gsl_rstat_kurtosis (const gsl_rstat_workspace * w)

   This function returns the kurtosis of all data added to the accumulator,
   defined as

   .. only:: not texinfo

      .. math::

         kurtosis = \left( {1 \over N} \sum 
          {\left(x_i - {\Hat\mu} \over {\Hat\sigma} \right)}^4 
          \right) 
          - 3

   .. only:: texinfo

      ::

         kurtosis = ((1/N) \sum ((x_i - \Hat\mu)/\Hat\sigma)^4)  - 3

.. function:: double gsl_rstat_median (const gsl_rstat_workspace * w)

   This function returns an estimate of the median of the data added to
   the accumulator. For an even number of data, the median is the average
   of the two middle values, as in :func:`gsl_stats_median`. The result is
   exact as long as at most :math:`k` data have been added.

.. function:: double gsl_rstat_quantile (const double p, const gsl_rstat_workspace * w)

   This function returns an estimate of the :data:`p`-quantile of the data
   added to the accumulator, where :data:`p` is between :math:`0` and :math:`1`.
   The quantile is interpolated linearly between order statistics, as in
   :func:`gsl_stats_quantile_from_sorted_data`. The values :math:`p = 0` and
   :math:`p = 1` give the exact minimum and maximum.

.. function:: double gsl_rstat_mad0 (const gsl_rstat_workspace * w)
              double gsl_rstat_mad (const gsl_rstat_workspace * w)

   These functions return estimates of the median absolute deviation
   (MAD) of the data added to the accumulator. The function
   :func:`gsl_rstat_mad0` returns the unscaled MAD,

   .. math:: MAD_0 = \textrm{median} \left\{ \left| x_i - \textrm{median}(x) \right| \right\}

   and :func:`gsl_rstat_mad` returns the scaled MAD, :math:`1.4826 \times MAD_0`,
   which is a consistent estimate of the standard deviation for Gaussian
   data, as in :func:`gsl_stats_mad0` and :func:`gsl_stats_mad`.

.. function:: double gsl_rstat_cdf (const double x, const gsl_rstat_workspace * w)

   This function returns an estimate of the empirical cumulative distribution
   function at :data:`x`, the fraction of the data added to the accumulator
   which are less than or equal to :data:`x`. This is the inverse of
   :func:`gsl_rstat_quantile`, and gives for example the percentile rank of
   :data:`x`. The result is exact as long as at most :math:`k` data have been
   added, and has an error of :math:`O(1/k)` beyond.

.. function:: double gsl_rstat_norm (const gsl_rstat_workspace * w)

   This function returns the Euclidean norm of all data added to the
   accumulator, defined as

   .. math:: \left|\left| x \right|\right| = \sqrt{\sum x_i^2}

Weighted Running Statistics
===========================

The functions in this section compute running statistics of weighted
data :math:`x_i` with weights :math:`w_i \ge 0`. They are declared in the
header file :file:`gsl_wrstat.h`.

.. type:: gsl_wrstat_workspace

   This workspace contains parameters used to calculate weighted statistics,
   which are updated after each data point is added to the accumulator.

.. function:: gsl_wrstat_workspace * gsl_wrstat_alloc (void)
              gsl_wrstat_workspace * gsl_wrstat_alloc2 (const size_t k)

   These functions allocate a workspace for computing weighted running
   statistics, with accuracy parameter :data:`k` for the weighted median,
   quantiles and median absolute deviation, as for :func:`gsl_rstat_alloc2`.
   The function :func:`gsl_wrstat_alloc` uses :math:`k = 200`.

.. function:: void gsl_wrstat_free (gsl_wrstat_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: int gsl_wrstat_reset (gsl_wrstat_workspace * w)

   This function resets the workspace :data:`w` to its initial state,
   so it can begin working on a new set of data.

.. function:: int gsl_wrstat_merge (const gsl_wrstat_workspace * src, gsl_wrstat_workspace * dest)

   This function merges the weighted data accumulated in :data:`src` into
   :data:`dest`, as described for :func:`gsl_rstat_merge`.

.. function:: int gsl_wrstat_add (const double wt, const double x, gsl_wrstat_workspace * w)

   This function adds the data point :data:`x` with weight :data:`wt` to the
   accumulator. The weight must be non-negative. Data with zero weight are
   counted, and included in the minimum and maximum, but do not affect the
   other statistics. Data with infinite or NaN weight, and NaN data, are
   excluded from the median, quantiles and median absolute deviation.

.. function:: size_t gsl_wrstat_n (const gsl_wrstat_workspace * w)

   This function returns the number of data so far added to the accumulator.

.. function:: double gsl_wrstat_min (const gsl_wrstat_workspace * w)
              double gsl_wrstat_max (const gsl_wrstat_workspace * w)

   These functions return the minimum and maximum value added to the accumulator.

.. function:: double gsl_wrstat_mean (const gsl_wrstat_workspace * w)

   This function returns the weighted mean,
   :math:`{\Hat\mu} = \sum w_i x_i / \sum w_i`.

.. function:: double gsl_wrstat_variance (const gsl_wrstat_workspace * w)
              double gsl_wrstat_sd (const gsl_wrstat_workspace * w)

   These functions return the weighted variance, as in :func:`gsl_stats_wvariance`,

   .. math:: {\Hat\sigma}^2 = {{\sum w_i} \over {(\sum w_i)^2 - \sum w_i^2}} \sum w_i (x_i - {\Hat\mu})^2

   and the weighted standard deviation, its square root. These are zero if
   fewer than two data with positive weight have been added.

.. function:: double gsl_wrstat_rms (const gsl_wrstat_workspace * w)

   This function returns the weighted root mean square,
   :math:`\sqrt{\sum w_i x_i^2 / \sum w_i}`.

.. function:: double gsl_wrstat_norm (const gsl_wrstat_workspace * w)

   This function returns the weighted Euclidean norm,
   :math:`\sqrt{\sum w_i x_i^2}`.

.. function:: double gsl_wrstat_neff (const gsl_wrstat_workspace * w)

   This function returns Kish's effective sample size (Kish, 1965),

   .. math:: n_{eff} = {{(\sum w_i)^2} \over {\sum w_i^2}}

   which equals the number of data with positive weight if all these
   weights are equal, and is smaller otherwise.

.. function:: double gsl_wrstat_sd_mean (const gsl_wrstat_workspace * w)

   This function returns the standard deviation of the weighted mean,
   :math:`{\Hat\sigma} / \sqrt{n_{eff}}`, which reduces to
   :func:`gsl_rstat_sd_mean` for equal weights.

.. function:: double gsl_wrstat_skew (const gsl_wrstat_workspace * w)
              double gsl_wrstat_kurtosis (const gsl_wrstat_workspace * w)

   These functions return the weighted skewness and kurtosis, as in
   :func:`gsl_stats_wskew` and :func:`gsl_stats_wkurtosis`,

   .. only:: not texinfo

      .. math::

         skew = {{\sum w_i (x_i - {\Hat\mu})^3} \over {{\Hat\sigma}^3 \sum w_i}}, \qquad
         kurtosis = {{\sum w_i (x_i - {\Hat\mu})^4} \over {{\Hat\sigma}^4 \sum w_i}} - 3

   .. only:: texinfo

      ::

         skew = (\sum w_i (x_i - \Hat\mu)^3) / (\Hat\sigma^3 \sum w_i)
         kurtosis = (\sum w_i (x_i - \Hat\mu)^4) / (\Hat\sigma^4 \sum w_i) - 3

   The central moments are updated with the numerically stable formulas
   of Pebay et al. (2016).

.. function:: double gsl_wrstat_median (const gsl_wrstat_workspace * w)
              double gsl_wrstat_quantile (const double p, const gsl_wrstat_workspace * w)

   These functions return estimates of the weighted median and the weighted
   :data:`p`-quantile, where :data:`p` is between :math:`0` and :math:`1`.
   With :math:`C(x)` the total weight of the data :math:`\le x` and :math:`W`
   the total weight, the weighted :math:`p`-quantile is defined as the average
   of the smallest data value :math:`x` with :math:`C(x) \ge pW` and the smallest
   data value with :math:`C(x) > pW`. The weighted median is the
   :math:`1/2`-quantile; for unit weights it equals the median of
   :func:`gsl_stats_median`.

.. function:: double gsl_wrstat_mad0 (const gsl_wrstat_workspace * w)
              double gsl_wrstat_mad (const gsl_wrstat_workspace * w)

   These functions return estimates of the weighted median absolute
   deviation, the weighted median of :math:`| x_i - \textrm{median} |` with
   weights :math:`w_i`, unscaled (:func:`gsl_wrstat_mad0`) and multiplied by
   :math:`1.4826` (:func:`gsl_wrstat_mad`).

.. function:: double gsl_wrstat_cdf (const double x, const gsl_wrstat_workspace * w)

   This function returns an estimate of the weighted empirical cumulative
   distribution function at :data:`x`, the fraction of the total weight of
   the data which are less than or equal to :data:`x`.

.. _sec_rstat-quantile-estimation:

Quantile Estimation
===================

The median, quantiles, median absolute deviation and cumulative distribution
function are computed with the KLL algorithm of Karnin, Lang and Liberty
(2016), which keeps a summary of about :math:`3k` data values, regardless of
how many data are added. The results are exact as long as at most :math:`k`
data have been added. Beyond that, the error is bounded in terms of rank:
the fraction of data below an estimated :math:`p`-quantile differs from
:math:`p` by :math:`O(1/k)`, whatever the order in which the data are added.
The default :math:`k = 200` gives a rank error below one percent. The
algorithm makes random choices internally, using a fixed seed, so results
are reproducible.

Single Quantile Estimation
==========================

The functions in this section estimate a single quantile dynamically
using the :math:`P^2` algorithm of Jain and Chlamtac, 1985, with
workspace of only a few bytes. Five points (markers) are stored which
represent the minimum and maximum of the data, as well as current
estimates of the :math:`p/2`-, :math:`p`-, and :math:`(1+p)/2`-quantiles.
Each time a new data point is added, the marker positions and heights are
updated by piecewise-parabolic interpolation. The :math:`P^2` algorithm is
a heuristic without error bounds; it is accurate for data which arrive in
random order, but can be inaccurate for sorted, drifting or otherwise
correlated data. For such data, :func:`gsl_rstat_quantile` should be used.

.. type:: gsl_rstat_quantile_workspace

   This workspace contains parameters for estimating quantiles of the
   current dataset

.. function:: gsl_rstat_quantile_workspace * gsl_rstat_quantile_alloc (const double p)

   This function allocates a workspace for the dynamic estimation of
   :data:`p`-quantiles, where :data:`p` is between :math:`0` and :math:`1`.
   The median corresponds to :math:`p = 0.5`. The size of the workspace
   is :math:`O(1)`.

.. function:: void gsl_rstat_quantile_free (gsl_rstat_quantile_workspace * w)

   This function frees the memory associated with the workspace :data:`w`.

.. function:: int gsl_rstat_quantile_reset (gsl_rstat_quantile_workspace * w)

   This function resets the workspace :data:`w` to its initial state,
   so it can begin working on a new set of data.

.. function:: int gsl_rstat_quantile_add (const double x, gsl_rstat_quantile_workspace * w)

   This function updates the estimate of the :math:`p`-quantile with
   the new data point :data:`x`.

.. function:: double gsl_rstat_quantile_get (gsl_rstat_quantile_workspace * w)

   This function returns the current estimate of the :math:`p`-quantile.

Examples
========

Here is a basic example of how to use the statistical functions:

.. include:: examples/rstat.c
   :code:

The program should produce the following output,

.. include:: examples/rstat.txt
   :code:

This next program estimates the lower quartile, median and upper
quartile of 10,000 samples of a random Rayleigh distribution with
:func:`gsl_rstat_quantile`, using the default accuracy parameter
:math:`k = 200`. The same data are added to one accumulator in the
order generated, and to a second accumulator in ascending order, to
show that the estimates do not depend on the order of the data.
The fraction of data less than or equal to 1 is also estimated with
:func:`gsl_rstat_cdf`. For comparison, the exact values are computed
from the sorted dataset.

.. include:: examples/rquantile.c
   :code:

The program should produce the following output,

.. include:: examples/rquantile.txt
   :code:

References and Further Reading
==============================

The KLL algorithm used for the median, quantiles and median absolute
deviation is described in the papers,

* Z. Karnin, K. Lang and E. Liberty.
  *Optimal Quantile Approximation in Streams*,
  IEEE 57th Annual Symposium on Foundations of Computer Science (FOCS), 2016,
  p. 71-78.

* N. Ivkin, E. Liberty, K. Lang, Z. Karnin and V. Braverman.
  *Streaming Quantiles Algorithms with Small Space and Update Time*,
  arXiv:1907.00236, 2019.

The :math:`P^2` algorithm used to dynamically estimate a single
:math:`p`-quantile is described in the paper,

* R. Jain and I. Chlamtac.
  *The P^2 algorithm for dynamic calculation of quantiles and histograms without storing observations*,
  Communications of the ACM, Volume 28 (October), Number 10, 1985,
  p. 1076-1085.

The algorithms used for the weighted mean and central moments, and for
merging accumulators, are described in,

* P. Pebay, T. Terriberry, H. Kolla and J. Bennett.
  *Numerically stable, scalable formulas for parallel and online computation
  of higher-order multivariate central moments with arbitrary weights*,
  Computational Statistics, Volume 31, 2016, p. 1305-1325.

The effective sample size of weighted data is described in,

* L. Kish.
  *Survey Sampling*,
  John Wiley & Sons, New York, 1965.
