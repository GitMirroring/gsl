/* linalg/lu_pivots.h
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

/*
 * The LU factorizations record their row interchanges as a pivot
 * sequence ipiv[0..npiv-1] (row i was swapped with row ipiv[i] >= i),
 * stored in the data array of the output permutation itself, so that
 * no workspace is needed.  lu_pivots_to_permutation() converts this
 * sequence in place into the permutation p of P A = L U.
 *
 * With tau_i the transposition (i ipiv[i]), the permutation is
 * P = tau_0 tau_1 ... tau_{npiv-1}, and swapping entries i and ipiv[i]
 * of an array q gives q tau_i.  Applied in increasing i this would
 * build P, but entry ipiv[i] > i still holds an unused pivot.  Applied
 * in decreasing i it builds Q = tau_{npiv-1} ... tau_0 = P^{-1} in
 * place: after step i, entries > i hold Q restricted to them, and the
 * entries < i, which Q_{i+1} = tau_{npiv-1} ... tau_{i+1} leaves fixed,
 * still hold their pivots.  P is then obtained by inverting Q in place,
 * reversing each cycle and using the top bit of size_t to mark the
 * entries done (Knuth, TAOCP vol. 1, 1.3.3, Algorithm I): a size_t
 * array cannot have 2^(bits-1) elements, so that bit is never part of
 * an index.  Both passes are O(size).
 */

#include <limits.h>
#include <gsl/gsl_permutation.h>

static void
lu_pivots_to_permutation (gsl_permutation * p, const size_t npiv, int * signum)
{
  const size_t M = p->size;
  const size_t mark = ((size_t) 1) << (sizeof (size_t) * CHAR_BIT - 1);
  size_t * d = p->data;
  size_t i;

  *signum = 1;

  for (i = npiv; i < M; ++i)
    d[i] = i;

  /* Q = P^{-1}, from the last interchange to the first */
  for (i = npiv; i-- > 0; )
    {
      const size_t j = d[i];

      if (j != i)
        {
          d[i] = d[j];
          d[j] = i;
          *signum = -(*signum);
        }
    }

  /* P = Q^{-1} */
  for (i = 0; i < M; ++i)
    {
      size_t prev, cur;

      if (d[i] & mark)
        continue;

      prev = i;
      cur = d[i];
      while (cur != i)
        {
          const size_t next = d[cur];
          d[cur] = prev | mark;
          prev = cur;
          cur = next;
        }
      d[i] = prev | mark;
    }

  for (i = 0; i < M; ++i)
    d[i] &= ~mark;
}
