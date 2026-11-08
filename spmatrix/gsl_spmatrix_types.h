/* spmatrix/gsl_spmatrix_types.h
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

/* definitions shared by the gsl_spmatrix_<type>.h headers */

#ifndef __GSL_SPMATRIX_TYPES_H__
#define __GSL_SPMATRIX_TYPES_H__

#include <stdlib.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

enum
{
  GSL_SPMATRIX_COO = 0, /* coordinate/triplet representation */
  GSL_SPMATRIX_CSC = 1, /* compressed sparse column */
  GSL_SPMATRIX_CSR = 2, /* compressed sparse row */
  GSL_SPMATRIX_TRIPLET = GSL_SPMATRIX_COO,
  GSL_SPMATRIX_CCS = GSL_SPMATRIX_CSC,
  GSL_SPMATRIX_CRS = GSL_SPMATRIX_CSR
};

/* memory pool for binary tree node allocation */
struct gsl_spmatrix_pool_node
{
  struct gsl_spmatrix_pool_node * next;
  void * block_ptr;          /* pointer to memory block, of size n*tree_node_size */
  unsigned char * free_slot; /* pointer to next available slot */
};

typedef struct gsl_spmatrix_pool_node gsl_spmatrix_pool;

#define GSL_SPMATRIX_ISCOO(m)         ((m)->sptype == GSL_SPMATRIX_COO)
#define GSL_SPMATRIX_ISCSC(m)         ((m)->sptype == GSL_SPMATRIX_CSC)
#define GSL_SPMATRIX_ISCSR(m)         ((m)->sptype == GSL_SPMATRIX_CSR)

#define GSL_SPMATRIX_ISTRIPLET(m)     GSL_SPMATRIX_ISCOO(m)
#define GSL_SPMATRIX_ISCCS(m)         GSL_SPMATRIX_ISCSC(m)
#define GSL_SPMATRIX_ISCRS(m)         GSL_SPMATRIX_ISCSR(m)

#define GSL_SPMATRIX_FLG_GROW         (1 << 0) /* allow size of matrix to grow as elements are added */
#define GSL_SPMATRIX_FLG_FIXED        (1 << 1) /* sparsity pattern is fixed */

/* compare matrix entries (ia,ja) and (ib,jb) - sort by rows first, then by columns */
#define GSL_SPMATRIX_COMPARE_ROWCOL(m,ia,ja,ib,jb)   ((ia) < (ib) ? -1 : ((ia) > (ib) ? 1 : ((ja) < (jb) ? -1 : ((ja) > (jb)))))

/* common/utility functions */

void gsl_spmatrix_cumsum(const size_t n, int * c);

__END_DECLS

#endif /* __GSL_SPMATRIX_TYPES_H__ */
