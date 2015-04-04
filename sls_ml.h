/***********************************************************************
 *                                                                     *
 *                   OCaml interface to Sundials                       *
 *                                                                     *
 *             Timothy Bourke, Jun Inoue, and Marc Pouzet              *
 *             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *
 *                                                                     *
 *  Copyright 2015 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a New BSD License, refer to the file LICENSE.                *
 *                                                                     *
 ***********************************************************************/

#ifndef _SLS_ML_H__
#define _SLS_ML_H__

#include <caml/mlvalues.h>
#include <sundials/sundials_sparse.h>

#define SLSMAT(v) (*(SlsMat *)Data_custom_val(v))

enum sls_densematrix_index {
  RECORD_SLS_SPARSEMATRIX_COLPTRS   = 0,
  RECORD_SLS_SPARSEMATRIX_ROWVALS,
  RECORD_SLS_SPARSEMATRIX_DATA,
  RECORD_SLS_SPARSEMATRIX_SLSMAT,
  RECORD_SLS_SPARSEMATRIX_VALID,
  RECORD_SLS_SPARSEMATRIX_SIZE /* This has to come last. */
};

CAMLprim value c_sls_invalidate(value);
CAMLprim value c_sls_sparse_wrap(SlsMat a, int finalize);

#endif

