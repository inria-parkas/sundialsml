/***********************************************************************
 *                                                                     *
 *                   OCaml interface to Sundials                       *
 *                                                                     *
 *             Timothy Bourke, Jun Inoue, and Marc Pouzet              *
 *             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *
 *                                                                     *
 *  Copyright 2014 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a New BSD License, refer to the file LICENSE.                *
 *                                                                     *
 ***********************************************************************/

#ifndef _DLS_ML_H__
#define _DLS_ML_H__

#include <caml/mlvalues.h>
#include <sundials/sundials_dense.h>

#define DLSMAT(v) (*(DlsMat *)Data_custom_val(v))

enum dls_densematrix_index {
  RECORD_DLS_DENSEMATRIX_PAYLOAD    = 0,
  RECORD_DLS_DENSEMATRIX_DLSMAT,
  RECORD_DLS_DENSEMATRIX_VALID,
  RECORD_DLS_DENSEMATRIX_SIZE /* This has to come last. */
};

enum dls_bandmatrix_index {
  RECORD_DLS_BANDMATRIX_PAYLOAD    = 0,
  RECORD_DLS_BANDMATRIX_DLSMAT,
  RECORD_DLS_BANDMATRIX_SMU,
  RECORD_DLS_BANDMATRIX_VALID,
  RECORD_DLS_BANDMATRIX_SIZE /* This has to come last. */
};

enum dls_bandmatrix_dims_index {
  RECORD_DLS_BANDMATRIX_DIMS_N    = 0,
  RECORD_DLS_BANDMATRIX_DIMS_MU,
  RECORD_DLS_BANDMATRIX_DIMS_SMU,
  RECORD_DLS_BANDMATRIX_DIMS_ML,
  RECORD_DLS_BANDMATRIX_DIMS_SIZE /* This has to come last. */
};

#endif
