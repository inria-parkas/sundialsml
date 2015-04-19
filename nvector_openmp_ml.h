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

#ifndef __NVECTOR_OPENMP_ML_H__
#define __NVECTOR_OPENMP_ML_H__

#include <sundials/sundials_nvector.h>
#include <caml/mlvalues.h>

/* OCaml interface to OpenMP NVectors.

   See the comments in nvector_ml.h for the interfacing principles common to
   all NVectors.

   OpenMP nvectors
   -----------------
   The payload is a Bigarray.

   The N_Vector content->data field points to the data in the C heap which
   underlies the payload Bigarray.

   The N_Vector ops are identical to those of a standard OpenMP N_Vector,
   except for nvclone, nvcloneempty, and nvdestroy which are functions,
   implemented in nvector_ml.c, to create the arrangement described here.
*/

// Creation functions
value ml_nvec_wrap_openmp(value nthreads, value payload, value checkfn);

#endif
