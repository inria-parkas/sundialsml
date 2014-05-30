/***********************************************************************
 *                                                                     *
 *               OCaml interface to (serial) Sundials                  *
 *                                                                     *
 *  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *
 *                                                                     *
 *  Copyright 2014 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a BSD 2-Clause License, refer to the file LICENSE.           *
 *                                                                     *
 ***********************************************************************/

/* The parts of the Sundials interface that distinguish between Serial
   NVectors (handled by Bigarrays) and generic NVectors (handled by a
   wrapper type). */

#include <cvodes/cvodes.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_types.h>

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/unixsupport.h>
#include <caml/bigarray.h>

/* linear solvers */
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_band.h>
#include <cvodes/cvodes_diag.h>
#include <cvodes/cvodes_spgmr.h>
#include <cvodes/cvodes_spbcgs.h>
#include <cvodes/cvodes_sptfqmr.h>
#include <cvodes/cvodes_bandpre.h>
#include <sundials/sundials_config.h>

#include "spils_ml.h"
#include "sundials_ml.h"
#include "cvodes_ml.h"
#include "nvector_ml.h"

#if SUNDIALS_BLAS_LAPACK == 1
#include <cvodes/cvodes_lapack.h>
#endif

#ifdef RESTRICT_INTERNAL_PRECISION
#ifdef __GNUC__
#include <fpu_control.h>
#endif
#endif

// Call with CVODE_ML_BIGARRAYS to compile for the Serial NVector to
// Bigarray interface code.

#ifdef CVODE_ML_BIGARRAYS

#define CVTYPE(fname) c_ba_cvodes_ ## fname
#include <nvector/nvector_serial.h>

#define WRAP_NVECTOR(v) caml_ba_alloc(BIGARRAY_FLOAT, 1, NV_DATA_S(v), &(NV_LENGTH_S(v)))
#define RELINQUISH_WRAPPEDNV(v_ba) Caml_ba_array_val(v_ba)->dim[0] = 0

#define NVECTORIZE_VAL(ba) N_VMake_Serial(Caml_ba_array_val(ba)->dim[0], (realtype *)Caml_ba_data_val(ba))
#define RELINQUISH_NVECTORIZEDVAL(nv) N_VDestroy(nv)

#else

#define CVTYPE(fname) c_nvec_cvodes_ ## fname
#include <sundials/sundials_nvector.h>

#define WRAP_NVECTOR(v) NVEC_DATA(v)
#define RELINQUISH_WRAPPEDNV(v) {}

#define NVECTORIZE_VAL(v) NVEC_VAL(v)
#define RELINQUISH_NVECTORIZEDVAL(nv) {}

#endif

#define DOQUOTE(text) #text
#define QUOTE(val) DOQUOTE(val)
#define CVTYPESTR(fname) QUOTE(CVTYPE(fname))

/* callbacks */

#define CAML_FN(name)					\
    static value *name;					\
    if (name == NULL)					\
	name = caml_named_value (CVTYPESTR (name));

