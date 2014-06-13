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

/*
 * This header defines all constants and functions which are common to CVODE
 * and IDA.  Mostly corresponds to the code in sundials.ml.
 */

#ifndef _SUNDIALS_ML_H__
#define _SUNDIALS_ML_H__

#include <caml/mlvalues.h>

#include "config.h"

value sundials_ml_big_real ();
value sundials_ml_unit_roundoff();

/* Interfacing with OCaml's bigarray infrastructure.  */
#define BIGARRAY_FLOAT (CAML_BA_FLOAT64 | CAML_BA_C_LAYOUT)
#define BIGARRAY_INT (CAML_BA_INT32 | CAML_BA_C_LAYOUT)

#define INT_ARRAY(v) ((int *)Caml_ba_data_val(v))
#define LONG_ARRAY(v) ((long int *)Caml_ba_data_val(v))
#define REAL_ARRAY(v) ((realtype *)Caml_ba_data_val(v))
#define REAL_ARRAY2(v) ((realtype **)Caml_ba_data_val(v))

#define ARRAY1_LEN(v) (Caml_ba_array_val((v))->dim[0])

#if HAVE_WEAK
CAMLprim value caml_weak_get (value ar, value n);
#define sundials_ml_weak_get caml_weak_get
#else
value sundials_ml_weak_get (value ar, value n);
#endif

#define WEAK_DEREF(dest, ptr)                                   \
  do {                                                          \
    dest = sundials_ml_weak_get ((ptr), Val_int (0));           \
    if (!Is_block (dest))                                       \
      caml_failwith ("Internal error: weak reference is dead"); \
    dest = Field (dest, 0);                                     \
  } while (0)

#define Val_none (Val_int(0))
#define Some_val(v) (Field((v), 0))

#define ARRAY2_DATA(v)  (Caml_ba_array_val(Field((v), 0)))
#define ARRAY2_ACOLS(v) ((realtype **) Data_custom_val(Field((v), 1)))
#define ARRAY2_NCOLS(v) (ARRAY2_DATA(v)->dim[0])
#define ARRAY2_NROWS(v) (ARRAY2_DATA(v)->dim[1])

enum sundials_error_details_index {
  RECORD_SUNDIALS_ERROR_DETAILS_ERROR_CODE    = 0,
  RECORD_SUNDIALS_ERROR_DETAILS_MODULE_NAME,
  RECORD_SUNDIALS_ERROR_DETAILS_FUNCTION_NAME,
  RECORD_SUNDIALS_ERROR_DETAILS_ERROR_MESSAGE,
  RECORD_SUNDIALS_ERROR_DETAILS_SIZE /* This has to come last. */
};

#endif /* _SUNDIALS_ML_H__ */
