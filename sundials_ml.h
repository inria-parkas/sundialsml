/***********************************************************************
 *                                                                     *
 *     OCaml interface to Sundials (serial) CVODE and IDA solvers      *
 *                                                                     *
 *  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *
 *                                                                     *
 *  Copyright 2013 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under the terms of the GNU Library General Public License, with    *
 *  the special exception on linking described in file LICENSE.        *
 *                                                                     *
 ***********************************************************************/

/*
 * This header defines all constants and functions which are common to CVODE
 * and IDA.  Mostly corresponds to the code in sundials.ml.
 */

#ifndef _SUNDIALS_ML_H__
#define _SUNDIALS_ML_H__

value sundials_ml_big_real ();
value sundials_ml_unit_roundoff();

/* Interfacing with OCaml's bigarray infrastructure.  */
#define BIGARRAY_FLOAT (CAML_BA_FLOAT64 | CAML_BA_C_LAYOUT)
#define BIGARRAY_INT (CAML_BA_INT32 | CAML_BA_C_LAYOUT)

#define INT_ARRAY(v) ((int *)Caml_ba_data_val(v))
#define LONG_ARRAY(v) ((long int *)Caml_ba_data_val(v))
#define REAL_ARRAY(v) ((realtype *)Caml_ba_data_val(v))
#define REAL_ARRAY2(v) ((realtype **)Caml_ba_data_val(v))

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

#endif /* _SUNDIALS_ML_H__ */
