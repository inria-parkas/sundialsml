/***********************************************************************
 *                                                                     *
 *                   OCaml interface to Sundials                       *
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

#if (NDEBUG+0) && (DEBUG+0)
#error "DEBUG and NDEBUG are both defined.  I'm confused."
#endif

#include <assert.h>

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
    assert (Is_block (dest));					\
    dest = Field (dest, 0);                                     \
  } while (0)

#define Val_none (Val_int(0))
#define Some_val(v) (Field((v), 0))
#define Store_some(dest, v) { dest = caml_alloc_small(1, 0); \
			      Field(dest, 0) = v; }

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

/* Exceptions are registered in one global variable,
   sundials_ml_exn_table.  This variable conceptually has the type
   type sundials_ml_exn_table =
     {
       sundials : exn array;
       cvode : exn array;
       cvodes : exn array;
       ...
     }
   and each field holds the exceptions used in each module.  The fields
   are all initialized to 0 by sundials_ml.c and then populated on a need-
   basis by the initialization code of each module.

   Note: C standard doesn't say 0 == (long)NULL, but this is a common
   assumption and is adopted by OCaml's runtime.  See e.g. caml_alloc.
 */
extern value sundials_ml_exn_table;

enum sundials_exn_set_index {
    SUNDIALS_EXN_SET = 0,
    CVODE_EXN_SET,
    CVODES_EXN_SET,
    IDA_EXN_SET,
    IDAS_EXN_SET,
    KINSOL_EXN_SET,
    DLS_EXN_SET,
    SPILS_EXN_SET,
    NVECTOR_PARALLEL_EXN_SET,
    SUNDIALS_NUM_EXN_SETS
};

/* Set a field in sundials_ml_exn_table */
void sundials_ml_register_exns (enum sundials_exn_set_index index, value exns);

/* This enum must list exceptions in the same order as the call to
 * c_register_exns in sundials.ml.  */
enum sundials_exn_index {
  SUNDIALS_EXN_RecoverableFailure = 0,
  SUNDIALS_EXN_NonPositiveEwt,
  SUNDIALS_EXN_SET_SIZE
};

#define SUNDIALS_EXN(name) (Field (Field (sundials_ml_exn_table,	\
					  SUNDIALS_EXN_SET),		\
				   SUNDIALS_EXN_ ## name))

#define REGISTER_EXNS(MODULE, exns)					\
    (assert (Wosize_val (exns) == MODULE ## _EXN_SET_SIZE),		\
     sundials_ml_register_exns (MODULE ## _EXN_SET, exns))

/* Callback functions are passed from OCaml to C by basically the same
 * mechanism as exceptions, but since callbacks are very frequently
 * accessed, we register them each as a separate generational global
 * root.
 *
 * Each module defines a file-local variable named callbacks, of type
 * value[num_callbacks], and an enum type for indexing this array with
 * prefix IX_.
 */
#define REGISTER_CALLBACKS(cbs)						\
    do {								\
	int _i;								\
	assert (Wosize_val (cbs) == NUM_CALLBACKS);			\
	for (_i = 0; _i < NUM_CALLBACKS; ++_i) {			\
	    callbacks[_i] = Field (Field (cbs, _i), 0);			\
	    caml_register_generational_global_root (&callbacks[_i]);	\
	}								\
    } while (0)

#define CAML_FN(fcn) (callbacks[IX_ ## fcn])

/* Generate trampolines needed for functions with >= 6 arguments.  */
#define COMMA ,
#define BYTE_STUB(fcn_name, extras)				\
    CAMLprim value fcn_name##_byte (value *args, int n)		\
    {								\
	return fcn_name (args[0], args[1], args[2],		\
			 args[3], args[4], args[5] extras);	\
    }
#define BYTE_STUB6(fcn_name)			\
    BYTE_STUB(fcn_name, /* empty */)
#define BYTE_STUB7(fcn_name)			\
    BYTE_STUB(fcn_name, COMMA args[6])
#define BYTE_STUB8(fcn_name)			\
    BYTE_STUB(fcn_name, COMMA args[6] COMMA args[7])
#define BYTE_STUB9(fcn_name)			\
    BYTE_STUB(fcn_name, COMMA args[6] COMMA args[7] COMMA args[8])
#define BYTE_STUB10(fcn_name)			\
    BYTE_STUB(fcn_name, COMMA args[6] COMMA args[7] COMMA args[8] COMMA args[9])

#endif /* _SUNDIALS_ML_H__ */
