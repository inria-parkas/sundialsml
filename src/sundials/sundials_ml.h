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

/*
 * This header defines all constants and functions which are common to CVODE
 * and IDA.  Mostly corresponds to the code in sundials.ml.
 */

#ifndef _SUNDIALS_ML_H__
#define _SUNDIALS_ML_H__

#include <sundials/sundials_types.h>

#include <caml/mlvalues.h>

#include "../config.h"
#include <stdio.h>
#include <assert.h>

#if SUNDIALS_LIB_VERSION < 260
#define SUNMAX(A, B) ((A) > (B) ? (A) : (B))
#define SUNMIN(A, B) ((A) < (B) ? (A) : (B))
#define SUNSQR(A) ((A)*(A))
#endif

#if SUNDIALS_LIB_VERSION < 600
#define sunrealtype realtype
#define sunbooleantype booleantype
#endif

#define CAML_NAME_SPACE

void sunml_warn_discarded_exn (value exn, const char *context);

/* Sundials and integers
 *
 * Sundials < 3.0.0
 *     nvector vspace op:         long int
 *     *_bbd callbacks:           long int
 *     sparse matrix indices:     int
 *     dense/band matrix indices: long int
 *
 * Sundials >= 3.0.0
 *    nvector vspace op:          sunindextype
 *    *_bbd callbacks:            sunindextype
 *    sparse matrix indices:      sunindextype
 *    dense/band matrix indices:  sunindextype
 */
#if SUNDIALS_LIB_VERSION >= 300
typedef sunindextype sundials_ml_index;
typedef sundials_ml_index sundials_ml_smat_index;
#define SmatIndex_val(x) Index_val(x)
#define Val_SmatIndex(x) Val_index(x)
#else
typedef long int sundials_ml_index;
typedef int sundials_ml_smat_index;
#define SmatIndex_val(x) Int_val(x)
#define Val_SmatIndex(x) Val_int(x)
#endif
#define INDEX_ARRAY(v) ((sundials_ml_index *)Caml_ba_data_val(v))

/* Interfacing with OCaml's bigarray infrastructure.  */
#define BIGARRAY_FLOAT (CAML_BA_FLOAT64 | CAML_BA_C_LAYOUT)
#define BIGARRAY_INDEX (CAML_BA_INDEX | CAML_BA_C_LAYOUT)

#define INT_ARRAY(v) ((int *)Caml_ba_data_val(v))
#define LONG_ARRAY(v) ((long int *)Caml_ba_data_val(v))
#define REAL_ARRAY(v) ((sunrealtype *)Caml_ba_data_val(v))

#define ARRAY1_LEN(v) (Caml_ba_array_val((v))->dim[0])

#define OCAMLSIZEOF(type) ((sizeof(type) + sizeof(value) - 1) / sizeof(value))

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

#if OCAML_VERSION < 41200
#define Val_none (Val_int(0))
#define Some_val(v) (Field((v), 0))
#define Is_none(v) ((v) == Val_none)
#endif

#define Store_some(dest, v)			\
    do {					\
      dest = caml_alloc_tuple(1);		\
      Store_field(dest, 0, v);			\
    } while (0)

#define ARRAY2_BA(v)    (Caml_ba_array_val(Field((v), 0)))
#define ARRAY2_DATA(v)  ((sunrealtype *)Caml_ba_data_val(Field((v), 0)))
#define ARRAY2_ACOLS(v) ((sunrealtype **) Data_custom_val(Field((v), 1)))
#define ARRAY2_NCOLS(v) (ARRAY2_DATA(v)->dim[0])
#define ARRAY2_NROWS(v) (ARRAY2_DATA(v)->dim[1])

// wrap a pointer to an integrator (cvode_mem, ida_mem. etc.)
value sunml_wrap_session_pointer(void *sun_mem);
#define SUNML_MEM(v) (*(void **)Data_custom_val(v))

// create a Sundials.RealArray2.t from C
CAMLprim value sunml_sundials_realarray2_create(int nc, int nr);

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
    NVEC_EXN_SET,
    CVODE_EXN_SET,
    CVODES_EXN_SET,
    IDA_EXN_SET,
    IDAS_EXN_SET,
    KINSOL_EXN_SET,
    ARKODE_EXN_SET,
    NVECTOR_PARALLEL_EXN_SET,
    MATRIX_EXN_SET,
    LSOLVER_EXN_SET,
    NLSOLVER_EXN_SET,
    SUNDIALS_NUM_EXN_SETS
};

enum sundials_version_number_index {
  RECORD_SUNDIALS_VERSION_NUMBER_MAJOR = 0,
  RECORD_SUNDIALS_VERSION_NUMBER_MINOR,
  RECORD_SUNDIALS_VERSION_NUMBER_PATCH,
  RECORD_SUNDIALS_VERSION_NUMBER_LABEL,
  RECORD_SUNDIALS_VERSION_NUMBER_SIZE
};

/* Set a field in sundials_ml_exn_table */
void sunml_register_exns (enum sundials_exn_set_index index, value exns);

#define REGISTER_EXNS(MODULE, exns)					\
    (assert (Wosize_val (exns) == MODULE ## _EXN_SET_SIZE),		\
     sunml_register_exns (MODULE ## _EXN_SET, exns))

#define REGISTERED_EXN_OBJ(MODULE, name)			\
    (Field (Field (sundials_ml_exn_table, MODULE ## _EXN_SET),	\
	    MODULE ## _EXN_ ## name))

#define REGISTERED_EXN_TAG(MODULE, name)	\
    (Field (REGISTERED_EXN_OBJ (MODULE, name), 0))

/* For constant exceptions only.  Exceptions with arguments should be
 * raised by extracting the tag with REGISTERED_EXN_TAG.  */
#define REGISTERED_EXN(MODULE, name) REGISTERED_EXN_OBJ(MODULE, name)

/* This enum must list exceptions in the same order as the call to
 * c_init_module in sundials.ml.  */
enum sundials_exn_index {
  SUNDIALS_EXN_RecoverableFailure = 0,
  SUNDIALS_EXN_NonPositiveEwt,
  SUNDIALS_EXN_NotImplementedBySundialsVersion,
  SUNDIALS_EXN_SET_SIZE
};

#define SUNDIALS_EXN(name)     REGISTERED_EXN(SUNDIALS, name)
#define SUNDIALS_EXN_TAG(name) REGISTERED_EXN_TAG(SUNDIALS, name)

/* Return the value that a vptr points to. See Sundials_impl.make_vptr. */
#define VPTRCROOT(x) (*(value **)Data_custom_val(Field(x, 1)))

/* Callback functions are passed from OCaml to C by basically the same
 * mechanism as exceptions. We register each as a separate generational global
 * root.
 *
 * Each module defines a file-local variable named ocaml_values, of type
 * value[NUM_OCAML_VALUES], and an enum type for indexing this array with
 * prefix IX_.
 */
#define REGISTER_OCAML_VALUES(cbs)					\
    do {								\
	int _i;								\
	assert (Wosize_val (cbs) == NUM_OCAML_VALUES);			\
	for (_i = 0; _i < NUM_OCAML_VALUES; ++_i) {			\
	    ocaml_values[_i] = Field (cbs, _i);				\
	    caml_register_generational_global_root (&ocaml_values[_i]);	\
	}								\
    } while (0)

#define OCAML_VALUE(fcn) (ocaml_values[IX_ ## fcn])

/* Accessing FILE* values */
#define ML_CFILE(v) (*(FILE **)Data_custom_val(v))

/* SUNOutputFormat values */
enum sundials_output_format_tag {
  VARIANT_SUNDIALS_OUTPUTFORMAT_TABLE = 0,
  VARIANT_SUNDIALS_OUTPUTFORMAT_CSV,
};

#if 620 <= SUNDIALS_LIB_VERSION
#define SUNML_OUTPUT_FORMAT(v) \
    ((Int_val(v) == VARIANT_SUNDIALS_OUTPUTFORMAT_CSV) \
	? SUN_OUTPUTFORMAT_CSV : SUN_OUTPUTFORMAT_TABLE)
#endif

/* Accessing SUNProfilter values */
#if 600 <= SUNDIALS_LIB_VERSION
#define ML_PROFILER(v) (*(SUNProfiler *)Data_custom_val(v))
#endif

/* Accessing SUNLogger values */
#if 620 <= SUNDIALS_LIB_VERSION
#define ML_LOGGER(v) (*(SUNLogger *)Data_custom_val(v))
#define SUNML_HAS_LOGGING 0 < SUNDIALS_LOGGING_LEVEL
#else
#define SUNML_HAS_LOGGING 0
#endif

enum sundials_logger_level_tag {
  VARIANT_SUNDIALS_LOGGER_LEVEL_ERROR = 0,
  VARIANT_SUNDIALS_LOGGER_LEVEL_WARNING,
  VARIANT_SUNDIALS_LOGGER_LEVEL_INFO,
  VARIANT_SUNDIALS_LOGGER_LEVEL_DEBUG,
};

/* Accessing SUNContext values */
#if 600 <= SUNDIALS_LIB_VERSION
#define ML_CCONTEXT(v) (*(SUNContext *)Data_custom_val(v))
#define ML_CONTEXT(v) (ML_CCONTEXT(Field(v, 0)))
#endif

CAMLprim value sunml_sundials_wrap_file(FILE* f);

/* Accessing callback values */
#define ML_CALLBACK_CPTR(v) (Field(v, 0))

/* Functions for storing OCaml values in the C heap, with registration of a
   global root and construction of a block header. */

value *sunml_sundials_malloc_value(value);
void sunml_sundials_free_value(value *heapref);

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
