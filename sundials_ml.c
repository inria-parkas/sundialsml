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

/* Sundials interface functions that are common to CVODE and IDA. */

#include <sundials/sundials_config.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_band.h>

#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/alloc.h>
#include <caml/callback.h>
#include <caml/bigarray.h>

#include <stdio.h>

#include "sundials_ml.h"

value sundials_ml_exn_table = 0;

void sundials_ml_register_exns(enum sundials_exn_set_index index, value exns)
{
    CAMLparam1 (exns);
    CAMLlocal1 (r);
    if (sundials_ml_exn_table == 0) {
	int i;
	r = caml_alloc_small (SUNDIALS_NUM_EXN_SETS, 0);
        for (i = 0; i < SUNDIALS_NUM_EXN_SETS; ++i)
	    Field (r, i) = 0;
	Field (r, index) = exns;
	sundials_ml_exn_table = r;
	caml_register_generational_global_root (&sundials_ml_exn_table);
    } else {
	Store_field (sundials_ml_exn_table, index, exns);
    }
    CAMLnoreturn;
}

/* Setting up access to Weak.get */

#if !HAVE_WEAK
static value weak_get = 0;
#endif

CAMLprim value c_sundials_init_module (value vweak_get, value exns)
{
    CAMLparam2 (vweak_get, exns);
    CAMLlocal1 (r);
    REGISTER_EXNS (SUNDIALS, exns);
#if !HAVE_WEAK
    weak_get = vweak_get;
    caml_register_generational_global_root (&weak_get);
#endif
    r = caml_alloc_tuple (3);
    Store_field (r, 0, Val_bool (SUNDIALS_BLAS_LAPACK));
    Store_field (r, 1, caml_copy_double(BIG_REAL));
    Store_field (r, 2, caml_copy_double(UNIT_ROUNDOFF));
    CAMLreturn (r);
}

#if !HAVE_WEAK
CAMLprim value sundials_ml_weak_get (value ar, value n)
{
    CAMLparam2 (ar, n);
    CAMLreturn (caml_callback2 (weak_get, ar, n));
}
#endif	/* !HAVE_WEAK */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * Two-dimensional real arrays based on BigArrays with a column-access table.
 *
 * We represent these as pairs:
 *	Field(ra, 0) = the underlying 2 dimensional big array
 *		       (for use from within OCaml)
 *	Field(ra, 1) = a custom value giving a realtype ** table pointing into
 *		       the columns of the big array (for passing to Sundials)
 */

CAMLprim value c_sundials_realarray2_wrap(value vba)
{
    CAMLparam1(vba);
    CAMLlocal1(r);

    struct caml_ba_array *ba = Caml_ba_array_val(vba);
    int nc = ba->dim[0];
    int nr = ba->dim[1];

    mlsize_t table_size = nc * sizeof(realtype *);
    value vtable = caml_alloc_final(1 + nc, NULL, table_size, table_size * 20);
    realtype **table = (realtype **)Data_custom_val(vtable);

    int j;

    table[0] = (realtype *)(ba->data);
    for (j = 1; j < nc; ++j) {
	table[j] = table[j - 1] + nr;
    }

    r = caml_alloc_tuple(2);
    Store_field(r, 0, vba);
    Store_field(r, 1, vtable);

    CAMLreturn(r);
}

CAMLprim void sundials_crash (value msg)
{
    CAMLparam1 (msg);
    fputs (String_val (msg), stderr);
    fflush (stderr);
    abort ();
    CAMLreturn0;
}
