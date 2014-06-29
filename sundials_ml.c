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

#include "sundials_ml.h"

CAMLprim value sundials_ml_blas_lapack_supported ()
{
    CAMLparam0();
    CAMLreturn(Val_bool (SUNDIALS_BLAS_LAPACK));
}

CAMLprim value sundials_ml_big_real()
{
    CAMLparam0();
    CAMLreturn(caml_copy_double(BIG_REAL));
}

CAMLprim value sundials_ml_unit_roundoff()
{
    CAMLparam0();
    CAMLreturn(caml_copy_double(UNIT_ROUNDOFF));
}


/* Setting up access to Weak.get */

#if !HAVE_WEAK
static value weak_get = 0;
#endif

CAMLprim void sundials_ml_register_weak_get (value vweak_get)
{
    CAMLparam1 (vweak_get);
#if !HAVE_WEAK
    caml_register_generational_global_root (&weak_get);
    weak_get = vweak_get;
#endif
    CAMLreturn0;
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

