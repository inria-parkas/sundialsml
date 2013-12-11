/***********************************************************************
 *                                                                     *
 *     OCaml interface to Sundials (serial) CVODE and IDA solvers      *
 *                                                                     *
 *  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *
 *                                                                     *
 *  Copyright 2013 Institut National de Recherche en Informatique et   *
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
