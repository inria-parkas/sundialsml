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

/* Sundials interface functions that are common to CVODE and IDA. */

#include <sundials/sundials_config.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_band.h>

#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/alloc.h>

#include "sundials_ml.h"

value sundials_ml_big_real()
{
    CAMLparam0();
    CAMLreturn(caml_copy_double(BIG_REAL));
}

value sundials_ml_unit_roundoff()
{
    CAMLparam0();
    CAMLreturn(caml_copy_double(UNIT_ROUNDOFF));
}

#if !HAVE_WEAK
value sundials_ml_weak_get (value ar, value n)
{
    CAMLparam2 (ar, n);
    value *weak_get;
    if (weak_get == NULL)
	weak_get = caml_named_value ("sundials_weak_get");

    CAMLreturn (caml_callback2 (*weak_get, ar, n));
}
#endif	/* !HAVE_WEAK */
