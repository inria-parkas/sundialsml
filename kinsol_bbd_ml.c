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

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/unixsupport.h>
#include <caml/bigarray.h>

#include <sundials/sundials_config.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>

/* linear solvers */
#include <kinsol/kinsol.h>
#include <kinsol/kinsol_bbdpre.h>

#include "dls_ml.h"
#include "spils_ml.h"
#include "kinsol_ml.h"

static int bbdlocal(long int nlocal, N_Vector u, N_Vector gval, void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 3);
    int r;
    value *backref = user_data;
    CAML_FN (call_bbdlocal);

    args[0] = *backref;
    args[1] = NVEC_BACKLINK(u);
    args[2] = NVEC_BACKLINK(gval);

    r = Int_val (caml_callbackN(*call_bbdlocal,
                                sizeof (args) / sizeof (*args),
                                args));

    CAMLreturnT(int, r);
}

static int bbdcomm(long int nlocal, N_Vector u, void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 2);
    int r;
    value *backref = user_data;
    CAML_FN (call_bbdcomm);

    args[0] = *backref;
    args[1] = NVEC_BACKLINK(u);

    r = Int_val (caml_callbackN(*call_bbdcomm,
                                sizeof (args) / sizeof (*args),
                                args));

    CAMLreturnT(int, r);
}

CAMLprim value c_kinsol_bbd_prec_init (value vkin_mem, value vlocaln,
				       value vbandwidths, value vdqrely,
				       value vhascomm)
{
    CAMLparam5(vkin_mem, vlocaln, vbandwidths, vdqrely, vhascomm);
    void *kin_mem = KINSOL_MEM_FROM_ML (vkin_mem);
    int flag;

    flag = KINBBDPrecInit (kin_mem,
	Long_val(vlocaln),
	Long_val(Field(vbandwidths, RECORD_KINSOL_BANDBLOCK_BANDWIDTHS_MUDQ)),
	Long_val(Field(vbandwidths, RECORD_KINSOL_BANDBLOCK_BANDWIDTHS_MLDQ)),
	Long_val(Field(vbandwidths, RECORD_KINSOL_BANDBLOCK_BANDWIDTHS_MUKEEP)),
	Long_val(Field(vbandwidths, RECORD_KINSOL_BANDBLOCK_BANDWIDTHS_MLKEEP)),
	Double_val(vdqrely),
	bbdlocal,
	Bool_val(vhascomm) ? bbdcomm : NULL);
    CHECK_FLAG ("KINBBDPrecInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_kinsol_bbd_get_work_space(value vkin_mem)
{
    CAMLparam1(vkin_mem);
    CAMLlocal1(r);

    int flag;
    long int lenrw;
    long int leniw;

    flag = KINBBDPrecGetWorkSpace(KINSOL_MEM_FROM_ML(vkin_mem), &lenrw, &leniw);
    CHECK_FLAG("KINBBDPrecGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrw));
    Store_field(r, 1, Val_long(leniw));

    CAMLreturn(r);
}

CAMLprim value c_kinsol_bbd_get_num_gfn_evals(value vkin_mem)
{
    CAMLparam1(vkin_mem);

    int flag;
    long int v;

    flag = KINBBDPrecGetNumGfnEvals(KINSOL_MEM_FROM_ML(vkin_mem), &v);
    CHECK_FLAG("KINBBDPrecGetNumGfnEvals", flag);

    CAMLreturn(Val_long(v));
}
  
