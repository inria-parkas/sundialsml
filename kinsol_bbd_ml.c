/***********************************************************************
 *                                                                     *
 *                   OCaml interface to Sundials                       *
 *                                                                     *
 *             Timothy Bourke, Jun Inoue, and Marc Pouzet              *
 *             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *
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

enum callback_index {
    IX_call_bbdlocal = 0,
    IX_call_bbdcomm,
    NUM_CALLBACKS
};

static value callbacks[NUM_CALLBACKS];

CAMLprim value c_kinsol_bbd_init_module (value cbs)
{
    CAMLparam1 (cbs);
    REGISTER_CALLBACKS (cbs);
    CAMLreturn (Val_unit);
}

/* Sundials 2.5.0 User's Guide incorrectly states that KINLocalFn
 * returns void.  The comment in kinsol_bbdpre.h says it should return
 * 0 for success, non-zero otherwise.  */
static int bbdlocal(long int nlocal, N_Vector u, N_Vector gval, void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 2);
    CAMLlocal2(session, cb);

    args[0] = NVEC_BACKLINK(u);
    args[1] = NVEC_BACKLINK(gval);

    WEAK_DEREF (session, *(value*)user_data);
    cb = KINSOL_LS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_KINSOL_BBD_CALLBACKS_LOCAL_FN);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (cb, args[0], args[1]);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, UNRECOVERABLE));
}

/* Sundials 2.5.0 User's Guide incorrectly states that KINCommFn
 * returns void.  The comment in kinsol_bbdpre.h says it should return
 * 0 for success, non-zero otherwise.  */
static int bbdcomm(long int nlocal, N_Vector u, void *user_data)
{
    CAMLparam0();
    CAMLlocal2(session, cb);

    cb = KINSOL_LS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_KINSOL_BBD_CALLBACKS_LOCAL_FN);
    cb = Field (cb, 0);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (cb, NVEC_BACKLINK(u));

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, UNRECOVERABLE));
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
  
