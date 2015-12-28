/***********************************************************************
 *                                                                     *
 *                   OCaml interface to Sundials                       *
 *                                                                     *
 *             Timothy Bourke, Jun Inoue, and Marc Pouzet              *
 *             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *
 *                                                                     *
 *  Copyright 2015 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a New BSD License, refer to the file LICENSE.                *
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

#include <arkode/arkode.h>
#include <arkode/arkode_bbdpre.h>

#include "sundials_ml.h"
#include "arkode_ml.h"
#include "nvector_ml.h"


/* Callbacks */

static int bbdlocal(long int nlocal, realtype t, N_Vector y, N_Vector glocal,
		    void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 3);
    CAMLlocal2(session, cb);

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK(y);
    args[2] = NVEC_BACKLINK(glocal);

    WEAK_DEREF (session, *(value*)user_data);
    cb = ARKODE_LS_PRECFNS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_ARKODE_BBD_PRECFNS_LOCAL_FN);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (cb, args[0], args[1], args[2]);

    CAMLreturnT(int, CHECK_EXCEPTION(session, r, RECOVERABLE));
}

static int bbdcomm(long int nlocal, realtype t, N_Vector y, void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 2);
    CAMLlocal2(session, cb);

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK(y);

    WEAK_DEREF (session, *(value*)user_data);
    cb = ARKODE_LS_PRECFNS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_ARKODE_BBD_PRECFNS_COMM_FN);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (cb, args[0], args[1]);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

CAMLprim value c_arkode_bbd_prec_init (value varkode_mem, value vlocaln,
				       value vbandwidths, value vdqrely,
				       value vhascomm)
{
    CAMLparam5(varkode_mem, vlocaln, vbandwidths, vdqrely, vhascomm);
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

    flag = ARKBBDPrecInit (arkode_mem,
	Long_val(vlocaln),
	Long_val(Field(vbandwidths, RECORD_ARKODE_BANDBLOCK_BANDWIDTHS_MUDQ)),
	Long_val(Field(vbandwidths, RECORD_ARKODE_BANDBLOCK_BANDWIDTHS_MLDQ)),
	Long_val(Field(vbandwidths, RECORD_ARKODE_BANDBLOCK_BANDWIDTHS_MUKEEP)),
	Long_val(Field(vbandwidths, RECORD_ARKODE_BANDBLOCK_BANDWIDTHS_MLKEEP)),
	Double_val(vdqrely),
	bbdlocal,
	Bool_val(vhascomm) ? bbdcomm : NULL);
    CHECK_FLAG ("ARKBBDPrecInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_arkode_bbd_prec_reinit (value varkode_mem, value vmudq,
					 value vmldq, value vdqrely)
{
    CAMLparam4(varkode_mem, vmudq, vmldq, vdqrely);
    void *arkode_mem = ARKODE_MEM_FROM_ML (varkode_mem);
    int flag;

    flag = ARKBBDPrecReInit (arkode_mem, Long_val(vmudq), Long_val(vmldq),
				      Double_val(vdqrely));
    CHECK_FLAG ("ARKBBDPrecReInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_arkode_bbd_get_work_space(value varkode_mem)
{
    CAMLparam1(varkode_mem);
    CAMLlocal1(r);

    int flag;
    long int lenrw;
    long int leniw;

    flag = ARKBBDPrecGetWorkSpace(ARKODE_MEM_FROM_ML(varkode_mem),
				  &lenrw, &leniw);
    CHECK_FLAG("ARKBBDPrecGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrw));
    Store_field(r, 1, Val_long(leniw));

    CAMLreturn(r);
}

CAMLprim value c_arkode_bbd_get_num_gfn_evals(value varkode_mem)
{
    CAMLparam1(varkode_mem);

    int flag;
    long int v;

    flag = ARKBBDPrecGetNumGfnEvals(ARKODE_MEM_FROM_ML(varkode_mem), &v);
    CHECK_FLAG("ARKBBDPrecGetNumGfnEvals", flag);

    CAMLreturn(Val_long(v));
}

