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

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/bigarray.h>

#include <sundials/sundials_config.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>

// Headers for IDA works for both IDA and IDAS, since this module uses
// a common subset of both.
#include <idas/idas.h>
#include <idas/idas_bbdpre.h>

#include "../sundials/sundials_ml.h"
#include "ida_ml.h"
#include "../nvectors/nvector_ml.h"

/* Callbacks */

static int bbdlocal(sundials_ml_index nlocal, realtype t, N_Vector y,
		    N_Vector yp, N_Vector gval, void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 4);
    CAMLlocal2(session, cb);

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK(y);
    args[2] = NVEC_BACKLINK(yp);
    args[3] = NVEC_BACKLINK(gval);

    WEAK_DEREF (session, *(value*)user_data);
    cb = IDA_LS_PRECFNS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_IDA_BBD_PRECFNS_LOCAL_FN);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (cb, 4, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

static int bbdcomm(sundials_ml_index nlocal, realtype t, N_Vector y,
		   N_Vector yp, void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 3);
    CAMLlocal2(session, cb);

    args[0] = caml_copy_double(t);
    args[1] = NVEC_BACKLINK(y);
    args[2] = NVEC_BACKLINK(yp);

    WEAK_DEREF (session, *(value*)user_data);
    cb = IDA_LS_PRECFNS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_IDA_BBD_PRECFNS_COMM_FN);
    cb = Some_val (cb);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn (cb, args[0], args[1], args[2]);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

CAMLprim value c_ida_bbd_prec_init (value vida_mem, value vlocaln,
				    value vbandwidths, value vdqrely,
				    value vhascomm)
{
    CAMLparam5(vida_mem, vlocaln, vbandwidths, vdqrely, vhascomm);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    int flag;

    flag = IDABBDPrecInit (ida_mem,
	Long_val(vlocaln),
	Long_val(Field(vbandwidths, RECORD_IDA_BANDWIDTHS_MUDQ)),
	Long_val(Field(vbandwidths, RECORD_IDA_BANDWIDTHS_MLDQ)),
	Long_val(Field(vbandwidths, RECORD_IDA_BANDWIDTHS_MUKEEP)),
	Long_val(Field(vbandwidths, RECORD_IDA_BANDWIDTHS_MLKEEP)),
	Double_val(vdqrely),
	bbdlocal,
	Bool_val(vhascomm) ? bbdcomm : NULL);
    CHECK_FLAG ("IDABBDPrecInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_bbd_prec_reinit (value vida_mem, value vmudq,
				      value vmldq, value vdqrely)
{
    CAMLparam4(vida_mem, vmudq, vmldq, vdqrely);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    int flag;

    flag = IDABBDPrecReInit (ida_mem, Long_val(vmudq), Long_val(vmldq),
				       Double_val(vdqrely));
    CHECK_FLAG ("IDABBDPrecReInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_ida_bbd_get_work_space(value vida_mem)
{
    CAMLparam1(vida_mem);
    CAMLlocal1(r);

    int flag;
    long int lenrw;
    long int leniw;

    flag = IDABBDPrecGetWorkSpace(IDA_MEM_FROM_ML(vida_mem), &lenrw, &leniw);
    CHECK_FLAG("IDABBDPrecGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrw));
    Store_field(r, 1, Val_long(leniw));

    CAMLreturn(r);
}

CAMLprim value c_ida_bbd_get_num_gfn_evals(value vida_mem)
{
    CAMLparam1(vida_mem);

    int flag;
    long int v;

    flag = IDABBDPrecGetNumGfnEvals(IDA_MEM_FROM_ML(vida_mem), &v);
    CHECK_FLAG("IDABBDPrecGetNumGfnEvals", flag);

    CAMLreturn(Val_long(v));
}
