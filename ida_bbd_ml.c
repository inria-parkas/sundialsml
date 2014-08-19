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

// Headers for IDA works for both IDA and IDAS, since this module uses
// a common subset of both.
#include <idas/idas.h>
#include <idas/idas_bbdpre.h>

#include "sundials_ml.h"
#include "ida_ml.h"
#include "nvector_ml.h"

static int bbdlocal(long int nlocal, realtype t, N_Vector y, N_Vector yp,
		    N_Vector gval, void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 5);
    int r;
    value *backref = user_data;
    CAML_FN (call_bbdlocal);

    args[0] = *backref;
    args[1] = caml_copy_double(t);
    args[2] = NVEC_BACKLINK(y);
    args[3] = NVEC_BACKLINK(yp);
    args[4] = NVEC_BACKLINK(gval);

    r = Int_val (caml_callbackN(*call_bbdlocal,
                                sizeof (args) / sizeof (*args),
                                args));

    CAMLreturnT(int, r);
}

static int bbdcomm(long int nlocal, realtype t, N_Vector y, N_Vector yp,
		   void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 4);
    int r;
    value *backref = user_data;
    CAML_FN (call_bbdcomm);

    args[0] = *backref;
    args[1] = caml_copy_double(t);
    args[2] = NVEC_BACKLINK(y);
    args[3] = NVEC_BACKLINK(yp);

    r = Int_val (caml_callbackN(*call_bbdcomm,
                                sizeof (args) / sizeof (*args),
                                args));

    CAMLreturnT(int, r);
}

CAMLprim void c_ida_bbd_prec_init (value vida_mem, value vlocaln,
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

    CAMLreturn0;
}

CAMLprim void c_ida_bbd_prec_reinit (value vida_mem, value vmudq,
				     value vmldq, value vdqrely)
{
    CAMLparam4(vida_mem, vmudq, vmldq, vdqrely);
    void *ida_mem = IDA_MEM_FROM_ML (vida_mem);
    int flag;

    flag = IDABBDPrecReInit (ida_mem, Long_val(vmudq), Long_val(vmldq),
				       Double_val(vdqrely));
    CHECK_FLAG ("IDABBDPrecReInit", flag);

    CAMLreturn0;
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
