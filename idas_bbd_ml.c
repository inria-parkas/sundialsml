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

#include <idas/idas.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_band.h>
#include <sundials/sundials_nvector.h>
#include <idas/idas_bbdpre.h>

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/unixsupport.h>
#include <caml/bigarray.h>

#include "dls_ml.h"
#include "spils_ml.h"
#include "sundials_ml.h"
#include "ida_ml.h"
#include "idas_ml.h"
#include "nvector_ml.h"

/* callbacks */

static int bbbdlocal(long int nlocal, realtype t, N_Vector yy, N_Vector yp,
		     N_Vector yyB, N_Vector ypB, N_Vector glocal,
		     void *user_data)
{
    CAMLparam0();
    CAMLlocal3(args, session, cb);

    args = caml_alloc_tuple (RECORD_IDAS_ADJ_BRESFN_ARGS_SIZE);
    Store_field (args, RECORD_IDAS_ADJ_BRESFN_ARGS_T, caml_copy_double (t));
    Store_field (args, RECORD_IDAS_ADJ_BRESFN_ARGS_Y, NVEC_BACKLINK (yy));
    Store_field (args, RECORD_IDAS_ADJ_BRESFN_ARGS_YP, NVEC_BACKLINK (yp));
    Store_field (args, RECORD_IDAS_ADJ_BRESFN_ARGS_YB, NVEC_BACKLINK (yyB));
    Store_field (args, RECORD_IDAS_ADJ_BRESFN_ARGS_YBP, NVEC_BACKLINK (ypB));

    WEAK_DEREF (session, *(value*)user_data);
    cb = IDA_LS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_IDAS_BBBD_CALLBACKS_LOCAL_FN);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (cb, args, NVEC_BACKLINK (glocal));

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

static int bbbdcomm(long int nlocal, realtype t,
		    N_Vector yy, N_Vector yp,
		    N_Vector yyB, N_Vector ypB,
		    void *user_data)
{
    CAMLparam0();
    CAMLlocal3(args, session, cb);

    args = caml_alloc_tuple (RECORD_IDAS_ADJ_BRESFN_ARGS_SIZE);
    Store_field (args, RECORD_IDAS_ADJ_BRESFN_ARGS_T, caml_copy_double (t));
    Store_field (args, RECORD_IDAS_ADJ_BRESFN_ARGS_Y, NVEC_BACKLINK (yy));
    Store_field (args, RECORD_IDAS_ADJ_BRESFN_ARGS_YP, NVEC_BACKLINK (yp));
    Store_field (args, RECORD_IDAS_ADJ_BRESFN_ARGS_YB, NVEC_BACKLINK (yyB));
    Store_field (args, RECORD_IDAS_ADJ_BRESFN_ARGS_YBP, NVEC_BACKLINK (ypB));

    WEAK_DEREF (session, *(value*)user_data);
    cb = IDA_LS_CALLBACKS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_IDAS_BBBD_CALLBACKS_COMM_FN);
    cb = Field (cb, 0);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (cb, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

CAMLprim value c_idas_bbd_prec_initb (value vparentwhich, value vlocaln,
				      value vbandwidths, value vdqrely,
				      value vhascomm)
{
    CAMLparam5(vparentwhich, vlocaln, vbandwidths, vdqrely, vhascomm);
    void *ida_mem = IDA_MEM_FROM_ML (Field(vparentwhich, 0));
    int flag;

    flag = IDABBDPrecInitB (ida_mem, Int_val(Field(vparentwhich, 1)),
	Long_val(vlocaln),
	Long_val(Field(vbandwidths, RECORD_IDA_BANDWIDTHS_MUDQ)),
	Long_val(Field(vbandwidths, RECORD_IDA_BANDWIDTHS_MLDQ)),
	Long_val(Field(vbandwidths, RECORD_IDA_BANDWIDTHS_MUKEEP)),
	Long_val(Field(vbandwidths, RECORD_IDA_BANDWIDTHS_MLKEEP)),
	Double_val(vdqrely),
	bbbdlocal,
	Bool_val(vhascomm) ? bbbdcomm : NULL);
    CHECK_FLAG ("IDABBDPrecInitB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_idas_bbd_prec_reinitb (value vparent, value vwhich,
					value vmudq, value vmldq,
					value vdqrely)
{
    CAMLparam5(vparent, vwhich, vmudq, vmldq, vdqrely);
    void *ida_mem = IDA_MEM_FROM_ML (vparent);
    int flag;

    flag = IDABBDPrecReInitB (ida_mem, Int_val(vwhich),
			      Long_val(vmudq), Long_val(vmldq),
			      Double_val(vdqrely));
    CHECK_FLAG ("IDABBDPrecReInitB", flag);

    CAMLreturn (Val_unit);
}

