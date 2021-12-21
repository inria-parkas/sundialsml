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

#include <cvodes/cvodes.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_band.h>
#include <sundials/sundials_nvector.h>
#include <cvodes/cvodes_bbdpre.h>

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/bigarray.h>

#include "../sundials/sundials_ml.h"
#include "../cvode/cvode_ml.h"
#include "cvodes_ml.h"
#include "../nvectors/nvector_ml.h"

/* callbacks */

static int bbbdlocal(sundials_ml_index nlocal, sunrealtype t, N_Vector y,
		     N_Vector yb, N_Vector glocal, void *user_data)
{
    CAMLparam0();
    CAMLlocal3(args, session, cb);

    args = caml_alloc_tuple (RECORD_CVODES_ADJ_BRHSFN_ARGS_SIZE);
    Store_field (args, RECORD_CVODES_ADJ_BRHSFN_ARGS_T, caml_copy_double (t));
    Store_field (args, RECORD_CVODES_ADJ_BRHSFN_ARGS_Y, NVEC_BACKLINK (y));
    Store_field (args, RECORD_CVODES_ADJ_BRHSFN_ARGS_YB, NVEC_BACKLINK (yb));

    WEAK_DEREF (session, *(value*)user_data);
    cb = CVODE_LS_PRECFNS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_CVODES_BBBD_PRECFNS_LOCAL_FN);
    assert (Tag_val (cb) == Closure_tag);


    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn (cb, args, NVEC_BACKLINK (glocal));

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

static int bbbdcomm(sundials_ml_index nlocal, sunrealtype t, N_Vector y,
		    N_Vector yb, void *user_data)
{
    CAMLparam0();
    CAMLlocal3(args, session, cb);

    args = caml_alloc_tuple (RECORD_CVODES_ADJ_BRHSFN_ARGS_SIZE);
    Store_field (args, RECORD_CVODES_ADJ_BRHSFN_ARGS_T, caml_copy_double (t));
    Store_field (args, RECORD_CVODES_ADJ_BRHSFN_ARGS_Y, NVEC_BACKLINK (y));
    Store_field (args, RECORD_CVODES_ADJ_BRHSFN_ARGS_YB, NVEC_BACKLINK (yb));

    WEAK_DEREF (session, *(value*)user_data);
    cb = CVODE_LS_PRECFNS_FROM_ML (session);
    cb = Field (cb, 0);
    cb = Field (cb, RECORD_CVODES_BBBD_PRECFNS_COMM_FN);
    cb = Some_val (cb);
    assert (Tag_val (cb) == Closure_tag);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn (cb, args);

    CAMLreturnT(int, CHECK_EXCEPTION (session, r, RECOVERABLE));
}

CAMLprim value sunml_cvodes_bbd_prec_initb (value vparentwhich, value vlocaln,
					value vbandwidths, value vdqrely,
					value vhascomm)
{
    CAMLparam5(vparentwhich, vlocaln, vbandwidths, vdqrely, vhascomm);
    void *cvode_mem = CVODE_MEM_FROM_ML (Field(vparentwhich, 0));
    int flag;

    flag = CVBBDPrecInitB (cvode_mem, Int_val(Field(vparentwhich, 1)),
	Long_val(vlocaln),
	Long_val(Field(vbandwidths, RECORD_CVODE_BANDBLOCK_BANDWIDTHS_MUDQ)),
	Long_val(Field(vbandwidths, RECORD_CVODE_BANDBLOCK_BANDWIDTHS_MLDQ)),
	Long_val(Field(vbandwidths, RECORD_CVODE_BANDBLOCK_BANDWIDTHS_MUKEEP)),
	Long_val(Field(vbandwidths, RECORD_CVODE_BANDBLOCK_BANDWIDTHS_MLKEEP)),
	Double_val(vdqrely),
	bbbdlocal,
	Bool_val(vhascomm) ? bbbdcomm : NULL);
    CHECK_FLAG ("CVBBDPrecInitB", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_cvodes_bbd_prec_reinitb (value vparent, value vwhich,
					  value vmudq, value vmldq,
					  value vdqrely)
{
    CAMLparam5(vparent, vwhich, vmudq, vmldq, vdqrely);
    void *cvode_mem = CVODE_MEM_FROM_ML (vparent);
    int flag;

    flag = CVBBDPrecReInitB (cvode_mem, Int_val(vwhich),
			     Long_val(vmudq), Long_val(vmldq),
			     Double_val(vdqrely));
    CHECK_FLAG ("CVBBDPrecReInitB", flag);

    CAMLreturn (Val_unit);
}

