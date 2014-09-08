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

// When we compile with sensitivity (CVODES), we are obliged to use the
// cvodes/cvodes_* header files. In fact, nearly everything functions
// correctly if the cvode/cvode_* header files are used instead (the
// function prototypes are identical), except for the function
// c_cvode_set_alternate which relies on the internal representation of
// CVodeMem (for cv_lsolve, etc.). In any case, it seems a better idea to
// use the appropriate header files even if this introduces a minor
// complication in the build system.

#ifdef SUNDIALSML_WITHSENS
/* CVODES (with sensitivity) */

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_bbdpre.h>

#else
/* CVODE (without sensitivity) */

#include <cvode/cvode.h>
#include <cvode/cvode_bbdpre.h>

#endif	/* SUNDIALSML_WITHSENS */

#include "sundials_ml.h"
#include "cvode_ml.h"
#include "nvector_ml.h"


/* Callbacks */

enum callback_index {
    IX_call_bbdlocal = 0,
    IX_call_bbdcomm,
    NUM_CALLBACKS
};

static value callbacks[NUM_CALLBACKS];

CAMLprim value c_cvode_bbd_init_module (value cbs)
{
    CAMLparam1 (cbs);
    REGISTER_CALLBACKS (cbs);
    CAMLreturn (Val_unit);
}

static int bbdlocal(long int nlocal, realtype t, N_Vector y, N_Vector glocal,
		    void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 4);
    int r;
    value *backref = user_data;

    args[0] = *backref;
    args[1] = caml_copy_double(t);
    args[2] = NVEC_BACKLINK(y);
    args[3] = NVEC_BACKLINK(glocal);

    r = Int_val (caml_callbackN(CAML_FN(call_bbdlocal),
                                sizeof (args) / sizeof (*args),
                                args));

    CAMLreturnT(int, r);
}

static int bbdcomm(long int nlocal, realtype t, N_Vector y, void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 3);
    int r;
    value *backref = user_data;

    args[0] = *backref;
    args[1] = caml_copy_double(t);
    args[2] = NVEC_BACKLINK(y);

    r = Int_val (caml_callbackN(CAML_FN(call_bbdcomm),
                                sizeof (args) / sizeof (*args),
                                args));

    CAMLreturnT(int, r);
}

CAMLprim value c_cvode_bbd_prec_init (value vcvode_mem, value vlocaln,
				      value vbandwidths, value vdqrely,
				      value vhascomm)
{
    CAMLparam5(vcvode_mem, vlocaln, vbandwidths, vdqrely, vhascomm);
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    int flag;

    flag = CVBBDPrecInit (cvode_mem,
	Long_val(vlocaln),
	Long_val(Field(vbandwidths, RECORD_CVODE_BANDBLOCK_BANDWIDTHS_MUDQ)),
	Long_val(Field(vbandwidths, RECORD_CVODE_BANDBLOCK_BANDWIDTHS_MLDQ)),
	Long_val(Field(vbandwidths, RECORD_CVODE_BANDBLOCK_BANDWIDTHS_MUKEEP)),
	Long_val(Field(vbandwidths, RECORD_CVODE_BANDBLOCK_BANDWIDTHS_MLKEEP)),
	Double_val(vdqrely),
	bbdlocal,
	Bool_val(vhascomm) ? bbdcomm : NULL);
    CHECK_FLAG ("CVBBDPrecInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_cvode_bbd_prec_reinit (value vcvode_mem, value vmudq,
					value vmldq, value vdqrely)
{
    CAMLparam4(vcvode_mem, vmudq, vmldq, vdqrely);
    void *cvode_mem = CVODE_MEM_FROM_ML (vcvode_mem);
    int flag;

    flag = CVBBDPrecReInit (cvode_mem, Long_val(vmudq), Long_val(vmldq),
				       Double_val(vdqrely));
    CHECK_FLAG ("CVBBDPrecReInit", flag);

    CAMLreturn (Val_unit);
}

CAMLprim value c_cvode_bbd_get_work_space(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);
    CAMLlocal1(r);

    int flag;
    long int lenrw;
    long int leniw;

    flag = CVBBDPrecGetWorkSpace(CVODE_MEM_FROM_ML(vcvode_mem), &lenrw, &leniw);
    CHECK_FLAG("CVBBDPrecGetWorkSpace", flag);

    r = caml_alloc_tuple(2);

    Store_field(r, 0, Val_long(lenrw));
    Store_field(r, 1, Val_long(leniw));

    CAMLreturn(r);
}

CAMLprim value c_cvode_bbd_get_num_gfn_evals(value vcvode_mem)
{
    CAMLparam1(vcvode_mem);

    int flag;
    long int v;

    flag = CVBBDPrecGetNumGfnEvals(CVODE_MEM_FROM_ML(vcvode_mem), &v);
    CHECK_FLAG("CVBBDPrecGetNumGfnEvals", flag);

    CAMLreturn(Val_long(v));
}
