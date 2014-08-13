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

static int bbbdlocal(long int nlocal, realtype t, N_Vector yy, N_Vector yp,
		     N_Vector yyB, N_Vector ypB, N_Vector glocal,
		     void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 7);
    int r;
    value *backref = user_data;
    CAML_FN (call_bbbdlocal);

    args[0] = *backref;
    args[1] = caml_copy_double(t);
    args[2] = NVEC_BACKLINK(yy);
    args[3] = NVEC_BACKLINK(yp);
    args[4] = NVEC_BACKLINK(yyB);
    args[5] = NVEC_BACKLINK(ypB);
    args[6] = NVEC_BACKLINK(glocal);

    r = Int_val (caml_callbackN(*call_bbbdlocal,
                                sizeof (args) / sizeof (*args),
                                args));

    CAMLreturnT(int, r);
}

static int bbbdcomm(long int nlocal, realtype t,
		    N_Vector yy, N_Vector yp,
		    N_Vector yyB, N_Vector ypB,
		    void *user_data)
{
    CAMLparam0();
    CAMLlocalN(args, 6);
    int r;
    value *backref = user_data;
    CAML_FN (call_bbbdcomm);

    args[0] = *backref;
    args[1] = caml_copy_double(t);
    args[2] = NVEC_BACKLINK(yy);
    args[3] = NVEC_BACKLINK(yp);
    args[4] = NVEC_BACKLINK(yyB);
    args[5] = NVEC_BACKLINK(ypB);

    r = Int_val (caml_callbackN(*call_bbbdcomm,
                                sizeof (args) / sizeof (*args),
                                args));

    CAMLreturnT(int, r);
}

CAMLprim void c_idas_bbd_prec_initb (value vparentwhich, value vlocaln,
				     value vbandwidths, value vdqrely,
				     value vhascomm)
{
    CAMLparam5(vparentwhich, vlocaln, vbandwidths, vdqrely, vhascomm);
    void *ida_mem = IDA_MEM_FROM_ML (Field(vparentwhich, 0));
    int flag;

    flag = IDABBDPrecInitB (ida_mem, Int_val(Field(vparentwhich, 1)),
	Long_val(vlocaln),
	Long_val(Field(vbandwidths, RECORD_IDA_BANDBLOCK_BANDWIDTHS_MUDQ)),
	Long_val(Field(vbandwidths, RECORD_IDA_BANDBLOCK_BANDWIDTHS_MLDQ)),
	Long_val(Field(vbandwidths, RECORD_IDA_BANDBLOCK_BANDWIDTHS_MUKEEP)),
	Long_val(Field(vbandwidths, RECORD_IDA_BANDBLOCK_BANDWIDTHS_MLKEEP)),
	Double_val(vdqrely),
	bbbdlocal,
	Bool_val(vhascomm) ? bbbdcomm : NULL);
    CHECK_FLAG ("IDABBDPrecInitB", flag);

    CAMLreturn0;
}

CAMLprim void c_idas_bbd_prec_reinitb (value vparent, value vwhich,
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

    CAMLreturn0;
}

