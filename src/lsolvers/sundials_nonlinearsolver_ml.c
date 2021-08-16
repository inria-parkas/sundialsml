/***********************************************************************
 *                                                                     *
 *                   OCaml interface to Sundials                       *
 *                                                                     *
 *             Timothy Bourke, Jun Inoue, and Marc Pouzet              *
 *             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *
 *                                                                     *
 *  Copyright 2020 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a New BSD License, refer to the file LICENSE.                *
 *                                                                     *
 ***********************************************************************/

#include "../config.h"

#include "../sundials/sundials_ml.h"
#include "../nvectors/nvector_ml.h"
#include "../lsolvers/sundials_nonlinearsolver_ml.h"

#if SUNDIALS_LIB_VERSION >= 400

#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

#define SYSFN_VAL(v) (*(SUNNonlinSolSysFn*)Data_custom_val(v))
#define LSETUPFN_VAL(v) (*(SUNNonlinSolLSetupFn*)Data_custom_val(v))
#define LSOLVEFN_VAL(v) (*(SUNNonlinSolLSolveFn*)Data_custom_val(v))
#define CONVTESTFN_VAL(v) (*(SUNNonlinSolConvTestFn*)Data_custom_val(v))

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Unfortunately sundials/sundials_nvector_senswrapper.h is not installed with
  the other header files, so the required internal definitions and macros are
  redeclared here. They must be kept synchronized with those in the source! */

struct _N_VectorContent_SensWrapper {
  N_Vector* vecs;
  int nvecs;
  booleantype own_vecs;
};

typedef struct _N_VectorContent_SensWrapper *N_VectorContent_SensWrapper;

#define NV_CONTENT_SW(v)  ( (N_VectorContent_SensWrapper)(v->content) )
#define NV_VECS_SW(v)     ( NV_CONTENT_SW(v)->vecs )
#define NV_NVECS_SW(v)    ( NV_CONTENT_SW(v)->nvecs )
#define NV_VEC_SW(v,i)    ( NV_VECS_SW(v)[i] )
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#endif

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>

#define MAX_ERRMSG_LEN 256

/* Nonlinear solvers (NLS) and callbacks from their "consumers".

   A consumer instantiates an NLS, configures callbacks (sysfn, lsetupfn,
   lsolvefn, convtestfn), and invokes init, setup, and solve.

   There are two types of "consumers": integrators (e.g., CVODE or IDA) and
   application code (written in OCaml).

   The interface code must consider four scenarios.

   1. C consumer (integrator) && C NLS && OCaml application (C/Cnls)
      'mem <> user (an internal value is provided by the integrator)

      - Callback configuration: sysfn, lsetupfn, lsolvefn
        Integrator only, occurs completely within C.
	Configuration from OCaml is not possible ('mem <> user) since the mem
	argument is controlled by the integrator and there is thus no way to
	link a C stub passed to the NLS with the OCaml callback.
	In any case, it is not clear that overriding the integrator's sys,
	setup, and solve function is very useful in practice.

      - Callback configuration: convtestfn
        Can be overridden from OCaml. In this case, the C stub passed to
	the NLS can use the SUNNonlinearSolver argument to resolve the
	OCaml callback (a "backlink" is guaranteed since the application
	has access to the NLS to configure the callback). The (opaque)
	'mem value is wrapped by the C stub before invocation of the OCaml
	code; it can later be used in calls to the wrapped C function
	retrieved through Newton/FixedPoint.get_sys_fn.

	NB: an NLS is not passed to the OCaml convtestfn since this would
	    require a cyclic "backlink" from the C-side NLS to its OCaml
	    counterpart. A closure should instead be used to provide the
	    NLS value necessary to invoke Newton/FixedPoint.get_sys_fn.
    
      - Callback invocation: sysfn, lsetupfn, lsolvefn
	Directly into C.

      - Callback invocation: convtestfn
        Either directly into C, or via the NLS "backlink" into OCaml.

      - Init/setup/solve invocation
        From C only. Not possible from OCaml ('mem <> C) since we would have
	to provide a suitable 'mem value for the callbacks into C.
	(This would not be impossible, but its utility is unclear.)

   2. C consumer (integrator) && OCaml (custom) NLS && OCaml application (C/Onls)
      'mem <> user (an internal value is provided by the integrator)

      - Callback configuration: sysfn, lsetupfn, lsolvefn
        The function pointers provided by the integrator are wrapped and
	passed to the NLS as OCaml functions. The wrappers take care of
	unwrapping the 'mem value passed from C to OCaml (and thus to C again).
	Technically, nothing would prevent an OCaml application from
	overriding these functions, but an extra phantom type argument would
	be required to distinguish this case from the previous one (where
	'mem comes from the C side, but where sysfn/lsetupfn/lsolvefn cannot
	be overridden by the application). The potential gain in utility does
	not seem to warrant the extra typing noise and complexity.

      - Callback configuration: convtestfn
        From C: wrapped as per sysfn, lsetupfn, lsolvefn.
        Can be overridden from OCaml, as per (1) but the "backlink" mechanism
	is unnecessary -- the callback occurs directly in OCaml without
	passing via C.

      - Callback invocation: 
	Occurs directly from OCaml.

      - Init/setup/solve invocation
        From C only. Not possible from OCaml ('mem <> C) since we would have
	to provide a suitable 'mem value for the callbacks into C.
	(This would not be impossible, but its utility is unclear.)

   3. OCaml consumer/application && C NLS (O/Cnls)
      'mem = user
      The internal value is provided by OCaml: we pass a pointer to the NLS's
      "backlink". The NLS must have been created from within OCaml (to be
      passed as an argument in callback configuration or init/setup/solve
      invocation) and thus must have a backlink. The backlink is stored on
      the C heap and already registered as a global root.

      - Callback configuration: sysfn, lsetupfn, lsolvefn, convtestfn
        From OCaml only: pass generic C stubs to the C NLS.

      - Callback invocation: sysfn, lsetupfn, lsolvefn
        The C stub uses the mem value (the backlink) to invoke the relevant
	OCaml closure.

      - Callback invocation: convtestfn
        The C stub uses the NLS value (with the same backlink) to invoke the
	relevant OCaml closure as per (1).

      - Init/setup/solve invocation
        From OCaml only: pass a pointer to the NLS backlink value through the
	C stub. The C NLS then transmits this value through to the callback
	C stubs.

   4. OCaml consumer/application && OCaml (custom) NLS (O/Onls)
      'mem = user
      The internal value is provided by OCaml but is unused.

      - Callback configuration: sysfn, lsetupfn, lsolvefn, convtestfn
        Store the OCaml closure directly.

      - Callback invocation: sysfn, lsetupfn, lsolvefn, convtestfn
        Invoke the stored OCaml closure directly.

      - Init/setup/solve invocation
        Call OCaml closure directly.

  Note that the only place where we really need the backlink is in the
  convtest callback (for the Cnls case), otherwise solve could just pass
  an OCaml value on the stack for the mem argument which is used in sysfn,
  lsetupfn, and solvefn (and potentially convtestfn) to callback into OCaml.

 */

/* - - - Exception handling - - - */

typedef enum {
    UNRECOVERABLE      = 0x0,
    RECOVERABLE        = 0x1,
    WITH_WARN	       = 0x2,
} recoverability;

#if 400 <= SUNDIALS_LIB_VERSION
static int translate_exception(value exn, recoverability recoverable)
{
    CAMLparam1(exn);

    if ((recoverable & 0x1)
	    && Field(exn, 0) == SUNDIALS_EXN_TAG (RecoverableFailure))
	CAMLreturnT (int, 1);

    if (recoverable & WITH_WARN)
	sunml_warn_discarded_exn (exn, "nonlinear solver");

    /* Unrecoverable error -1. Unfortunately we lose the exception. */
    CAMLreturnT (int, -1);
}

#define CHECK_NLS_EXCEPTION(result, recoverable)			      \
    (Is_exception_result (result)					      \
     ? translate_exception (result = Extract_exception (result), recoverable) \
     : 0)
#endif

CAMLprim void sunml_nlsolver_init_module (value exns)
{
    CAMLparam1 (exns);
    REGISTER_EXNS (NLSOLVER, exns);
    CAMLreturn0;
}

#if 400 <= SUNDIALS_LIB_VERSION
void sunml_nlsolver_check_flag(const char *call, int flag)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == SUN_NLS_SUCCESS
	    || flag == SUN_NLS_CONTINUE
	    || flag == SUN_NLS_CONV_RECVR) return;

    switch (flag) {
	case SUN_NLS_ILL_INPUT:
	    caml_invalid_argument(call);

	case SUN_NLS_VECTOROP_ERR:
	    caml_raise_constant(NLSOLVER_EXN(VectorOpError));

	case SUN_NLS_MEM_NULL:
	    caml_raise_constant(NLSOLVER_EXN(IncorrectUse));

	case SUN_NLS_MEM_FAIL:
	    caml_raise_out_of_memory();

#if 500 <= SUNDIALS_LIB_VERSION
	case SUN_NLS_EXT_FAIL:
	    caml_raise_constant(NLSOLVER_EXN(ExtFail));
#endif

	default:
	    if (flag > 0) {
		caml_raise_constant(SUNDIALS_EXN(RecoverableFailure));
	    } else {
		snprintf(exmsg, MAX_ERRMSG_LEN, "%s: unexpected error", call);
		caml_failwith(exmsg);
	    }
    }
}
#endif

/* - - - Senswrappers - - - */

#if 400 <= SUNDIALS_LIB_VERSION

#define SENSWRAPPER(v) (*(N_Vector *)Data_custom_val(v))

static intnat senswrapper_hash(value vsw)
{
    // treat the senswrapper pointer directly as a hash value
    return ((intnat)SENSWRAPPER(vsw));
}

static int senswrapper_compare(value vsw1, value vsw2)
{
    // just compare the pointer addresses
    return ((long)SENSWRAPPER(vsw1) - (long)SENSWRAPPER(vsw2));
}

static struct custom_operations senswrapper_custom_ops = {
    identifier: "sunml_senswrapper",
    finalize:    custom_finalize_default,
    compare:     senswrapper_compare,
    hash:        senswrapper_hash,
    serialize:   custom_serialize_default,
    deserialize: custom_deserialize_default
};

CAMLprim value sunml_senswrapper_wrap(N_Vector sw)
{
    CAMLparam0();
    CAMLlocal1(vsw);

    vsw = caml_alloc_custom(&senswrapper_custom_ops, 1, 0, 1);
    SENSWRAPPER(vsw) = sw;
    
    CAMLreturn(vsw);
}

static void invalidate_senswrapper(value vsw)
{
    SENSWRAPPER(vsw) = NULL;
}
#endif

CAMLprim value sunml_senswrapper_data(value vsw)
{
    CAMLparam1(vsw);
    CAMLlocal2(vr, vnv);
#if 400 <= SUNDIALS_LIB_VERSION
    N_Vector sw = SENSWRAPPER(vsw);
    int nvecs, i;

    if (sw == NULL) {
	caml_raise_constant(NLSOLVER_EXN(IncorrectUse));
    } else {
	nvecs = NV_NVECS_SW(sw);
	vr = caml_alloc_tuple(nvecs);
	for (i = 0; i < nvecs; ++i) {
	    vnv = NVEC_BACKLINK(NV_VEC_SW(sw, i));
	    Store_field(vr, i, vnv);
	}
    }
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(vr);
}


/* - - - O/Cnls functions (OCaml invoking init/setup/solve in a C NLS) - - - */

CAMLprim value sunml_nlsolver_init(value vnls)
{
    CAMLparam1(vnls);
#if SUNDIALS_LIB_VERSION >= 400
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);

    int flag = SUNNonlinSolInitialize(nls);
    NLS_CHECK_FLAG("SUNNonlinSolInitialize", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_setup(value vnls, value vy)
{
    CAMLparam2(vnls, vy);
#if SUNDIALS_LIB_VERSION >= 400
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    N_Vector y  = NVEC_VAL(vy);
    void *mem = &(NLS_CALLBACKS(nls));

    int flag = SUNNonlinSolSetup(nls, y, mem);
    NLS_CHECK_FLAG("SUNNonlinSolSetup", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_solve(value vnls, value vargs, value vmem)
{
    CAMLparam3(vnls, vargs, vmem);
#if SUNDIALS_LIB_VERSION >= 400
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    N_Vector y0 = NVEC_VAL(Field(vargs, 0));
    N_Vector y  = NVEC_VAL(Field(vargs, 1));
    N_Vector w  = NVEC_VAL(Field(vargs, 2));
    double tol  = Double_val(Field(vargs, 3));
    int vcallsetup = Bool_val(Field(vargs, 4));
    void *mem = &(NLS_CALLBACKS(nls));

    int flag = SUNNonlinSolSolve(nls, y0, y, w, tol, vcallsetup, mem);
    NLS_CHECK_FLAG("SUNNonlinSolSolve", flag);
#endif
    CAMLreturn (Val_unit);
}

/* - - - O/Cnls: callback invocation (also C/Cnls for convtestfn) - - - */

#if 400 <= SUNDIALS_LIB_VERSION
static int sys_callback(N_Vector y, N_Vector F, void* mem)
{
    CAMLparam0();
    CAMLlocalN(args, 3);
    CAMLlocal1(callbacks);

    args[0] = NVEC_BACKLINK(y);
    args[1] = NVEC_BACKLINK(F);
    args[2] = Val_unit;
    callbacks = *((value *)mem);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (
	Field(callbacks, RECORD_NLSOLVER_CALLBACKS_SYSFN), 3, args);
    if (!Is_exception_result (r)) CAMLreturnT(int, 0);

    r = Extract_exception (r);
    CAMLreturnT(int, translate_exception (r, RECOVERABLE));
}
#endif

#if 400 <= SUNDIALS_LIB_VERSION
#if 500 <= SUNDIALS_LIB_VERSION
static int lsetup_callback(booleantype jbad, booleantype* jcur, void* mem)
#else
static int lsetup_callback(N_Vector y, N_Vector F, booleantype jbad,
			   booleantype* jcur, void* mem)
#endif
{
    CAMLparam0();
    CAMLlocalN(args, 2);
    CAMLlocal1(callbacks);

    args[0] = Val_bool(jbad);
    args[1] = Val_unit;
    callbacks = *((value *)mem);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (
	Field(callbacks, RECORD_NLSOLVER_CALLBACKS_LSETUPFN), 2, args);

    /* Update jcur; leave it unchanged if an error occurred. */
    if (!Is_exception_result (r)) {
	*jcur = Bool_val (r);
	CAMLreturnT(int, 0);
    }

    r = Extract_exception (r);
    CAMLreturnT(int, translate_exception (r, RECOVERABLE));
}
#endif

#if 400 <= SUNDIALS_LIB_VERSION
#if 500 <= SUNDIALS_LIB_VERSION
static int lsolve_callback(N_Vector b, void* mem)
#else
static int lsolve_callback(N_Vector y, N_Vector b, void* mem)
#endif
{
    CAMLparam0();
    CAMLlocalN(args, 2);
    CAMLlocal1(callbacks);

    args[0] = NVEC_BACKLINK(b);
    args[1] = Val_unit;
    callbacks = *((value *)mem);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (
	Field(callbacks, RECORD_NLSOLVER_CALLBACKS_LSOLVEFN), 2, args);

    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}
#endif

#if 400 <= SUNDIALS_LIB_VERSION
static int convtest_callback(SUNNonlinearSolver nls, N_Vector y, N_Vector del,
			     realtype tol, N_Vector ewt, void* mem)
{
    CAMLparam0();
    CAMLlocalN(args, 5);
    CAMLlocal1(vcallbacks);

    vcallbacks = NLS_CALLBACKS(nls);

    args[0] = NVEC_BACKLINK(y);
    args[1] = NVEC_BACKLINK(del);
    args[2] = caml_copy_double(tol);
    args[3] = NVEC_BACKLINK(ewt);

    /* - For C/Cnls with OCaml convtest override ('mem <> user): we wrap the
         value passed from the C integrator before passing it into OCaml. The
	 only way this value can be used is in a later call to a sysfn
	 retrieved via FixedPoint/Newton.get_sys_fn; the original value is
	 recovered at that point
	 (This is true in Sundials >= 500 since all internal calls to
	  SUNNonlinSolSetConvTestFn pass a pointer to the integrator data
	  structure.)

	-For O/Cnls, mem is NULL in Sundials >=500 (see
	 sunml_nlsolver_set_convtest_fn), and otherwise a pointer to the
	 integrator data structure.
	 There is thus no need to wrap it, but we have no way to
	 distinguish the two cases at this point. In this case, the given
	 value (of type "user") cannot be manipulated in any way from OCaml,
	 so wrapping it again is not problematic (although it is slightly
	 inefficient). */
    args[4] = caml_alloc_final(1, NULL, 0, 1);
    NLMEM_VAL(args[4]) = mem;

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (
	Field(vcallbacks, RECORD_NLSOLVER_CALLBACKS_CONVTESTFN), 5, args);

    if (!Is_exception_result (r)) {
	switch (Int_val(r)) {
	case VARIANT_NLSOLVER_CONVTEST_SUCCESS:
	    CAMLreturnT(int, SUN_NLS_SUCCESS);

	case VARIANT_NLSOLVER_CONVTEST_CONTINUE:
	    CAMLreturnT(int, SUN_NLS_CONTINUE);

	case VARIANT_NLSOLVER_CONVTEST_RECOVER:
	    CAMLreturnT(int, SUN_NLS_CONV_RECVR);

	default:
	    CAMLreturnT(int, -1);
	}
    }

    r = Extract_exception (r);
    CAMLreturnT(int, translate_exception (r, UNRECOVERABLE));
}
#endif

/* - - - O/Cnls: callback configuration - - - */

CAMLprim value sunml_nlsolver_set_sys_fn(value vnls)
{
    CAMLparam1(vnls);
#if SUNDIALS_LIB_VERSION >= 400
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);

    int flag = SUNNonlinSolSetSysFn(nls, sys_callback);
    NLS_CHECK_FLAG("SUNNonlinSolSetSysFn", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_set_lsetup_fn(value vnls)
{
    CAMLparam1(vnls);
#if SUNDIALS_LIB_VERSION >= 400
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);

    int flag = SUNNonlinSolSetLSetupFn(nls, lsetup_callback);
    NLS_CHECK_FLAG("SUNNonlinSolSetLSetupFn", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_set_lsolve_fn(value vnls)
{
    CAMLparam1(vnls);
#if SUNDIALS_LIB_VERSION >= 400
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);

    int flag = SUNNonlinSolSetLSolveFn(nls, lsolve_callback);
    NLS_CHECK_FLAG("SUNNonlinSolSetLSolveFn", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_set_convtest_fn(value vnls)
{
    CAMLparam1(vnls);
#if SUNDIALS_LIB_VERSION >= 400
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);

#if SUNDIALS_LIB_VERSION >= 500
    int flag = SUNNonlinSolSetConvTestFn(nls, convtest_callback, NULL);
#else
    int flag = SUNNonlinSolSetConvTestFn(nls, convtest_callback);
#endif
    NLS_CHECK_FLAG("SUNNonlinSolSetConvTestFn", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_set_max_iters(value vnls, value vi)
{
    CAMLparam2(vnls, vi);
#if SUNDIALS_LIB_VERSION >= 400
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);

    int flag = SUNNonlinSolSetMaxIters(nls, Int_val(vi));
    NLS_CHECK_FLAG("SUNNonlinSolSetMaxIters", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_get_num_iters(value vnls)
{
    CAMLparam1(vnls);
    long int niters = 0;
#if SUNDIALS_LIB_VERSION >= 400
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    int flag = SUNNonlinSolGetNumIters(nls, &niters);
    NLS_CHECK_FLAG("SUNNonlinSolGetNumIters", flag);
#endif
    CAMLreturn (Val_long(niters));
}

CAMLprim value sunml_nlsolver_get_cur_iter(value vnls)
{
    CAMLparam1(vnls);
    int iter = 0;
#if SUNDIALS_LIB_VERSION >= 400
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    int flag = SUNNonlinSolGetCurIter(nls, &iter);
    NLS_CHECK_FLAG("SUNNonlinSolGetCurIter", flag);
#endif
    CAMLreturn (Val_int(iter));
}

CAMLprim value sunml_nlsolver_get_num_conv_fails(value vnls)
{
    CAMLparam1(vnls);
    long int nconvfails = 0;
#if SUNDIALS_LIB_VERSION >= 400
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    int flag = SUNNonlinSolGetNumConvFails(nls, &nconvfails);
    NLS_CHECK_FLAG("SUNNonlinSolGetNumConvFails", flag);
#endif
    CAMLreturn (Val_long(nconvfails));
}

CAMLprim void sunml_nlsolver_set_info_file_newton(value vnls, value vfile)
{
    CAMLparam2(vnls, vfile);
#if 530 <= SUNDIALS_LIB_VERSION
    int flag = SUNNonlinSolSetInfoFile_Newton(NLSOLVER_VAL(vnls),
					      ML_CFILE(vfile));
    NLS_CHECK_FLAG("SUNNonlinSolSetInfofile_Newton", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn0;
}

CAMLprim void sunml_nlsolver_set_info_file_fixedpoint(value vnls, value vfile)
{
    CAMLparam2(vnls, vfile);
#if 530 <= SUNDIALS_LIB_VERSION
    int flag = SUNNonlinSolSetInfoFile_FixedPoint(NLSOLVER_VAL(vnls),
						  ML_CFILE(vfile));
    NLS_CHECK_FLAG("SUNNonlinSolSetInfofile_FixedPoint", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn0;
}

CAMLprim void sunml_nlsolver_set_print_level_newton(value vnls, value vlevel)
{
    CAMLparam2(vnls, vlevel);
#if 530 <= SUNDIALS_LIB_VERSION
    int flag = SUNNonlinSolSetPrintLevel_Newton(NLSOLVER_VAL(vnls),
					        Int_val(vlevel));
    NLS_CHECK_FLAG("SUNNonlinSolSetPrintLevel_Newton", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn0;
}

CAMLprim void sunml_nlsolver_set_print_level_fixedpoint(value vnls, value vlevel)
{
    CAMLparam2(vnls, vlevel);
#if 530 <= SUNDIALS_LIB_VERSION
    int flag = SUNNonlinSolSetPrintLevel_FixedPoint(NLSOLVER_VAL(vnls),
						    Int_val(vlevel));
    NLS_CHECK_FLAG("SUNNonlinSolSetPrintLevel_FixedPoint", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn0;
}

#if SUNDIALS_LIB_VERSION >= 400
static void finalize_nls(value vnls)
{
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    caml_remove_generational_global_root(&NLS_CALLBACKS(nls));
    SUNNonlinSolFree(nls);
}
#endif

CAMLprim value sunml_nlsolver_call_sys_fn(value vsysfn, value vy,
					  value vfg, value vmem)
{
    CAMLparam4(vsysfn, vy, vfg, vmem);
#if SUNDIALS_LIB_VERSION >= 400
    SUNNonlinSolSysFn sysfn = SYSFN_VAL(vsysfn);
    N_Vector y  = NVEC_VAL(vy);
    N_Vector fg = NVEC_VAL(vfg);
    void *mem = NLMEM_VAL(vmem);

    int flag = (*sysfn)(y, fg, mem);
    NLS_CHECK_FLAG("SUNNonlinSolSysFn", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_call_sys_fn_sens(value vsysfn, value vy,
					       value vfg, value vmem)
{
    CAMLparam4(vsysfn, vy, vfg, vmem);
#if SUNDIALS_LIB_VERSION >= 400
    SUNNonlinSolSysFn sysfn = SYSFN_VAL(vsysfn);
    N_Vector y  = SENSWRAPPER(vy);
    N_Vector fg = SENSWRAPPER(vfg);
    void *mem = NLMEM_VAL(vmem);

#if SUNDIALS_ML_SAFE == 1
    if ((y == NULL) || (fg == NULL))
	caml_raise_constant(NLSOLVER_EXN(IncorrectUse));
#endif

    int flag = (*sysfn)(y, fg, mem);
    NLS_CHECK_FLAG("SUNNonlinSolSysFn (sens)", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_call_lsetup_fn(value vlsetupfn,
					     value vjbad, value vmem)
{
    CAMLparam3(vlsetupfn, vjbad, vmem);
#if SUNDIALS_LIB_VERSION >= 400
    SUNNonlinSolLSetupFn lsetupfn = LSETUPFN_VAL(vlsetupfn);
    void *mem = NLMEM_VAL(vmem);
    int jcur = 0;

#if SUNDIALS_LIB_VERSION >= 500
    int flag = (*lsetupfn)(Bool_val(vjbad), &jcur, mem);
#else
    int flag = (*lsetupfn)(NULL, NULL, Bool_val(vjbad), &jcur, mem);
#endif
    NLS_CHECK_FLAG("SUNNonlinSolLSetupFn", flag);

    CAMLreturn (Val_bool(jcur));
#endif
    CAMLreturn (Val_bool(0));
}

CAMLprim value sunml_nlsolver_call_lsetup_fn_sens(value vlsetupfn,
					          value vjbad, value vmem)
{
    CAMLparam3(vlsetupfn, vjbad, vmem);
#if SUNDIALS_LIB_VERSION >= 400
    SUNNonlinSolLSetupFn lsetupfn = LSETUPFN_VAL(vlsetupfn);
    void *mem = NLMEM_VAL(vmem);
    int jcur = 0;

#if SUNDIALS_LIB_VERSION >= 500
    int flag = (*lsetupfn)(Bool_val(vjbad), &jcur, mem);
#else
    int flag = (*lsetupfn)(NULL, NULL, Bool_val(vjbad), &jcur, mem);
#endif
    NLS_CHECK_FLAG("SUNNonlinSolLSetupFn (sens)", flag);

    CAMLreturn (Val_bool(jcur));
#endif
    CAMLreturn (Val_bool(0));
}

CAMLprim value sunml_nlsolver_call_lsolve_fn(value vlsolvefn,
					     value vb, value vmem)
{
    CAMLparam3(vlsolvefn, vb, vmem);
#if SUNDIALS_LIB_VERSION >= 400
    SUNNonlinSolLSolveFn lsolvefn = LSOLVEFN_VAL(vlsolvefn);
    N_Vector b = NVEC_VAL(vb);
    void *mem = NLMEM_VAL(vmem);

#if SUNDIALS_LIB_VERSION >= 500
    int flag = (*lsolvefn)(b, mem);
#else
    int flag = (*lsolvefn)(NULL, b, mem);
#endif
    NLS_CHECK_FLAG("SUNNonlinSolLSolveFn", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_call_lsolve_fn_sens(value vlsolvefn,
						  value vb, value vmem)
{
    CAMLparam3(vlsolvefn, vb, vmem);
#if SUNDIALS_LIB_VERSION >= 400
    SUNNonlinSolLSolveFn lsolvefn = LSOLVEFN_VAL(vlsolvefn);
    N_Vector b = SENSWRAPPER(vb);
    void *mem = NLMEM_VAL(vmem);

#if SUNDIALS_ML_SAFE == 1
    if (b == NULL) caml_raise_constant(NLSOLVER_EXN(IncorrectUse));
#endif

#if SUNDIALS_LIB_VERSION >= 500
    int flag = (*lsolvefn)(b, mem);
#else
    int flag = (*lsolvefn)(NULL, b, mem);
#endif
    NLS_CHECK_FLAG("SUNNonlinSolLSolveFn (sens)", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_call_convtest_fn(value vnls, value vconvtestfn,
					       value vargs)
{
    CAMLparam3(vnls, vconvtestfn, vargs);
#if SUNDIALS_LIB_VERSION >= 400
    SUNNonlinSolConvTestFn convtestfn = CONVTESTFN_VAL(vconvtestfn);

    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    N_Vector y   = NVEC_VAL(Field(vargs, 0));
    N_Vector del = NVEC_VAL(Field(vargs, 1));
    realtype tol = Double_val(Field(vargs, 2));
    N_Vector ewt = NVEC_VAL(Field(vargs, 3));
    void *mem    = NLMEM_VAL(Field(vargs, 4));

    int flag = (*convtestfn)(nls, y, del, tol, ewt, mem);
    switch (flag) {
    case SUN_NLS_SUCCESS:
	CAMLreturn (VARIANT_NLSOLVER_CONVTEST_SUCCESS);

    case SUN_NLS_CONTINUE:
	CAMLreturn (VARIANT_NLSOLVER_CONVTEST_CONTINUE);

    case SUN_NLS_CONV_RECVR:
	CAMLreturn (VARIANT_NLSOLVER_CONVTEST_RECOVER);

    default:
	sunml_nlsolver_check_flag("SUNNonlinSolConvTestFn", flag);
    }
#endif
    CAMLreturn (VARIANT_NLSOLVER_CONVTEST_SUCCESS);
}

CAMLprim value sunml_nlsolver_call_convtest_fn_sens(value vnls,
						value vconvtestfn, value vargs)
{
    CAMLparam3(vnls, vconvtestfn, vargs);
#if SUNDIALS_LIB_VERSION >= 400
    SUNNonlinSolConvTestFn convtestfn = CONVTESTFN_VAL(vconvtestfn);

    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    N_Vector y   = SENSWRAPPER(Field(vargs, 0));
    N_Vector del = SENSWRAPPER(Field(vargs, 1));
    realtype tol =  Double_val(Field(vargs, 2));
    N_Vector ewt = SENSWRAPPER(Field(vargs, 3));
    void *mem    =   NLMEM_VAL(Field(vargs, 4));

#if SUNDIALS_ML_SAFE == 1
    if ((y == NULL) || (del == NULL) || (ewt == NULL))
	caml_raise_constant(NLSOLVER_EXN(IncorrectUse));
#endif

    int flag = (*convtestfn)(nls, y, del, tol, ewt, mem);
    switch (flag) {
    case SUN_NLS_SUCCESS:
	CAMLreturn (VARIANT_NLSOLVER_CONVTEST_SUCCESS);

    case SUN_NLS_CONTINUE:
	CAMLreturn (VARIANT_NLSOLVER_CONVTEST_CONTINUE);

    case SUN_NLS_CONV_RECVR:
	CAMLreturn (VARIANT_NLSOLVER_CONVTEST_RECOVER);

    default:
	sunml_nlsolver_check_flag("SUNNonlinSolConvTestFn (sens)", flag);
    }
#endif
    CAMLreturn (VARIANT_NLSOLVER_CONVTEST_SUCCESS);
}

#if SUNDIALS_LIB_VERSION >= 400
static value rewrap_nlsolver(SUNNonlinearSolver nls0, value vcallbacks)
{
    CAMLparam1(vcallbacks);
    CAMLlocal1(vr);
    SUNNonlinearSolver nls;

    if (nls0 == NULL) caml_raise_out_of_memory();

    nls = (SUNNonlinearSolver)malloc(sizeof(struct sunml_nls));
    if (nls == NULL) {
	SUNNonlinSolFree(nls0);
	caml_raise_out_of_memory();
    }

    nls->content = nls0->content;
    nls->ops = nls0->ops;
    free(nls0);

    NLS_CALLBACKS(nls) = vcallbacks;
    caml_register_generational_global_root(&NLS_CALLBACKS(nls));

    // Setup the OCaml-side
    vr = caml_alloc_final(1, &finalize_nls, 1, 20);
    NLSOLVER_VAL(vr) = nls;

    CAMLreturn (vr);
}
#endif

CAMLprim value sunml_nlsolver_newton_make(value vy, value vcallbacks)
{
    CAMLparam2(vy, vcallbacks);
#if SUNDIALS_LIB_VERSION >= 400
    CAMLlocal1(vr);
    N_Vector y = NVEC_VAL(vy);

    vr = rewrap_nlsolver(
	    SUNNonlinSol_Newton(y),
	    vcallbacks);
    CAMLreturn (vr);
#else
    CAMLreturn (Val_unit);
#endif
}

CAMLprim value sunml_nlsolver_newton_make_sens(value vcount, value vy,
					       value vcallbacks)
{
    CAMLparam3(vcount, vy, vcallbacks);
#if SUNDIALS_LIB_VERSION >= 400
    CAMLlocal1(vr);
    N_Vector y = NVEC_VAL(vy);

    vr = rewrap_nlsolver(
	    SUNNonlinSol_NewtonSens(Int_val(vcount), y),
	    vcallbacks);
    CAMLreturn (vr);
#else
    CAMLreturn (Val_unit);
#endif
}

CAMLprim value sunml_nlsolver_newton_get_sys_fn(value vnls)
{
    CAMLparam1(vnls);
    CAMLlocal2(vr, vsysfn);
#if SUNDIALS_LIB_VERSION >= 400
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    SUNNonlinSolSysFn sysfn;

    int flag = SUNNonlinSolGetSysFn_Newton(nls, &sysfn);
    NLS_CHECK_FLAG("SUNNonlinSolGetSysFn_Newton", flag);

    if (sysfn == NULL) {
	vr = Val_none;
    } else {
	vsysfn = caml_alloc_final(1, NULL, 0, 1);
	SYSFN_VAL(vsysfn) = sysfn;

	Store_some(vr, vsysfn);
    }

    CAMLreturn (vr);
#else
    CAMLreturn (Val_unit);
#endif
}

CAMLprim value sunml_nlsolver_fixedpoint_make(value vy, value vm,
					      value vcallbacks)
{
    CAMLparam3(vy, vm, vcallbacks);
#if SUNDIALS_LIB_VERSION >= 400
    CAMLlocal1(vr);
    N_Vector y = NVEC_VAL(vy);

    vr = rewrap_nlsolver(
	    SUNNonlinSol_FixedPoint(y, Int_val(vm)),
	    vcallbacks);
    CAMLreturn (vr);
#else
    CAMLreturn (Val_unit);
#endif
}

CAMLprim value sunml_nlsolver_fixedpoint_make_sens(value vcount, value vy,
						   value vm, value vcallbacks)
{
    CAMLparam4(vcount, vy, vm, vcallbacks);
#if SUNDIALS_LIB_VERSION >= 400
    CAMLlocal1(vr);
    N_Vector y = NVEC_VAL(vy);

    vr = rewrap_nlsolver(
	    SUNNonlinSol_FixedPointSens(Int_val(vcount), y, Int_val(vm)),
	    vcallbacks);
    CAMLreturn (vr);
#else
    CAMLreturn (Val_unit);
#endif
}

CAMLprim value sunml_nlsolver_fixedpoint_set_damping(value vnls, value vbeta)
{
    CAMLparam2(vnls, vbeta);

#if SUNDIALS_LIB_VERSION >= 510
    int flag = SUNNonlinSolSetDamping_FixedPoint(NLSOLVER_VAL(vnls),
						 Double_val(vbeta));
    NLS_CHECK_FLAG("SUNNonlinSolSetDamping_FixedPoint", flag);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif

    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_fixedpoint_get_sys_fn(value vnls)
{
    CAMLparam1(vnls);
    CAMLlocal2(vr, vsysfn);
#if SUNDIALS_LIB_VERSION >= 400
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    SUNNonlinSolSysFn sysfn;

    int flag = SUNNonlinSolGetSysFn_FixedPoint(nls, &sysfn);
    NLS_CHECK_FLAG("SUNNonlinSolGetSysFn_FixedPoint", flag);

    if (sysfn == NULL) {
	vr = Val_none;
    } else {
	vsysfn = caml_alloc_final(1, NULL, 0, 1);
	SYSFN_VAL(vsysfn) = sysfn;

	Store_some(vr, vsysfn);
    }

    CAMLreturn (vr);
#else
    CAMLreturn (Val_unit);
#endif
}

/* - - - O/Onls && C/Onls: custom nonlinear solvers - - - */

#if SUNDIALS_LIB_VERSION >= 400

#define NLSOLV_OP_TABLE(nls)  ((value)((nls)->content))
#define GET_OP(nls, x) (Field((value)NLSOLV_OP_TABLE(nls), RECORD_NLSOLVER_OPS_ ## x))
#define GET_SOME_OP(nls, x) (Some_val (GET_OP(nls, x)))

static SUNNonlinearSolver_Type callml_custom_gettype(SUNNonlinearSolver nls)
{
    CAMLparam0();
    SUNNonlinearSolver_Type r = SUNNONLINEARSOLVER_ROOTFIND;

    switch (Int_val(GET_OP(nls, TYPE))) {
    case VARIANT_NLSOLVER_TYPE_ROOT_FIND:
	break;

    case VARIANT_NLSOLVER_TYPE_FIXED_POINT:
	r = SUNNONLINEARSOLVER_FIXEDPOINT;
	break;
    }

    CAMLreturnT(SUNNonlinearSolver_Type, r);
}

static int callml_custom_initialize(SUNNonlinearSolver nls)
{
    CAMLparam0();
    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn(GET_SOME_OP(nls, INIT), Val_unit);
    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static int callml_custom_setup(SUNNonlinearSolver nls, N_Vector y, void* mem)
{
    CAMLparam0();
    CAMLlocalN(args, 2);

    args[0] = NVEC_BACKLINK(y);
    args[1] = caml_alloc_final(1, NULL, 0, 1);
    NLMEM_VAL(args[1]) = mem;

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(GET_SOME_OP(nls, SETUP), 2, args);
    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static int callml_custom_setup_sens(SUNNonlinearSolver nls, N_Vector y, void* mem)
{
    CAMLparam0();
    CAMLlocalN(args, 2);

    args[0] = sunml_senswrapper_wrap(y);
    args[1] = caml_alloc_final(1, NULL, 0, 1);
    NLMEM_VAL(args[1]) = mem;

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(GET_SOME_OP(nls, SETUP), 2, args);

    invalidate_senswrapper(args[0]);

    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static int callml_custom_solve(SUNNonlinearSolver nls,
			       N_Vector y0, N_Vector y, N_Vector w,
			       realtype tol, booleantype callLSetup, void* mem)
{
    CAMLparam0();
    CAMLlocalN(args, 6);

    args[0] = NVEC_BACKLINK(y0);
    args[1] = NVEC_BACKLINK(y);
    args[2] = NVEC_BACKLINK(w);
    args[3] = caml_copy_double(tol);
    args[4] = Val_bool(callLSetup);
    args[5] = caml_alloc_final(1, NULL, 0, 1);
    NLMEM_VAL(args[5]) = mem;

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(GET_OP(nls, SOLVE), 6, args);
    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, RECOVERABLE));
}

static int callml_custom_solve_sens(SUNNonlinearSolver nls,
				    N_Vector y0, N_Vector y, N_Vector w,
				    realtype tol, booleantype callLSetup,
				    void* mem)
{
    CAMLparam0();
    CAMLlocalN(args, 6);

    args[0] = sunml_senswrapper_wrap(y0);
    args[1] = sunml_senswrapper_wrap(y);
    args[2] = sunml_senswrapper_wrap(w);
    args[3] = caml_copy_double(tol);
    args[4] = Val_bool(callLSetup);
    args[5] = caml_alloc_final(1, NULL, 0, 1);
    NLMEM_VAL(args[5]) = mem;

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(GET_OP(nls, SOLVE), 6, args);

    invalidate_senswrapper(args[0]);
    invalidate_senswrapper(args[1]);
    invalidate_senswrapper(args[2]);

    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, RECOVERABLE));
}

static int callml_custom_free(SUNNonlinearSolver nls)
{
    if (nls == NULL) return(SUN_NLS_SUCCESS);

    nls->content = NULL;
    caml_remove_generational_global_root((value *)&(nls->content));

    if (nls->ops) {
	free(nls->ops);
	nls->ops = NULL;
    }

    free(nls);

    return(SUN_NLS_SUCCESS);
}

static int callml_custom_setsysfn(SUNNonlinearSolver nls,
				  SUNNonlinSolSysFn sysfn)
{
    CAMLparam0();
    CAMLlocal1(vsysfn);
    vsysfn = caml_alloc_final(1, NULL, 0, 1);
    SYSFN_VAL(vsysfn) = sysfn;

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn(
		*caml_named_value("Sundials_NonlinearSolver.set_c_sys_fn"),
		NLSOLV_OP_TABLE(nls),
		vsysfn);
    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static int callml_custom_setsysfn_sens(SUNNonlinearSolver nls,
				       SUNNonlinSolSysFn sysfn)
{
    CAMLparam0();
    CAMLlocal1(vsysfn);
    vsysfn = caml_alloc_final(1, NULL, 0, 1);
    SYSFN_VAL(vsysfn) = sysfn;

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn(
		*caml_named_value("Sundials_NonlinearSolver.set_c_sys_fn_sens"),
		NLSOLV_OP_TABLE(nls),
		vsysfn);
    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static int callml_custom_setlsetupfn(SUNNonlinearSolver nls,
				     SUNNonlinSolLSetupFn lsetupfn)
{
    CAMLparam0();
    CAMLlocal1(vlsetupfn);
    vlsetupfn = caml_alloc_final(1, NULL, 0, 1);
    LSETUPFN_VAL(vlsetupfn) = lsetupfn;

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn(
		*caml_named_value("Sundials_NonlinearSolver.set_c_lsetup_fn"),
		NLSOLV_OP_TABLE(nls),
		vlsetupfn);
    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static int callml_custom_setlsetupfn_sens(SUNNonlinearSolver nls,
					  SUNNonlinSolLSetupFn lsetupfn)
{
    CAMLparam0();
    CAMLlocal1(vlsetupfn);
    vlsetupfn = caml_alloc_final(1, NULL, 0, 1);
    LSETUPFN_VAL(vlsetupfn) = lsetupfn;

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn(
		*caml_named_value("Sundials_NonlinearSolver.set_c_lsetup_fn_sens"),
		NLSOLV_OP_TABLE(nls),
		vlsetupfn);
    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static int callml_custom_setlsolvefn(SUNNonlinearSolver nls,
				     SUNNonlinSolLSolveFn lsolvefn)
{
    CAMLparam0();
    CAMLlocal1(vlsolvefn);
    vlsolvefn = caml_alloc_final(1, NULL, 0, 1);
    LSOLVEFN_VAL(vlsolvefn) = lsolvefn;

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn(
		*caml_named_value("Sundials_NonlinearSolver.set_c_lsolve_fn"),
		NLSOLV_OP_TABLE(nls),
		vlsolvefn);
    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static int callml_custom_setlsolvefn_sens(SUNNonlinearSolver nls,
				          SUNNonlinSolLSolveFn lsolvefn)
{
    CAMLparam0();
    CAMLlocal1(vlsolvefn);
    vlsolvefn = caml_alloc_final(1, NULL, 0, 1);
    LSOLVEFN_VAL(vlsolvefn) = lsolvefn;

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn(
		*caml_named_value("Sundials_NonlinearSolver.set_c_lsolve_fn_sens"),
		NLSOLV_OP_TABLE(nls),
		vlsolvefn);
    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static int callml_custom_setconvtestfn(SUNNonlinearSolver nls,
				       SUNNonlinSolConvTestFn convtestfn,
				       void *cdata)
{
    CAMLparam0();
    CAMLlocal2(vconvtestfn, vnls);

    vconvtestfn = caml_alloc_final(1, NULL, 0, 1);
    CONVTESTFN_VAL(vconvtestfn) = convtestfn;

    vnls = caml_alloc_final(1, NULL, 0, 1); // no finalizer, just pass pointer
    NLSOLVER_VAL(vnls) = nls;

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn(
		*caml_named_value("Sundials_NonlinearSolver.set_c_convtest_fn"),
		vnls,
		NLSOLV_OP_TABLE(nls),
		vconvtestfn);
    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static int callml_custom_setconvtestfn_sens(SUNNonlinearSolver nls,
				            SUNNonlinSolConvTestFn convtestfn,
					    void *cdata)
{
    CAMLparam0();
    CAMLlocal2(vconvtestfn, vnls);

    vconvtestfn = caml_alloc_final(1, NULL, 0, 1);
    CONVTESTFN_VAL(vconvtestfn) = convtestfn;

    vnls = caml_alloc_final(1, NULL, 0, 1); // no finalizer, just pass pointer
    NLSOLVER_VAL(vnls) = nls;

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn(
		*caml_named_value("Sundials_NonlinearSolver.set_c_convtest_fn_sens"),
		vnls,
		NLSOLV_OP_TABLE(nls),
		vconvtestfn);
    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static int callml_custom_setmaxiters(SUNNonlinearSolver nls, int maxiters)
{
    CAMLparam0();
    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn(GET_SOME_OP(nls, SET_MAX_ITERS),
				Val_int(maxiters));
    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static int callml_custom_getnumiters(SUNNonlinearSolver nls, long int *niters)
{
    CAMLparam0();
    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn(GET_SOME_OP(nls, GET_NUM_ITERS), Val_unit);

    /* Update niters; leave it unchanged if an error occurred. */
    if (!Is_exception_result (r)) {
	*niters = Long_val (r);
	CAMLreturnT(int, 0);
    }

    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static int callml_custom_getcuriter(SUNNonlinearSolver nls, int *iter)
{
    CAMLparam0();
    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn(GET_SOME_OP(nls, GET_CUR_ITER), Val_unit);

    /* Update niters; leave it unchanged if an error occurred. */
    if (!Is_exception_result (r)) {
	*iter = Int_val (r);
	CAMLreturnT(int, 0);
    }

    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static int callml_custom_getnumconvfails(SUNNonlinearSolver nls,
				 	 long int *nconvfails)
{
    CAMLparam0();
    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn(GET_SOME_OP(nls, GET_NUM_CONV_FAILS), Val_unit);

    /* Update niters; leave it unchanged if an error occurred. */
    if (!Is_exception_result (r)) {
	*nconvfails = Long_val (r);
	CAMLreturnT(int, 0);
    }

    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static CAMLprim value custom_make(int sens, value vcallbacks, value vops)
{
    CAMLparam2(vcallbacks, vops);
    CAMLlocal1(vnls);
    SUNNonlinearSolver nls;
    SUNNonlinearSolver_Ops ops;

    // Create the cptr
    nls = (SUNNonlinearSolver) malloc(sizeof(struct sunml_nls));
    if (nls == NULL) caml_raise_out_of_memory();

    ops = (SUNNonlinearSolver_Ops) malloc(sizeof *ops);
    if (ops == NULL) {
	free(nls);
	caml_raise_out_of_memory();
    }

    nls->content = (void *) vops;
    caml_register_generational_global_root((value *)&(nls->content));
    nls->ops = ops;

    /* Attach operations */
    ops->gettype = callml_custom_gettype;

    ops->initialize      = (Field(vops, RECORD_NLSOLVER_OPS_INIT)
	    == Val_none ? NULL : callml_custom_initialize);

    ops->setup           = (Field(vops, RECORD_NLSOLVER_OPS_SETUP)
	    == Val_none ? NULL : (sens ? callml_custom_setup_sens
				       : callml_custom_setup));

    ops->solve           = (sens ? callml_custom_solve_sens
			         : callml_custom_solve);
    ops->free            = callml_custom_free;
    ops->setsysfn        = (sens ? callml_custom_setsysfn_sens
				 : callml_custom_setsysfn);

    ops->setlsetupfn     = (Field(vops, RECORD_NLSOLVER_OPS_SET_LSETUP_FN)
	    == Val_none ? NULL : (sens ? callml_custom_setlsetupfn_sens
				       : callml_custom_setlsetupfn));

    ops->setlsolvefn     = (Field(vops, RECORD_NLSOLVER_OPS_SET_LSOLVE_FN)
	    == Val_none ? NULL : (sens ? callml_custom_setlsolvefn_sens
				       : callml_custom_setlsolvefn));

    ops->setctestfn      = (Field(vops, RECORD_NLSOLVER_OPS_SET_CONVTEST_FN)
	    == Val_none ? NULL : (sens ? callml_custom_setconvtestfn_sens
				       : callml_custom_setconvtestfn));

    ops->setmaxiters     = (Field(vops, RECORD_NLSOLVER_OPS_SET_MAX_ITERS)
	    == Val_none ? NULL : callml_custom_setmaxiters);

    ops->getnumiters     = (Field(vops, RECORD_NLSOLVER_OPS_GET_NUM_ITERS)
	    == Val_none ? NULL : callml_custom_getnumiters);

    ops->getcuriter      = (Field(vops, RECORD_NLSOLVER_OPS_GET_CUR_ITER)
	    == Val_none ? NULL : callml_custom_getcuriter);

    ops->getnumconvfails = (Field(vops, RECORD_NLSOLVER_OPS_GET_NUM_CONV_FAILS)
	    == Val_none ? NULL : callml_custom_getnumconvfails);

    NLS_CALLBACKS(nls) = vcallbacks;
    caml_register_generational_global_root(&NLS_CALLBACKS(nls));

    vnls = caml_alloc_final(1, &finalize_nls, 1, 20);
    NLSOLVER_VAL(vnls) = nls;

    CAMLreturn (vnls);
}
#endif

CAMLprim value sunml_nlsolver_custom_make(value vcallbacks, value vops)
{
    CAMLparam2(vcallbacks, vops);
#if SUNDIALS_LIB_VERSION >= 400
    CAMLreturn (custom_make(0, vcallbacks, vops));
#else
    CAMLreturn (Val_unit);
#endif
}

CAMLprim value sunml_nlsolver_custom_make_sens(value vcallbacks, value vops)
{
    CAMLparam2(vcallbacks, vops);
#if SUNDIALS_LIB_VERSION >= 400
    CAMLreturn (custom_make(1, vcallbacks, vops));
#else
    CAMLreturn (Val_unit);
#endif
}

