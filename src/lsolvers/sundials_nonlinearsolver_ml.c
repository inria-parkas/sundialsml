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

#if 400 <= SUNDIALS_LIB_VERSION

#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

typedef void* (*FromValFn)(value);

#define SYSFN_VAL(v) (*(SUNNonlinSolSysFn*)Data_custom_val(v))
#define LSETUPFN_VAL(v) (*(SUNNonlinSolLSetupFn*)Data_custom_val(v))
#define LSOLVEFN_VAL(v) (*(SUNNonlinSolLSolveFn*)Data_custom_val(v))
#define FROMVAL_VAL(v) (*(FromValFn*)Data_custom_val(v))

struct convtestfn_data {
    SUNNonlinSolConvTestFn convtestfn;
    FromValFn fromvalfn;
};

#define CONVTESTFN_VAL(v) (((struct convtestfn_data*)Data_custom_val(v))->convtestfn)
#define CONVTESTFN_FROMVALFN(v) ((( struct convtestfn_data*)Data_custom_val(v))->fromvalfn)

struct nls_callback_value {
    void *mem;
    value *pvcallbacks;				// holds Weak pointer
    struct sunml_nls_callbacks *pccallbacks;
    value *pvmem;
};
typedef struct nls_callback_value* p_nls_cbv;

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

/* Nonlinear solvers (NLS) - use cases

   A consumer instantiates an NLS, configures callbacks (sysfn, lsetupfn,
   lsolvefn, convtestfn), and invokes init, setup, and solve.
   There are two types of "consumers": integrators (e.g., CVODE or IDA) and
   application code (written in OCaml).

   There are three types of invocations:

   - setfn (set_sysfn, set_lsetupfn, set_lsolvefn, set_convtestfn)
     tells the NLS which callback to use

   - setup/solve
     call the NLS to solve a problem

   - callbackfn (sysfn, lsetupfn, lsolvefn, convtestfn)
     during setup/solve, the NLS calls back into the given functions

   Two types of data must be passed: nvectors, state.
   C passes N_Vector, OCaml receives payload
   OCaml passes Nvector.t, C receives N_Vector

   An nls_callback_value is used passed as the state. It includes both C and
   OCaml versions of the integrator state (conversions are done via to_value
   and from_value functions provided by integrators), and pointers to OCaml
   and C callbacks.

   1. C-integrator using (Newton/FixedPoint) C-NLS

      The integrator and NLS interact directly in C, but this interaction
      passes through wrapper functions. The wrapper functions support
      the use of nls_callback_value, allowing largely unrestricted use of the
      nls from C and OCaml. A more complicated implementation thus allows
      simpler types in the OCaml interface.

      setup/solve: the original function is cached in sunml_nls, Sundials
      instead calls into a wrapper function that creates an nls_callback_value
      and applies special treatment for the convtestfn.
      E.g.,
	  1. SUNNonlinSolSolve (arg: N_Vector)
	  2. sunml_nlsolver_wrapped_solve
	  3. callback to original C solve function

      setfn: the original callback is cached in sunml_nls, Sundials is given
      a wrapper callback.
      E.g.,
	  1. SUNNonlinSolSetSysFn (arg : SUNNonlinSolSysFn)
	  2. sunml_nlsolver_wrapped_setsysfn
	  3. original C setsysfn with callback sunml_nlsolver_wrapped_sysfn

      callbackfn: the wrapper function is invoked, it extracts the relevant
      part of nls_callback_value and calls the cached original callbac,.
      E.g.,
	  1. sunml_nlsolver_wrapped_sysfn
	  2. original sysfn

   2. C-integrator using (Custom) OCaml-NLS

      The rawptr is passed to the integrator.
      Ops point to the callml_custom_* or callml_custom_*_sens functions.

      setup/solve:
      The callml_custom_* functions pass into OCaml via the nls-solver-ops
      structure.
      E.g.,
	1. SUNNonlinSolSolve (arg: N_Vector)
	2. callml_custom_solve (via C ops table, uses NVEC_BACKLINK)
	3. OCaml ops.solve is invoked (arg: payload)

      setfn:
      The callml_custom_* functions pass into OCaml via registered callbacks
      "Sundials_NonlinearSolver.set_c_*" that call into
      Custom.set_c_*(_sens). The Custom.set_c_* functions (a) wrap the given
      function using sunml_nlsolver_call_* and (b) call the nls-solver-ops
      set_* function on the result.
      E.g.,
	1. SUNNonlinSolSetSysFn (arg : SUNNonlinSolSysFn)
	2. callml_custom_setsysfn (via C ops table)
	3. Sundials.NonlinearSolver.Custom.set_c_sys_fn (via caml_named_value)
	4. (apply sunml_nlsolver_call_sys_fn to arg giving an OCaml closure)
	5. OCaml ops.set_sys_fn to store closure

      callbackfn:
      The OCaml solver invokes the closure created by setfn, the
      sunml_nlsolver_call_* function calls the C function passed earlier.
      E.g.,
	1. sunml_nlsolver_call_sys_fn (arg: Nvector.t)
	2. calls *sysfn (via pointer in closure, uses NVEC_VAL)
	3. sysfn (arg: N_Vector)

      setup/solve is C -> OCaml so N_Vectors become payloads ('data)
      callbackfn is OCaml -> C so Nvector.t becomes N_Vector
      (This is unfortunate: it means that the OCaml solver must rewrap the
       data in an Nvector.t when using the callback.)

      with sensitivities: senswrappers are passed from C -> OCaml -> C

   3. OCaml-app using C-NLS

      setup/solve:
      The function is called via Sundials.NonlinearSolver and uses the
      FixedPointSolver/NewtonSolver case, which invokes sunml_nlsolver_*,
      which in turn calls into Sundials.
      E.g.,
	1. call to Sundials.NonlinearSolver.solve (arg: Nvector.t)
	2. sunml_nlsolver_solve (uses NVEC_VAL)
	3. callback to original C solve function (arg: N_Vector)

      setfn:
      The function is called via Sundials.NonlinearSolver using the
      FixedPointSolver/NewtonSolver case, which invokes sunml_nlsolver_set_*,
      which in turn calls into Sundials passing a generic callback function.
      E.g.,
	1. call to Sundials.NonlinearSolver.set_sys_fn (stores closure
	   in the OCaml nonlinear_solver-callbacks structure)
	2. sunml_nlsolver_set_sys_fn
	3. SUNNonlinSolSetSysFn specifying the generic sys_callback

      callbackfn:
      The generic callback function invokes the OCaml closure from the
      callback structure.
      E.g.,
	1. sys_callback (arg: N_Vector, uses NVEC_BACKLINK)
	2. OCaml callbacks.sysfn (arg: payload)

      setup/solve is OCaml -> C so an Nvector.t is required
      callbackfn is C to OCaml so an N_Vector becomes a payload

   4. OCaml-app using OCaml-NLS

      setup/solve:
      The function is called via Sundials.NonlinearSolver and uses the
      CustomSolver case, which directly invokes the OCaml function in
      nls-solver-ops.
      E.g.,
	1. call to Sundials.NonlinearSolver.solve (arg: Nvector.t)
	2. OCaml ops.solve is invoked (arg: payload)
      The Nvector.t is unwrapped for compatibility with case 2.
      
      setfn:
      The function is called via Sundials.NonlinearSolver using the
      CustomSolver case, which directly invokes the OCaml function in
      nls-solver-ops.
      E.g.,
	1. call to Sundials.NonlinearSolver.set_sys_fn (arg: OCaml closure)
	2. OCaml ops.set_sys_fn is invoked

      callbackfn:
      The closure passed to setfn is simply invoked directly.
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
    .identifier = "sunml_senswrapper",
    .finalize =    custom_finalize_default,
    .compare =     senswrapper_compare,
    .hash =        senswrapper_hash,
    .serialize =   custom_serialize_default,
    .deserialize = custom_deserialize_default
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

/* - - - O/Cnls: callback invocation (also C/Cnls for convtestfn) - - - */

#if 400 <= SUNDIALS_LIB_VERSION
static int sys_callback(N_Vector y, N_Vector F, void* mem)
{
    CAMLparam0();
    CAMLlocalN(args, 3);
    CAMLlocal1(vcallbacks);

    p_nls_cbv cbv = (p_nls_cbv)mem;

    args[0] = NVEC_BACKLINK(y);
    args[1] = NVEC_BACKLINK(F);
    args[2] = *(cbv->pvmem);
    WEAK_DEREF(vcallbacks, *(cbv->pvcallbacks));

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (
	Field(vcallbacks, RECORD_NLSOLVER_CALLBACKS_SYSFN), 3, args);
    if (!Is_exception_result (r)) CAMLreturnT(int, 0);

    r = Extract_exception (r);
    CAMLreturnT(int, translate_exception (r, RECOVERABLE));
}

static int sys_callback_sens(N_Vector y, N_Vector F, void* mem)
{
    CAMLparam0();
    CAMLlocalN(args, 3);
    CAMLlocal1(vcallbacks);

    p_nls_cbv cbv = (p_nls_cbv)mem;

    args[0] = sunml_senswrapper_wrap(y);
    args[1] = sunml_senswrapper_wrap(F);
    args[2] = *(cbv->pvmem);
    WEAK_DEREF(vcallbacks, *(cbv->pvcallbacks));

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (
	Field(vcallbacks, RECORD_NLSOLVER_CALLBACKS_SYSFN), 3, args);
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
    CAMLlocal1(vcallbacks);

    p_nls_cbv cbv = (p_nls_cbv)mem;

    args[0] = Val_bool(jbad);
    args[1] = *(cbv->pvmem);
    WEAK_DEREF(vcallbacks, *(cbv->pvcallbacks));

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (
	Field(vcallbacks, RECORD_NLSOLVER_CALLBACKS_LSETUPFN), 2, args);

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
    CAMLlocal1(vcallbacks);

    p_nls_cbv cbv = (p_nls_cbv)mem;

    args[0] = NVEC_BACKLINK(b);
    args[1] = *(cbv->pvmem);
    WEAK_DEREF(vcallbacks, *(cbv->pvcallbacks));

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (
	Field(vcallbacks, RECORD_NLSOLVER_CALLBACKS_LSOLVEFN), 2, args);

    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}
#endif

#if 400 <= SUNDIALS_LIB_VERSION
#if 500 <= SUNDIALS_LIB_VERSION
static int lsolve_callback_sens(N_Vector b, void* mem)
#else
static int lsolve_callback_sens(N_Vector y, N_Vector b, void* mem)
#endif
{
    CAMLparam0();
    CAMLlocalN(args, 2);
    CAMLlocal1(vcallbacks);

    p_nls_cbv cbv = (p_nls_cbv)mem;

    args[0] = sunml_senswrapper_wrap(b);
    args[1] = *(cbv->pvmem);
    WEAK_DEREF(vcallbacks, *(cbv->pvcallbacks));

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn (
	Field(vcallbacks, RECORD_NLSOLVER_CALLBACKS_LSOLVEFN), 2, args);

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

    p_nls_cbv cbv = (p_nls_cbv)mem;

    WEAK_DEREF(vcallbacks, *(cbv->pvcallbacks));

    args[0] = NVEC_BACKLINK(y);
    args[1] = NVEC_BACKLINK(del);
    args[2] = caml_copy_double(tol);
    args[3] = NVEC_BACKLINK(ewt);
    args[4] = *(cbv->pvmem);

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

#if 400 <= SUNDIALS_LIB_VERSION
static int convtest_callback_sens(SUNNonlinearSolver nls, N_Vector y, N_Vector del,
				  realtype tol, N_Vector ewt, void* mem)
{
    CAMLparam0();
    CAMLlocalN(args, 5);
    CAMLlocal1(vcallbacks);

    p_nls_cbv cbv = (p_nls_cbv)mem;

    WEAK_DEREF(vcallbacks, *(cbv->pvcallbacks));

    args[0] = sunml_senswrapper_wrap(y);
    args[1] = sunml_senswrapper_wrap(del);
    args[2] = caml_copy_double(tol);
    args[3] = sunml_senswrapper_wrap(ewt);
    args[4] = *(cbv->pvmem);

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


/* - - - O/Cnls functions (OCaml invoking init/setup/solve in a C NLS) - - - */

CAMLprim value sunml_nlsolver_init(value vnls)
{
    CAMLparam1(vnls);
#if 400 <= SUNDIALS_LIB_VERSION
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);

    int flag = SUNNonlinSolInitialize(nls);
    NLS_CHECK_FLAG("SUNNonlinSolInitialize", flag);
#endif
    CAMLreturn (Val_unit);
}

#if 500 <= SUNDIALS_LIB_VERSION
static int sunml_nlsolver_failed_ctestfn(SUNNonlinearSolver nls,
					 N_Vector y,
					 N_Vector del,
					 realtype tol,
					 N_Vector ewt,
					 void* mem)
{
    // this callback should never be invoked
    assert(0);
    return 0;
}
#endif

#if 400 <= SUNDIALS_LIB_VERSION
static int sunml_nlsolver_wrapped_ctestfn(SUNNonlinearSolver nls,
					  N_Vector y,
					  N_Vector del,
					  realtype tol,
					  N_Vector ewt,
					  void* mem)
{
    p_nls_cbv cbv = (p_nls_cbv)mem;
    return (cbv->pccallbacks->ctestfn(nls, y, del, tol, ewt, cbv->mem));
}
#endif

#if 500 <= SUNDIALS_LIB_VERSION
static int sunml_nlsolver_install_ctestfn(p_sunml_nls snls,
					  p_nls_cbv pcbmem)
{
    // we set the convtest function before calling setup, so that the ctdata
    // is set to &cbmem (which lives in our stack frame).
    SUNNonlinSolConvTestFn ctestfn = snls->c_callbacks.ctestfn;
    if (ctestfn == convtest_callback || ctestfn == convtest_callback_sens) {
	// callback into OCaml - no need for wrapping
	return (snls->orig_setctestfn(&(snls->nls), ctestfn, pcbmem));
    } else {
	// callback into C - need to unwrap the mem argument
	return (snls->orig_setctestfn(&(snls->nls), sunml_nlsolver_wrapped_ctestfn,
		pcbmem));
    }
}
#endif

CAMLprim value sunml_nlsolver_setup(value vnls, value vy, value vmem)
{
    CAMLparam3(vnls, vy, vmem);
#if 400 <= SUNDIALS_LIB_VERSION
    int flag;
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    N_Vector y  = NVEC_VAL(vy);

    p_sunml_nls snls = NLS_EXTENDED(nls);
    struct nls_callback_value cbmem = {
	.mem = (snls->from_value == NULL) ? NULL : snls->from_value(vmem),
	.pvcallbacks = &(snls->callbacks),
	.pccallbacks = &(snls->c_callbacks),
	.pvmem = &vmem
    };

#if 500 <= SUNDIALS_LIB_VERSION
    flag = sunml_nlsolver_install_ctestfn(snls, &cbmem);
    NLS_CHECK_FLAG("SUNNonlinSolSetConvTestFn", flag);
#endif

    // call setup (which may invoke ctestfn callback
    flag = snls->orig_setup(nls, y, &cbmem);
    NLS_CHECK_FLAG("SUNNonlinSolSetup", flag);

#if 500 <= SUNDIALS_LIB_VERSION && defined(SUNDIALS_ML_DEBUG)
    snls->orig_setctestfn(nls, sunml_nlsolver_failed_ctestfn, (void *)0xdeadbeef);
#endif

#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_solve(value vnls, value vargs, value vmem)
{
    CAMLparam3(vnls, vargs, vmem);
#if 400 <= SUNDIALS_LIB_VERSION
    int flag;
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    N_Vector y0 = NVEC_VAL(Field(vargs, 0));
    N_Vector y  = NVEC_VAL(Field(vargs, 1));
    N_Vector w  = NVEC_VAL(Field(vargs, 2));
    double tol  = Double_val(Field(vargs, 3));
    int callLSetup = Bool_val(Field(vargs, 4));

    p_sunml_nls snls = NLS_EXTENDED(nls);
    struct nls_callback_value cbmem = {
	.mem = (snls->from_value == NULL) ? NULL : snls->from_value(vmem),
	.pvcallbacks = &(snls->callbacks),
	.pccallbacks = &(snls->c_callbacks),
	.pvmem = &vmem
    };

#if 500 <= SUNDIALS_LIB_VERSION
    flag = sunml_nlsolver_install_ctestfn(snls, &cbmem);
    NLS_CHECK_FLAG("SUNNonlinSolSetConvTestFn", flag);
#endif

    flag = snls->orig_solve(nls, y0, y, w, tol, callLSetup, &cbmem);
    NLS_CHECK_FLAG("SUNNonlinSolSolve", flag);

#if 500 <= SUNDIALS_LIB_VERSION && defined(SUNDIALS_ML_DEBUG)
    snls->orig_setctestfn(nls, sunml_nlsolver_failed_ctestfn, (void *)0xdeadbeef);
#endif

#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_setup_sens(value vnls, value vy, value vmem)
{
    CAMLparam3(vnls, vy, vmem);
#if 400 <= SUNDIALS_LIB_VERSION
    int flag;
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    N_Vector y  = SENSWRAPPER(vy);

    p_sunml_nls snls = NLS_EXTENDED(nls);
    struct nls_callback_value cbmem = {
	.mem = (snls->from_value == NULL) ? NULL : snls->from_value(vmem),
	.pvcallbacks = &(snls->callbacks),
	.pccallbacks = &(snls->c_callbacks),
	.pvmem = &vmem
    };

#if 500 <= SUNDIALS_LIB_VERSION
    flag = sunml_nlsolver_install_ctestfn(snls, &cbmem);
    NLS_CHECK_FLAG("SUNNonlinSolSetConvTestFn", flag);
#endif

    // call setup (which may invoke ctestfn callback
    flag = snls->orig_setup(nls, y, &cbmem);
    NLS_CHECK_FLAG("SUNNonlinSolSetup", flag);

#if 500 <= SUNDIALS_LIB_VERSION && defined(SUNDIALS_ML_DEBUG)
    snls->orig_setctestfn(nls, sunml_nlsolver_failed_ctestfn, (void *)0xdeadbeef);
#endif
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_solve_sens(value vnls, value vargs, value vmem)
{
    CAMLparam3(vnls, vargs, vmem);
#if 400 <= SUNDIALS_LIB_VERSION
    int flag;
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    N_Vector y0 = SENSWRAPPER(Field(vargs, 0));
    N_Vector y  = SENSWRAPPER(Field(vargs, 1));
    N_Vector w  = SENSWRAPPER(Field(vargs, 2));
    double tol  = Double_val(Field(vargs, 3));
    int callLSetup = Bool_val(Field(vargs, 4));

    p_sunml_nls snls = NLS_EXTENDED(nls);
    struct nls_callback_value cbmem = {
	.mem = (snls->from_value == NULL) ? NULL : snls->from_value(vmem),
	.pvcallbacks = &(snls->callbacks),
	.pccallbacks = &(snls->c_callbacks),
	.pvmem = &vmem
    };

#if 500 <= SUNDIALS_LIB_VERSION
    flag = sunml_nlsolver_install_ctestfn(snls, &cbmem);
    NLS_CHECK_FLAG("SUNNonlinSolSetConvTestFn", flag);
#endif

    flag = snls->orig_solve(nls, y0, y, w, tol, callLSetup, &cbmem);
    NLS_CHECK_FLAG("SUNNonlinSolSolve", flag);

#if 500 <= SUNDIALS_LIB_VERSION && defined(SUNDIALS_ML_DEBUG)
    snls->orig_setctestfn(nls, sunml_nlsolver_failed_ctestfn, (void *)0xdeadbeef);
#endif

#endif
    CAMLreturn (Val_unit);
}

/* - - - O/Cnls: callback configuration - - - */

CAMLprim value sunml_nlsolver_set_sys_fn(value vnls)
{
    CAMLparam1(vnls);
#if 400 <= SUNDIALS_LIB_VERSION
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    p_sunml_nls snls = NLS_EXTENDED(nls);

    int flag = snls->orig_setsysfn(nls, sys_callback);
    NLS_CHECK_FLAG("SUNNonlinSolSetSysFn", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_set_lsetup_fn(value vnls)
{
    CAMLparam1(vnls);
#if 400 <= SUNDIALS_LIB_VERSION
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    p_sunml_nls snls = NLS_EXTENDED(nls);

    int flag = snls->orig_setlsetupfn(nls, lsetup_callback);
    NLS_CHECK_FLAG("SUNNonlinSolSetLSetupFn", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_set_lsolve_fn(value vnls)
{
    CAMLparam1(vnls);
#if 400 <= SUNDIALS_LIB_VERSION
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    p_sunml_nls snls = NLS_EXTENDED(nls);

    int flag = snls->orig_setlsolvefn(nls, lsolve_callback);
    NLS_CHECK_FLAG("SUNNonlinSolSetLSolveFn", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_set_convtest_fn(value vnls, value vcfun)
{
    CAMLparam2(vnls, vcfun);
#if 400 <= SUNDIALS_LIB_VERSION
    CAMLlocal1(vcptr);
    int flag;
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    p_sunml_nls snls = NLS_EXTENDED(nls);
    vcptr = ML_CALLBACK_CPTR(vcfun);

    snls->c_callbacks.ctestfn = CONVTESTFN_VAL(vcptr);
    snls->from_value = CONVTESTFN_FROMVALFN(vcptr);

#if 500 <= SUNDIALS_LIB_VERSION
    // the ctestfn is set just before invoking setup or solve so that
    // the integrator pointer can be wrapped and passed as ctestdata.
    // we must set a value here because some nls check that ctestfn <> NULL
    flag = snls->orig_setctestfn(nls, sunml_nlsolver_failed_ctestfn,
	    (void *)0xdeadbeef);
#else
    flag = snls->orig_setctestfn(nls, snls->c_callbacks.ctestfn);
#endif
    NLS_CHECK_FLAG("SUNNonlinSolSetConvTestFn", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_set_convtest_fn_callback(value vnls)
{
    CAMLparam1(vnls);
#if 400 <= SUNDIALS_LIB_VERSION
    int flag;
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    p_sunml_nls snls = NLS_EXTENDED(nls);

    snls->c_callbacks.ctestfn = convtest_callback;
#if 500 <= SUNDIALS_LIB_VERSION
    // the ctestfn is set just before invoking setup or solve so that
    // the integrator pointer can be wrapped and passed as ctestdata.
    // we must set a value here because some nls check that ctestfn <> NULL
    flag = snls->orig_setctestfn(nls, sunml_nlsolver_failed_ctestfn,
	    (void *)0xdeadbeef);
#else
    flag = snls->orig_setctestfn(nls, convtest_callback);
#endif
    NLS_CHECK_FLAG("SUNNonlinSolSetConvTestFn", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_set_sys_fn_sens(value vnls)
{
    CAMLparam1(vnls);
#if 400 <= SUNDIALS_LIB_VERSION
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    p_sunml_nls snls = NLS_EXTENDED(nls);

    int flag = snls->orig_setsysfn(nls, sys_callback_sens);
    NLS_CHECK_FLAG("SUNNonlinSolSetSysFn", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_set_lsolve_fn_sens(value vnls)
{
    CAMLparam1(vnls);
#if 400 <= SUNDIALS_LIB_VERSION
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    p_sunml_nls snls = NLS_EXTENDED(nls);

    int flag = snls->orig_setlsolvefn(nls, lsolve_callback_sens);
    NLS_CHECK_FLAG("SUNNonlinSolSetLSolveFn", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_set_convtest_fn_callback_sens(value vnls)
{
    CAMLparam1(vnls);
#if 400 <= SUNDIALS_LIB_VERSION
    int flag;
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    p_sunml_nls snls = NLS_EXTENDED(nls);

    snls->c_callbacks.ctestfn = convtest_callback_sens;
#if 500 <= SUNDIALS_LIB_VERSION
    // the ctestfn is set just before invoking setup or solve so that
    // the integrator pointer can be wrapped and passed as ctestdata.
    // we must set a value here because some nls check that ctestfn <> NULL
    flag = snls->orig_setctestfn(nls, sunml_nlsolver_failed_ctestfn,
	    (void *)0xdeadbeef);
#else
    flag = snls->orig_setctestfn(nls, convtest_callback_sens);
#endif
    NLS_CHECK_FLAG("SUNNonlinSolSetConvTestFn", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_set_max_iters(value vnls, value vi)
{
    CAMLparam2(vnls, vi);
#if 400 <= SUNDIALS_LIB_VERSION
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
#if 400 <= SUNDIALS_LIB_VERSION
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
#if 400 <= SUNDIALS_LIB_VERSION
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
#if 400 <= SUNDIALS_LIB_VERSION
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

#if 400 <= SUNDIALS_LIB_VERSION
static void finalize_nls(value vnls)
{
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    caml_remove_generational_global_root(&NLS_CALLBACKS(nls));
    SUNNonlinSolFree(nls);
}
#endif

CAMLprim value sunml_nlsolver_call_sys_fn(value vfns, value vy,
					  value vfg, value vmem)
{
    CAMLparam4(vfns, vy, vfg, vmem);
#if 400 <= SUNDIALS_LIB_VERSION
    SUNNonlinSolSysFn sysfn = SYSFN_VAL(Field(vfns, 0));
    void *(*from_value)(value) = FROMVAL_VAL(Field(vfns, 1));
    N_Vector y  = NVEC_VAL(vy);
    N_Vector fg = NVEC_VAL(vfg);

    int flag = (*sysfn)(y, fg, from_value(vmem));
    NLS_CHECK_FLAG("SUNNonlinSolSysFn", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_call_sys_fn_sens(value vfns, value vy,
					       value vfg, value vmem)
{
    CAMLparam4(vfns, vy, vfg, vmem);
#if 400 <= SUNDIALS_LIB_VERSION
    SUNNonlinSolSysFn sysfn = SYSFN_VAL(Field(vfns, 0));
    void *(*from_value)(value) = FROMVAL_VAL(Field(vfns, 1));
    N_Vector y  = SENSWRAPPER(vy);
    N_Vector fg = SENSWRAPPER(vfg);

#if SUNDIALS_ML_SAFE == 1
    if ((y == NULL) || (fg == NULL))
	caml_raise_constant(NLSOLVER_EXN(IncorrectUse));
#endif

    int flag = (*sysfn)(y, fg, from_value(vmem));
    NLS_CHECK_FLAG("SUNNonlinSolSysFn (sens)", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_call_lsetup_fn(value vfns,
					     value vjbad, value vmem)
{
    CAMLparam3(vfns, vjbad, vmem);
#if 400 <= SUNDIALS_LIB_VERSION
    SUNNonlinSolLSetupFn lsetupfn = LSETUPFN_VAL(Field(vfns, 0));
    void *(*from_value)(value) = FROMVAL_VAL(Field(vfns, 1));
    int jcur = 0;

#if 500 <= SUNDIALS_LIB_VERSION
    int flag = (*lsetupfn)(Bool_val(vjbad), &jcur, from_value(vmem));
#else
    int flag = (*lsetupfn)(NULL, NULL, Bool_val(vjbad), &jcur, from_value(vmem));
#endif
    NLS_CHECK_FLAG("SUNNonlinSolLSetupFn", flag);

    CAMLreturn (Val_bool(jcur));
#endif
    CAMLreturn (Val_bool(0));
}

CAMLprim value sunml_nlsolver_call_lsetup_fn_sens(value vfns,
					          value vjbad, value vmem)
{
    CAMLparam3(vfns, vjbad, vmem);
#if 400 <= SUNDIALS_LIB_VERSION
    SUNNonlinSolLSetupFn lsetupfn = LSETUPFN_VAL(Field(vfns, 0));
    void *(*from_value)(value) = FROMVAL_VAL(Field(vfns, 1));
    int jcur = 0;

#if 500 <= SUNDIALS_LIB_VERSION
    int flag = (*lsetupfn)(Bool_val(vjbad), &jcur, from_value(vmem));
#else
    int flag = (*lsetupfn)(NULL, NULL, Bool_val(vjbad), &jcur, from_value(vmem));
#endif
    NLS_CHECK_FLAG("SUNNonlinSolLSetupFn (sens)", flag);

    CAMLreturn (Val_bool(jcur));
#endif
    CAMLreturn (Val_bool(0));
}

CAMLprim value sunml_nlsolver_call_lsolve_fn(value vfns, value vb, value vmem)
{
    CAMLparam3(vfns, vb, vmem);
#if 400 <= SUNDIALS_LIB_VERSION
    SUNNonlinSolLSolveFn lsolvefn = LSOLVEFN_VAL(Field(vfns, 0));
    void *(*from_value)(value) = FROMVAL_VAL(Field(vfns, 1));
    N_Vector b = NVEC_VAL(vb);

#if 500 <= SUNDIALS_LIB_VERSION
    int flag = (*lsolvefn)(b, from_value(vmem));
#else
    int flag = (*lsolvefn)(NULL, b, from_value(vmem));
#endif
    NLS_CHECK_FLAG("SUNNonlinSolLSolveFn", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_call_lsolve_fn_sens(value vfns,
						  value vb, value vmem)
{
    CAMLparam3(vfns, vb, vmem);
#if 400 <= SUNDIALS_LIB_VERSION
    SUNNonlinSolLSolveFn lsolvefn = LSOLVEFN_VAL(Field(vfns, 0));
    void *(*from_value)(value) = FROMVAL_VAL(Field(vfns, 1));
    N_Vector b = SENSWRAPPER(vb);

#if SUNDIALS_ML_SAFE == 1
    if (b == NULL) caml_raise_constant(NLSOLVER_EXN(IncorrectUse));
#endif

#if 500 <= SUNDIALS_LIB_VERSION
    int flag = (*lsolvefn)(b, from_value(vmem));
#else
    int flag = (*lsolvefn)(NULL, b, from_value(vmem));
#endif
    NLS_CHECK_FLAG("SUNNonlinSolLSolveFn (sens)", flag);
#endif
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_nlsolver_call_convtest_fn(value vnls, value vconvtestfn,
					       value vargs)
{
    CAMLparam3(vnls, vconvtestfn, vargs);
#if 400 <= SUNDIALS_LIB_VERSION
    SUNNonlinSolConvTestFn convtestfn = CONVTESTFN_VAL(vconvtestfn);
    FromValFn from_value = CONVTESTFN_FROMVALFN(vconvtestfn);

    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    N_Vector y   = NVEC_VAL(Field(vargs, 0));
    N_Vector del = NVEC_VAL(Field(vargs, 1));
    realtype tol = Double_val(Field(vargs, 2));
    N_Vector ewt = NVEC_VAL(Field(vargs, 3));
    void *mem    = from_value(Field(vargs, 4));

    int flag = (*convtestfn)(nls, y, del, tol, ewt, mem);
    switch (flag) {
    case SUN_NLS_SUCCESS:
	CAMLreturn (Val_int(VARIANT_NLSOLVER_CONVTEST_SUCCESS));

    case SUN_NLS_CONTINUE:
	CAMLreturn (Val_int(VARIANT_NLSOLVER_CONVTEST_CONTINUE));

    case SUN_NLS_CONV_RECVR:
	CAMLreturn (Val_int(VARIANT_NLSOLVER_CONVTEST_RECOVER));

    default:
	sunml_nlsolver_check_flag("SUNNonlinSolConvTestFn", flag);
    }
#endif
    CAMLreturn (Val_int(VARIANT_NLSOLVER_CONVTEST_SUCCESS));
}

CAMLprim value sunml_nlsolver_call_convtest_fn_sens(value vnls,
						    value vconvtestfn,
						    value vargs)
{
    CAMLparam3(vnls, vconvtestfn, vargs);
#if 400 <= SUNDIALS_LIB_VERSION
    SUNNonlinSolConvTestFn convtestfn = CONVTESTFN_VAL(vconvtestfn);
    FromValFn from_value = CONVTESTFN_FROMVALFN(vconvtestfn);

    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    N_Vector y   = SENSWRAPPER(Field(vargs, 0));
    N_Vector del = SENSWRAPPER(Field(vargs, 1));
    realtype tol = Double_val(Field(vargs, 2));
    N_Vector ewt = SENSWRAPPER(Field(vargs, 3));
    void *mem    = from_value(Field(vargs, 4));

#if SUNDIALS_ML_SAFE == 1
    if ((y == NULL) || (del == NULL) || (ewt == NULL))
	caml_raise_constant(NLSOLVER_EXN(IncorrectUse));
#endif

    int flag = (*convtestfn)(nls, y, del, tol, ewt, mem);
    switch (flag) {
    case SUN_NLS_SUCCESS:
	CAMLreturn (Val_int(VARIANT_NLSOLVER_CONVTEST_SUCCESS));

    case SUN_NLS_CONTINUE:
	CAMLreturn (Val_int(VARIANT_NLSOLVER_CONVTEST_CONTINUE));

    case SUN_NLS_CONV_RECVR:
	CAMLreturn (Val_int(VARIANT_NLSOLVER_CONVTEST_RECOVER));

    default:
	sunml_nlsolver_check_flag("SUNNonlinSolConvTestFn (sens)", flag);
    }
#endif
    CAMLreturn (Val_int(VARIANT_NLSOLVER_CONVTEST_SUCCESS));
}

#if 400 <= SUNDIALS_LIB_VERSION
void sunml_nlsolver_set_to_from_mem(SUNNonlinearSolver nls,
				    value (*to_value)(void *),
				    void* (*from_value)(value))
{
    p_sunml_nls snls = NLS_EXTENDED(nls);
    snls->to_value = to_value;
    snls->from_value = from_value;
}
#endif

#if 400 <= SUNDIALS_LIB_VERSION
static int sunml_nlsolver_wrapped_setup(SUNNonlinearSolver nls,
					N_Vector y, void* mem)
{
    CAMLparam0();
    CAMLlocal1(vmem);
    int r = 0;

    p_sunml_nls snls = NLS_EXTENDED(nls);
    if (snls->to_value != NULL) {
	vmem = snls->to_value(mem);
    } else {
	vmem = Val_unit;
    }

    struct nls_callback_value cbmem = {
	.mem = mem,
	.pvcallbacks = &(snls->callbacks),
	.pccallbacks = &(snls->c_callbacks),
	.pvmem = &vmem,
    };

#if 500 <= SUNDIALS_LIB_VERSION
    if (sunml_nlsolver_install_ctestfn(snls, &cbmem) != SUN_NLS_SUCCESS) goto done;
#endif

    // call setup (which may invoke ctestfn callback
    r = snls->orig_setup(nls, y, &cbmem);

#if 500 <= SUNDIALS_LIB_VERSION && defined(SUNDIALS_ML_DEBUG)
    snls->orig_setctestfn(nls, sunml_nlsolver_failed_ctestfn, (void *)0xdeadbeef);
#endif

#if 500 <= SUNDIALS_LIB_VERSION
done:
#endif
    CAMLreturnT(int, r);
}

static int sunml_nlsolver_wrapped_solve(SUNNonlinearSolver nls,
					N_Vector y0, N_Vector y,
					N_Vector w, realtype tol,
					booleantype callLSetup, void *mem)
{
    CAMLparam0();
    CAMLlocal1(vmem);
    int r = 0;

    p_sunml_nls snls = NLS_EXTENDED(nls);
    if (snls->to_value != NULL) {
	vmem = snls->to_value(mem);
    } else {
	vmem = Val_unit;
    }

    struct nls_callback_value cbmem = {
	.mem = mem,
	.pvcallbacks = &(snls->callbacks),
	.pccallbacks = &(snls->c_callbacks),
	.pvmem = &vmem,
    };

#if 500 <= SUNDIALS_LIB_VERSION
    if (sunml_nlsolver_install_ctestfn(snls, &cbmem) != SUN_NLS_SUCCESS) goto done;
#endif

    r = snls->orig_solve(nls, y0, y, w, tol, callLSetup, &cbmem);

#if 500 <= SUNDIALS_LIB_VERSION && defined(SUNDIALS_ML_DEBUG)
    snls->orig_setctestfn(nls, sunml_nlsolver_failed_ctestfn, (void *)0xdeadbeef);
#endif

#if 500 <= SUNDIALS_LIB_VERSION
done:
#endif
    CAMLreturnT(int, r);
}

static int sunml_nlsolver_wrapped_sysfn(N_Vector y, N_Vector F, void* mem)
{
    p_nls_cbv cbv = (p_nls_cbv)mem;
    return (cbv->pccallbacks->sysfn(y, F, cbv->mem));
}

#if 500 <= SUNDIALS_LIB_VERSION
static int sunml_nlsolver_wrapped_lsetupfn(booleantype jbad, booleantype* jcur,
					   void* mem)
{
    p_nls_cbv cbv = (p_nls_cbv)mem;
    return (cbv->pccallbacks->lsetupfn(jbad, jcur, cbv->mem));
}
#else
static int sunml_nlsolver_wrapped_lsetupfn(N_Vector y, N_Vector F,
					   booleantype jbad, booleantype* jcur,
					   void* mem)
{
    p_nls_cbv cbv = (p_nls_cbv)mem;
    return (cbv->pccallbacks->lsetupfn(y, F, jbad, jcur, cbv->mem));
}
#endif

#if 500 <= SUNDIALS_LIB_VERSION
static int sunml_nlsolver_wrapped_lsolvefn(N_Vector b, void* mem)
{
    p_nls_cbv cbv = (p_nls_cbv)mem;
    return (cbv->pccallbacks->lsolvefn(b, cbv->mem));
}
#else
static int sunml_nlsolver_wrapped_lsolvefn(N_Vector y, N_Vector b, void* mem)
{
    p_nls_cbv cbv = (p_nls_cbv)mem;
    return (cbv->pccallbacks->lsolvefn(y, b, cbv->mem));
}
#endif

static int sunml_nlsolver_wrapped_setsysfn(SUNNonlinearSolver nls,
					   SUNNonlinSolSysFn sysfn)
{
    p_sunml_nls snls = NLS_EXTENDED(nls);
    snls->c_callbacks.sysfn = sysfn;
    return (snls->orig_setsysfn(nls, sunml_nlsolver_wrapped_sysfn));
}

static int sunml_nlsolver_wrapped_setlsetupfn(SUNNonlinearSolver nls,
					      SUNNonlinSolLSetupFn lsetupfn)
{
    p_sunml_nls snls = NLS_EXTENDED(nls);
    snls->c_callbacks.lsetupfn = lsetupfn;
#if 500 <= SUNDIALS_LIB_VERSION
    return (snls->orig_setlsetupfn(nls, sunml_nlsolver_wrapped_lsetupfn));
#else
    return (snls->orig_setlsetupfn(nls, sunml_nlsolver_wrapped_lsetupfn));
#endif
}

static int sunml_nlsolver_wrapped_setlsolvefn(SUNNonlinearSolver nls,
					      SUNNonlinSolLSolveFn lsolvefn)
{
    p_sunml_nls snls = NLS_EXTENDED(nls);
    snls->c_callbacks.lsolvefn = lsolvefn;
    return (snls->orig_setlsolvefn(nls, sunml_nlsolver_wrapped_lsolvefn));
}

#if 500 <= SUNDIALS_LIB_VERSION
static int sunml_nlsolver_wrapped_setctestfn(SUNNonlinearSolver nls,
					     SUNNonlinSolConvTestFn ctestfn,
					     void* ctdata)
#else
static int sunml_nlsolver_wrapped_setctestfn(SUNNonlinearSolver nls,
					     SUNNonlinSolConvTestFn ctestfn)
#endif
{
    p_sunml_nls snls = NLS_EXTENDED(nls);
    snls->c_callbacks.ctestfn = ctestfn;

#if SUNDIALS_LIB_VERSION < 500
    return (snls->orig_setctestfn(nls, sunml_nlsolver_wrapped_ctestfn));
#else
    // We do not call orig_setctestfn here.
    // Rather this is done by sunml_nlsolver_wrapped_setup
    // and sunml_nlsolver_wrapped_solve so that they can set
    // ctdata to the session passed by the integrator. The Sundials
    // integrators pass the same value to solve/setup as they do to
    // SUNNonlinSolSetConvTestFn.
    // We cannot pass an nls_callback_value here without having
    // allocated it on the heap (but then who will later free it?)
    // or in the sunml_nls structure (but then it holds OCaml values, and
    // that potentially with a cycle).

    return (snls->orig_setctestfn(nls, sunml_nlsolver_failed_ctestfn,
				  (void *)0xdeadbeef));
#endif
}

static value rewrap_nlsolver(SUNNonlinearSolver nls0, value vcallbacks)
{
    CAMLparam1(vcallbacks);
    CAMLlocal1(vr);
    SUNNonlinearSolver nls;
    p_sunml_nls snls;

    if (nls0 == NULL) caml_raise_out_of_memory();

    nls = (SUNNonlinearSolver)calloc(1, sizeof(struct sunml_nls));
    if (nls == NULL) {
	SUNNonlinSolFree(nls0);
	caml_raise_out_of_memory();
    }

    nls->content = nls0->content;
    nls->ops = nls0->ops;
    free(nls0);

    snls = NLS_EXTENDED(nls);
    snls->callbacks = vcallbacks;
    caml_register_generational_global_root(&(snls->callbacks));

    /* add wrapping to allow callbacks into OCaml */
    snls->orig_setup = nls->ops->setup;
    nls->ops->setup =
	(snls->orig_setup == NULL) ? NULL : sunml_nlsolver_wrapped_setup;
    snls->orig_solve = nls->ops->solve;
    nls->ops->solve = sunml_nlsolver_wrapped_solve;

    if (nls->ops->setsysfn != NULL) {
	snls->orig_setsysfn = nls->ops->setsysfn;
	nls->ops->setsysfn = sunml_nlsolver_wrapped_setsysfn;
    }
    if (nls->ops->setlsetupfn != NULL) {
	snls->orig_setlsetupfn = nls->ops->setlsetupfn;
	nls->ops->setlsetupfn = sunml_nlsolver_wrapped_setlsetupfn;
    }
    if (nls->ops->setlsolvefn != NULL) {
	snls->orig_setlsolvefn = nls->ops->setlsolvefn;
	nls->ops->setlsolvefn = sunml_nlsolver_wrapped_setlsolvefn;
    }
    if (nls->ops->setctestfn != NULL) {
	snls->orig_setctestfn = nls->ops->setctestfn;
	nls->ops->setctestfn = sunml_nlsolver_wrapped_setctestfn;
    }

    // Setup the OCaml-side
    vr = caml_alloc_final(1, &finalize_nls, 1, 20);
    NLSOLVER_VAL(vr) = nls;

    CAMLreturn (vr);
}

static value sunml_nlsolver_wrap_from_value(SUNNonlinearSolver nls)
{
    CAMLparam0();
    CAMLlocal1(vfromvalue);

    p_sunml_nls snls = NLS_EXTENDED(nls);

    // an integrator cannot call callml_custom_set*fn without having
    // installed to_value/from_value
    assert(snls->from_value != NULL);

    vfromvalue = caml_alloc_final(1, NULL, 0, 1);
    FROMVAL_VAL(vfromvalue) = snls->from_value;

    CAMLreturn(vfromvalue);
}

#endif

CAMLprim value sunml_nlsolver_newton_make(value vy, value vcallbacks)
{
    CAMLparam2(vy, vcallbacks);
#if 400 <= SUNDIALS_LIB_VERSION
    CAMLlocal1(vr);
    N_Vector y = NVEC_VAL(vy);

    vr = rewrap_nlsolver(SUNNonlinSol_Newton(y), vcallbacks);
    CAMLreturn (vr);
#else
    CAMLreturn (Val_unit);
#endif
}

CAMLprim value sunml_nlsolver_newton_make_sens(value vcount, value vy,
					       value vcallbacks)
{
    CAMLparam3(vcount, vy, vcallbacks);
#if 400 <= SUNDIALS_LIB_VERSION
    CAMLlocal1(vr);
    N_Vector y = NVEC_VAL(vy);

    vr = rewrap_nlsolver(SUNNonlinSol_NewtonSens(Int_val(vcount), y), vcallbacks);
    CAMLreturn (vr);
#else
    CAMLreturn (Val_unit);
#endif
}

CAMLprim value sunml_nlsolver_newton_get_sys_fn(value vnls)
{
    CAMLparam1(vnls);
    CAMLlocal3(vr, vfns, vsysfn);
#if 400 <= SUNDIALS_LIB_VERSION
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    p_sunml_nls snls = NLS_EXTENDED(nls);
    SUNNonlinSolSysFn sysfn;

    int flag = SUNNonlinSolGetSysFn_Newton(nls, &sysfn);
    NLS_CHECK_FLAG("SUNNonlinSolGetSysFn_Newton", flag);

    vr = caml_alloc_tuple(2);

    if (sysfn == NULL) {
	Store_field(vr, 0, Val_false);
	Store_field(vr, 1, Val_none);

    } else if (sysfn == sys_callback) {
	Store_field(vr, 0, Val_true);
	Store_field(vr, 1, Val_none);

    } else {
	assert (sysfn == sunml_nlsolver_wrapped_sysfn);
	vsysfn = caml_alloc_final(1, NULL, 0, 1);
	SYSFN_VAL(vsysfn) = snls->c_callbacks.sysfn;

	vfns = caml_alloc_tuple(2);
	Store_field(vfns, 0, vsysfn);
	Store_field(vfns, 1, sunml_nlsolver_wrap_from_value(nls));
	Store_some(vfns, vfns);

	Store_field(vr, 0, Val_false);
	Store_field(vr, 1, vfns);
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
#if 400 <= SUNDIALS_LIB_VERSION
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
#if 400 <= SUNDIALS_LIB_VERSION
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

#if 510 <= SUNDIALS_LIB_VERSION
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
    CAMLlocal3(vr, vfns, vsysfn);
#if 400 <= SUNDIALS_LIB_VERSION
    SUNNonlinearSolver nls = NLSOLVER_VAL(vnls);
    p_sunml_nls snls = NLS_EXTENDED(nls);
    SUNNonlinSolSysFn sysfn;

    int flag = SUNNonlinSolGetSysFn_FixedPoint(nls, &sysfn);
    NLS_CHECK_FLAG("SUNNonlinSolGetSysFn_FixedPoint", flag);

    vr = caml_alloc_tuple(2);

    if (sysfn == NULL) {
	Store_field(vr, 0, Val_false);
	Store_field(vr, 1, Val_none);

    } else if (sysfn == sys_callback) {
	Store_field(vr, 0, Val_true);
	Store_field(vr, 1, Val_none);

    } else {
	assert (sysfn == sunml_nlsolver_wrapped_sysfn);
	vsysfn = caml_alloc_final(1, NULL, 0, 1);
	SYSFN_VAL(vsysfn) = snls->c_callbacks.sysfn;

	vfns = caml_alloc_tuple(2);
	Store_field(vfns, 0, vsysfn);
	Store_field(vfns, 1, sunml_nlsolver_wrap_from_value(nls));
	Store_some(vfns, vfns);

	Store_field(vr, 0, Val_false);
	Store_field(vr, 1, vfns);
    }

    CAMLreturn (vr);
#else
    CAMLreturn (Val_unit);
#endif
}

/* - - - O/Onls && C/Onls: custom nonlinear solvers - - - */

#if 400 <= SUNDIALS_LIB_VERSION

#define NLSOLV_OP_TABLE(nls)  ((value)((nls)->content))
#define GET_OP(vops, x) (Field((vops), RECORD_NLSOLVER_OPS_ ## x))
#define GET_SOME_OP(vops, x) (Some_val (GET_OP(vops, x)))

static SUNNonlinearSolver_Type callml_custom_gettype(SUNNonlinearSolver nls)
{
    CAMLparam0();
    CAMLlocal1(vops);
    WEAK_DEREF(vops, NLSOLV_OP_TABLE(nls));
    SUNNonlinearSolver_Type r = SUNNONLINEARSOLVER_ROOTFIND;

    switch (Int_val(GET_OP(vops, TYPE))) {
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
    CAMLlocal1(vops);
    WEAK_DEREF(vops, NLSOLV_OP_TABLE(nls));
    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn(GET_SOME_OP(vops, INIT), Val_unit);
    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static int callml_custom_setup(SUNNonlinearSolver nls, N_Vector y, void* mem)
{
    CAMLparam0();
    CAMLlocal1(vops);
    WEAK_DEREF(vops, NLSOLV_OP_TABLE(nls));
    CAMLlocalN(args, 2);

    p_sunml_nls snls = NLS_EXTENDED(nls);
    // an integrator cannot call callml_custom_setup without having
    // installed to_value/from_value
    assert(snls->to_value != NULL);

    args[0] = NVEC_BACKLINK(y);
    args[1] = snls->to_value(mem);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(GET_SOME_OP(vops, SETUP), 2, args);
    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static int callml_custom_setup_sens(SUNNonlinearSolver nls, N_Vector y, void* mem)
{
    CAMLparam0();
    CAMLlocal1(vops);
    WEAK_DEREF(vops, NLSOLV_OP_TABLE(nls));
    CAMLlocalN(args, 2);

    p_sunml_nls snls = NLS_EXTENDED(nls);
    // an integrator cannot call callml_custom_setup_sens without having
    // installed to_value/from_value
    assert(snls->to_value != NULL);

    args[0] = sunml_senswrapper_wrap(y);
    args[1] = snls->to_value(mem);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(GET_SOME_OP(vops, SETUP), 2, args);

    invalidate_senswrapper(args[0]);

    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static int callml_custom_solve(SUNNonlinearSolver nls,
			       N_Vector y0, N_Vector y, N_Vector w,
			       realtype tol, booleantype callLSetup, void* mem)
{
    CAMLparam0();
    CAMLlocal1(vops);
    WEAK_DEREF(vops, NLSOLV_OP_TABLE(nls));
    CAMLlocalN(args, 6);

    p_sunml_nls snls = NLS_EXTENDED(nls);
    // an integrator cannot call callml_custom_solve without having
    // installed to_value/from_value
    assert(snls->to_value != NULL);

    args[0] = NVEC_BACKLINK(y0);
    args[1] = NVEC_BACKLINK(y);
    args[2] = NVEC_BACKLINK(w);
    args[3] = caml_copy_double(tol);
    args[4] = Val_bool(callLSetup);
    args[5] = snls->to_value(mem);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(GET_OP(vops, SOLVE), 6, args);
    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, RECOVERABLE));
}

static int callml_custom_solve_sens(SUNNonlinearSolver nls,
				    N_Vector y0, N_Vector y, N_Vector w,
				    realtype tol, booleantype callLSetup,
				    void* mem)
{
    CAMLparam0();
    CAMLlocal1(vops);
    WEAK_DEREF(vops, NLSOLV_OP_TABLE(nls));
    CAMLlocalN(args, 6);

    p_sunml_nls snls = NLS_EXTENDED(nls);
    // an integrator cannot call callml_custom_solve_sens without having
    // installed to_value/from_value
    assert(snls->to_value != NULL);

    args[0] = sunml_senswrapper_wrap(y0);
    args[1] = sunml_senswrapper_wrap(y);
    args[2] = sunml_senswrapper_wrap(w);
    args[3] = caml_copy_double(tol);
    args[4] = Val_bool(callLSetup);
    args[5] = snls->to_value(mem);

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callbackN_exn(GET_OP(vops, SOLVE), 6, args);

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
    caml_remove_generational_global_root(&NLS_CALLBACKS(nls));

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
    CAMLlocal1(vops);
    WEAK_DEREF(vops, NLSOLV_OP_TABLE(nls));
    vsysfn = caml_alloc_final(1, NULL, 0, 1);
    SYSFN_VAL(vsysfn) = sysfn;

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn(
		*caml_named_value("Sundials_NonlinearSolver.set_c_sys_fn"),
		vops, vsysfn, sunml_nlsolver_wrap_from_value(nls));
    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static int callml_custom_setsysfn_sens(SUNNonlinearSolver nls,
				       SUNNonlinSolSysFn sysfn)
{
    CAMLparam0();
    CAMLlocal1(vsysfn);
    CAMLlocal1(vops);
    WEAK_DEREF(vops, NLSOLV_OP_TABLE(nls));
    vsysfn = caml_alloc_final(1, NULL, 0, 1);
    SYSFN_VAL(vsysfn) = sysfn;

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn(
		*caml_named_value("Sundials_NonlinearSolver.set_c_sys_fn_sens"),
		vops, vsysfn, sunml_nlsolver_wrap_from_value(nls));
    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static int callml_custom_setlsetupfn(SUNNonlinearSolver nls,
				     SUNNonlinSolLSetupFn lsetupfn)
{
    CAMLparam0();
    CAMLlocal1(vlsetupfn);
    CAMLlocal1(vops);
    WEAK_DEREF(vops, NLSOLV_OP_TABLE(nls));
    vlsetupfn = caml_alloc_final(1, NULL, 0, 1);
    LSETUPFN_VAL(vlsetupfn) = lsetupfn;

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn(
		*caml_named_value("Sundials_NonlinearSolver.set_c_lsetup_fn"),
		vops, vlsetupfn, sunml_nlsolver_wrap_from_value(nls));
    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static int callml_custom_setlsetupfn_sens(SUNNonlinearSolver nls,
					  SUNNonlinSolLSetupFn lsetupfn)
{
    CAMLparam0();
    CAMLlocal1(vlsetupfn);
    CAMLlocal1(vops);
    WEAK_DEREF(vops, NLSOLV_OP_TABLE(nls));
    vlsetupfn = caml_alloc_final(1, NULL, 0, 1);
    LSETUPFN_VAL(vlsetupfn) = lsetupfn;

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn(
		*caml_named_value("Sundials_NonlinearSolver.set_c_lsetup_fn_sens"),
		vops, vlsetupfn, sunml_nlsolver_wrap_from_value(nls));
    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static int callml_custom_setlsolvefn(SUNNonlinearSolver nls,
				     SUNNonlinSolLSolveFn lsolvefn)
{
    CAMLparam0();
    CAMLlocal1(vlsolvefn);
    CAMLlocal1(vops);
    WEAK_DEREF(vops, NLSOLV_OP_TABLE(nls));
    vlsolvefn = caml_alloc_final(1, NULL, 0, 1);
    LSOLVEFN_VAL(vlsolvefn) = lsolvefn;

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn(
		*caml_named_value("Sundials_NonlinearSolver.set_c_lsolve_fn"),
		vops, vlsolvefn, sunml_nlsolver_wrap_from_value(nls));
    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static int callml_custom_setlsolvefn_sens(SUNNonlinearSolver nls,
				          SUNNonlinSolLSolveFn lsolvefn)
{
    CAMLparam0();
    CAMLlocal1(vlsolvefn);
    CAMLlocal1(vops);
    WEAK_DEREF(vops, NLSOLV_OP_TABLE(nls));
    vlsolvefn = caml_alloc_final(1, NULL, 0, 1);
    LSOLVEFN_VAL(vlsolvefn) = lsolvefn;

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback3_exn(
		*caml_named_value("Sundials_NonlinearSolver.set_c_lsolve_fn_sens"),
		vops, vlsolvefn, sunml_nlsolver_wrap_from_value(nls));
    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

#if 500 <= SUNDIALS_LIB_VERSION
static int callml_custom_setconvtestfn(SUNNonlinearSolver nls,
				       SUNNonlinSolConvTestFn convtestfn,
				       void *cdata)
#else
static int callml_custom_setconvtestfn(SUNNonlinearSolver nls,
				       SUNNonlinSolConvTestFn convtestfn)
#endif
{
    CAMLparam0();
    CAMLlocal1(vconvtestfn);
    CAMLlocal1(vops);
    WEAK_DEREF(vops, NLSOLV_OP_TABLE(nls));
    p_sunml_nls snls = NLS_EXTENDED(nls);

    vconvtestfn = caml_alloc_final(
		    OCAMLSIZEOF(struct convtestfn_data), NULL, 0, 1);
    CONVTESTFN_VAL(vconvtestfn) = convtestfn;
    CONVTESTFN_FROMVALFN(vconvtestfn) = snls->from_value;
    
    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn(
		*caml_named_value("Sundials_NonlinearSolver.set_c_convtest_fn"),
		vops, vconvtestfn);
    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

#if 500 <= SUNDIALS_LIB_VERSION
static int callml_custom_setconvtestfn_sens(SUNNonlinearSolver nls,
				            SUNNonlinSolConvTestFn convtestfn,
					    void *cdata)
#else
static int callml_custom_setconvtestfn_sens(SUNNonlinearSolver nls,
				            SUNNonlinSolConvTestFn convtestfn)
#endif
{
    CAMLparam0();
    CAMLlocal1(vconvtestfn);
    CAMLlocal1(vops);
    WEAK_DEREF(vops, NLSOLV_OP_TABLE(nls));
    p_sunml_nls snls = NLS_EXTENDED(nls);

    vconvtestfn = caml_alloc_final(
		    OCAMLSIZEOF(struct convtestfn_data), NULL, 0, 1);
    CONVTESTFN_VAL(vconvtestfn) = convtestfn;
    CONVTESTFN_FROMVALFN(vconvtestfn) = snls->from_value;

    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback2_exn(
		*caml_named_value("Sundials_NonlinearSolver.set_c_convtest_fn_sens"),
		vops, vconvtestfn);
    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static int callml_custom_setmaxiters(SUNNonlinearSolver nls, int maxiters)
{
    CAMLparam0();
    CAMLlocal1(vops);
    WEAK_DEREF(vops, NLSOLV_OP_TABLE(nls));
    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn(GET_SOME_OP(vops, SET_MAX_ITERS),
				Val_int(maxiters));
    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static int callml_custom_getnumiters(SUNNonlinearSolver nls, long int *niters)
{
    CAMLparam0();
    CAMLlocal1(vops);
    WEAK_DEREF(vops, NLSOLV_OP_TABLE(nls));
    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn(GET_SOME_OP(vops, GET_NUM_ITERS), Val_unit);

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
    CAMLlocal1(vops);
    WEAK_DEREF(vops, NLSOLV_OP_TABLE(nls));
    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn(GET_SOME_OP(vops, GET_CUR_ITER), Val_unit);

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
    CAMLlocal1(vops);
    WEAK_DEREF(vops, NLSOLV_OP_TABLE(nls));
    /* NB: Don't trigger GC while processing this return value!  */
    value r = caml_callback_exn(GET_SOME_OP(vops, GET_NUM_CONV_FAILS), Val_unit);

    /* Update niters; leave it unchanged if an error occurred. */
    if (!Is_exception_result (r)) {
	*nconvfails = Long_val (r);
	CAMLreturnT(int, 0);
    }

    CAMLreturnT(int, CHECK_NLS_EXCEPTION (r, UNRECOVERABLE));
}

static CAMLprim value custom_make(int sens, value vcallbacks, value vweakops)
{
    CAMLparam2(vcallbacks, vweakops);
    CAMLlocal2(vnls, vops);
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

    nls->content = (void *) vweakops;
    caml_register_generational_global_root((value *)&(nls->content));
    nls->ops = ops;
    WEAK_DEREF(vops, vweakops);

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
#if 400 <= SUNDIALS_LIB_VERSION
    CAMLreturn (custom_make(0, vcallbacks, vops));
#else
    CAMLreturn (Val_unit);
#endif
}

CAMLprim value sunml_nlsolver_custom_make_sens(value vcallbacks, value vops)
{
    CAMLparam2(vcallbacks, vops);
#if 400 <= SUNDIALS_LIB_VERSION
    CAMLreturn (custom_make(1, vcallbacks, vops));
#else
    CAMLreturn (Val_unit);
#endif
}

