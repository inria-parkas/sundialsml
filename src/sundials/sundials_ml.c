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

/* Sundials interface functions that are common to CVODE and IDA. */

#include <stdio.h>
#include <errno.h>
#include <string.h>

#include <sundials/sundials_config.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_band.h>

#include <caml/mlvalues.h>
#include <caml/gc.h>
#include <caml/fail.h>
#include <caml/memory.h>
#include <caml/alloc.h>
#include <caml/custom.h>
#include <caml/callback.h>
#include <caml/bigarray.h>

#include "sundials_ml.h"

#if 580 <= SUNDIALS_LIB_VERSION
#include <sundials/sundials_math.h>
#endif

#if 600 <= SUNDIALS_LIB_VERSION
#include <sundials/sundials_context.h>
#endif

value sundials_ml_exn_table = 0;

void sunml_register_exns(enum sundials_exn_set_index index, value exns)
{
    CAMLparam1 (exns);
    CAMLlocal1 (r);
    if (sundials_ml_exn_table == 0) {
	int i;
	r = caml_alloc_small (SUNDIALS_NUM_EXN_SETS, 0);
        for (i = 0; i < SUNDIALS_NUM_EXN_SETS; ++i)
	    Field (r, i) = 0;
	Store_field (r, index, exns);
	sundials_ml_exn_table = r;
	caml_register_generational_global_root (&sundials_ml_exn_table);
    } else {
	Store_field (sundials_ml_exn_table, index, exns);
    }
    CAMLreturn0;
}

static value warn_discarded_exn = 0;

void sunml_warn_discarded_exn (value exn, const char *context)
{
    CAMLparam1 (exn);
    CAMLlocal1 (vcontext);
    vcontext = caml_copy_string (context);
    caml_callback2_exn (warn_discarded_exn, exn, vcontext);
    CAMLreturn0;
}

/* Setting up access to Weak.get */

#if !HAVE_WEAK
static value weak_get = 0;
#endif

CAMLprim void sunml_sundials_init_module (value vwarn_discarded_exn,
				      value vweak_get, value exns)
{
    CAMLparam2 (vweak_get, exns);
    REGISTER_EXNS (SUNDIALS, exns);
#if !HAVE_WEAK
    weak_get = vweak_get;
    caml_register_generational_global_root (&weak_get);
#endif
    warn_discarded_exn = vwarn_discarded_exn;
    caml_register_generational_global_root (&warn_discarded_exn);

    CAMLreturn0;
}

CAMLprim value sunml_sundials_get_constants (void)
{
    CAMLparam0 ();
    CAMLlocal1 (r);

    r = caml_alloc_tuple (3);
#if 600 <= SUNDIALS_LIB_VERSION
    Store_field (r, 0, caml_copy_double(SUN_BIG_REAL));
    Store_field (r, 1, caml_copy_double(SUN_SMALL_REAL));
    Store_field (r, 2, caml_copy_double(SUN_UNIT_ROUNDOFF));
#else
    Store_field (r, 0, caml_copy_double(BIG_REAL));
    Store_field (r, 1, caml_copy_double(SMALL_REAL));
    Store_field (r, 2, caml_copy_double(UNIT_ROUNDOFF));
#endif

    CAMLreturn (r);
}

#if !HAVE_WEAK
CAMLprim value sundials_ml_weak_get (value ar, value n)
{
    CAMLparam2 (ar, n);
    CAMLreturn (caml_callback2_exn (weak_get, ar, n));
}
#endif	/* !HAVE_WEAK */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * Two-dimensional real arrays based on BigArrays with a column-access table.
 *
 * We represent these as pairs:
 *	Field(ra, 0) = the underlying 2 dimensional big array
 *		       (for use from within OCaml)
 *	Field(ra, 1) = a custom value giving a sunrealtype ** table pointing into
 *		       the columns of the big array (for passing to Sundials)
 */

CAMLprim value sunml_sundials_realarray2_wrap(value vba)
{
    CAMLparam1(vba);
    CAMLlocal2(r, vtable);

    struct caml_ba_array *ba = Caml_ba_array_val(vba);
    sunrealtype *ba_data = ba->data; // ba is invalid after caml_alloc_final
    int nc = ba->dim[0];
    int nr = ba->dim[1];

    vtable = caml_alloc_final(nc, NULL, 1, 20);
    sunrealtype **table = (sunrealtype **)Data_custom_val(vtable);

    if (nc > 0) {
	int j;
	table[0] = ba_data;
	for (j = 1; j < nc; ++j) {
	    table[j] = table[j - 1] + nr;
	}
    }

    r = caml_alloc_tuple(2);
    Store_field(r, 0, vba);
    Store_field(r, 1, vtable);

    CAMLreturn(r);
}

CAMLprim value sunml_sundials_realarray2_create(int nc, int nr)
{
    CAMLparam0();
    CAMLlocal1(vba);

    vba = caml_ba_alloc_dims(BIGARRAY_FLOAT, 2, NULL, nc, nr);
    CAMLreturn(sunml_sundials_realarray2_wrap(vba));
}

CAMLprim void sunml_crash (value msg)
{
    CAMLparam1 (msg);
    fputs (String_val (msg), stderr);
    fflush (stderr);
    abort ();
    CAMLreturn0;
}

CAMLprim int sunml_sundials_compare_tol(value va, value vb, value vtol)
{
    CAMLparam3(va, vb, vtol);
    CAMLlocal1(vr);
#if 580 <= SUNDIALS_LIB_VERSION
    int r = SUNRCompareTol(Double_val(va), Double_val(vb), Double_val(vtol));
    vr = Int_val(r);
#else
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
#endif
    CAMLreturn(vr);
}

/* Functions for storing OCaml values in the C heap. */

value *sunml_sundials_malloc_value(value v)
{
    header_t *block;
    block = (header_t *)malloc(Bhsize_wosize(1));
    if (block == NULL) return NULL;
#if OCAML_VERSION < 41200
    *block = Make_header(1, 0, Caml_black);
#else
    // see notes in ocaml/runtime/caml/address_class.h
    *block = Caml_out_of_heap_header(1, 0);
#endif
    Field(Val_hp(block), 0) = v;
    caml_register_generational_global_root (Op_hp(block));
    return Op_hp(block);
}

void sunml_sundials_free_value(value *pv)
{
    caml_remove_generational_global_root (pv);
    free (Hp_op(pv));
}

/* Functions for storing pointers to integrators (cvode_mem, ida_mem, etc.) */

static void finalize_session_pointer(value vmem) {
    // there is nothing to finalize, but a distinct function is necessary
    // because caml_final_custom_operations uses the address of this
    // function as a key to find the custom operations table.
}

static int compare_session_pointers(value vmem1, value vmem2)
{
    void* mem1 = SUNML_MEM(vmem1);
    void* mem2 = SUNML_MEM(vmem2);

    // only works for equality, other pointer comparisons are undefined.
    // the comparison is negated to return 0 if equal and 1 otherwise
    return (mem1 != mem2);
}

static struct custom_operations session_pointer_ops = {
    .identifier   = "sunml_session_pointer",
    .finalize     = finalize_session_pointer,
    .compare      = compare_session_pointers,
    .hash         = custom_hash_default,
    .serialize    = custom_serialize_default,
    .deserialize  = custom_deserialize_default,
    .compare_ext  = custom_compare_ext_default,
#if 40800 <= OCAML_VERSION
    .fixed_length = custom_fixed_length_default,
#endif
};

value sunml_wrap_session_pointer(void *sun_mem)
{
    CAMLparam0();
    CAMLlocal1(vmem);

    vmem = caml_alloc_custom(&session_pointer_ops, sizeof(value), 1, 15);
    SUNML_MEM(vmem) = sun_mem;

    CAMLreturn(vmem);
}

/* Functions for sharing OCaml values with C. */

static void sunml_finalize_vptr(value cptr)
{
    value *croot = (*(value **)Data_custom_val(cptr));
    if (croot != NULL) {
	caml_remove_generational_global_root(croot);
	free(croot);
    }
}

CAMLprim value sunml_make_vptr(value v)
{
    CAMLparam1(v);
    CAMLlocal2(cptr, vptr);
    value *croot;

    // create the croot and wrap it
    cptr = caml_alloc_final(1, &sunml_finalize_vptr, 1, 20);
    croot = malloc(sizeof(value));
    (*(value **)Data_custom_val(cptr)) = croot;
    if (croot == NULL) caml_raise_out_of_memory();
    *croot = v;
    caml_register_generational_global_root(croot);

    // create the vptr pair
    vptr = caml_alloc_tuple(2);
    Store_field(vptr, 0, v);
    Store_field(vptr, 1, cptr);

    CAMLreturn(vptr);
}

/* Functions for manipulating FILE pointers. */

static void finalize_cfile(value vf)
{
    FILE *file = ML_CFILE(vf);

    if (file != NULL) {
	fclose(file);
    }
}

CAMLprim value sunml_sundials_fopen(value vpath, value vtrunc)
{
    CAMLparam2(vpath, vtrunc);
    CAMLlocal1(vr);

    char *mode = Bool_val(vtrunc) ? "w" : "a";
    FILE *file = fopen(String_val(vpath), mode);

    if (file == NULL) {
	// uerror("fopen", vpath); /* depends on unix.cma */
	caml_failwith(strerror(errno));
    }

    vr = caml_alloc_final(1, &finalize_cfile, 1, 10);
    ML_CFILE(vr) = file;

    CAMLreturn (vr);
}

CAMLprim void sunml_sundials_write(value vfile, value vdata)
{
    CAMLparam2(vfile, vdata);
    FILE *file = ML_CFILE(vfile);
    size_t len = caml_string_length(vdata);
    size_t w;
#if 40600 <= OCAML_VERSION
    w = fwrite(Bytes_val(vdata), 1, len, file);
#else
    w = fwrite(Bp_val(vdata), 1, len, file);
#endif
    if (w < len) caml_failwith(strerror(errno));
    CAMLreturn0;
}

CAMLprim value sunml_sundials_fflush(value vfile)
{
    CAMLparam1(vfile);
    FILE *file = ML_CFILE(vfile);
    fflush(file);
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_sundials_close(value vfile)
{
    CAMLparam1(vfile);
    FILE *file = ML_CFILE(vfile);
    fclose(file);
    ML_CFILE(vfile) = NULL;
    CAMLreturn (Val_unit);
}

CAMLprim value sunml_sundials_stderr(value vunit)
{
    CAMLparam1(vunit);
    CAMLlocal1(vr);

    vr = caml_alloc_final(1, NULL, 0, 1);
    ML_CFILE(vr) = stderr;

    CAMLreturn (vr);
}

CAMLprim value sunml_sundials_stdout(value vunit)
{
    CAMLparam1(vunit);
    CAMLlocal1(vr);

    vr = caml_alloc_final(1, NULL, 0, 1);
    ML_CFILE(vr) = stdout;

    CAMLreturn (vr);
}

CAMLprim value sunml_sundials_wrap_file(FILE* f)
{
    CAMLparam0();
    CAMLlocal1(vr);

    if (f == NULL) {
	vr = Val_none;
    } else {
	vr = caml_alloc_final(1, NULL, 0, 1);
	ML_CFILE(vr) = f;
	Store_some(vr, vr);
    }

    CAMLreturn (vr);
}

/* Functions for profiling */

#if 600 <= SUNDIALS_LIB_VERSION && defined(SUNDIALS_BUILD_WITH_PROFILING)
static void finalize_profiler(value vprofiler)
{
    SUNProfiler profiler = ML_PROFILER(vprofiler);
    if (profiler != NULL) SUNProfiler_Free(&profiler);
}
#endif

CAMLprim value sunml_profiler_make(value vname)
{
    CAMLparam1(vname);
    CAMLlocal1(vprofiler);

#if 600 <= SUNDIALS_LIB_VERSION && defined(SUNDIALS_BUILD_WITH_PROFILING)
    SUNProfiler profiler = NULL;

    SUNProfiler_Create(NULL, String_val(vname), &profiler);
    if (profiler == NULL) caml_raise_out_of_memory();

    vprofiler = caml_alloc_final(1, &finalize_profiler, 1, 10);
    ML_PROFILER(vprofiler) = profiler;
#else
    vprofiler = Val_unit;
#endif

    CAMLreturn(vprofiler);
}

CAMLprim void sunml_profiler_begin(value vprofiler, value vname)
{
#if 600 <= SUNDIALS_LIB_VERSION && defined(SUNDIALS_BUILD_WITH_PROFILING)
    SUNProfiler_Begin(ML_PROFILER(vprofiler), String_val(vname));
#endif
}

CAMLprim void sunml_profiler_end(value vprofiler, value vname)
{
#if 600 <= SUNDIALS_LIB_VERSION && defined(SUNDIALS_BUILD_WITH_PROFILING)
    SUNProfiler_End(ML_PROFILER(vprofiler), String_val(vname));
#endif
}

CAMLprim void sunml_profiler_print(value vprofiler, value vfile)
{
    CAMLparam2(vprofiler, vfile);
#if 600 <= SUNDIALS_LIB_VERSION && defined(SUNDIALS_BUILD_WITH_PROFILING)
    FILE *file = ML_CFILE(vfile);
    SUNProfiler_Print(ML_PROFILER(vprofiler), file);
#endif
    CAMLreturn0;
}

/* Functions for manipulating contexts */

#if 600 <= SUNDIALS_LIB_VERSION
static void finalize_context(value vctx)
{
    SUNContext ctx = ML_CCONTEXT(vctx);
    if (ctx != NULL) SUNContext_Free(&ctx);
}
#endif

CAMLprim value sunml_context_make(void)
{
    CAMLparam0();
    CAMLlocal1(vctx);

#if 600 <= SUNDIALS_LIB_VERSION
    SUNContext ctx;

    SUNContext_Create(NULL, &ctx);
    if (ctx == NULL) caml_raise_out_of_memory();

    vctx = caml_alloc_final(1, &finalize_context, 1, 10);
    ML_CCONTEXT(vctx) = ctx;
#else
    vctx = Val_unit;
#endif

    CAMLreturn (vctx);
}

CAMLprim void sunml_context_set_profiler(value vctx, value vprofiler)
{
    CAMLparam2(vctx, vprofiler);
#if 600 <= SUNDIALS_LIB_VERSION && defined(SUNDIALS_BUILD_WITH_PROFILING)
    SUNContext_SetProfiler(ML_CCONTEXT(vctx), ML_PROFILER(vprofiler));
#endif
    CAMLreturn0;
}

#ifdef MPI_ENABLED

/* Must correspond with camlmpi.h */
#define Comm_val(comm) (*((MPI_Comm *) &Field(comm, 1)))

CAMLprim value sunml_profiler_make(value vcomm, value vname)
{
    CAMLparam2(vcomm, vname);
    CAMLlocal1(vprofiler);

#if 600 <= SUNDIALS_LIB_VERSION && defined(SUNDIALS_BUILD_WITH_PROFILING)
    SUNProfiler profiler = NULL;
    MPI_Comm comm = Comm_val(vcomm);

    SUNProfiler_Create(comm, String_val(vname), &profiler);
    if (profiler == NULL) caml_raise_out_of_memory();

    vprofiler = caml_alloc_final(1, &finalize_profiler, 1, 10);
    ML_PROFILER(vprofiler) = profiler;
#else
    vprofiler = Val_unit;
#endif

    CAMLreturn(vprofiler);
}

CAMLprim value sunml_context_make_parallel(value vcomm)
{
    CAMLparam1(vcomm);
    CAMLlocal1(vctx);

#if 600 <= SUNDIALS_LIB_VERSION
    SUNContext ctx;
    MPI_Comm comm = Comm_val(vcomm);

    SUNContext_Create(comm, &ctx);
    if (ctx == NULL) caml_raise_out_of_memory();

    vctx = caml_alloc_final(1, &finalize_context, 1, 10);
    ML_CCONTEXT(vctx) = ctx;
#else
    vctx = Val_unit;
#endif

    CAMLreturn (vctx);
}

#endif

