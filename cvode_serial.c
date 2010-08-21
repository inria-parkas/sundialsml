/* Aug 2010, Timothy Bourke (INRIA) */

#include "caml/mlvalues.h"
#include "caml/bigarray.h"
#include "caml/memory.h"
#include "caml/callback.h"
#include "caml/custom.h"
#include "caml/fail.h"

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_types.h> /* definition of type realtype */

#include <stdio.h>

/* XXX:
 * - realtype must equal double
 */

/*
 * TODO:
 * - call Gc.full_major () from the f and roots routines to see if we get
 *   any segmentation fault problems.
 */

/*
 * Possible extensions:
 * - allow different solvers to be specified (in init)
 * - add a closure for returning the Jacobian.
 * - add handling for flag == CV_TSTOP_RETURN to solver
 * - allow T0 to be configured
 * - add an interface for CVodeSVtolerances(cvode_mem, reltol, abstol);
 *
 */


#define T0    RCONST(0.0)      /* initial time */
#define BIGARRAY_FLAGS (CAML_BA_FLOAT64 | CAML_BA_C_LAYOUT)
#define MAX_ERRMSG_LEN 256

static void check_flag(const char *call, int flag, void *to_free);

// TODO: Is there any risk that the Ocaml GC will try to free the two
//	 closures? Do we have to somehow record that we're using them,
//	 and then release them again in the free routine?
//	 SEE: ml_cvode_data_alloc and finalize
struct ml_cvode_data {
    void *cvode_mem;
    intnat num_roots;
    value *closure_f;
    value *closure_roots;
};

typedef struct ml_cvode_data* ml_cvode_data_p;

static void finalize(value vdata)
{
    ml_cvode_data_p data = (ml_cvode_data_p)Data_custom_val(vdata);

    // TODO:
    // The Ocaml Manual (18.9.1) says:
    // ``Note: the finalize, compare, hash, serialize and deserialize
    // functions attached to custom block descriptors must never trigger a
    // garbage collection. Within these functions, do not call any of the
    // Caml allocation functions, and do not perform a callback into Caml
    // code. Do not use CAMLparam to register the parameters to these
    // functions, and do not use CAMLreturn to return the result.''
    //
    // But, obviously, we're calling two caml functions. We need to find out
    // if this is ok.
    caml_remove_generational_global_root(data->closure_f);
    caml_remove_generational_global_root(data->closure_roots);

    if (data->cvode_mem != NULL) {
	CVodeFree(&(data->cvode_mem));
    }
}

// TODO:
// The Ocaml Manual (18.9.3) says:
// ``The contents of custom blocks are not scanned by the garbage collector,
// and must therefore not contain any pointer inside the Caml heap. In other
// terms, never store a Caml value in a custom block, and do not use Field,
// Store_field nor caml_modify to access the data part of a custom block.
// Conversely, any C data structure (not containing heap pointers) can be
// stored in a custom block.''
//
// But, obviously, we're storing two closure values in the struct. We need
// to find out if and when this is ok.
//
static value ml_cvode_data_alloc(mlsize_t approx_size)
{
    return caml_alloc_final(sizeof(struct ml_cvode_data), &finalize,
			    approx_size, 10);
}


static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    value *closure_f = ((ml_cvode_data_p)user_data)->closure_f;

    intnat y_l = NV_LENGTH_S(y);
    intnat ydot_l = NV_LENGTH_S(ydot);

    value y_ba = caml_ba_alloc(BIGARRAY_FLAGS, 1, NV_DATA_S(y), &y_l);
    value ydot_ba = caml_ba_alloc(BIGARRAY_FLAGS, 1, NV_DATA_S(ydot), &ydot_l);

    // XXX: the data payloads inside y_ba and ydot_ba are only valid during
    //      this call, afterward that memory goes back to cvode. These
    //      bigarrays must not be retained by closure_f! If it wants a
    //      permanent copy, then it has to make it manually.
    //
    //      Eventually y_ba and ydot_ba will be reclaimed by the ocaml gc,
    //      which should not, however, free the attached payload.
    value r = caml_callback3(*closure_f, caml_copy_double(t), y_ba, ydot_ba);
    
    return(Int_val(r));
}

static int roots(realtype t, N_Vector y, realtype *gout, void *user_data)
{
    ml_cvode_data_p data = (ml_cvode_data_p)user_data;

    intnat y_l = NV_LENGTH_S(y);
    value y_ba = caml_ba_alloc(BIGARRAY_FLAGS, 1, NV_DATA_S(y), &y_l);

    value gout_ba = caml_ba_alloc(BIGARRAY_FLAGS, 1, gout, &(data->num_roots));

    // XXX: see notes for f()
    value r = caml_callback3(*(data->closure_roots), caml_copy_double(t),
			     y_ba, gout_ba);
    
    return(Int_val(r));
}

static mlsize_t approx_size_cvode_mem(void *cvode_mem)
{
    mlsize_t used = 0;
    long int lenrw = 0;
    long int leniw = 0;
    int flag = CVodeGetWorkSpace(cvode_mem, &lenrw, &leniw);

    if (flag == CV_SUCCESS) {
    	used = lenrw * sizeof(realtype) + leniw * sizeof(long int);
    }

    return used;
}

 
CAMLprim value c_init(value initial, value num_roots)
{
    CAMLparam2(initial, num_roots);

    int flag;

    int initial_l = Caml_ba_array_val(initial)->dim[0];
    realtype *initial_d = Caml_ba_data_val(initial);
    N_Vector initial_nv = N_VMake_Serial(initial_l, initial_d);

    void *cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);

    value datav = ml_cvode_data_alloc(approx_size_cvode_mem(cvode_mem));
    ml_cvode_data_p data = (ml_cvode_data_p)Data_custom_val(datav);

    data->cvode_mem = cvode_mem;
    data->closure_f = caml_named_value("cvode_serial_callback_f");
    data->closure_roots = caml_named_value("cvode_serial_callback_roots");
    // TODO: check if these two calls are necessary and ok:
    caml_register_generational_global_root(data->closure_f);
    caml_register_generational_global_root(data->closure_roots);

    data->num_roots = Int_val(num_roots);

    if (data->cvode_mem == NULL) {
	free(data);
	caml_failwith("CVodeCreate returned NULL");
	CAMLreturn0;
    }

    flag = CVodeInit(data->cvode_mem, f, T0, initial_nv);
    N_VDestroy(initial_nv);
    check_flag("CVodeInit", flag, data);

    flag = CVodeRootInit(data->cvode_mem, data->num_roots, roots);
    check_flag("CVodeRootInit", flag, data);

    CVodeSetUserData(data->cvode_mem, (void *)data);

    // TODO: where does this really belong? How should it work?
    //	     Is it even necessary?
    N_Vector abstol = N_VNew_Serial(initial_l); 
    int i;
    for (i=0; i < initial_l; ++i) {
	NV_Ith_S(abstol, i) = RCONST(1.0e-8);
    }
    flag = CVodeSVtolerances(data->cvode_mem, RCONST(1.0e-4), abstol);
    N_VDestroy_Serial(abstol);

    CAMLreturn(datav);
}

CAMLprim value c_free(value vdata)
{
    CAMLparam1(vdata);
    finalize(vdata);
    Store_field(vdata, 1, (value)NULL);
    CAMLreturn0;
}

static void check_flag(const char *call, int flag, void *to_free)
{
    static char exmsg[MAX_ERRMSG_LEN] = "";

    if (flag == CV_SUCCESS
	|| flag == CV_ROOT_RETURN
	|| flag == CV_TSTOP_RETURN) return;

    if (to_free != NULL) free(to_free);

    snprintf(exmsg, MAX_ERRMSG_LEN, "%s: %s", call,
	     CVodeGetReturnFlagName(flag));
    caml_failwith(exmsg);
}

static value solver(value vdata, value nextt, value y, int onestep)
{
    CAMLparam2(vdata, nextt);

    realtype t = 0.0;
    ml_cvode_data_p data = (ml_cvode_data_p)Data_custom_val(vdata);
    if (data->cvode_mem == NULL) {
	caml_failwith("This session has been freed");
    }

    int leny = Bigarray_val(y)->dim[0];

    N_Vector y_nv = N_VMake_Serial(leny, Caml_ba_data_val(y));

    // TODO:
    // The payload of y (a big array) must not be shifted by the Ocaml GC
    // during this function call, even though Caml will be reentered
    // through the callback f. Is this guaranteed?
    int flag = CVode(data->cvode_mem, Double_val(nextt), y_nv, &t,
		     onestep ? CV_ONE_STEP : CV_NORMAL);
    N_VDestroy(y_nv);
    check_flag("CVode", flag, NULL);

    value r = caml_alloc_tuple(2);
    Store_field(r, 0, caml_copy_double(t));
    Store_field(r, 1, flag == CV_ROOT_RETURN ? Val_true : Val_false);

    CAMLreturn(r);
}

CAMLprim value c_advance(value vdata, value nextt, value y)
{
    CAMLparam0();
    CAMLreturn(solver(vdata, nextt, y, 0));
}

CAMLprim value c_step(value vdata, value nextt, value y)
{
    CAMLparam0();
    CAMLreturn(solver(vdata, nextt, y, 1));
}

