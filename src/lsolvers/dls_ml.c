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

#include <stdio.h>

#include <sundials/sundials_types.h>
#include <sundials/sundials_band.h>

#include <sundials/sundials_dense.h>

#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/fail.h>
#include <caml/bigarray.h>

#include "../sundials/sundials_ml.h"
#include "../lsolvers/lsolver_ml.h"

/* Array dense matrix functions */

CAMLprim value c_arraydensematrix_scale(value vc, value va)
{
    CAMLparam2(vc, va);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

    denseScale(Double_val(vc), ARRAY2_ACOLS(va), m, n);
    CAMLreturn (Val_unit);
}

CAMLprim value c_arraydensematrix_add_identity(value va)
{
    CAMLparam1(va);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[0];

    if (m != n)
	caml_invalid_argument("matrix not square.");
#endif

    denseAddIdentity(ARRAY2_ACOLS(va), m);
    CAMLreturn (Val_unit);
}

#if SUNDIALS_LIB_VERSION >= 260
CAMLprim value c_arraydensematrix_matvec(value va, value vx, value vy)
{
    CAMLparam3(va, vx, vy);
    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

#if SUNDIALS_ML_SAFE == 1
    if (ARRAY1_LEN(vx) < n)
	caml_invalid_argument("x array too small.");
    if (ARRAY1_LEN(vy) < m)
	caml_invalid_argument("y array too small.");
#endif
    denseMatvec(ARRAY2_ACOLS(va), REAL_ARRAY(vx), REAL_ARRAY(vy), m, n);
    CAMLreturn (Val_unit);
}
#else
CAMLprim value c_arraydensematrix_matvec(value va, value vx, value vy)
{
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
}
#endif

CAMLprim value c_arraydensematrix_getrf(value va, value vp)
{
    CAMLparam2(va, vp);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

#if SUNDIALS_ML_SAFE == 1
    if (ARRAY1_LEN(vp) < n)
	caml_invalid_argument("pivot array too small.");
#endif

    int r = denseGETRF(ARRAY2_ACOLS(va), m, n, INDEX_ARRAY(vp));

    if (r != 0) {
	caml_raise_with_arg(LSOLVER_EXN_TAG(ZeroDiagonalElement),
			    Val_long(r));
    }
    CAMLreturn (Val_unit);
}

CAMLprim value c_arraydensematrix_getrs(value va, value vp, value vb)
{
    CAMLparam3(va, vp, vb);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[0];
    if (m != n)
	caml_invalid_argument("matrix not square.");
    if (ARRAY1_LEN(vb) < n)
	caml_invalid_argument("solution vector too small.");
    if (ARRAY1_LEN(vp) < n)
	caml_invalid_argument("pivot array too small.");
#endif

    denseGETRS(ARRAY2_ACOLS(va), m, INDEX_ARRAY(vp), REAL_ARRAY(vb));
    CAMLreturn (Val_unit);
}

CAMLprim value c_arraydensematrix_getrs_off(value va, value vp,
					    value vb, value vboff)
{
    CAMLparam4(va, vp, vb, vboff);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat boff = Int_val(vboff);

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[0];
    if (m != n)
	caml_invalid_argument("matrix not square.");
    if (ARRAY1_LEN(vb) - boff < n)
	caml_invalid_argument("b is too small.");
    if (ARRAY1_LEN(vp) < n)
	caml_invalid_argument("p is too small.");
#endif

    denseGETRS(ARRAY2_ACOLS(va), m, INDEX_ARRAY(vp), REAL_ARRAY(vb) + boff);
    CAMLreturn (Val_unit);
}

CAMLprim value c_arraydensematrix_potrf(value va)
{
    CAMLparam1(va);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[0];
    if (m != n)
	caml_invalid_argument("matrix not square");
#endif

    densePOTRF(ARRAY2_ACOLS(va), m);
    CAMLreturn (Val_unit);
}

CAMLprim value c_arraydensematrix_potrs(value va, value vb)
{
    CAMLparam2(va, vb);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[0];
    if (m != n)
	caml_invalid_argument("matrix not square.");
    if (ARRAY1_LEN(vb) < m)
	caml_invalid_argument("b is too small.");
#endif

    densePOTRS(ARRAY2_ACOLS(va), m, REAL_ARRAY(vb));
    CAMLreturn (Val_unit);
}

CAMLprim value c_arraydensematrix_geqrf(value va, value vbeta, value vv)
{
    CAMLparam3(va, vbeta, vv);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

#if SUNDIALS_ML_SAFE == 1
    if (m < n)
	caml_invalid_argument("fewer rows than columns.");
    if (ARRAY1_LEN(vbeta) < n)
	caml_invalid_argument("beta is too small.");
    if (ARRAY1_LEN(vv) < m)
	caml_invalid_argument("work is too small.");
#endif

    denseGEQRF(ARRAY2_ACOLS(va), m, n, REAL_ARRAY(vbeta), REAL_ARRAY(vv));
    CAMLreturn (Val_unit);
}

CAMLprim value c_arraydensematrix_ormqr(value va, value vormqr)
{
    CAMLparam2(va, vormqr);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[1];
    intnat n = ba->dim[0];

    realtype *beta = REAL_ARRAY(Field(vormqr, 0));
    realtype *vv   = REAL_ARRAY(Field(vormqr, 1));
    realtype *vw   = REAL_ARRAY(Field(vormqr, 2));
    realtype *work = REAL_ARRAY(Field(vormqr, 3));

#if SUNDIALS_ML_SAFE == 1
    if (m < n)
	caml_invalid_argument("fewer rows than columns.");
    if (ARRAY1_LEN(Field(vormqr, 0)) < n)
	caml_invalid_argument("beta is too small.");
    if (ARRAY1_LEN(Field(vormqr, 1)) < n)
	caml_invalid_argument("v is too small.");
    if (ARRAY1_LEN(Field(vormqr, 2)) < m)
	caml_invalid_argument("w is too small.");
    if (ARRAY1_LEN(Field(vormqr, 3)) < m)
	caml_invalid_argument("work is too small.");
#endif

    denseORMQR(ARRAY2_ACOLS(va), m, n, beta, vv, vw, work);
    CAMLreturn (Val_unit);
}

/* Array Band matrix functions */

CAMLprim value c_arraybandmatrix_copy(value va, value vb, value vsizes)
{
    CAMLparam3(va, vb, vsizes);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat am = ba->dim[0];

    long int a_smu  = Long_val(Field(vsizes, 0));
    long int b_smu  = Long_val(Field(vsizes, 1));
    long int copymu = Long_val(Field(vsizes, 2));
    long int copyml = Long_val(Field(vsizes, 3));

#if SUNDIALS_ML_SAFE == 1
    intnat an = ba->dim[1];
    struct caml_ba_array *bb = ARRAY2_DATA(vb);

    intnat bm = bb->dim[0];
    intnat bn = bb->dim[1];

    if (an < copymu + copyml + 1)
	caml_invalid_argument("source matrix too small.");
    if (bn < copymu + copyml + 1)
	caml_invalid_argument("destination matrix too small.");
    if ((am != bm) || (bm != bn))
	caml_invalid_argument("matrix sizes differ.");
#endif

    bandCopy(ARRAY2_ACOLS(va), ARRAY2_ACOLS(vb), am, a_smu, b_smu,
	     copymu, copyml);
    CAMLreturn (Val_unit);
}

CAMLprim value c_arraybandmatrix_scale(value vc, value va, value vsizes)
{
    CAMLparam3(vc, va, vsizes);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[0];

    long int mu  = Long_val(Field(vsizes, 0));
    long int ml  = Long_val(Field(vsizes, 1));
    long int smu = Long_val(Field(vsizes, 2));

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[1];

    if (n < mu + ml + 1)
	caml_invalid_argument("matrix badly sized.");
#endif

    bandScale(Double_val(vc), ARRAY2_ACOLS(va), m, mu, ml, smu);
    CAMLreturn (Val_unit);
}

CAMLprim value c_arraybandmatrix_add_identity(value va, value vsmu)
{
    CAMLparam2(va, vsmu);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[0];
    intnat smu = Long_val(vsmu);

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[1];

    if (n <= smu)
	caml_invalid_argument("matrix badly sized.");
#endif

    bandAddIdentity(ARRAY2_ACOLS(va), m, smu);
    CAMLreturn (Val_unit);
}

#if SUNDIALS_LIB_VERSION >= 260
CAMLprim value c_arraybandmatrix_matvec(value va, value vsizes,
					value vx, value vy)
{
    CAMLparam4(va, vsizes, vx, vy);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[0];

    long int mu  = Long_val(Field(vsizes, 0));
    long int ml  = Long_val(Field(vsizes, 1));
    long int smu = Long_val(Field(vsizes, 2));

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[1];

    if (n < mu + ml + 1)
	caml_invalid_argument("matrix badly sized.");

    if (ARRAY1_LEN(vx) < n)
	caml_invalid_argument("x array too small.");
    if (ARRAY1_LEN(vy) < m)
	caml_invalid_argument("y array too small.");
#endif
    bandMatvec(ARRAY2_ACOLS(va), REAL_ARRAY(vx), REAL_ARRAY(vy),
	       m, mu, ml, smu);
    CAMLreturn (Val_unit);
}
#else
CAMLprim value c_arraybandmatrix_matvec(value va, value vx, value vy,
					value vsizes)
{
    caml_raise_constant(SUNDIALS_EXN(NotImplementedBySundialsVersion));
}
#endif

CAMLprim value c_arraybandmatrix_gbtrf(value va, value vsizes, value vp)
{
    CAMLparam3(va, vsizes, vp);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[0];

    long int mu  = Long_val(Field(vsizes, 0));
    long int ml  = Long_val(Field(vsizes, 1));
    long int smu = Long_val(Field(vsizes, 2));

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[1];

    if (n < mu + ml + 1)
	caml_invalid_argument("matrix badly sized.");
    if (ARRAY1_LEN(vp) < m)
	caml_invalid_argument("p is too small.");
#endif

    bandGBTRF(ARRAY2_ACOLS(va), m, mu, ml, smu, INDEX_ARRAY(vp));
    CAMLreturn (Val_unit);
}

CAMLprim value c_arraybandmatrix_gbtrs(value va, value vsizes, value vp, value vb)
{
    CAMLparam4(va, vsizes, vp, vb);

    struct caml_ba_array *ba = ARRAY2_DATA(va);
    intnat m = ba->dim[0];

    long int smu = Long_val(Field(vsizes, 0));
    long int ml  = Long_val(Field(vsizes, 1));

#if SUNDIALS_ML_SAFE == 1
    intnat n = ba->dim[1];

    if (n < smu + ml + 1)
	caml_invalid_argument("matrix badly sized.");
    if (ARRAY1_LEN(vp) < m)
	caml_invalid_argument("p is too small.");
    if (ARRAY1_LEN(vb) < m)
	caml_invalid_argument("b is too small.");
#endif

    bandGBTRS(ARRAY2_ACOLS(va), m, smu, ml, INDEX_ARRAY(vp), REAL_ARRAY(vb));
    CAMLreturn (Val_unit);
}

