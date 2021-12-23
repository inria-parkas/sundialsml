#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_direct.h>
#include <sundials/sundials_dense.h>
#include "../../../src/config.h"
#include "../../../src/sundials/sundials_ml.h"

#define NROWS 3
#define NCOLS 3

sunrealtype a_init[NROWS][NCOLS] = {
    {  1.0,  2.0,  3.0},
    {  2.0, -4.0,  6.0},
    {  3.0, -9.0, -3.0}
};

void print_mat(DlsMat m, sundials_ml_index nr, sundials_ml_index nc) {
    int i, j;

    for (i=0; i < nr; ++i) {
	for (j=0; j < nc; ++j) {
	    printf(" % e", DENSE_ELEM(m, i, j));
	}
	printf("\n");
    }
}

void print_vec(sunrealtype* m, sundials_ml_index nr) {
    int i;

    for (i=0; i < nr; ++i) {
	printf(" % e", m[i]);
    }
    printf("\n");
}

void print_pivots(sundials_ml_index* m, sundials_ml_index nr) {
    int i;

    for (i=0; i < nr; ++i) {
	printf(" % lld", (long long)m[i]);
    }
    printf("\n");
}

int main(int argc, char** argv)
{
#if 600 <= SUNDIALS_LIB_VERSION
    DlsMat a = SUNDlsMat_NewDenseMat(NROWS, NCOLS);
    DlsMat b = SUNDlsMat_NewDenseMat(NROWS, NCOLS);
#else
    DlsMat a = NewDenseMat(NROWS, NCOLS);
    DlsMat b = NewDenseMat(NROWS, NCOLS);
#endif
    sundials_ml_index p[NROWS] = { 0.0 };
    sunrealtype s[NROWS] = { 5.0, 18.0, 6.0 };
    int i, j;

    for (i=0; i < NROWS; ++i) {
	for (j=0; j < NCOLS; ++j) {
	    DENSE_ELEM(a, i, j) = a_init[i][j];
	}
    }

    printf("initially: a=\n");
    print_mat(a, NROWS, NCOLS);
    printf("\n");

#if SUNDIALS_LIB_VERSION >= 260
    {
	sunrealtype x[NCOLS] = { 1.0,  2.0, 3.0 };
	sunrealtype y[NROWS] = { 0.0 };
	printf("matvec: y=\n");
#if 600 <= SUNDIALS_LIB_VERSION
	SUNDlsMat_DenseMatvec(a, x, y);
#else
	DenseMatvec(a, x, y);
#endif
	print_vec(y, NROWS);
	printf("\n");
    }
#endif

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_DenseCopy(a, b);
#else
    DenseCopy(a, b);
#endif

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_DenseScale(2.0, b);
#else
    DenseScale(2.0, b);
#endif
    printf("scale copy x2: b=\n");
    print_mat(b, NROWS, NCOLS);
    printf("\n");

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_denseAddIdentity(b->cols, NROWS);
#else
    denseAddIdentity(b->cols, NROWS);
#endif
    printf("add identity: b=\n");
    print_mat(b, NROWS, NCOLS);
    printf("\n");

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_DenseGETRF(a, p);
#else
    DenseGETRF(a, p);
#endif
    printf("getrf: a=\n");
    print_mat(a, NROWS, NCOLS);
    printf("\n       p=\n");
    print_pivots(p, NROWS);
    printf("\n");

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_DenseGETRS(a, p, s);
#else
    DenseGETRS(a, p, s);
#endif
    printf("getrs: s=\n");
    print_vec(s, NROWS);

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_DestroyMat(a);
    SUNDlsMat_DestroyMat(b);
#else
    DestroyMat(a);
    DestroyMat(b);
#endif

    return 0;
}

