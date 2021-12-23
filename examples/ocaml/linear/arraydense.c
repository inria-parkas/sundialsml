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

void print_mat(sunrealtype** m, sundials_ml_index nr, sundials_ml_index nc) {
    int i, j;

    for (i=0; i < nr; ++i) {
	for (j=0; j < nc; ++j) {
	    printf(" % e", m[j][i]);
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
    sunrealtype **a = SUNDlsMat_newDenseMat(NROWS, NCOLS);
    sunrealtype **b = SUNDlsMat_newDenseMat(NROWS, NCOLS);
#else
    sunrealtype **a = newDenseMat(NROWS, NCOLS);
    sunrealtype **b = newDenseMat(NROWS, NCOLS);
#endif
    sundials_ml_index p[NROWS] = { 0.0 };
    sunrealtype s[NROWS] = { 5.0, 18.0, 6.0 };
    int i, j;

    for (i=0; i < NROWS; ++i) {
	for (j=0; j < NCOLS; ++j) {
	    a[j][i] = a_init[i][j];
	}
    }

    printf("initially: a=\n");
    print_mat(a, NROWS, NCOLS);
    printf("\n");

#if 260 <= SUNDIALS_LIB_VERSION
    {
	sunrealtype x[NCOLS] = { 1.0,  2.0, 3.0 };
	sunrealtype y[NROWS] = { 0.0 };
	printf("matvec: y=\n");
#if 600 <= SUNDIALS_LIB_VERSION
	SUNDlsMat_denseMatvec(a, x, y, NROWS, NCOLS);
#else
	denseMatvec(a, x, y, NROWS, NCOLS);
#endif
	print_vec(y, NROWS);
	printf("\n");
    }
#endif

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_denseCopy(a, b, NROWS, NCOLS);
    SUNDlsMat_denseScale(2.0, b, NROWS, NCOLS);
#else
    denseCopy(a, b, NROWS, NCOLS);
    denseScale(2.0, b, NROWS, NCOLS);
#endif
    printf("scale copy x2: b=\n");
    print_mat(b, NROWS, NCOLS);
    printf("\n");

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_denseAddIdentity(b, NROWS);
#else
    denseAddIdentity(b, NROWS);
#endif
    printf("add identity: b=\n");
    print_mat(b, NROWS, NCOLS);
    printf("\n");

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_denseGETRF(a, NROWS, NCOLS, p);
#else
    denseGETRF(a, NROWS, NCOLS, p);
#endif
    printf("getrf: a=\n");
    print_mat(a, NROWS, NCOLS);
    printf("\n       p=\n");
    print_pivots(p, NROWS);
    printf("\n");

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_denseGETRS(a, NROWS, p, s);
#else
    denseGETRS(a, NROWS, p, s);
#endif
    printf("getrs: s=\n");
    print_vec(s, NROWS);

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_destroyMat(a);
    SUNDlsMat_destroyMat(b);
#else
    destroyMat(a);
    destroyMat(b);
#endif

    return 0;
}

