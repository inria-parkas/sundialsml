#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_direct.h>
#include <sundials/sundials_dense.h>
#include "../../../src/config.h"
#include "../../../src/sundials/sundials_ml.h"

#define NROWS 3
#define NCOLS 3

realtype a_init[NROWS][NCOLS] = {
    {  1.0,  2.0,  3.0},
    {  2.0, -4.0,  6.0},
    {  3.0, -9.0, -3.0}
};

void print_mat(realtype** m, sundials_ml_index nr, sundials_ml_index nc) {
    int i, j;

    for (i=0; i < nr; ++i) {
	for (j=0; j < nc; ++j) {
	    printf(" % e", m[j][i]);
	}
	printf("\n");
    }
}

void print_vec(realtype* m, sundials_ml_index nr) {
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
    realtype **a = newDenseMat(NROWS, NCOLS);
    realtype **b = newDenseMat(NROWS, NCOLS);
    sundials_ml_index p[NROWS] = { 0.0 };
    realtype s[NROWS] = { 5.0, 18.0, 6.0 };
    int i, j;

    for (i=0; i < NROWS; ++i) {
	for (j=0; j < NCOLS; ++j) {
	    a[j][i] = a_init[i][j];
	}
    }

    printf("initially: a=\n");
    print_mat(a, NROWS, NCOLS);
    printf("\n");

#if SUNDIALS_LIB_VERSION >= 260
    {
	realtype x[NCOLS] = { 1.0,  2.0, 3.0 };
	realtype y[NROWS] = { 0.0 };
	printf("matvec: y=\n");
	denseMatvec(a, x, y, NROWS, NCOLS);
	print_vec(y, NROWS);
	printf("\n");
    }
#endif

    denseCopy(a, b, NROWS, NCOLS);
    denseScale(2.0, b, NROWS, NCOLS);
    printf("scale copy x2: b=\n");
    print_mat(b, NROWS, NCOLS);
    printf("\n");

    denseAddIdentity(b, NROWS);
    printf("add identity: b=\n");
    print_mat(b, NROWS, NCOLS);
    printf("\n");

    denseGETRF(a, NROWS, NCOLS, p);
    printf("getrf: a=\n");
    print_mat(a, NROWS, NCOLS);
    printf("\n       p=\n");
    print_pivots(p, NROWS);
    printf("\n");

    denseGETRS(a, NROWS, p, s);
    printf("getrs: s=\n");
    print_vec(s, NROWS);

    destroyMat(a);
    destroyMat(b);

    return 0;
}

