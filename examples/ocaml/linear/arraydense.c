#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_direct.h>
#include <sundials/sundials_dense.h>

#define NROWS 3
#define NCOLS 3

realtype a_init[NROWS][NCOLS] = {
    {  1.0,  2.0,  3.0},
    {  2.0, -4.0,  6.0},
    {  3.0, -9.0, -3.0}
};

void print_mat(realtype** m, long int nr, long int nc) {
    int i, j;

    for (i=0; i < nr; ++i) {
	for (j=0; j < nc; ++j) {
	    printf(" % e", m[j][i]);
	}
	printf("\n");
    }
}

void print_vec(realtype* m, long int nr) {
    int i;

    for (i=0; i < nr; ++i) {
	printf(" % e", m[i]);
    }
    printf("\n");
}

void print_pivots(long int* m, long int nr) {
    int i;

    for (i=0; i < nr; ++i) {
	printf(" % ld", m[i]);
    }
    printf("\n");
}

int main(int argc, char** argv)
{
    realtype **a = newDenseMat(NROWS, NCOLS);
    realtype **b = newDenseMat(NROWS, NCOLS);
    long int p[NROWS] = { 0.0 };
    realtype s[NROWS] = { 5.0, 18.0, 6.0 };
    int i, j;

    for (i=0; i < NROWS; ++i) {
	for (j=0; j < NCOLS; ++j) {
	    a[j][i] = a_init[i][j];
	}
    }

    printf("initially a=\n");
    print_mat(a, NROWS, NCOLS);
    printf("\n");

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

