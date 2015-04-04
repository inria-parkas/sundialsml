#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>

#define NROWS  3
#define NCOLS  3
#define NZEROS 5

realtype a_init[NROWS][NCOLS] = {
    {  0.0,  2.0,  0.0},
    {  0.0, -4.0,  0.0},
    {  0.0, -9.0, -3.0}
};

void print_vec(realtype* m, long int nr)
{
    int i;

    for (i=0; i < nr; ++i) {
	printf(" % e", m[i]);
    }
    printf("\n");
}

void print_pivots(long int* m, long int nr)
{
    int i;

    for (i=0; i < nr; ++i) {
	printf(" % ld", m[i]);
    }
    printf("\n");
}

void print_mat(SlsMat m)
{
    int i, j;

    printf("matrix (M=%d, N=%d, NNZ=%d):\n", m->M, m->N, m->NNZ);
    for (i = 0; i < m->N; ++i) {
	printf("  col %d:", i);
	for (j = m->colptrs[i]; j < m->colptrs[i+1]; ++j) {
	    printf(" (%d: % .02e)", m->rowvals[j], m->data[j]);
	}
	printf("\n");
    }
}

int main(int argc, char** argv)
{
    DlsMat da = NewDenseMat(NROWS, NCOLS);
    SlsMat a;
    SlsMat b  = NewSparseMat(NROWS, NCOLS, 2);
    realtype x[NCOLS] = { 2.0, 3.0, 4.0 };
    realtype y[NROWS];
    int i, j;

    for (i=0; i < NROWS; ++i) {
	for (j=0; j < NCOLS; ++j) {
	    DENSE_ELEM(da, i, j) = a_init[i][j];
	}
    }

    printf("initially da=\n");
    PrintMat(da);
    printf("\n");

    a = SlsConvertDls(da);
    printf("initially a=\n");
    PrintSparseMat(a);
    print_mat(a);
    printf("\n");

    AddIdentitySparseMat(a);
    printf("a + 1=\n");
    print_mat(a);
    printf("\n");

    CopySparseMat(a, b);
    ScaleSparseMat(2.0, b);
    printf("scale copy x2: b=\n");
    print_mat(b);
    printf("\n");

    SlsSetToZero(b);
    printf("set to zero: b=\n");
    print_mat(b);
    printf("\n");

    b->colptrs[0] = 0;
    b->colptrs[1] = 0;
    b->colptrs[2] = 0;
    b->colptrs[3] = 1;

    b->rowvals[0] = 1;
    b->data[0] = 7.0;

    printf("b with element [r=1, c=2] set to 7.0:\n");
    print_mat(b);
    printf("\n");

    SlsAddMat(a, b);
    printf("a = a + b: a=\n");
    print_mat(a);
    printf("\n");

    SlsMatvec(a, x, y);
    printf("y = A*x: x=\n");
    print_vec(x, NCOLS);
    printf("y=\n");
    print_vec(y, NROWS);

    DestroyMat(da);
    DestroySparseMat(a);
    DestroySparseMat(b);

    return 0;
}

