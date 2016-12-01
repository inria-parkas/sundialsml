#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>
#include "../../../src/sundials/sundials_ml.h"

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
#if SUNDIALS_LIB_VERSION >= 270
	for (j = m->indexptrs[i]; j < m->indexptrs[i+1]; ++j) {
	    printf(" (%d: % .02e)", m->indexvals[j], m->data[j]);
	}
#else
	for (j = m->colptrs[i]; j < m->colptrs[i+1]; ++j) {
	    printf(" (%d: % .02e)", m->rowvals[j], m->data[j]);
	}
#endif
	printf("\n");
    }
}

int main(int argc, char** argv)
{
    DlsMat da = NewDenseMat(NROWS, NCOLS);
    SlsMat a;
#if SUNDIALS_LIB_VERSION >= 270
    SlsMat b  = SparseNewMat(NROWS, NCOLS, 2, CSC_MAT);
#else
    SlsMat b  = NewSparseMat(NROWS, NCOLS, 2);
#endif
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

#if SUNDIALS_LIB_VERSION >= 270
    a = SparseFromDenseMat(da, CSC_MAT);
    printf("initially a=\n");
    SparsePrintMat(a, stdout);
#else
    a = SlsConvertDls(da);
    printf("initially a=\n");
    PrintSparseMat(a);
#endif
    print_mat(a);
    printf("\n");

#if SUNDIALS_LIB_VERSION >= 270
    SparseAddIdentityMat(a);
#else
    AddIdentitySparseMat(a);
#endif
    printf("a + 1=\n");
    print_mat(a);
    printf("\n");

#if SUNDIALS_LIB_VERSION >= 270
    SparseCopyMat(a, b);
    SparseScaleMat(2.0, b);
#else
    CopySparseMat(a, b);
    ScaleSparseMat(2.0, b);
#endif
    printf("scale copy x2: b=\n");
    print_mat(b);
    printf("\n");

#if SUNDIALS_LIB_VERSION >= 270
    SparseSetMatToZero(b);
#else
    SlsSetToZero(b);
#endif
    b->NNZ = 9;
#if SUNDIALS_LIB_VERSION >= 270
    b->indexvals = realloc(b->indexvals, b->NNZ*sizeof(int));
#else
    b->rowvals = realloc(b->rowvals, b->NNZ*sizeof(int));
#endif
    b->data    = realloc(b->data, b->NNZ*sizeof(realtype));
    printf("set to zero (NNZ=9): b=\n");
    print_mat(b);
    printf("\n");

#if SUNDIALS_LIB_VERSION >= 270
    b->indexptrs[0] = 0;
    b->indexptrs[1] = 0;
    b->indexptrs[2] = 0;
    b->indexptrs[3] = 1;

    b->indexvals[0] = 1;
#else
    b->colptrs[0] = 0;
    b->colptrs[1] = 0;
    b->colptrs[2] = 0;
    b->colptrs[3] = 1;

    b->rowvals[0] = 1;
#endif
    b->data[0] = 7.0;

    printf("b with element [r=1, c=2] set to 7.0:\n");
    print_mat(b);
    printf("\n");

#if SUNDIALS_LIB_VERSION >= 270
    SparseAddMat(a, b);
#else
    SlsAddMat(a, b);
#endif
    printf("a = a + b: a=\n");
    print_mat(a);
    printf("\n");

#if SUNDIALS_LIB_VERSION >= 270
    SparseMatvec(a, x, y);
#else
    SlsMatvec(a, x, y);
#endif
    printf("y = A*x: x=\n");
    print_vec(x, NCOLS);
    printf("y=\n");
    print_vec(y, NROWS);

    DestroyMat(da);
#if SUNDIALS_LIB_VERSION >= 270
    SparseDestroyMat(a);
    SparseDestroyMat(b);
#else
    DestroySparseMat(a);
    DestroySparseMat(b);
#endif

    return 0;
}

