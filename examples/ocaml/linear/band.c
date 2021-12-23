#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_direct.h>
#include <sundials/sundials_band.h>
#include "../../../src/config.h"
#include "../../../src/sundials/sundials_ml.h"

#define SIZE 5

#define ML 1
#define MU 1
#define SMU 2	/* min(SIZE -1, MU + ML) */

void print_mat(DlsMat m) {
    int i, j;

    for (i=0; i < m->M; ++i) {
	for (j=0; j < m->N; ++j) {
	    if ((i > j + m->ml) || (j > i + m->mu)) {
		printf("       --     ");
	    } else {
		printf(" % e", BAND_ELEM(m, i, j));
	    }
	}
	printf("\n");
    }
}

void print_factored_mat(DlsMat m) {
    int i, j;

    for (i=0; i < m->M; ++i) {
	for (j=0; j < m->N; ++j) {
	    if ((j > i + m->mu) && (j <= i + m->s_mu)) {
		printf(" (% e)", BAND_ELEM(m, i, j));
	    } else if ((i > j + m->ml) || (j > i + m->s_mu)) {
		printf("        --      ");
	    } else {
		printf("  % e ", BAND_ELEM(m, i, j));
	    }
	}
	printf("\n");
    }
}

void print_vec(sunrealtype* m, sundials_ml_index nr)
{
    int i;

    for (i=0; i < nr; ++i) {
	printf(" % e", m[i]);
    }
    printf("\n");
}

void print_pivots(sundials_ml_index* m, sundials_ml_index nr)
{
    int i;

    for (i=0; i < nr; ++i) {
	printf(" % lld", (long long)m[i]);
    }
    printf("\n");
}

int main(int argc, char** argv)
{
#if 600 <= SUNDIALS_LIB_VERSION
    DlsMat a = SUNDlsMat_NewBandMat(SIZE, MU, ML, SMU);
    DlsMat b = SUNDlsMat_NewBandMat(SIZE, MU, ML, SMU);
#else
    DlsMat a = NewBandMat(SIZE, MU, ML, SMU);
    DlsMat b = NewBandMat(SIZE, MU, ML, SMU);
#endif
    sundials_ml_index p[SIZE] = { 0.0 };
    sunrealtype s[SIZE] = { 5.0, 15.0, 31.0, 53.0, 45.0 };

    BAND_ELEM(a,0,0) = 1.0;
    BAND_ELEM(a,0,1) = 2.0;

    BAND_ELEM(a,1,0) = 2.0;
    BAND_ELEM(a,1,1) = 2.0;
    BAND_ELEM(a,1,2) = 3.0;

    BAND_ELEM(a,2,1) = 3.0;
    BAND_ELEM(a,2,2) = 3.0;
    BAND_ELEM(a,2,3) = 4.0;

    BAND_ELEM(a,3,2) = 4.0;
    BAND_ELEM(a,3,3) = 4.0;
    BAND_ELEM(a,3,4) = 5.0;

    BAND_ELEM(a,4,3) = 5.0;
    BAND_ELEM(a,4,4) = 5.0;

    printf("initially: a=\n");
    print_mat(a);
    printf("\n");

#if SUNDIALS_LIB_VERSION >= 260
    {
	sunrealtype x[SIZE] = { 1.0,  2.0, 3.0, 4.0, 5.0 };
	sunrealtype y[SIZE] = { 0.0 };
	printf("matvec: y=\n");
#if 600 <= SUNDIALS_LIB_VERSION
	SUNDlsMat_BandMatvec(a, x, y);
#else
	BandMatvec(a, x, y);
#endif
	print_vec(y, SIZE);
	printf("\n");
    }
#endif

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_BandCopy(a, b, MU, ML);
    SUNDlsMat_BandScale(2.0, b);
#else
    BandCopy(a, b, MU, ML);
    BandScale(2.0, b);
#endif
    printf("scale copy x2: b=\n");
    print_mat(b);
    printf("\n");

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_bandAddIdentity(b->cols, b->N, b->s_mu);
#else
    bandAddIdentity(b->cols, b->N, b->s_mu);
#endif
    printf("add identity: b=\n");
    print_mat(b);
    printf("\n");

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_BandGBTRF(a, p);
#else
    BandGBTRF(a, p);
#endif
    printf("getrf: a=\n");
    print_factored_mat(a);
    printf("\n       p=\n");
    print_pivots(p, SIZE);
    printf("\n");

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_BandGBTRS(a, p, s);
#else
    BandGBTRS(a, p, s);
#endif
    printf("getrs: s=\n");
    print_vec(s, SIZE);

#if 600 <= SUNDIALS_LIB_VERSION
    SUNDlsMat_DestroyMat(a);
    SUNDlsMat_DestroyMat(b);
#else
    DestroyMat(a);
    DestroyMat(b);
#endif

    return 0;
}

