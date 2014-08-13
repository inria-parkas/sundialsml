#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_direct.h>
#include <sundials/sundials_band.h>

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
    DlsMat a = NewBandMat(SIZE, MU, ML, SMU);
    DlsMat b = NewBandMat(SIZE, MU, ML, SMU);
    long int p[SIZE] = { 0.0 };
    realtype s[SIZE] = { 5.0, 15.0, 31.0, 53.0, 45.0 };

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

    printf("initially a=\n");
    print_mat(a);
    printf("\n");

    BandCopy(a, b, MU, ML);
    BandScale(2.0, b);
    printf("scale copy x2: b=\n");
    print_mat(b);
    printf("\n");

    bandAddIdentity(b->cols, b->N, b->s_mu);
    printf("add identity: b=\n");
    print_mat(b);
    printf("\n");

    BandGBTRF(a, p);
    printf("getrf: a=\n");
    print_factored_mat(a);
    printf("\n       p=\n");
    print_pivots(p, SIZE);
    printf("\n");

    BandGBTRS(a, p, s);
    printf("getrs: s=\n");
    print_vec(s, SIZE);

    return 0;
}

