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

#define BANDELEM(A,smu,i,j) (A[j][(i)-(j)+(smu)])

void print_mat_data(realtype** m, sundials_ml_index nr, sundials_ml_index nc) {
    int i, j;

    for (j=0; j < nr; ++j) {
	for (i=0; i < nc; ++i) {
	    printf(" % e", m[j][i]);
	}
	printf("\n");
    }
}

void zero_mat_data(realtype** m, sundials_ml_index nr, sundials_ml_index nc) {
    int i, j;

    for (j=0; j < nr; ++j) {
	for (i=0; i < nc; ++i) {
	    m[j][i] = 0;
	}
    }
}

void print_mat(realtype** m,
	       sundials_ml_index n,
	       sundials_ml_index mu,
	       sundials_ml_index ml,
	       sundials_ml_index smu)
{
    int i, j;

    for (i=0; i < n; ++i) {
	for (j=0; j < n; ++j) {
	    if ((i > j + ml) || (j > i + mu)) {
		printf("       --     ");
	    } else {
		printf(" % e", BANDELEM(m, smu, i, j));
	    }
	}
	printf("\n");
    }
}

void print_factored_mat(realtype** m,
			sundials_ml_index n,
			sundials_ml_index mu,
			sundials_ml_index ml,
			sundials_ml_index smu)
{
    int i, j;

    for (i=0; i < n; ++i) {
	for (j=0; j < n; ++j) {
	    if ((j > i + mu) && (j <= i + smu)) {
		printf(" (% e)", BANDELEM(m, smu, i, j));
	    } else if ((i > j + ml) || (j > i + mu)) {
		printf("        --      ");
	    } else {
		printf("  % e ", BANDELEM(m, smu, i, j));
	    }
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
    realtype **a = newBandMat(SIZE, SMU, ML);
    realtype **b = newBandMat(SIZE, SMU, ML);
    sundials_ml_index p[SIZE] = { 0.0 };
    realtype s[SIZE] = { 5.0, 15.0, 31.0, 53.0, 45.0 };

    zero_mat_data(a, SIZE, SMU + ML + 1);
    zero_mat_data(b, SIZE, SMU + ML + 1);

    a[0][0] = 0.0;
    a[0][1] = 0.0;
    BANDELEM(a,SMU,0,0) = 1.0;
    BANDELEM(a,SMU,0,1) = 2.0;

    BANDELEM(a,SMU,1,0) = 2.0;
    BANDELEM(a,SMU,1,1) = 2.0;
    BANDELEM(a,SMU,1,2) = 3.0;

    BANDELEM(a,SMU,2,1) = 3.0;
    BANDELEM(a,SMU,2,2) = 3.0;
    BANDELEM(a,SMU,2,3) = 4.0;

    BANDELEM(a,SMU,3,2) = 4.0;
    BANDELEM(a,SMU,3,3) = 4.0;
    BANDELEM(a,SMU,3,4) = 5.0;

    BANDELEM(a,SMU,4,3) = 5.0;
    BANDELEM(a,SMU,4,4) = 5.0;

    printf("initially: a.data=\n");
    print_mat_data(a, SIZE, SMU + ML + 1);
    printf("\n");

    printf("initially: a=\n");
    print_mat(a, SIZE, MU, ML, SMU);
    printf("\n");

#if SUNDIALS_LIB_VERSION >= 260
    {
	realtype x[SIZE] = { 1.0,  2.0, 3.0, 4.0, 5.0 };
	realtype y[SIZE] = { 0.0 };
	printf("matvec: y=\n");
	bandMatvec(a, x, y, SIZE, MU, ML, SMU);
	print_vec(y, SIZE);
	printf("\n");
    }
#endif

    bandCopy(a, b, SIZE, SMU, SMU, MU, ML);
    bandScale(2.0, b, SIZE, MU, ML, SMU);
    printf("scale copy x2: b=\n");
    print_mat(b, SIZE, MU, ML, SMU);
    printf("\n");

    bandAddIdentity(b, SIZE, SMU);
    printf("add identity: b=\n");
    print_mat(b, SIZE, MU, ML, SMU);
    printf("\n");

    bandGBTRF(a, SIZE, MU, ML, SMU, p);
    printf("getrf: a=\n");
    print_factored_mat(a, SIZE, MU, ML, SMU);
    printf("\n       p=\n");
    print_pivots(p, SIZE);
    printf("\n");

    bandGBTRS(a, SIZE, SMU, ML, p, s);
    printf("getrs: s=\n");
    print_vec(s, SIZE);

    destroyMat(a);
    destroyMat(b);

    return 0;
}

