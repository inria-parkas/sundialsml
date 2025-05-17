// cc -isystem SUNDIALS_DIR_INCLUDE -o compare_inc_with_lib compare_inc_with_lib.c -lsundials_generic

#include <stdio.h>
#include <string.h>
#include "sundials/sundials_version.h"
#include "sundials/sundials_config.h"

#define LEN_VERSION_STR 10

static char get_version[LEN_VERSION_STR];

int main(int argc, char* argv[])
{
    int verbose = (argc == 2) && (strcmp(argv[1], "-v") == 0);

    if (verbose) printf("SUNDIALS_VERSION=%s\n", SUNDIALS_VERSION);

    if (SUNDIALSGetVersion(get_version, LEN_VERSION_STR) != 0) {
	printf("SUNDIALSGetVersion failed\n");
	return 1;
    }

    if (verbose) printf("SUNDIALSGetVersion=%s\n", get_version);

    if (strncmp(get_version, SUNDIALS_VERSION, LEN_VERSION_STR) != 0) {
	printf("SUNDIALSGetVersion=%s != SUNDIALS_VERSION=%s\n",
	       get_version, SUNDIALS_VERSION);
    }

    return 0;
}

