/* Nov 2010, Timothy Bourke (INRIA)
 *
 * Ocaml interface for custom Sundials 2.4.0 CVode NVectors.
 *
 */

#ifndef __NVECTOR_ML_H__
#define __NVECTOR_ML_H__

#include <sundials/sundials_nvector.h>
#include <caml/mlvalues.h>

struct _ml_nvec_content {
    value data;
    value callbacks;
};

typedef struct _ml_nvec_content *ml_nvec_content;

value ml_nvec_new(value mlops, value data);
N_Vector ml_nvec_val(value nvec);

N_Vector callml_vclone(N_Vector w);
N_Vector callml_vcloneempty(N_Vector w);
void callml_vdestroy(N_Vector v);
void callml_vspace(N_Vector v, long int *lrw, long int *liw);
void callml_vlinearsum(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);
void callml_vconst(realtype c, N_Vector z);
void callml_vprod(N_Vector x, N_Vector y, N_Vector z);
void callml_vdiv(N_Vector x, N_Vector y, N_Vector z);
void callml_vscale(realtype c, N_Vector x, N_Vector z);
void callml_vabs(N_Vector x, N_Vector z);
void callml_vinv(N_Vector x, N_Vector z);
void callml_vaddconst(N_Vector x, realtype b, N_Vector z);
realtype callml_vdotprod(N_Vector x, N_Vector y);
realtype callml_vmaxnorm(N_Vector x);
realtype callml_vwrmsnorm(N_Vector x, N_Vector w);
realtype callml_vwrmsnormmask(N_Vector x, N_Vector w, N_Vector id);
realtype callml_vmin(N_Vector x);
realtype callml_vwl2norm(N_Vector x, N_Vector w);
realtype callml_vl1norm(N_Vector x);
void callml_vcompare(realtype c, N_Vector x, N_Vector z);
booleantype callml_vinvtest(N_Vector x, N_Vector z);
booleantype callml_vconstrmask(N_Vector c, N_Vector x, N_Vector m);
realtype callml_vminquotient(N_Vector num, N_Vector denom);

#define NVEC_DATA(nvec) (((ml_nvec_content)nvec->content)->data)
#define NVEC_VAL(v) ((N_Vector)Data_custom_val(v))

#endif

