Sundials/ML
===========

Sundials/ML is an interface to the Sundials suite of numerical solvers.

[Full documentation](http://inria-parkas.github.io/sundialsml/) is online.

[Sundials](http://computation.llnl.gov/casc/sundials/main.html) is a 
collection of six numerical solvers:
CVODE, CVODES, IDA, IDAS, ARKODE, and KINSOL.
This interface provides access to all features of the underlying library 
except the Hypre, PETSC, CUDA, RAJA, and OpenMPdev nvectors (since these 
require other libraries that are not yet interfaced for OCaml).

Installation
------------
See our [detailed notes](http://inria-parkas.github.io/sundialsml/#running), 
or try:

1. [Download Sundials](http://computation.llnl.gov/casc/sundials/download/download.php), extract, and install it:
    1. `mkdir sundials-build`
    2. `cd sundials-build`
    3. `cmake -Wno-dev ../sundials-3.1.1`, optionally adding:
        - `-DOPENMP_ENABLE=1` for OpenMP nvectors,
        - `-DPTHREAD_ENABLE=1` for Pthreads nvectors,
        - `-DMPI_ENABLE=1` for parallel nvectors,
        - `-DLAPACK_ENABLE=1` to use LAPACK routines,
        - `-DSUPERLUMT_ENABLE=1` for 
          [SuperLUMT](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/#superlu_mt) 
          solvers, and
        - `-DKLU_ENABLE=1` for 
          [KLU](http://faculty.cse.tamu.edu/davis/suitesparse.html) solvers.
    4. `make install`
2. Either install from OPAM: `opam install sundialsml`, or
    1. [Download Sundials/ML](https://github.com/inria-parkas/sundialsml/releases), extract, and install it:
    2. `./configure`
    3. `make`
    4. `make install` or `make install-findlib`
3. [Start coding!](http://inria-parkas.github.io/sundialsml/#running)

Contact
-------
* [Support (public OCaml list)](mailto:caml-list@inria.fr?subject=Sundials/ML:)
* [Bug reports/Feature requests](https://github.com/inria-parkas/sundialsml/issues/new)
* [Pull requests](https://github.com/inria-parkas/sundialsml/compare)

Known Limitations
-----------------
Our goal is to provide access from OCaml to as much of Sundials as we can. 
We do not, however, expose features for which there is no corresponding 
OCaml library (e.g., some of the nvector modules) or whose implementation 
would be overly complicated (e.g., in terms of types or mixed C/OCaml data 
structures) or inefficient.

* Unsupported nvector modules: CUDA, Hypre ParVector, PETSC, RAJA, and 
  Trilinos (if you provide the required OCaml interface, we'll add the 
  nvector).
* No support for CUDA (including SUNMATRIX_CUSPARSE and 
  SUNLinearSolver_cuSolverSp_batchQR linear solver).
* No support for fused kernels.
* No support for the SUNLinSol_KLUGetSymbolic, SUNLinSol_KLUGetNumeric,
  and SUNLinSol_KLUGetCommon functions (patches welcome).
* No support for SUNMatrix_SLUNRloc and SUNLinearSolver_SuperLUDIST (if you 
  provide a solid OCaml interface to SuperLU_DIST, we'll add support for 
  this matrix type and linear solver).

