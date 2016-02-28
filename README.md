Sundials/ML
===========

Sundials/ML is an interface to the Sundials suite of numerical solvers.
[Full documentation](http://inria-parkas.github.io/sundialsml/) is online.

[Sundials](http://computation.llnl.gov/casc/sundials/main.html) is a 
collection of six numerical solvers:
CVODE, CVODES, IDA, IDAS, ARKODE, and KINSOL.
This interface provides access to all features of the underlying library.

Installation
------------
1. [Download Sundials](http://computation.llnl.gov/casc/sundials/download/download.php), extract, and install it:
    1. `mkdir sundials-build`
    2. `cd sundials-build`
    3. `cmake -Wno-dev ../sundials-2.6.3`, optionally adding:
        - `-DOPENMP_ENABLE=1` for OpenMP nvectors,
        - `-DPTHREAD_ENABLE=1` for Pthreads nvectors,
        - `-DMPI_ENABLE=1` for parallel nvectors,
        - `-DLAPACK_ENABLE=1` to use LAPACK routines,
        - `-DSUPERLUMT_ENABLE=1` for 
          [SuperLUMT](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/#superlu_mt) 
          solvers (with `-DSUPERLUMT_THREAD_TYPE=Pthread` or `OpenMP`,
           `-DSUPERLUMT_LIBRARY_DIR=...`, `-DSUPERLUMT_INCLUDE_DIR=...`), 
           and
        - `-DKLU_ENABLE=1` for 
          [KLU](http://faculty.cse.tamu.edu/davis/suitesparse.html) solvers
          (with `-DKLU_LIBRARY_DIR=...`, `-DKLU_INCLUDE_DIR=...`)
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
