Data structures shared between OCaml and C
==========================================

Session types
-------------

Four session types are used to interface with the six solvers:
`Cvode_impl.session`, `Ida_impl.session`, `Kinsol_impl.session`,
and `Arkode_impl.session`.

They each contain three types of pointers into the C heap:

* `cvode`/`ida`/`kinsol`/`arkode`: pointer to session record. It is created
  within the `c_*_init` functions by a call to the appropriate Sundials
  create function (e.g., `CVodeCreate`), which return a pointer to
  `malloc`ed memory. These fields are accessed frequently on the C side
  via the macros `*_MEM_FROM_ML`. They are never accessed on the OCaml side.

* `backref`: pointer to a global root containing a weak pointer back to the
  session value. It is `malloc`ed and set by the `c_*_init` functions. It is
  used by the `c_*_session_finalize` functions to remove the GC root before
  its memory is freed. This field is never used from OCaml. It is also used
  by the `c_*_set_err_handler_fn` functions to retrieve the weak session
  pointer. Access is via the `*_BACKREF_FROM_ML` macro.

* `err_file`, `info_file`, `diag_file`: pointer to file handle (`FILE *`)
  returned by `fopen` (called within our code). Accessed from C in the
  functions `c_*_init` (to set to NULL), `c_*_session_finalize` (to close
  the handle) and `c_*_set_error_file` (to open the handle).

Nvectors
--------

An `Nvector.t` (used by all nvector types: serial, parallel, Pthreads,
OpenMP) contains a `cnvec` value, which is realised as a custom block with a
finalizer around an `N_Vector`. The `NVEC_CVAL` macro is used to access the
payload of such a value (or `NVEC_VAL` to access it via a field within a
value of type `Nvector.t`).

Matrices
--------

There are two types of matrix pointer: `dlsmat` (`DenseMatrix.t` and
`BandMatrix.t`) and `slsmat` (`SparseMatrix.t`). Both are created as custom
blocks using `caml_alloc_final` (by `c_dls_dense_wrap` and
`c_sls_sparse_wrap`).

SPILS Solvers
-------------

Pointers can be created to sessions of the various SPILS solvers (e.g., the
field of type `Spils.SPGMR.memrec` within `Spils.SPGMR.t`). These pointers
are created as custom blocks using `caml_alloc_final` (by `c_spils_*_make`).

