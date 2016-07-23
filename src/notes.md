Data structures shared between OCaml and C
==========================================

Session types
-------------

Four session types are used to interface with the six solvers:
`Cvode_impl.session`, `Ida_impl.session`, `Kinsol_impl.session`,
and `Arkode_impl.session`.

They each contain two types of pointers into the C heap:

* `cvode`/`ida`/`kinsol`/`arkode`: pointer to session record. It is created
  within the `c_*_init` functions by a call to the appropriate Sundials
  create function (e.g., `CVodeCreate`). The create functions return a
  pointer to `malloc`ed memory, which is then wrapped as a custom block.
  This block has no finalizer since the finalizer of the enclosing OCaml
  session value frees the associated memory (e.g., `CVodeFree`). The two
  objects cannot exist independently (due to user data and callbacks
  configure across them both). These fields are accessed frequently on the
  C side via the macros `*_MEM_FROM_ML`. They are never accessed on the
  OCaml side.

* `backref`: pointer to a global root containing a weak pointer back to the
  session value. The `c_sundials_malloc_value` function `malloc`s and
  configures an OCaml block on the C-heap and registers its only value as a
  global root (the header is valid and marked black to satisfy
  -no-naked-pointers). This function is called by the `c_*_init` functions.
  A `c_sundials_Free_value` function deregisters the global root and `free`s
  the memory. It is called by the `c_*_session_finalize` functions. This
  field is also used by the `c_*_set_err_handler_fn` functions to retrieve
  the weak session pointer. Access is via the `*_BACKREF_FROM_ML` macro. It
  is never used from OCaml.

Note that the sensitivity solvers also contain pointers to the session
record (from, e.g., `CVodeSetUserDataB`), to a `backref`, and to file
handles; see the `c_cvodes_adj_init_backward` and c_idas_adj_init_backward`
functions.

File handles
------------
The Sundials.Logfile module provides limited access to libC file handles 
(`FILE *`) for use in setting error, informational, and diagnostic log 
files. These handles (created by calls to `fopen`) are wrapped in custom 
blocks.

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

