Sundials/ML 2.6.2p0 (March 2016)
--------------------------------
Sundials/ML v2.6.2p0 adds support for v2.6.x of the Sundials Suite of
numerical solvers, including:
* the new ARKODE solver,
* sparse matrices and the KLU and SuperLU/MT linear solvers,
* OpenMP and Pthreads nvectors, and
* various new functions and linear solvers in existing solvers.

OCaml 3.12.1 or greater is required, and optionally OCamlMPI 1.01.

We continue to provide support for the Sundials 2.5.x series which is still
found in many packaging systems (like Homebrew and Debian/Ubuntu).

Notes:
* The source files have been reorganized into subdirectories.
* Sensitivity features are now disabled via a findlib predicate.
* The Spils jac_times_vec function is no longer associated with individual
  preconditioners, but rather with the linear solvers directly.
* Adjoint linear solver callbacks in CVODES and IDAS may now depend on
  forward sensitivities, the types dense_jac_fn, band_jac_fn, and
  jac_times_vec_fn become variants, new preconditioner functions are
  provided.
* The Kinsol interface changes for new features (new strategies and Anderson
  iteration).
* Incompatibility: The {Cvode,Ida,...}.serial_session type synonyms gain a 
  polymorphic variable to admit OpenMP and Pthreads nvectors.

Sundials/ML 2.5.0p0 (November 2014)
-----------------------------------
Sundials/ML v2.5.0p0 is an OCaml interface to v2.5.0 of the Sundials suite
of numerical solvers (CVODE, CVODES, IDA, IDAS, KINSOL).

It requires OCaml 3.12.1 or greater, Sundials 2.5.0, and optionally
OCamlMPI 1.01.

* When building Sundials manually, we recommend applying the
  `sundials-2.5.0.patch` file and building with examples and shared library
  support:

      patch -p1 < path/to/sundials-2.5.0.patch
      ./configure --enable-examples --enable-shared

  Sundials/ML will function correctly if the patch is not applied, but some
  examples will fail with incorrect results.

* The backward preconditioner, banded and dense jacobian, and jacobian
  times vector callbacks in Cvodes.Adjoint and Idas.Adjoint do not function
  correctly due to an issue in the underlying C library.

