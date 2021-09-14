Sundials/ML 5.7.0p0 (September 2021)
------------------------------------
Sundials/ML v5.7.0p0 adds support for v5.x of the Sundials Suite of
numerical solvers.

OCaml 4.03.0 or greater is now required and optionally OCamlMPI 1.05.

Notes:
* Reworking of the Arkode.MRIStep interface, which now requires the 
  underlying library to be v5.x or above.
* Changes to the callback interfaces NonlinearSolver.{lsetupfn, lsolvefn} 
  and the name (and use) of the second argument of the function 
  NonlinearSolver.solve.
* Support for most new 5.x features in Cvode, Cvodes, Ida, Idas, and Kinsol
  including the new ManyVector, MPIManyVector, and MPIPlusX nvectors.

Compatibility
* The naming scheme for nvector operations has changed. Operations no longer 
  commence with the three characters `n_v`. For instance, `n_vdiv` becomes 
  simply `div`.

Sundials/ML 4.1.0p0 (September 2020)
------------------------------------
Sundials/ML v4.1.0p0 adds support for v4.x of the Sundials Suite of
numerical solvers.

Notes:
* New Sundials.Nonlinear module and corresponding changes to the Cvode, 
  Cvodes, Ida, Idas, Kinsol, and ARKode integrators.
* Support for the new fused and array nvector operations.
* Testing of all OCaml nvector operations and several bug fixes.
* Removal of the _Solver_.Alternate modules for attaching custom linear
  solvers. This interface was superseded since Sundials 3.x by the new 
  linear solver interface which provides similar functionality.
* ARKode: reworked interface to new ARKStep, ERKStep, and MRIStep modules.
* LinearSolver: use a generic linear solver type, replacing the previous 
  distinction between Dls and Iterative linear solvers.
* It is best to avoid specifying Nonlinear solvers with 4.0.0 <= Sundials <= 
  4.0.2 as there is a bug that can lead to OCaml-owned data structures being 
  freed (e.g., examples/cvode/serial/cvDirectDemo_ls.ml sometimes crashes 
  with a segmentation error)
* The arkode/C_openmp/ark_head1D_omp.c example gives inconsistent results; 
  they change depending on the number of OpenMP threads.

Compatibility:
* When initializing Cvode, the new non-linear solver interface replaces the 
  previous iteration argument, so
    Cvode.init lmm Cvode.Functional tol ...
  becomes
    Cvode.init lmm tol ~nlsolver:(NonlinearSolver.FixedPoint.make y0 0) ...
  Otherwise, not passing an nlsolver argument specifies a default nonlinear 
  solver based on Newton iteration. The linear solver must now be specified 
  with the `~lsolver` label. Note that the tolerances must now be passed as 
  the second argument.

* When initializing Ida, the tolerance argument is now given first (for 
  consistency with Cvode) and the linear solver arguments has a label (to 
  avoid confusion with the nonlinear solver and for consistency with other
  integrators where this argument is optional), so
    let ida = Ida.init solver (Ida.SStolerances (1e-9, 1e-9)) ...
  becomes
    let ida = Ida.init (Ida.SStolerances (1e-9, 1e-9)) ~lsolver:solver ...

* For Kinsol.init, the ?linsolv argument has been renamed to ?lsolver for 
  consistency with the other solvers.

* The sensitivity method for the Cvodes/Idas.Sensitivity solvers now takes a 
  nonlinear solver as an optional argument.

* When initializing Arkode(.ARKStep), the call
    let arkode_mem = Arkode.(
      init
        (Implicit (f udata,
                   Newton Spils.(solver (pcg ~maxl:n_mesh y)
                                        ~jac_times_vec:(None, jac)
                                        prec_none),
                   linearity))
        (SStolerances (rtol, atol))
        t0
        y
    ) ...
  becomes
    let arkode_mem = Arkode.ARKStep.(
      init
        (implicit
          ~lsolver:Spils.(solver (pcg ~maxl:n_mesh y)
                                 ~jac_times_vec:(None, jac)
                                 prec_none)
          ~linearity
          (f udata))
        (SStolerances (rtol, atol))
        t0
        y
    ) ...

* *.Dls.get_num_rhs_evals -> *.Dls.get_num_lin_rhs_evals
* *.Dls.get_num_func_evals -> *.Dls.get_num_lin_func_evals
* *.Spils.get_num_rhs_evals -> *.Spils.get_num_lin_rhs_evals
* *.Spils.get_num_conv_fails -> *.Spils.get_num_lin_conv_fails
* *.Spils.get_num_func_evals -> *.Spils.get_num_lin_func_evals

* *.blit now requires ~src and ~dst labels.
* *.blit_some becomes *.blitn and requires labels.

* Sundials.LinearSolver.linear_solver -> Sundials.LinearSolver.t
* *.session_linear_solver -> *.linear_solver

Sundials/ML 3.1.1p0 (July 2018)
------------------------------------
Sundials/ML v3.1.1p0 adds support for v3.1.x of the Sundials Suite of
numerical solvers.

Notably this release adds support for the new generic matrix and linear 
solver interfaces. The OCaml interface changes but the library is backward 
compatible with Sundials 2.7.0.

OCaml 4.02.3 or greater is now required and optionally OCamlMPI 1.03.

Notes:
* New Sundials.Matrix and Sundials.LinearSolver modules.
* Better treatment of integer type used for matrix indexing.
* Refactor Dls and Sls modules into Sundials.Matrix.
* Add confidence intervals to performance graph.
* Miscellaneous improvements to configure script.
* Potential incompatibility: changes to some label names: comm_fn -> comm;
  iter_type -> iter.
* Untangle the ARKODE mass-solver interface from the Jacobian interface.

Sundials/ML 2.7.0p0 (December 2016)
------------------------------------
Sundials/ML v2.7.0p0 adds support for v2.7.x of the Sundials Suite of
numerical solvers.

Notes:
* Arkode: the interfaces to the Butcher tables have changed.
* The sparse matrix interface has changed:
  Sls.SparseMatrix:
    make       -> make_csc
    create     -> create_csc
    from_dense -> csc_from_dense
    from_band  -> csc_from_band
* The Klu and Superlumt linear solver interfaces have changed.
    *.Klu.solver -> Klu.solver_csc
    *.Superlumt.solver -> Superlumt.solver_csc

Sundials/ML 2.6.2p1 (September 2016)
------------------------------------
Sundials/ML v2.6.2p1 includes several bug fixes and minor 
additions/improvements:
* Add pretty printers with automatic installation
  (thanks to Nils Becker for the suggestion).
* Improve Opam integration allowing pin from source code
  (thanks to Gabriel Scherer for the suggestion).
* Ensure compatibility with OCaml no-naked-pointers mode.
* Fix segfaulting on exceptions in newer versions of OCaml.
* Fix bug in RealArray2.size.
* Update the set_err_file/set_info_file/set_diagnostics interface
  (minor incompatibility).
* Miscellaneous improvements to the build system.
* Remove the Kinsol.set_linear_solver function due to
  [memory leak issues] [1].

[1]: http://sundials.2283335.n4.nabble.com/KINSOL-documentation-td4653693.html

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

