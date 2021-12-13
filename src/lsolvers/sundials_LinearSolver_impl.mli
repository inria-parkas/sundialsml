val e : exn
type linear_solver_type =
    Direct
  | Iterative
  | MatrixIterative
  | MatrixEmbedded
type linear_solver_id =
    Band
  | Dense
  | Klu
  | LapackBand
  | LapackDense
  | Pcg
  | Spbcgs
  | Spfgmr
  | Spgmr
  | Sptfqmr
  | Superludist
  | Superlumt
  | Custom
module Klu :
  sig
    type ordering = Amd | ColAmd | Natural
    type info = {
      mutable ordering : ordering option;
      mutable reinit : int -> int option -> unit;
      mutable set_ordering : ordering -> unit;
    }
    val info : unit -> info
  end
module Superlumt :
  sig
    type ordering = Natural | MinDegreeProd | MinDegreeSum | ColAmd
    type info = {
      mutable ordering : ordering option;
      mutable set_ordering : ordering -> unit;
      num_threads : int;
    }
    val info : int -> info
  end
module Iterative :
  sig
    type gramschmidt_type = ModifiedGS | ClassicalGS
    type preconditioning_type = PrecNone | PrecLeft | PrecRight | PrecBoth
    type info = {
      mutable maxl : int;
      mutable gs_type : gramschmidt_type option;
      mutable max_restarts : int option;
      mutable set_maxl : int -> unit;
      mutable set_gs_type : gramschmidt_type -> unit;
      mutable set_max_restarts : int -> unit;
      mutable set_prec_type : preconditioning_type -> unit;
    }
    val info : info
  end
module Custom :
  sig
    type ('data, 'kind) atimes_with_data
    external call_atimes :
      ('data, 'kind) atimes_with_data ->
      ('data, 'kind) Nvector.t -> ('data, 'kind) Nvector.t -> unit
      = "sunml_lsolver_call_atimes"
    type ('data, 'kind) precond_with_data
    external call_psetup : ('data, 'kind) precond_with_data -> unit
      = "sunml_lsolver_call_psetup"
    external call_psolve :
      ('data, 'kind) precond_with_data ->
      ('data, 'kind) Nvector.t ->
      ('data, 'kind) Nvector.t -> float -> bool -> unit
      = "sunml_lsolver_call_psolve"
    type ('matrix, 'data, 'kind, 't) ops = {
      init : 't -> unit;
      setup : 't -> 'matrix -> unit;
      solve : 't -> 'matrix -> 'data -> 'data -> float -> unit;
      set_atimes : 't -> ('data, 'kind) atimes_with_data -> unit;
      set_preconditioner :
        't -> ('data, 'kind) precond_with_data -> bool -> bool -> unit;
      set_scaling_vectors : 't -> 'data option -> 'data option -> unit;
      set_zero_guess : 't -> bool -> unit;
      get_id : 't -> linear_solver_id;
      get_num_iters : 't -> int;
      get_res_norm : 't -> float;
      get_res_id : 't -> ('data, 'kind) Nvector.t;
      get_last_flag : 't -> int;
      get_work_space : 't -> int * int;
      set_prec_type : 't -> Iterative.preconditioning_type -> unit;
    }
    type has_ops = {
      has_init : bool;
      has_setup : bool;
      has_set_atimes : bool;
      has_set_preconditioner : bool;
      has_set_scaling_vectors : bool;
      has_set_zero_guess : bool;
      has_get_num_iters : bool;
      has_get_res_norm : bool;
      has_get_res_id : bool;
      has_get_last_flag : bool;
      has_get_work_space : bool;
    }
  end
exception LinearSolverInUse
type ('m, 'nd, 'nk) cptr
type (_, 'nd, 'nk, _) solver_data =
    Spbcgs : ('m, 'nd, 'nk, [> `Spbcgs ]) solver_data
  | Spfgmr : ('m, 'nd, 'nk, [> `Spfgmr ]) solver_data
  | Spgmr : ('m, 'nd, 'nk, [> `Spgmr ]) solver_data
  | Sptfqmr : ('m, 'nd, 'nk, [> `Sptfqmr ]) solver_data
  | Pcg : ('m, 'nd, 'nk, [> `Pcg ]) solver_data
  | Dense : (Sundials.Matrix.Dense.t, 'nd, 'nk, [> `Dls ]) solver_data
  | LapackDense : (Sundials.Matrix.Dense.t, 'nd, 'nk, [> `Dls ]) solver_data
  | Band : (Sundials.Matrix.Band.t, 'nd, 'nk, [> `Dls ]) solver_data
  | LapackBand : (Sundials.Matrix.Band.t, 'nd, 'nk, [> `Dls ]) solver_data
  | Klu :
      Klu.info -> ('s Sundials.Matrix.Sparse.t, 'nd, 'nk, [> `Klu ])
                  solver_data
  | Superlumt :
      Superlumt.info -> ('s Sundials.Matrix.Sparse.t, 'nd, 'nk, [> `Slu ])
                        solver_data
  | Custom :
      ('t * ('m, 'nd, 'nk, 't) Custom.ops) -> ('m, 'nd, 'nk,
                                               [> `Custom of 't ])
                                              solver_data
type 'd atimesfn = 'd -> 'd -> unit
type psetupfn = unit -> unit
type 'd psolvefn = 'd -> 'd -> float -> bool -> unit
type ('d, 'k) ocaml_callbacks = {
  mutable ocaml_atimes : 'd atimesfn;
  mutable ocaml_psetup : psetupfn;
  mutable ocaml_psolve : 'd psolvefn;
  mutable scaling_vector1 : ('d, 'k) Nvector.t option;
  mutable scaling_vector2 : ('d, 'k) Nvector.t option;
}
type ('m, 'mk, 'nd, 'nk, 't) linear_solver_data = {
  rawptr : ('m, 'nd, 'nk) cptr;
  solver : ('m, 'nd, 'nk, 't) solver_data;
  matrix : ('mk, 'm, 'nd, 'nk) Sundials.Matrix.t option;
  compat : Iterative.info;
  mutable check_prec_type : Iterative.preconditioning_type -> bool;
  ocaml_callbacks : ('nd, 'nk) ocaml_callbacks Sundials_impl.Vptr.vptr;
  mutable info_file : Sundials.Logfile.t option;
  mutable attached : bool;
}
val empty_ocaml_callbacks :
  unit -> ('a, 'b) ocaml_callbacks Sundials_impl.Vptr.vptr
type ('m, 'nd, 'nk, 't) linear_solver =
    LS :
      ('m, 'mk, 'nd, 'nk, 't) linear_solver_data -> ('m, 'nd, 'nk, 't)
                                                    linear_solver
type held_linear_solver =
    HLS : ('m, 'mk, 'nd, 'nk, 't) linear_solver_data -> held_linear_solver
  | NoHLS : held_linear_solver
type ('nd, 'nk) solver =
    S : ('m, 'nd, 'nk, 't) solver_data -> ('nd, 'nk) solver
val attach : ('a, 'b, 'c, 'd) linear_solver -> unit
external c_set_prec_type :
  ('m, 'nd, 'nk) cptr ->
  ('m, 'nd, 'nk, 't) solver_data ->
  Iterative.preconditioning_type -> bool -> unit
  = "sunml_lsolver_set_prec_type"
val impl_set_prec_type :
  ('m, 'nd, 'nk) cptr ->
  ('m, 'nd, 'nk, 't) solver_data ->
  Iterative.preconditioning_type -> bool -> unit
external c_make_custom :
  linear_solver_type ->
  ('t * ('m, 'nd, 'nk, 't) Custom.ops) Weak.t ->
  Custom.has_ops -> ('m, 'nd, 'nk) cptr = "sunml_lsolver_make_custom"
