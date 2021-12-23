val e : exn
module LSI = Sundials_LinearSolver_impl
module NLSI = Sundials_NonlinearSolver_impl
type 'a triple = 'a * 'a * 'a
type ('t, 'a) jacobian_arg = {
  jac_t : float;
  jac_y : 'a;
  jac_fy : 'a;
  jac_tmp : 't;
}
module DirectTypes :
  sig
    type 'm jac_fn =
        (Sundials.RealArray.t triple, Sundials.RealArray.t) jacobian_arg ->
        'm -> unit
    type 'm jac_callback = { jacfn : 'm jac_fn; mutable jmat : 'm option; }
    type 'm linsys_fn =
        (Sundials.RealArray.t triple, Sundials.RealArray.t) jacobian_arg ->
        'm -> 'm option -> bool -> float -> bool
    val no_callback : 'a -> 'b -> 'c
    val no_linsysfn : 'a -> 'b -> 'c -> 'd -> 'e -> 'f
  end
module SpilsCommonTypes :
  sig
    type 'a prec_solve_arg = {
      rhs : 'a;
      gamma : float;
      delta : float;
      left : bool;
    }
  end
module SpilsTypes' :
  sig
    type 'a prec_solve_arg =
      'a SpilsCommonTypes.prec_solve_arg = {
      rhs : 'a;
      gamma : float;
      delta : float;
      left : bool;
    }
    type 'a prec_solve_fn =
        (unit, 'a) jacobian_arg -> 'a prec_solve_arg -> 'a -> unit
    type 'a prec_setup_fn = (unit, 'a) jacobian_arg -> bool -> float -> bool
    type 'd jac_times_setup_fn = (unit, 'd) jacobian_arg -> unit
    type 'a jac_times_vec_fn = ('a, 'a) jacobian_arg -> 'a -> 'a -> unit
    type 'a precfns = {
      prec_solve_fn : 'a prec_solve_fn;
      prec_setup_fn : 'a prec_setup_fn option;
    }
  end
module MassTypes' :
  sig
    module Direct' :
      sig
        type 'm mass_fn = float -> Sundials.RealArray.t triple -> 'm -> unit
        type 'm mass_callback = {
          massfn : 'm mass_fn;
          mutable mmat : 'm option;
        }
        val no_mass_callback : 'a -> 'b -> 'c
      end
    module Iterative' :
      sig
        type 'd prec_solve_arg = { rhs : 'd; delta : float; left : bool; }
        type 'd prec_solve_fn = float -> 'd prec_solve_arg -> 'd -> unit
        type 'd prec_setup_fn = float -> unit
        type mass_times_setup_fn = float -> unit
        type 'd mass_times_vec_fn = float -> 'd -> 'd -> unit
        type 'a precfns = {
          prec_solve_fn : 'a prec_solve_fn;
          prec_setup_fn : (float -> unit) option;
        }
      end
  end
module ArkodeBbdParamTypes :
  sig
    type 'a local_fn = float -> 'a -> 'a -> unit
    type 'a comm_fn = float -> 'a -> unit
    type 'a precfns = {
      local_fn : 'a local_fn;
      comm_fn : 'a comm_fn option;
    }
  end
module ArkodeBbdTypes :
  sig
    type bandwidths = { mudq : int; mldq : int; mukeep : int; mlkeep : int; }
  end
type 'step arkode_mem
type c_weak_ref
type arkstep = [ `ARKStep ]
type erkstep = [ `ERKStep ]
type mristep = [ `MRIStep ]
module Global :
  sig
    type 'a rhsfn = float -> 'a -> 'a -> unit
    type 'a rootsfn = float -> 'a -> Sundials.RealArray.t -> unit
    type error_handler = Sundials.Util.error_details -> unit
    type 'a error_weight_fun = 'a -> 'a -> unit
    type adaptivity_args = {
      h1 : float;
      h2 : float;
      h3 : float;
      e1 : float;
      e2 : float;
      e3 : float;
      q : int;
      p : int;
    }
    type 'd adaptivity_fn = float -> 'd -> adaptivity_args -> float
    type 'd stability_fn = float -> 'd -> float
    type 'd resize_fn = 'd -> 'd -> unit
    type 'd postprocess_step_fn = float -> 'd -> unit
    type 'd stage_predict_fn = float -> 'd -> unit
    type 'd pre_inner_fn = float -> 'd array -> unit
    type 'd post_inner_fn = float -> 'd -> unit
  end
type ('d, 'k) inner_stepper_cptr
type fullrhs_mode = Start | End | Other
type 'd inner_stepper_callbacks = {
  evolve_fn : float -> float -> 'd -> unit;
  full_rhs_fn : float -> 'd -> 'd -> fullrhs_mode -> unit;
  reset_fn : float -> 'd -> unit;
}
type 'a res_weight_fun = 'a -> 'a -> unit
type ('a, 'kind, 'step) session = {
  arkode : 'step arkode_mem;
  backref : c_weak_ref;
  nroots : int;
  mutable checkvec : ('a, 'kind) Nvector.t -> unit;
  mutable uses_resv : bool;
  context : Sundials.Context.t;
  mutable exn_temp : exn option;
  mutable problem : problem_type;
  mutable rhsfn1 : 'a Global.rhsfn;
  mutable rhsfn2 : 'a Global.rhsfn;
  mutable rootsfn : 'a Global.rootsfn;
  mutable errh : Global.error_handler;
  mutable errw : 'a Global.error_weight_fun;
  mutable resw : 'a res_weight_fun;
  mutable error_file : Sundials.Logfile.t option;
  mutable diag_file : Sundials.Logfile.t option;
  mutable adaptfn : 'a Global.adaptivity_fn;
  mutable stabfn : 'a Global.stability_fn;
  mutable resizefn : 'a Global.resize_fn;
  mutable poststepfn : 'a Global.postprocess_step_fn;
  mutable stagepredictfn : 'a Global.stage_predict_fn;
  mutable preinnerfn : 'a Global.pre_inner_fn;
  mutable postinnerfn : 'a Global.post_inner_fn;
  mutable preinnerarray : 'a array;
  mutable linsolver : ('a, 'kind, 'step) lin_solver option;
  mutable ls_solver : LSI.held_linear_solver;
  mutable ls_callbacks : ('a, 'kind) linsolv_callbacks;
  mutable ls_precfns : 'a linsolv_precfns;
  mutable mass_solver : LSI.held_linear_solver;
  mutable mass_callbacks : ('a, 'kind) mass_callbacks;
  mutable mass_precfns : 'a mass_precfns;
  mutable nls_solver :
    ('a, 'kind, ('a, 'kind, 'step) session, [ `Nvec ]) NLSI.nonlinear_solver
    option;
  mutable nls_rhsfn : 'a Global.rhsfn;
  mutable inner_session : ('a, 'kind) inner_stepper option;
}
and problem_type = ImplicitOnly | ExplicitOnly | ImplicitAndExplicit
and ('data, 'kind, 'step) lin_solver =
    ('data, 'kind, 'step) session -> ('data, 'kind) Nvector.t -> unit
and ('a, 'kind) linsolv_callbacks =
    NoCallbacks
  | DlsDenseCallback of Sundials.Matrix.Dense.t DirectTypes.jac_callback *
      Sundials.Matrix.Dense.t DirectTypes.linsys_fn
  | DlsBandCallback of Sundials.Matrix.Band.t DirectTypes.jac_callback *
      Sundials.Matrix.Band.t DirectTypes.linsys_fn
  | SlsKluCallback : 's Sundials.Matrix.Sparse.t DirectTypes.jac_callback *
      's Sundials.Matrix.Sparse.t DirectTypes.linsys_fn -> ('a, 'kind)
                                                           linsolv_callbacks
  | SlsSuperlumtCallback :
      's Sundials.Matrix.Sparse.t DirectTypes.jac_callback *
      's Sundials.Matrix.Sparse.t DirectTypes.linsys_fn -> ('a, 'kind)
                                                           linsolv_callbacks
  | DirectCustomCallback : 'm DirectTypes.jac_callback *
      'm DirectTypes.linsys_fn -> ('a, 'kind) linsolv_callbacks
  | SpilsCallback1 of 'a SpilsTypes'.jac_times_vec_fn option *
      'a SpilsTypes'.jac_times_setup_fn option
  | SpilsCallback2 of 'a Global.rhsfn
and 'a linsolv_precfns =
    NoPrecFns
  | PrecFns of 'a SpilsTypes'.precfns
  | BandedPrecFns
  | BBDPrecFns of 'a ArkodeBbdParamTypes.precfns
and ('a, 'kind) mass_callbacks =
    NoMassCallbacks
  | DlsDenseMassCallback of
      Sundials.Matrix.Dense.t MassTypes'.Direct'.mass_callback *
      Sundials.Matrix.Dense.t
  | DlsBandMassCallback of
      Sundials.Matrix.Band.t MassTypes'.Direct'.mass_callback *
      Sundials.Matrix.Band.t
  | SlsKluMassCallback :
      's Sundials.Matrix.Sparse.t MassTypes'.Direct'.mass_callback *
      's Sundials.Matrix.Sparse.t -> ('a, 'kind) mass_callbacks
  | SlsSuperlumtMassCallback :
      's Sundials.Matrix.Sparse.t MassTypes'.Direct'.mass_callback *
      's Sundials.Matrix.Sparse.t -> ('a, 'kind) mass_callbacks
  | DirectCustomMassCallback : 'm MassTypes'.Direct'.mass_callback *
      'm -> ('a, 'kind) mass_callbacks
  | SpilsMassCallback of 'a MassTypes'.Iterative'.mass_times_vec_fn *
      MassTypes'.Iterative'.mass_times_setup_fn option
and 'a mass_precfns =
    NoMassPrecFns
  | MassPrecFns of 'a MassTypes'.Iterative'.precfns
and ('d, 'k) istepper =
    ARKStepInnerStepper of ('d, 'k, arkstep) session
  | SundialsInnerStepper
  | CustomInnerStepper of 'd inner_stepper_callbacks
and ('d, 'k) inner_stepper = {
  rawptr : ('d, 'k) inner_stepper_cptr;
  istepper : ('d, 'k) istepper;
  mutable icheckvec : (('d, 'k) Nvector.t -> unit) option;
}
val ls_check_direct : ('a, 'b, 'c) session -> unit
val ls_check_spils : ('a, 'b, 'c) session -> unit
val ls_check_spils_band : ('a, 'b, 'c) session -> unit
val ls_check_spils_bbd : ('a, 'b, 'c) session -> unit
val mass_check_direct : ('a, 'b, 'c) session -> unit
val mass_check_spils : ('a, 'b, 'c) session -> unit
type ('a, 'step) serial_session = (Nvector_serial.data, 'a, 'step) session
  constraint 'a = [> Nvector_serial.kind ]
module SpilsTypes :
  sig
    type 'a prec_solve_arg =
      'a SpilsCommonTypes.prec_solve_arg = {
      rhs : 'a;
      gamma : float;
      delta : float;
      left : bool;
    }
    type 'a prec_solve_fn =
        (unit, 'a) jacobian_arg -> 'a prec_solve_arg -> 'a -> unit
    type 'a prec_setup_fn = (unit, 'a) jacobian_arg -> bool -> float -> bool
    type 'd jac_times_setup_fn = (unit, 'd) jacobian_arg -> unit
    type 'a jac_times_vec_fn = ('a, 'a) jacobian_arg -> 'a -> 'a -> unit
    type 'a precfns =
      'a SpilsTypes'.precfns = {
      prec_solve_fn : 'a prec_solve_fn;
      prec_setup_fn : 'a prec_setup_fn option;
    }
    type ('a, 'k, 's) set_preconditioner =
        ('a, 'k, 's) session -> ('a, 'k) Nvector.t -> unit
    type ('a, 'k, 's) preconditioner =
        LSI.Iterative.preconditioning_type * ('a, 'k, 's) set_preconditioner
    type ('a, 's) serial_preconditioner =
        (Nvector_serial.data, 'a, 's) preconditioner
      constraint 'a = [> Nvector_serial.kind ]
  end
module MassTypes :
  sig
    type ('data, 'kind) solver =
        ('data, 'kind, arkstep) session -> ('data, 'kind) Nvector.t -> unit
    type 'a serial_solver = (Nvector_serial.data, 'a) solver
      constraint 'a = [> Nvector_serial.kind ]
    module Direct' :
      sig
        type 'm mass_fn = float -> Sundials.RealArray.t triple -> 'm -> unit
        type 'm mass_callback =
          'm MassTypes'.Direct'.mass_callback = {
          massfn : 'm mass_fn;
          mutable mmat : 'm option;
        }
        val no_mass_callback : 'a -> 'b -> 'c
      end
    module Iterative' :
      sig
        type 'd prec_solve_arg =
          'd MassTypes'.Iterative'.prec_solve_arg = {
          rhs : 'd;
          delta : float;
          left : bool;
        }
        type 'd prec_solve_fn = float -> 'd prec_solve_arg -> 'd -> unit
        type 'd prec_setup_fn = float -> unit
        type mass_times_setup_fn = float -> unit
        type 'd mass_times_vec_fn = float -> 'd -> 'd -> unit
        type 'a precfns =
          'a MassTypes'.Iterative'.precfns = {
          prec_solve_fn : 'a prec_solve_fn;
          prec_setup_fn : (float -> unit) option;
        }
        type ('a, 'k) set_preconditioner =
            ('a, 'k, arkstep) session -> ('a, 'k) Nvector.t -> unit
        type ('a, 'k) preconditioner =
            LSI.Iterative.preconditioning_type * ('a, 'k) set_preconditioner
        type 'a serial_preconditioner =
            (Nvector_serial.data, 'a) preconditioner
          constraint 'a = [> Nvector_serial.kind ]
      end
  end
val read_weak_ref :
  ('a, 'kind, 'step) session Weak.t -> ('a, 'kind, 'step) session
val empty_preinnerarray : unit -> 'a array
val dummy_rhsfn1 : 'a -> 'b -> 'c -> 'd
val dummy_rhsfn2 : 'a -> 'b -> 'c -> 'd
val dummy_nlsrhsfn : 'a -> 'b -> 'c -> 'd
val dummy_rootsfn : 'a -> 'b -> 'c -> 'd
val dummy_errh : 'a -> 'b
val dummy_errw : 'a -> 'b -> 'c
val dummy_resw : 'a -> 'b -> 'c
val dummy_adaptfn : 'a -> 'b -> 'c -> 'd
val dummy_stabfn : 'a -> 'b -> 'c
val dummy_resizefn : 'a -> 'b -> 'c
val dummy_poststepfn : 'a -> 'b -> 'c
val dummy_stagepredictfn : 'a -> 'b -> 'c
val dummy_preinnerfn : 'a -> 'b -> 'c
val dummy_postinnerfn : 'a -> 'b -> 'c
