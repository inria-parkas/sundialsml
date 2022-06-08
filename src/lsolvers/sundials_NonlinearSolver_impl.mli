val e : exn
exception VectorOpError
exception IncorrectUse
exception ExtFail
exception NonlinearSolverInUse
module Senswrapper :
  sig
    type ('d, 'k) t
    external data : ('d, 'k) t -> 'd array = "sunml_senswrapper_data"
  end
type 'a integrator = Integrator of 'a
type ('nv, 's) sysfn = 'nv -> 'nv -> 's -> unit
type 's lsetupfn = bool -> 's -> bool
type ('nv, 's) lsolvefn = 'nv -> 's -> unit
type convtest = Success | Continue | Recover
type ('nv, 's) convtestfn' = 'nv -> 'nv -> float -> 'nv -> 's -> convtest
type ('d, 'k, 's, 'v) cptr
type ('nv, 's) callbacks = {
  mutable sysfn : ('nv, 's) sysfn;
  mutable lsetupfn : 's lsetupfn;
  mutable lsolvefn : ('nv, 's) lsolvefn;
  mutable convtestfn : ('nv, 's) convtestfn';
}
type nonlinear_solver_type = RootFind | FixedPoint
type ('d, 'k, 's, _) solver =
    NewtonSolver : ('d, 's) callbacks -> ('d, 'k, 's, [ `Nvec ]) solver
  | NewtonSolverSens :
      (('d, 'k) Senswrapper.t, 's) callbacks -> ('d, 'k, 's, [ `Sens ])
                                                solver
  | FixedPointSolver : ('d, 's) callbacks *
      int -> ('d, 'k, 's, [ `Nvec ]) solver
  | FixedPointSolverSens : (('d, 'k) Senswrapper.t, 's) callbacks *
      int -> ('d, 'k, 's, [ `Sens ]) solver
  | CustomSolver : (('d, 'k) Nvector.t, 's) callbacks *
      (('d, 'k) Nvector.t, 'd, 's, [ `Nvec ]) ops -> ('d, 'k, 's, [ `Nvec ])
                                                     solver
  | CustomSolverSens : (('d, 'k) Senswrapper.t, 's) callbacks *
      (('d, 'k) Senswrapper.t, ('d, 'k) Senswrapper.t, 's, [ `Sens ]) ops -> 
      ('d, 'k, 's, [ `Sens ]) solver
and ('d, 'k, 's, 'v) nonlinear_solver = {
  rawptr : ('d, 'k, 's, 'v) cptr;
  solver : ('d, 'k, 's, 'v) solver;
  context : Sundials.Context.t;
  mutable info_file : Sundials.Logfile.t option;
  mutable attached : bool;
}
and 's convtest_callback = {
  f :
    'd1 'k1 't2 'd2 'k2.
      ('d1, 'k1, 't2, [ `Nvec ]) nonlinear_solver ->
      (('d2, 'k2) Nvector.t, 's) convtestfn';
} [@@unboxed]
and 's convtest_callback_sens = {
  f :
    'd1 'k1 't2 'd2 'k2.
      ('d1, 'k1, 't2, [ `Sens ]) nonlinear_solver ->
      (('d2, 'k2) Senswrapper.t, 's) convtestfn';
} [@@unboxed]
and ('nv, 's, 'v) convtestfn =
    CConvTest :
      's convtest_callback Sundials.cfun -> ('nv, 's, [ `Nvec ]) convtestfn
  | CSensConvTest :
      's convtest_callback_sens Sundials.cfun -> ('nv, 's, [ `Sens ])
                                                 convtestfn
  | OConvTest of ('nv, 's) convtestfn'
and ('nv, 'd, 's, 'v) ops = {
  nls_type : nonlinear_solver_type;
  init : (unit -> unit) option;
  setup : ('d -> 's -> unit) option;
  solve : 'd -> 'd -> 'd -> float -> bool -> 's -> unit;
  set_sys_fn : ('nv, 's) sysfn -> unit;
  set_lsetup_fn : ('s lsetupfn -> unit) option;
  set_lsolve_fn : (('nv, 's) lsolvefn -> unit) option;
  set_convtest_fn : (('d, 's, 'v) convtestfn -> unit) option;
  set_max_iters : (int -> unit) option;
  get_num_iters : (unit -> int) option;
  get_cur_iter : (unit -> int) option;
  get_num_conv_fails : (unit -> int) option;
}
type ('d, 'k, 's) nonlinear_solver_hold =
    NoNLS
  | NLS of ('d, 'k, 's, [ `Nvec ]) nonlinear_solver
  | NLS_sens of ('d, 'k, 's, [ `Sens ]) nonlinear_solver
val attach : ('a, 'b, 'c, 'd) nonlinear_solver -> unit
val detach : ('a, 'b, 'c, 'd) nonlinear_solver -> unit
val get_type : ('d, 'k, 's, 'v) nonlinear_solver -> nonlinear_solver_type
val empty_sysfn : 'a -> 'b -> 'c -> 'd
val empty_lsetupfn : 'a -> 'b -> 'c
val empty_lsolvefn : 'a -> 'b -> 'c
val empty_convtestfn : 'a -> 'b -> 'c -> 'd -> 'e
val empty_callbacks : unit -> ('a, 'b) callbacks
