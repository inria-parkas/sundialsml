val e : exn
module LSI = Sundials_LinearSolver_impl
type ('a, 'k) nvector = ('a, 'k) Nvector.t
type 'a double = 'a * 'a
type ('t, 'a) jacobian_arg = { jac_u : 'a; jac_fu : 'a; jac_tmp : 't; }
type real_array = Sundials.RealArray.t
type bandrange = { mupper : int; mlower : int; }
module DirectTypes :
  sig
    type 'm jac_fn =
        (real_array double, real_array) jacobian_arg -> 'm -> unit
    type 'm jac_callback = { jacfn : 'm jac_fn; mutable jmat : 'm option; }
    val no_callback : 'a -> 'b -> 'c
  end
module SpilsTypes' :
  sig
    type 'a solve_arg = { uscale : 'a; fscale : 'a; }
    type 'a prec_solve_fn =
        (unit, 'a) jacobian_arg -> 'a solve_arg -> 'a -> unit
    type 'a prec_setup_fn = (unit, 'a) jacobian_arg -> 'a solve_arg -> unit
    type 'a jac_times_vec_fn = 'a -> 'a -> 'a -> bool -> bool
    type 'a precfns = {
      prec_solve_fn : 'a prec_solve_fn;
      prec_setup_fn : 'a prec_setup_fn option;
    }
  end
module KinsolBbdParamTypes :
  sig
    type 'a local_fn = 'a -> 'a -> unit
    type 'a comm_fn = 'a -> unit
    type 'a precfns = {
      local_fn : 'a local_fn;
      comm_fn : 'a comm_fn option;
    }
  end
module KinsolBbdTypes :
  sig
    type bandwidths = { mudq : int; mldq : int; mukeep : int; mlkeep : int; }
  end
type kin_mem
type c_weak_ref
type 'a sysfn = 'a -> 'a -> unit
type errh = Sundials.Util.error_details -> unit
type infoh = Sundials.Util.error_details -> unit
type ('a, 'k) session = {
  kinsol : kin_mem;
  backref : c_weak_ref;
  initvec : ('a, 'k) Nvector.t;
  checkvec : ('a, 'k) Nvector.t -> unit;
  context : Sundials.Context.t;
  mutable neqs : int;
  mutable exn_temp : exn option;
  mutable sysfn : 'a sysfn;
  mutable errh : errh;
  mutable infoh : infoh;
  mutable error_file : Sundials.Logfile.t option;
  mutable info_file : Sundials.Logfile.t option;
  mutable ls_solver : LSI.held_linear_solver;
  mutable ls_callbacks : ('a, 'k) linsolv_callbacks;
  mutable ls_precfns : 'a linsolv_precfns;
}
and ('a, 'kind) linsolv_callbacks =
    NoCallbacks
  | DlsDenseCallback of Sundials.Matrix.Dense.t DirectTypes.jac_callback
  | DlsBandCallback of Sundials.Matrix.Band.t DirectTypes.jac_callback
  | SlsKluCallback :
      's Sundials.Matrix.Sparse.t DirectTypes.jac_callback -> ('a, 'kind)
                                                              linsolv_callbacks
  | SlsSuperlumtCallback :
      's Sundials.Matrix.Sparse.t DirectTypes.jac_callback -> ('a, 'kind)
                                                              linsolv_callbacks
  | DirectCustomCallback :
      'm DirectTypes.jac_callback -> ('a, 'kind) linsolv_callbacks
  | SpilsCallback1 of 'a SpilsTypes'.jac_times_vec_fn option
  | SpilsCallback2 of 'a sysfn
and 'a linsolv_precfns =
    NoPrecFns
  | PrecFns of 'a SpilsTypes'.precfns
  | BBDPrecFns of 'a KinsolBbdParamTypes.precfns
val ls_check_direct : ('a, 'b) session -> unit
val ls_check_spils : ('a, 'b) session -> unit
val ls_check_spils_bbd : ('a, 'b) session -> unit
type 'a serial_session = (Nvector_serial.data, 'a) session
  constraint 'a = [> Nvector_serial.kind ]
type ('data, 'kind) linear_solver =
    ('data, 'kind) session -> ('data, 'kind) nvector -> unit
module SpilsTypes :
  sig
    type 'a solve_arg =
      'a SpilsTypes'.solve_arg = {
      uscale : 'a;
      fscale : 'a;
    }
    type 'a prec_solve_fn =
        (unit, 'a) jacobian_arg -> 'a solve_arg -> 'a -> unit
    type 'a prec_setup_fn = (unit, 'a) jacobian_arg -> 'a solve_arg -> unit
    type 'a jac_times_vec_fn = 'a -> 'a -> 'a -> bool -> bool
    type 'a precfns =
      'a SpilsTypes'.precfns = {
      prec_solve_fn : 'a prec_solve_fn;
      prec_setup_fn : 'a prec_setup_fn option;
    }
    type ('a, 'k) set_preconditioner =
        ('a, 'k) session -> ('a, 'k) nvector -> unit
    type ('a, 'k) preconditioner =
        LSI.Iterative.preconditioning_type * ('a, 'k) set_preconditioner
  end
val read_weak_ref : ('a, 'k) session Weak.t -> ('a, 'k) session
val dummy_sysfn : 'a -> 'b -> 'c
val dummy_errh : 'a -> 'b
val dummy_infoh : 'a -> 'b
