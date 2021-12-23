val e : exn
type 'a double = 'a * 'a
type 'a triple = 'a * 'a * 'a
type ('data, 'kind) nvector = ('data, 'kind) Nvector.t
module LSI = Sundials_LinearSolver_impl
module NLSI = Sundials_NonlinearSolver_impl
type ('t, 'a) jacobian_arg = {
  jac_t : float;
  jac_y : 'a;
  jac_y' : 'a;
  jac_res : 'a;
  jac_coef : float;
  jac_tmp : 't;
}
module DirectTypes :
  sig
    type 'm jac_fn =
        (Sundials.RealArray.t triple, Sundials.RealArray.t) jacobian_arg ->
        'm -> unit
    type 'm jac_callback = { jacfn : 'm jac_fn; mutable jmat : 'm option; }
    val no_callback : 'a -> 'b -> 'c
  end
module SpilsTypes' :
  sig
    type 'a prec_solve_fn =
        (unit, 'a) jacobian_arg -> 'a -> 'a -> float -> unit
    type 'a prec_setup_fn = (unit, 'a) jacobian_arg -> unit
    type 'd jac_times_setup_fn = (unit, 'd) jacobian_arg -> unit
    type 'a jac_times_vec_fn =
        ('a double, 'a) jacobian_arg -> 'a -> 'a -> unit
    type 'a precfns = {
      prec_solve_fn : 'a prec_solve_fn;
      prec_setup_fn : 'a prec_setup_fn option;
    }
  end
module IdaBbdParamTypes :
  sig
    type 'a local_fn = float -> 'a -> 'a -> 'a -> unit
    type 'a comm_fn = float -> 'a -> 'a -> unit
    type 'a precfns = {
      local_fn : 'a local_fn;
      comm_fn : 'a comm_fn option;
    }
  end
module IdaBbdTypes :
  sig
    type bandwidths = { mudq : int; mldq : int; mukeep : int; mlkeep : int; }
  end
module QuadratureTypes :
  sig type 'a quadrhsfn = float -> 'a -> 'a -> 'a -> unit end
module SensitivityTypes :
  sig
    type 'd sensresfn_args = {
      t : float;
      y : 'd;
      y' : 'd;
      res : 'd;
      s : 'd array;
      s' : 'd array;
      tmp : 'd triple;
    }
    type 'd sensresfn = 'd sensresfn_args -> 'd array -> unit
    module QuadratureTypes :
      sig
        type 'd quadsensrhsfn = 'd quadsensrhsfn_args -> 'd array -> unit
        and 'd quadsensrhsfn_args = {
          t : float;
          y : 'd;
          y' : 'd;
          s : 'd array;
          s' : 'd array;
          q : 'd;
          tmp : 'd triple;
        }
      end
  end
module AdjointTypes' :
  sig
    type 'd bresfn_args = { t : float; y : 'd; y' : 'd; yb : 'd; yb' : 'd; }
    type 'a bresfn_no_sens = 'a bresfn_args -> 'a -> unit
    and 'a bresfn_with_sens =
        'a bresfn_args -> 'a array -> 'a array -> 'a -> unit
    type 'a bresfn =
        NoSens of 'a bresfn_no_sens
      | WithSens of 'a bresfn_with_sens
    module QuadratureTypes :
      sig
        type 'd bquadrhsfn_args = {
          t : float;
          y : 'd;
          y' : 'd;
          yb : 'd;
          yb' : 'd;
        }
        type 'a bquadrhsfn =
            NoSens of 'a bquadrhsfn_no_sens
          | WithSens of 'a bquadrhsfn_with_sens
        and 'a bquadrhsfn_no_sens = 'a bquadrhsfn_args -> 'a -> unit
        and 'a bquadrhsfn_with_sens =
            'a bquadrhsfn_args -> 'a array -> 'a array -> 'a -> unit
      end
    type ('t, 'a) jacobian_arg = {
      jac_t : float;
      jac_y : 'a;
      jac_y' : 'a;
      jac_yb : 'a;
      jac_yb' : 'a;
      jac_resb : 'a;
      jac_coef : float;
      jac_tmp : 't;
    }
    module DirectTypes :
      sig
        type 'm jac_fn_no_sens =
            (Sundials.RealArray.t triple, Sundials.RealArray.t) jacobian_arg ->
            'm -> unit
        type 'm jac_fn_with_sens =
            (Sundials.RealArray.t triple, Sundials.RealArray.t) jacobian_arg ->
            Sundials.RealArray.t array ->
            Sundials.RealArray.t array -> 'm -> unit
        type 'm jac_fn =
            NoSens of 'm jac_fn_no_sens
          | WithSens of 'm jac_fn_with_sens
        type 'm jac_callback_no_sens = {
          jacfn : 'm jac_fn_no_sens;
          mutable jmat : 'm option;
        }
        val no_callback : 'a -> 'b -> 'c
        type 'm jac_callback_with_sens = {
          jacfn_sens : 'm jac_fn_with_sens;
          mutable jmat : 'm option;
        }
      end
    module SpilsTypes' :
      sig
        type 'a prec_solve_fn =
            (unit, 'a) jacobian_arg -> 'a -> 'a -> float -> unit
        type 'a prec_setup_fn = (unit, 'a) jacobian_arg -> unit
        type 'a precfns_no_sens = {
          prec_solve_fn : 'a prec_solve_fn;
          prec_setup_fn : 'a prec_setup_fn option;
        }
        type 'd jac_times_setup_fn_no_sens = (unit, 'd) jacobian_arg -> unit
        type 'a jac_times_vec_fn_no_sens =
            ('a, 'a) jacobian_arg -> 'a -> 'a -> unit
        type 'd prec_solve_fn_with_sens =
            (unit, 'd) jacobian_arg ->
            'd array -> 'd array -> 'd -> 'd -> float -> unit
        type 'd prec_setup_fn_with_sens =
            (unit, 'd) jacobian_arg -> 'd array -> 'd array -> unit
        type 'a precfns_with_sens = {
          prec_solve_fn_sens : 'a prec_solve_fn_with_sens;
          prec_setup_fn_sens : 'a prec_setup_fn_with_sens option;
        }
        type 'd jac_times_setup_fn_with_sens =
            (unit, 'd) jacobian_arg -> 'd array -> unit
        type 'd jac_times_vec_fn_with_sens =
            ('d, 'd) jacobian_arg -> 'd array -> 'd array -> 'd -> 'd -> unit
      end
  end
module IdasBbdParamTypes :
  sig
    type 'a local_fn = 'a AdjointTypes'.bresfn_args -> 'a -> unit
    type 'a comm_fn = 'a AdjointTypes'.bresfn_args -> unit
    type 'a precfns = {
      local_fn : 'a local_fn;
      comm_fn : 'a comm_fn option;
    }
  end
module IdasBbdTypes = IdaBbdTypes
type ida_mem
type c_weak_ref
type 'a resfn = float -> 'a -> 'a -> 'a -> unit
type 'a rootsfn = float -> 'a -> 'a -> Sundials.RealArray.t -> unit
type error_handler = Sundials.Util.error_details -> unit
type 'a error_weight_fun = 'a -> 'a -> unit
type ('a, 'kind) session = {
  ida : ida_mem;
  backref : c_weak_ref;
  nroots : int;
  checkvec : ('a, 'kind) Nvector.t -> unit;
  context : Sundials.Context.t;
  mutable exn_temp : exn option;
  mutable id_set : bool;
  mutable resfn : 'a resfn;
  mutable rootsfn : 'a rootsfn;
  mutable errh : error_handler;
  mutable errw : 'a error_weight_fun;
  mutable error_file : Sundials.Logfile.t option;
  mutable ls_solver : LSI.held_linear_solver;
  mutable ls_callbacks : ('a, 'kind) linsolv_callbacks;
  mutable ls_precfns : 'a linsolv_precfns;
  mutable nls_solver :
    ('a, 'kind, ('a, 'kind) session, [ `Nvec ]) NLSI.nonlinear_solver option;
  mutable nls_resfn : 'a resfn;
  mutable sensext : ('a, 'kind) sensext;
}
and ('a, 'kind) sensext =
    NoSensExt
  | FwdSensExt of ('a, 'kind) fsensext
  | BwdSensExt of ('a, 'kind) bsensext
and ('a, 'kind) fsensext = {
  mutable quadrhsfn : 'a QuadratureTypes.quadrhsfn;
  mutable checkquadvec : ('a, 'kind) Nvector.t -> unit;
  mutable has_quad : bool;
  mutable num_sensitivities : int;
  mutable sensarray1 : 'a array;
  mutable sensarray2 : 'a array;
  mutable sensarray3 : 'a array;
  mutable senspvals : Sundials.RealArray.t option;
  mutable sensresfn : 'a SensitivityTypes.sensresfn;
  mutable quadsensrhsfn : 'a SensitivityTypes.QuadratureTypes.quadsensrhsfn;
  mutable fnls_solver :
    ('a, 'kind, ('a, 'kind) session, [ `Sens ]) NLSI.nonlinear_solver option;
  mutable bsessions : ('a, 'kind) session list;
}
and ('a, 'kind) bsensext = {
  parent : ('a, 'kind) session;
  which : int;
  bnum_sensitivities : int;
  bsensarray1 : 'a array;
  bsensarray2 : 'a array;
  mutable bresfn : 'a AdjointTypes'.bresfn_no_sens;
  mutable bresfn_sens : 'a AdjointTypes'.bresfn_with_sens;
  mutable bquadrhsfn : 'a AdjointTypes'.QuadratureTypes.bquadrhsfn_no_sens;
  mutable bquadrhsfn_sens :
    'a AdjointTypes'.QuadratureTypes.bquadrhsfn_with_sens;
  mutable checkbquadvec : ('a, 'kind) Nvector.t -> unit;
}
and ('a, 'kind) linsolv_callbacks =
    NoCallbacks
  | DlsDenseCallback of Sundials.Matrix.Dense.t DirectTypes.jac_callback
  | DlsBandCallback of Sundials.Matrix.Band.t DirectTypes.jac_callback
  | BDlsDenseCallback of
      Sundials.Matrix.Dense.t AdjointTypes'.DirectTypes.jac_callback_no_sens
  | BDlsDenseCallbackSens of
      Sundials.Matrix.Dense.t
      AdjointTypes'.DirectTypes.jac_callback_with_sens
  | BDlsBandCallback of
      Sundials.Matrix.Band.t AdjointTypes'.DirectTypes.jac_callback_no_sens
  | BDlsBandCallbackSens of
      Sundials.Matrix.Band.t AdjointTypes'.DirectTypes.jac_callback_with_sens
  | SlsKluCallback :
      's Sundials.Matrix.Sparse.t DirectTypes.jac_callback -> ('a, 'kind)
                                                              linsolv_callbacks
  | BSlsKluCallback :
      's Sundials.Matrix.Sparse.t
      AdjointTypes'.DirectTypes.jac_callback_no_sens -> ('a, 'kind)
                                                        linsolv_callbacks
  | BSlsKluCallbackSens :
      's Sundials.Matrix.Sparse.t
      AdjointTypes'.DirectTypes.jac_callback_with_sens -> ('a, 'kind)
                                                          linsolv_callbacks
  | SlsSuperlumtCallback :
      's Sundials.Matrix.Sparse.t DirectTypes.jac_callback -> ('a, 'kind)
                                                              linsolv_callbacks
  | BSlsSuperlumtCallback :
      's Sundials.Matrix.Sparse.t
      AdjointTypes'.DirectTypes.jac_callback_no_sens -> ('a, 'kind)
                                                        linsolv_callbacks
  | BSlsSuperlumtCallbackSens :
      's Sundials.Matrix.Sparse.t
      AdjointTypes'.DirectTypes.jac_callback_with_sens -> ('a, 'kind)
                                                          linsolv_callbacks
  | DirectCustomCallback :
      'm DirectTypes.jac_callback -> ('a, 'kind) linsolv_callbacks
  | BDirectCustomCallback :
      'm AdjointTypes'.DirectTypes.jac_callback_no_sens -> ('a, 'kind)
                                                           linsolv_callbacks
  | BDirectCustomCallbackSens :
      'm AdjointTypes'.DirectTypes.jac_callback_with_sens -> ('a, 'kind)
                                                             linsolv_callbacks
  | SpilsCallback1 of 'a SpilsTypes'.jac_times_vec_fn option *
      'a SpilsTypes'.jac_times_setup_fn option
  | SpilsCallback2 of 'a resfn
  | BSpilsCallback of
      'a AdjointTypes'.SpilsTypes'.jac_times_vec_fn_no_sens option *
      'a AdjointTypes'.SpilsTypes'.jac_times_setup_fn_no_sens option
  | BSpilsCallbackSens of
      'a AdjointTypes'.SpilsTypes'.jac_times_vec_fn_with_sens option *
      'a AdjointTypes'.SpilsTypes'.jac_times_setup_fn_with_sens option
  | BSpilsCallbackJTRhsfn of 'a resfn
and 'a linsolv_precfns =
    NoPrecFns
  | PrecFns of 'a SpilsTypes'.precfns
  | BPrecFns of 'a AdjointTypes'.SpilsTypes'.precfns_no_sens
  | BPrecFnsSens of 'a AdjointTypes'.SpilsTypes'.precfns_with_sens
  | BandedPrecFns
  | BBDPrecFns of 'a IdaBbdParamTypes.precfns
  | BBBDPrecFns of 'a IdasBbdParamTypes.precfns
val revlookup_bsession :
  ('d, 'k) session -> ida_mem -> ('d, 'k) session option
val ls_check_direct : ('a, 'b) session -> unit
val ls_check_spils : ('a, 'b) session -> unit
val ls_check_spils_bbd : ('a, 'b) session -> unit
type 'a serial_session = (Nvector_serial.data, 'a) session
  constraint 'a = [> Nvector_serial.kind ]
type ('data, 'kind) linear_solver =
    ('data, 'kind) session -> ('data, 'kind) Nvector.t -> unit
type 'a serial_linear_solver = (Nvector_serial.data, 'a) linear_solver
  constraint 'a = [> Nvector_serial.kind ]
module SpilsTypes :
  sig
    type 'a prec_solve_fn =
        (unit, 'a) jacobian_arg -> 'a -> 'a -> float -> unit
    type 'a prec_setup_fn = (unit, 'a) jacobian_arg -> unit
    type 'd jac_times_setup_fn = (unit, 'd) jacobian_arg -> unit
    type 'a jac_times_vec_fn =
        ('a double, 'a) jacobian_arg -> 'a -> 'a -> unit
    type 'a precfns =
      'a SpilsTypes'.precfns = {
      prec_solve_fn : 'a prec_solve_fn;
      prec_setup_fn : 'a prec_setup_fn option;
    }
    type ('a, 'k) set_preconditioner =
        ('a, 'k) session -> ('a, 'k) Nvector.t -> unit
    type ('a, 'k) preconditioner =
        LSI.Iterative.preconditioning_type * ('a, 'k) set_preconditioner
    type 'a serial_preconditioner = (Nvector_serial.data, 'a) preconditioner
      constraint 'a = [> Nvector_serial.kind ]
  end
module AdjointTypes :
  sig
    type 'd bresfn_args =
      'd AdjointTypes'.bresfn_args = {
      t : float;
      y : 'd;
      y' : 'd;
      yb : 'd;
      yb' : 'd;
    }
    type 'a bresfn_no_sens = 'a bresfn_args -> 'a -> unit
    and 'a bresfn_with_sens =
        'a bresfn_args -> 'a array -> 'a array -> 'a -> unit
    type 'a bresfn =
      'a AdjointTypes'.bresfn =
        NoSens of 'a bresfn_no_sens
      | WithSens of 'a bresfn_with_sens
    module QuadratureTypes = AdjointTypes'.QuadratureTypes
    type ('t, 'a) jacobian_arg =
      ('t, 'a) AdjointTypes'.jacobian_arg = {
      jac_t : float;
      jac_y : 'a;
      jac_y' : 'a;
      jac_yb : 'a;
      jac_yb' : 'a;
      jac_resb : 'a;
      jac_coef : float;
      jac_tmp : 't;
    }
    module DirectTypes = AdjointTypes'.DirectTypes
    module SpilsTypes' = AdjointTypes'.SpilsTypes'
    type ('a, 'k) bsession = Bsession of ('a, 'k) session
    type 'a serial_bsession = (Nvector_serial.data, 'a) bsession
      constraint 'a = [> Nvector_serial.kind ]
    val tosession : ('a, 'b) bsession -> ('a, 'b) session
    val parent_and_which : ('a, 'b) bsession -> ('a, 'b) session * int
    type ('data, 'kind) linear_solver =
        ('data, 'kind) bsession -> ('data, 'kind) Nvector.t -> unit
    type 'a serial_linear_solver = (Nvector_serial.data, 'a) linear_solver
      constraint 'a = [> Nvector_serial.kind ]
    module SpilsTypes :
      sig
        type 'a prec_solve_fn =
            (unit, 'a) AdjointTypes'.jacobian_arg ->
            'a -> 'a -> float -> unit
        type 'a prec_setup_fn = (unit, 'a) AdjointTypes'.jacobian_arg -> unit
        type 'a precfns_no_sens =
          'a AdjointTypes'.SpilsTypes'.precfns_no_sens = {
          prec_solve_fn : 'a prec_solve_fn;
          prec_setup_fn : 'a prec_setup_fn option;
        }
        type 'd jac_times_setup_fn_no_sens =
            (unit, 'd) AdjointTypes'.jacobian_arg -> unit
        type 'a jac_times_vec_fn_no_sens =
            ('a, 'a) AdjointTypes'.jacobian_arg -> 'a -> 'a -> unit
        type 'd prec_solve_fn_with_sens =
            (unit, 'd) AdjointTypes'.jacobian_arg ->
            'd array -> 'd array -> 'd -> 'd -> float -> unit
        type 'd prec_setup_fn_with_sens =
            (unit, 'd) AdjointTypes'.jacobian_arg ->
            'd array -> 'd array -> unit
        type 'a precfns_with_sens =
          'a AdjointTypes'.SpilsTypes'.precfns_with_sens = {
          prec_solve_fn_sens : 'a prec_solve_fn_with_sens;
          prec_setup_fn_sens : 'a prec_setup_fn_with_sens option;
        }
        type 'd jac_times_setup_fn_with_sens =
            (unit, 'd) AdjointTypes'.jacobian_arg -> 'd array -> unit
        type 'd jac_times_vec_fn_with_sens =
            ('d, 'd) AdjointTypes'.jacobian_arg ->
            'd array -> 'd array -> 'd -> 'd -> unit
        type ('a, 'k) set_preconditioner =
            ('a, 'k) bsession ->
            ('a, 'k) session -> int -> ('a, 'k) Nvector.t -> unit
        type ('a, 'k) preconditioner =
            LSI.Iterative.preconditioning_type * ('a, 'k) set_preconditioner
        type 'a serial_preconditioner =
            (Nvector_serial.data, 'a) preconditioner
          constraint 'a = [> Nvector_serial.kind ]
      end
  end
val read_weak_ref : ('a, 'kind) session Weak.t -> ('a, 'kind) session
val dummy_resfn : 'a -> 'b -> 'c -> 'd -> 'e
val dummy_nlsresfn : 'a -> 'b -> 'c -> 'd -> 'e
val dummy_rootsfn : 'a -> 'b -> 'c -> 'd -> 'e
val dummy_errh : 'a -> 'b
val dummy_errw : 'a -> 'b -> 'c
val dummy_bresfn_no_sens : 'a -> 'b -> 'c
val dummy_bresfn_with_sens : 'a -> 'b -> 'c -> 'd -> 'e
val dummy_bquadrhsfn_no_sens : 'a -> 'b -> 'c
val dummy_bquadrhsfn_with_sens : 'a -> 'b -> 'c -> 'd -> 'e
val dummy_quadrhsfn : 'a -> 'b -> 'c -> 'd -> 'e
val dummy_sensresfn : 'a -> 'b -> 'c
val dummy_quadsensrhsfn : 'a -> 'b -> 'c
