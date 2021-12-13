val e : exn
module LSI = Sundials_LinearSolver_impl
module NLSI = Sundials_NonlinearSolver_impl
type ('data, 'kind) nvector = ('data, 'kind) Nvector.t
type 'a double = 'a * 'a
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
        'm -> bool -> float -> bool
    val no_callback : 'a -> 'b -> 'c
    val no_linsysfn : 'a -> 'b -> 'c -> 'd -> 'e
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
module CvodeBbdParamTypes :
  sig
    type 'a local_fn = float -> 'a -> 'a -> unit
    type 'a comm_fn = float -> 'a -> unit
    type 'a precfns = {
      local_fn : 'a local_fn;
      comm_fn : 'a comm_fn option;
    }
  end
module CvodeBbdTypes :
  sig
    type bandwidths = { mudq : int; mldq : int; mukeep : int; mlkeep : int; }
  end
module QuadratureTypes :
  sig type 'a quadrhsfn = float -> 'a -> 'a -> unit end
module SensitivityTypes :
  sig
    type 'd sensrhsfn_args = { t : float; y : 'd; y' : 'd; tmp : 'd double; }
    type 'a sensrhsfn_all = 'a sensrhsfn_args -> 'a array -> 'a array -> unit
    type 'a sensrhsfn1 = int -> 'a sensrhsfn_args -> 'a -> 'a -> unit
    type 'a sensrhsfn =
        AllAtOnce of 'a sensrhsfn_all option
      | OneByOne of 'a sensrhsfn1 option
    module QuadratureTypes :
      sig
        type 'd quadsensrhsfn_args = {
          t : float;
          y : 'd;
          s : 'd array;
          yq' : 'd;
          tmp : 'd double;
        }
        type 'a quadsensrhsfn = 'a quadsensrhsfn_args -> 'a array -> unit
      end
  end
module AdjointTypes' :
  sig
    type 'd brhsfn_args = { t : float; y : 'd; yb : 'd; }
    type 'a brhsfn_no_sens = 'a brhsfn_args -> 'a -> unit
    type 'a brhsfn_with_sens = 'a brhsfn_args -> 'a array -> 'a -> unit
    type 'a brhsfn =
        NoSens of 'a brhsfn_no_sens
      | WithSens of 'a brhsfn_with_sens
    module QuadratureTypes :
      sig
        type 'd bquadrhsfn_args = { t : float; y : 'd; yb : 'd; }
        type 'a bquadrhsfn_no_sens = 'a bquadrhsfn_args -> 'a -> unit
        type 'a bquadrhsfn_with_sens =
            'a bquadrhsfn_args -> 'a array -> 'a -> unit
        type 'a bquadrhsfn =
            NoSens of 'a bquadrhsfn_no_sens
          | WithSens of 'a bquadrhsfn_with_sens
      end
    type ('t, 'a) jacobian_arg = {
      jac_t : float;
      jac_y : 'a;
      jac_yb : 'a;
      jac_fyb : 'a;
      jac_tmp : 't;
    }
    module DirectTypes :
      sig
        type 'm jac_fn_no_sens =
            (Sundials.RealArray.t triple, Sundials.RealArray.t) jacobian_arg ->
            'm -> unit
        type 'm jac_fn_with_sens =
            (Sundials.RealArray.t triple, Sundials.RealArray.t) jacobian_arg ->
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
        type 'm linsys_fn_no_sens =
            (Sundials.RealArray.t triple, Sundials.RealArray.t) jacobian_arg ->
            'm -> bool -> float -> bool
        type 'm linsys_fn_with_sens =
            (Sundials.RealArray.t triple, Sundials.RealArray.t) jacobian_arg ->
            Sundials.RealArray.t array -> 'm -> bool -> float -> bool
        type 'm linsys_fn =
            LNoSens of 'm linsys_fn_no_sens
          | LWithSens of 'm linsys_fn_with_sens
        val no_linsysfn : 'a linsys_fn
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
        type 'a prec_setup_fn =
            (unit, 'a) jacobian_arg -> bool -> float -> bool
        type 'a precfns_no_sens = {
          prec_solve_fn : 'a prec_solve_fn;
          prec_setup_fn : 'a prec_setup_fn option;
        }
        type 'd jac_times_setup_fn_no_sens = (unit, 'd) jacobian_arg -> unit
        type 'a jac_times_vec_fn_no_sens =
            ('a, 'a) jacobian_arg -> 'a -> 'a -> unit
        type 'd prec_solve_fn_with_sens =
            (unit, 'd) jacobian_arg ->
            'd prec_solve_arg -> 'd array -> 'd -> unit
        type 'd prec_setup_fn_with_sens =
            (unit, 'd) jacobian_arg -> 'd array -> bool -> float -> bool
        type 'a precfns_with_sens = {
          prec_solve_fn_sens : 'a prec_solve_fn_with_sens;
          prec_setup_fn_sens : 'a prec_setup_fn_with_sens option;
        }
        type 'd jac_times_setup_fn_with_sens =
            (unit, 'd) jacobian_arg -> 'd array -> unit
        type 'd jac_times_vec_fn_with_sens =
            ('d, 'd) jacobian_arg -> 'd array -> 'd -> 'd -> unit
      end
  end
module CvodesBbdParamTypes :
  sig
    type 'a local_fn = 'a AdjointTypes'.brhsfn_args -> 'a -> unit
    type 'a comm_fn = 'a AdjointTypes'.brhsfn_args -> unit
    type 'a precfns = {
      local_fn : 'a local_fn;
      comm_fn : 'a comm_fn option;
    }
  end
module CvodesBbdTypes = CvodeBbdTypes
type cvode_mem
type c_weak_ref
type 'a rhsfn = float -> 'a -> 'a -> unit
type 'a rootsfn = float -> 'a -> Sundials.RealArray.t -> unit
type error_handler = Sundials.Util.error_details -> unit
type 'a error_weight_fun = 'a -> 'a -> unit
type 'd proj_fn = float -> 'd -> 'd -> float -> 'd option -> unit
val no_rhsfn : 'a -> 'b -> 'c -> 'd
type ('a, 'kind) session = {
  cvode : cvode_mem;
  backref : c_weak_ref;
  nroots : int;
  checkvec : ('a, 'kind) Nvector.t -> unit;
  mutable exn_temp : exn option;
  mutable rhsfn : 'a rhsfn;
  mutable rootsfn : 'a rootsfn;
  mutable errh : error_handler;
  mutable errw : 'a error_weight_fun;
  mutable error_file : Sundials.Logfile.t option;
  mutable projfn : 'a proj_fn;
  mutable monitorfn : ('a, 'kind) session -> unit;
  mutable ls_solver : LSI.held_linear_solver;
  mutable ls_callbacks : ('a, 'kind) linsolv_callbacks;
  mutable ls_precfns : 'a linsolv_precfns;
  mutable nls_solver :
    ('a, 'kind, ('a, 'kind) session, [ `Nvec ]) NLSI.nonlinear_solver option;
  mutable nls_rhsfn : 'a rhsfn;
  mutable sensext : ('a, 'kind) sensext;
}
and ('a, 'kind) linsolv_callbacks =
    NoCallbacks
  | DiagNoCallbacks
  | DlsDenseCallback of Sundials.Matrix.Dense.t DirectTypes.jac_callback *
      Sundials.Matrix.Dense.t DirectTypes.linsys_fn
  | DlsBandCallback of Sundials.Matrix.Band.t DirectTypes.jac_callback *
      Sundials.Matrix.Band.t DirectTypes.linsys_fn
  | BDlsDenseCallback of
      Sundials.Matrix.Dense.t AdjointTypes'.DirectTypes.jac_callback_no_sens *
      Sundials.Matrix.Dense.t AdjointTypes'.DirectTypes.linsys_fn
  | BDlsDenseCallbackSens of
      Sundials.Matrix.Dense.t
      AdjointTypes'.DirectTypes.jac_callback_with_sens *
      Sundials.Matrix.Dense.t AdjointTypes'.DirectTypes.linsys_fn
  | BDlsBandCallback of
      Sundials.Matrix.Band.t AdjointTypes'.DirectTypes.jac_callback_no_sens *
      Sundials.Matrix.Band.t AdjointTypes'.DirectTypes.linsys_fn
  | BDlsBandCallbackSens of
      Sundials.Matrix.Band.t AdjointTypes'.DirectTypes.jac_callback_with_sens *
      Sundials.Matrix.Band.t AdjointTypes'.DirectTypes.linsys_fn
  | SlsKluCallback : 's Sundials.Matrix.Sparse.t DirectTypes.jac_callback *
      's Sundials.Matrix.Sparse.t DirectTypes.linsys_fn -> ('a, 'kind)
                                                           linsolv_callbacks
  | BSlsKluCallback :
      's Sundials.Matrix.Sparse.t
      AdjointTypes'.DirectTypes.jac_callback_no_sens *
      's Sundials.Matrix.Sparse.t AdjointTypes'.DirectTypes.linsys_fn -> 
      ('a, 'kind) linsolv_callbacks
  | BSlsKluCallbackSens :
      's Sundials.Matrix.Sparse.t
      AdjointTypes'.DirectTypes.jac_callback_with_sens *
      's Sundials.Matrix.Sparse.t AdjointTypes'.DirectTypes.linsys_fn -> 
      ('a, 'kind) linsolv_callbacks
  | SlsSuperlumtCallback :
      's Sundials.Matrix.Sparse.t DirectTypes.jac_callback *
      's Sundials.Matrix.Sparse.t DirectTypes.linsys_fn -> ('a, 'kind)
                                                           linsolv_callbacks
  | BSlsSuperlumtCallback :
      's Sundials.Matrix.Sparse.t
      AdjointTypes'.DirectTypes.jac_callback_no_sens *
      's Sundials.Matrix.Sparse.t AdjointTypes'.DirectTypes.linsys_fn -> 
      ('a, 'kind) linsolv_callbacks
  | BSlsSuperlumtCallbackSens :
      's Sundials.Matrix.Sparse.t
      AdjointTypes'.DirectTypes.jac_callback_with_sens *
      's Sundials.Matrix.Sparse.t AdjointTypes'.DirectTypes.linsys_fn -> 
      ('a, 'kind) linsolv_callbacks
  | DirectCustomCallback : 'm DirectTypes.jac_callback *
      'm DirectTypes.linsys_fn -> ('a, 'kind) linsolv_callbacks
  | BDirectCustomCallback :
      'm AdjointTypes'.DirectTypes.jac_callback_no_sens *
      'm AdjointTypes'.DirectTypes.linsys_fn -> ('a, 'kind) linsolv_callbacks
  | BDirectCustomCallbackSens :
      'm AdjointTypes'.DirectTypes.jac_callback_with_sens *
      'm AdjointTypes'.DirectTypes.linsys_fn -> ('a, 'kind) linsolv_callbacks
  | SpilsCallback1 of 'a SpilsTypes'.jac_times_vec_fn option *
      'a SpilsTypes'.jac_times_setup_fn option
  | SpilsCallback2 of 'a rhsfn
  | BSpilsCallbackJTRhsfn of 'a rhsfn
  | BSpilsCallbackNoSens of
      'a AdjointTypes'.SpilsTypes'.jac_times_vec_fn_no_sens option *
      'a AdjointTypes'.SpilsTypes'.jac_times_setup_fn_no_sens option
  | BSpilsCallbackWithSens of
      'a AdjointTypes'.SpilsTypes'.jac_times_vec_fn_with_sens option *
      'a AdjointTypes'.SpilsTypes'.jac_times_setup_fn_with_sens option
and 'a linsolv_precfns =
    NoPrecFns
  | PrecFns of 'a SpilsTypes'.precfns
  | BPrecFns of 'a AdjointTypes'.SpilsTypes'.precfns_no_sens
  | BPrecFnsSens of 'a AdjointTypes'.SpilsTypes'.precfns_with_sens
  | BandedPrecFns
  | BBDPrecFns of 'a CvodeBbdParamTypes.precfns
  | BBBDPrecFns of 'a CvodesBbdParamTypes.precfns
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
  mutable senspvals : Sundials.RealArray.t option;
  mutable sensrhsfn : 'a SensitivityTypes.sensrhsfn_all;
  mutable sensrhsfn1 : 'a SensitivityTypes.sensrhsfn1;
  mutable quadsensrhsfn : 'a SensitivityTypes.QuadratureTypes.quadsensrhsfn;
  mutable fnls_solver :
    ('a, 'kind, ('a, 'kind) session) NLSI.nonlinear_solver_hold;
  mutable bsessions : ('a, 'kind) session list;
}
and ('a, 'kind) bsensext = {
  parent : ('a, 'kind) session;
  which : int;
  bnum_sensitivities : int;
  bsensarray : 'a array;
  mutable brhsfn : 'a AdjointTypes'.brhsfn_no_sens;
  mutable brhsfn_sens : 'a AdjointTypes'.brhsfn_with_sens;
  mutable bquadrhsfn : 'a AdjointTypes'.QuadratureTypes.bquadrhsfn_no_sens;
  mutable bquadrhsfn_sens :
    'a AdjointTypes'.QuadratureTypes.bquadrhsfn_with_sens;
  mutable checkbquadvec : ('a, 'kind) Nvector.t -> unit;
}
val revlookup_bsession :
  ('d, 'k) session -> cvode_mem -> ('d, 'k) session option
val ls_check_diag : ('a, 'b) session -> unit
val ls_check_direct : ('a, 'b) session -> unit
val ls_check_spils : ('a, 'b) session -> unit
val ls_check_spils_band : ('a, 'b) session -> unit
val ls_check_spils_bbd : ('a, 'b) session -> unit
type 'a serial_session = (Nvector_serial.data, 'a) session
  constraint 'a = [> Nvector_serial.kind ]
type ('data, 'kind) linear_solver =
    ('data, 'kind) session -> ('data, 'kind) Nvector.t -> unit
type 'a serial_linear_solver = (Nvector_serial.data, 'a) linear_solver
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
    type ('a, 'k) set_preconditioner =
        ('a, 'k) session -> ('a, 'k) Nvector.t -> unit
    type ('a, 'k) preconditioner =
        LSI.Iterative.preconditioning_type * ('a, 'k) set_preconditioner
    type 'a serial_preconditioner = (Nvector_serial.data, 'a) preconditioner
      constraint 'a = [> Nvector_serial.kind ]
  end
module AdjointTypes :
  sig
    type 'd brhsfn_args =
      'd AdjointTypes'.brhsfn_args = {
      t : float;
      y : 'd;
      yb : 'd;
    }
    type 'a brhsfn_no_sens = 'a brhsfn_args -> 'a -> unit
    type 'a brhsfn_with_sens = 'a brhsfn_args -> 'a array -> 'a -> unit
    type 'a brhsfn =
      'a AdjointTypes'.brhsfn =
        NoSens of 'a brhsfn_no_sens
      | WithSens of 'a brhsfn_with_sens
    module QuadratureTypes = AdjointTypes'.QuadratureTypes
    type ('t, 'a) jacobian_arg =
      ('t, 'a) AdjointTypes'.jacobian_arg = {
      jac_t : float;
      jac_y : 'a;
      jac_yb : 'a;
      jac_fyb : 'a;
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
        ('data, 'kind) bsession -> ('data, 'kind) nvector -> unit
    type 'a serial_linear_solver = (Nvector_serial.data, 'a) linear_solver
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
            (unit, 'a) AdjointTypes'.jacobian_arg ->
            'a prec_solve_arg -> 'a -> unit
        type 'a prec_setup_fn =
            (unit, 'a) AdjointTypes'.jacobian_arg -> bool -> float -> bool
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
            'd prec_solve_arg -> 'd array -> 'd -> unit
        type 'd prec_setup_fn_with_sens =
            (unit, 'd) AdjointTypes'.jacobian_arg ->
            'd array -> bool -> float -> bool
        type 'a precfns_with_sens =
          'a AdjointTypes'.SpilsTypes'.precfns_with_sens = {
          prec_solve_fn_sens : 'a prec_solve_fn_with_sens;
          prec_setup_fn_sens : 'a prec_setup_fn_with_sens option;
        }
        type 'd jac_times_setup_fn_with_sens =
            (unit, 'd) AdjointTypes'.jacobian_arg -> 'd array -> unit
        type 'd jac_times_vec_fn_with_sens =
            ('d, 'd) AdjointTypes'.jacobian_arg ->
            'd array -> 'd -> 'd -> unit
        type ('a, 'k) set_preconditioner =
            ('a, 'k) bsession ->
            ('a, 'k) session -> int -> ('a, 'k) nvector -> unit
        type ('a, 'k) preconditioner =
            LSI.Iterative.preconditioning_type * ('a, 'k) set_preconditioner
        type 'a serial_preconditioner =
            (Nvector_serial.data, 'a) preconditioner
          constraint 'a = [> Nvector_serial.kind ]
      end
  end
val read_weak_ref : ('a, 'kind) session Weak.t -> ('a, 'kind) session
val dummy_rhsfn : 'a -> 'b -> 'c -> 'd
val dummy_nlsrhsfn : 'a -> 'b -> 'c -> 'd
val dummy_rootsfn : 'a -> 'b -> 'c -> 'd
val dummy_errh : 'a -> 'b
val dummy_errw : 'a -> 'b -> 'c
val dummy_projfn : 'a -> 'b -> 'c -> 'd -> 'e -> 'f
val dummy_monitorfn : 'a -> 'b
val dummy_brhsfn_no_sens : 'a -> 'b -> 'c
val dummy_brhsfn_with_sens : 'a -> 'b -> 'c -> 'd
val dummy_bquadrhsfn_no_sens : 'a -> 'b -> 'c
val dummy_bquadrhsfn_with_sens : 'a -> 'b -> 'c -> 'd
val dummy_quadrhsfn : 'a -> 'b -> 'c -> 'd
val dummy_sensrhsfn : 'a -> 'b -> 'c -> 'd
val dummy_sensrhsfn1 : 'a -> 'b -> 'c -> 'd -> 'e
val dummy_quadsensrhsfn : 'a -> 'b -> 'c
