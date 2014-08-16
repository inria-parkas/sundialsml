type ('data, 'kind) nvector = ('data, 'kind) Sundials.nvector
type real_array = Sundials.RealArray.t
type root_array = Sundials.Roots.t
type root_val_array = Sundials.Roots.val_array
type 'a single_tmp = 'a
type 'a double_tmp = 'a * 'a
type 'a triple_tmp = 'a * 'a * 'a
type ('t, 'a) jacobian_arg =
  ('t, 'a) Ida_impl.jacobian_arg = {
  jac_t : float;
  jac_y : 'a;
  jac_y' : 'a;
  jac_res : 'a;
  jac_coef : float;
  jac_tmp : 't;
}
and 'a spils_callbacks =
  'a Ida_impl.spils_callbacks = {
  prec_solve_fn :
    (('a single_tmp, 'a) jacobian_arg -> 'a -> 'a -> float -> unit) option;
  prec_setup_fn : (('a triple_tmp, 'a) jacobian_arg -> unit) option;
  jac_times_vec_fn :
    (('a double_tmp, 'a) jacobian_arg -> 'a -> 'a -> unit) option;
}
type bandrange = Ida_impl.bandrange = { mupper : int; mlower : int; }
type dense_jac_fn =
    (real_array triple_tmp, real_array) jacobian_arg ->
    Dls.DenseMatrix.t -> unit
type band_jac_fn =
    bandrange ->
    (real_array triple_tmp, real_array) jacobian_arg ->
    Dls.BandMatrix.t -> unit
type 'a quadrhsfn = 'a Ida_impl.quadrhsfn (* = float -> 'a -> 'a -> 'a -> unit*)
type 'a sensresfn =
    float ->
    'a ->
    'a -> 'a -> 'a array -> 'a array -> 'a array -> 'a -> 'a -> 'a -> unit
type 'a quadsensrhsfn =
    float ->
    'a ->
    'a -> 'a array -> 'a array -> 'a -> 'a array -> 'a -> 'a -> 'a -> unit
module Bbd :
  sig
    type 'data callbacks =
      'data Ida_impl.Bbd.callbacks = {
      local_fn : float -> 'data -> 'data -> 'data -> unit;
      comm_fn : (float -> 'data -> 'data -> unit) option;
    }
  end
module B :
  sig
    type 'a resfnb = float -> 'a -> 'a -> 'a -> 'a -> 'a -> unit
    and 'a resfnbs =
        float -> 'a -> 'a -> 'a array -> 'a array -> 'a -> 'a -> 'a -> unit
    type 'a bresfn =
      'a Ida_impl.B.bresfn =
        Basic of 'a resfnb
      | WithSens of 'a resfnbs
    type 'a bquadrhsfn =
      'a Ida_impl.B.bquadrhsfn =
        Basic of 'a bquadrhsfn_basic
      | WithSens of 'a bquadrhsfn_withsens
    and 'a bquadrhsfn_basic = float -> 'a -> 'a -> 'a -> 'a -> 'a -> unit
    and 'a bquadrhsfn_withsens =
        float -> 'a -> 'a -> 'a array -> 'a array -> 'a -> 'a -> 'a -> unit
    type ('t, 'a) jacobian_arg =
      ('t, 'a) Ida_impl.B.jacobian_arg = {
      jac_t : float;
      jac_y : 'a;
      jac_y' : 'a;
      jac_yb : 'a;
      jac_y'b : 'a;
      jac_resb : 'a;
      jac_coef : float;
      jac_tmp : 't;
    }
    type 'a spils_callbacks =
      'a Ida_impl.B.spils_callbacks = {
      prec_solve_fn :
        (('a single_tmp, 'a) jacobian_arg -> 'a -> 'a -> float -> unit)
        option;
      prec_setup_fn : (('a triple_tmp, 'a) jacobian_arg -> unit) option;
      jac_times_vec_fn :
        (('a single_tmp, 'a) jacobian_arg -> 'a -> 'a -> unit) option;
    }
    type dense_jac_fn =
        (real_array triple_tmp, real_array) jacobian_arg ->
        Dls.DenseMatrix.t -> unit
    type band_jac_fn =
        bandrange ->
        (real_array triple_tmp, real_array) jacobian_arg ->
        Dls.BandMatrix.t -> unit
    module Bbd :
      sig
        type 'data callbacks =
          'data Ida_impl.B.Bbd.callbacks =
          {
            local_fn : float -> 'data -> 'data -> 'data
                       -> 'data -> 'data -> unit;
            comm_fn  : (float -> 'data -> 'data -> 'data -> 'data -> unit)
                       option;
          }
      end
  end
type conv_fail = Ida_impl.conv_fail = NoFailures | FailBadJ | FailOther
type 'a alternate_linsolv =
  'a Ida_impl.alternate_linsolv = {
  linit : (unit -> bool) option;
  lsetup : (conv_fail -> 'a -> 'a -> 'a triple_tmp -> bool) option;
  lsolve : 'a -> 'a -> 'a -> 'a -> unit;
  lfree : (unit -> unit) option;
}
type ida_mem = Ida_impl.ida_mem
type c_weak_ref = Ida_impl.c_weak_ref
type ida_file = Ida_impl.ida_file
type ('a, 'kind) linsolv_callbacks =
  ('a, 'kind) Ida_impl.linsolv_callbacks =
    NoCallbacks
  | DenseCallback of dense_jac_fn
  | BandCallback of band_jac_fn
  | SpilsCallback of 'a spils_callbacks
  | BBDCallback of 'a Bbd.callbacks
  | AlternateCallback of 'a alternate_linsolv
  | BDenseCallback of B.dense_jac_fn
  | BBandCallback of B.band_jac_fn
  | BSpilsCallback of 'a B.spils_callbacks
  | BBBDCallback of 'a B.Bbd.callbacks
type ('a, 'kind) session =
  ('a, 'kind) Ida_impl.session = {
  ida : ida_mem;
  backref : c_weak_ref;
  nroots : int;
  err_file : ida_file;
  mutable exn_temp : exn option;
  mutable resfn : float -> 'a -> 'a -> 'a -> unit;
  mutable rootsfn : float -> 'a -> 'a -> root_val_array -> unit;
  mutable errh : Sundials.error_details -> unit;
  mutable errw : 'a -> 'a -> unit;
  mutable ls_callbacks : ('a, 'kind) linsolv_callbacks;
  mutable sensext : ('a, 'kind) sensext;
  mutable safety_check_flags : int;
}
and ('a, 'kind) sensext =
  ('a, 'kind) Ida_impl.sensext =
    NoSensExt
  | FwdSensExt of ('a, 'kind) fsensext
  | BwdSensExt of ('a, 'kind) bsensext
and ('a, 'kind) fsensext =
  ('a, 'kind) Ida_impl.fsensext = {
  mutable quadrhsfn : 'a quadrhsfn;
  mutable num_sensitivities : int;
  mutable sensarray1 : 'a array;
  mutable sensarray2 : 'a array;
  mutable sensarray3 : 'a array;
  mutable senspvals : Sundials.RealArray.t option;
  mutable sensresfn : 'a sensresfn;
  mutable quadsensrhsfn : 'a quadsensrhsfn;
  mutable bsessions : ('a, 'kind) session list;
}
and ('a, 'kind) bsensext =
  ('a, 'kind) Ida_impl.bsensext = {
  parent : ('a, 'kind) session;
  which : int;
  bnum_sensitivities : int;
  bsensarray1 : 'a array;
  bsensarray2 : 'a array;
  mutable resfnb : 'a B.resfnb;
  mutable resfnbs : 'a B.resfnbs;
  mutable bquadrhsfn : 'a B.bquadrhsfn_basic;
  mutable bquadrhsfn1 : 'a B.bquadrhsfn_withsens;
}
type ('a, 'k) bsession =
  ('a, 'k) Ida_impl.bsession =
    Bsession of ('a, 'k) session
val tosession : ('a, 'b) bsession -> ('a, 'b) session
type ('data, 'kind) linear_solver =
    ('data, 'kind) session ->
    ('data, 'kind) nvector -> ('data, 'kind) nvector -> unit
type ('data, 'kind) blinear_solver =
    ('data, 'kind) bsession ->
    ('data, 'kind) nvector -> ('data, 'kind) nvector -> unit
external crash : string -> unit = "sundials_crash"
val dummy_resfn : 'a -> 'b -> 'c -> 'd -> unit
val dummy_rootsfn : 'a -> 'b -> 'c -> 'd -> unit
val dummy_errh : 'a -> unit
val dummy_errw : 'a -> 'b -> unit
val dummy_resfnb : 'a -> 'b -> 'c -> 'd -> 'e -> 'f -> unit
val dummy_resfnbs : 'a -> 'b -> 'c -> 'd -> 'e -> 'f -> 'g -> 'h -> unit
val dummy_bquadrhsfn : 'a -> 'b -> 'c -> 'd -> 'e -> 'f -> unit
val dummy_bquadrhsfn1 : 'a -> 'b -> 'c -> 'd -> 'e -> 'f -> 'g -> 'h -> unit
val dummy_quadrhsfn : 'a -> 'b -> 'c -> 'd -> unit
val dummy_sensresfn :
  'a -> 'b -> 'c -> 'd -> 'e -> 'f -> 'g -> 'h -> 'i -> 'j -> unit
val dummy_quadsensrhsfn :
  'a -> 'b -> 'c -> 'd -> 'e -> 'f -> 'g -> 'h -> 'i -> 'j -> unit
type serial_session = (real_array, Nvector_serial.kind) session
external c_alloc_nvector_array : int -> 'a array
  = "c_idas_alloc_nvector_array"
val add_fwdsensext : ('a, 'b) session -> unit
val num_sensitivities : ('a, 'b) session -> int
val read_weak_ref : ('a, 'kind) session Weak.t -> ('a, 'kind) session
val read_weak_fwd_ref :
  ('a, 'b) session Weak.t -> ('a, 'b) session * ('a, 'b) fsensext
val read_weak_bwd_ref :
  ('a, 'b) session Weak.t -> ('a, 'b) session * ('a, 'b) bsensext
val adjust_retcode : ('a, 'b) session -> ('c -> 'd) -> 'c -> int
val adjust_retcode_and_bool :
  ('a, 'b) session -> ('c -> bool) -> 'c -> bool * int
val call_quadrhsfn :
  ('a, 'b) session Weak.t -> float -> 'a -> 'a -> 'a -> int
val call_bquadrhsfn :
  ('a, 'b) session Weak.t -> float -> 'a -> 'a -> 'a -> 'a -> 'a -> int
val call_bprecsetupfn :
  ('a, 'b) session Weak.t -> ('a triple_tmp, 'a) B.jacobian_arg -> int
val call_bprecsolvefn :
  ('a, 'b) session Weak.t ->
  ('a single_tmp, 'a) B.jacobian_arg -> 'a -> 'a -> float -> int
val call_bjactimesfn :
  ('a, 'b) session Weak.t ->
  ('a single_tmp, 'a) B.jacobian_arg -> 'a -> 'a -> int
val call_bjacfn :
  ('a, 'b) session Weak.t ->
  (real_array triple_tmp, real_array) B.jacobian_arg ->
  Dls.DenseMatrix.t -> int
val call_bbandjacfn :
  ('a, 'b) session Weak.t ->
  bandrange ->
  (real_array triple_tmp, real_array) B.jacobian_arg ->
  Dls.BandMatrix.t -> int
module Quadrature :
  sig
    exception QuadNotInitialized
    exception QuadRhsFuncFailure
    exception FirstQuadRhsFuncErr
    exception RepeatedQuadRhsFuncErr
    exception UnrecoverableQuadRhsFuncErr
    val fwdsensext : ('a, 'b) session -> ('a, 'b) fsensext
    type 'a quadrhsfn = float -> 'a -> 'a -> 'a -> unit
    external c_quad_init : ('a, 'k) session -> ('a, 'k) nvector -> unit
      = "c_idas_quad_init"
    val init : ('a, 'b) session -> 'a quadrhsfn -> ('a, 'b) nvector -> unit
    external reinit : ('a, 'k) session -> ('a, 'k) nvector -> unit
      = "c_idas_quad_reinit"
    external set_err_con : ('a, 'k) session -> bool -> unit
      = "c_idas_quad_set_err_con"
    external sv_tolerances :
      ('a, 'k) session -> float -> ('a, 'k) nvector -> unit
      = "c_idas_quad_sv_tolerances"
    external ss_tolerances : ('a, 'k) session -> float -> float -> unit
      = "c_idas_quad_ss_tolerances"
    type ('a, 'k) tolerance =
        NoStepSizeControl
      | SStolerances of float * float
      | SVtolerances of float * ('a, 'k) nvector
    val set_tolerances : ('a, 'b) session -> ('a, 'b) tolerance -> unit
    external get : ('a, 'k) session -> ('a, 'k) nvector -> float
      = "c_idas_quad_get"
    external get_dky :
      ('a, 'k) session -> float -> int -> ('a, 'k) nvector -> unit
      = "c_idas_quad_get_dky"
    external get_num_rhs_evals : ('a, 'k) session -> int
      = "c_idas_quad_get_num_rhs_evals"
    external get_num_err_test_fails : ('a, 'k) session -> int
      = "c_idas_quad_get_num_err_test_fails"
    external get_err_weights : ('a, 'k) session -> ('a, 'k) nvector -> unit
      = "c_idas_quad_get_err_weights"
    external get_stats : ('a, 'k) session -> int * int
      = "c_idas_quad_get_stats"
  end
module Sensitivity :
  sig
    type ('a, 'k) tolerance =
        SStolerances of float * Sundials.RealArray.t
      | SVtolerances of float * ('a, 'k) nvector array
      | EEtolerances
    external set_err_con : ('a, 'k) session -> bool -> unit
      = "c_idas_sens_set_err_con"
    external ss_tolerances :
      ('a, 'k) session -> float -> Sundials.RealArray.t -> unit
      = "c_idas_sens_ss_tolerances"
    external ee_tolerances : ('a, 'k) session -> unit
      = "c_idas_sens_ee_tolerances"
    external sv_tolerances :
      ('a, 'k) session -> float -> ('a, 'k) nvector array -> unit
      = "c_idas_sens_sv_tolerances"
    val set_tolerances : ('a, 'b) session -> ('a, 'b) tolerance -> unit
    exception SensNotInitialized
    exception SensResFuncErr
    exception FirstSensResFuncErr
    exception RepeatedSensRhsFuncErr
    exception UnrecoverableSensRhsFuncErr
    exception BadIS
    val fwdsensext : ('a, 'b) session -> ('a, 'b) fsensext
    type sens_method = Simultaneous | Staggered
    type 'a sensresfn = 'a Ida_impl.sensresfn
    type sens_params = {
      pvals : Sundials.RealArray.t option;
      pbar : Sundials.RealArray.t option;
      plist : int array option;
    }
    val no_sens_params : sens_params
    external c_sens_init :
      ('a, 'k) session ->
      sens_method ->
      bool -> ('a, 'k) nvector array -> ('a, 'k) nvector array -> unit
      = "c_idas_sens_init"
    external c_set_params : ('a, 'k) session -> sens_params -> unit
      = "c_idas_sens_set_params"
    val set_params : ('a, 'b) session -> sens_params -> unit
    val init :
      ('a, 'b) session ->
      ('a, 'b) tolerance ->
      sens_method ->
      sens_params ->
      'a sensresfn option ->
      ('a, 'b) nvector array -> ('a, 'b) nvector array -> unit
    external c_reinit :
      ('a, 'k) session ->
      sens_method -> ('a, 'k) nvector array -> ('a, 'k) nvector array -> unit
      = "c_idas_sens_reinit"
    val reinit :
      ('a, 'b) session ->
      sens_method -> ('a, 'b) nvector array -> ('a, 'b) nvector array -> unit
    external toggle_off : ('a, 'k) session -> unit = "c_idas_sens_toggle_off"
    external c_get : ('a, 'k) session -> ('a, 'k) nvector array -> float
      = "c_idas_sens_get"
    val get : ('a, 'b) session -> ('a, 'b) nvector array -> float
    external c_get_dky :
      ('a, 'k) session -> float -> int -> ('a, 'k) nvector array -> unit
      = "c_idas_sens_get_dky"
    val get_dky :
      ('a, 'b) session -> float -> int -> ('a, 'b) nvector array -> unit
    external get1 : ('a, 'k) session -> int -> ('a, 'k) nvector -> float
      = "c_idas_sens_get1"
    external get_dky1 :
      ('a, 'k) session -> float -> int -> int -> ('a, 'k) nvector -> unit
      = "c_idas_sens_get_dky1"
    type dq_method = DQCentered | DQForward
    external set_dq_method : ('a, 'k) session -> dq_method -> float -> unit
      = "c_idas_sens_set_dq_method"
    external set_max_nonlin_iters : ('a, 'k) session -> int -> unit
      = "c_idas_sens_set_max_nonlin_iters"
    external get_num_res_evals : ('a, 'k) session -> int
      = "c_idas_sens_get_num_res_evals"
    external get_num_res_evals_sens : ('a, 'k) session -> int
      = "c_idas_sens_get_num_res_evals_sens"
    external get_num_err_test_fails : ('a, 'k) session -> int
      = "c_idas_sens_get_num_err_test_fails"
    external get_num_lin_solv_setups : ('a, 'k) session -> int
      = "c_idas_sens_get_num_lin_solv_setups"
    type sensitivity_stats = {
      num_res_evals : int;
      num_sens_evals : int;
      num_err_test_fails : int;
      num_lin_solv_setups : int;
    }
    external get_stats : ('a, 'k) session -> sensitivity_stats
      = "c_idas_sens_get_stats"
    external c_get_err_weights :
      ('a, 'k) session -> ('a, 'k) nvector array -> unit
      = "c_idas_sens_get_err_weights"
    val get_err_weights : ('a, 'b) session -> ('a, 'b) nvector array -> unit
    external c_sens_calc_ic_ya_yd' :
      ('a, 'k) session ->
      ('a, 'k) nvector option ->
      ('a, 'k) nvector option ->
      ('a, 'k) nvector array option ->
      ('a, 'k) nvector array option -> ('a, 'k) nvector -> float -> unit
      = "c_ida_sens_calc_ic_ya_ydp_byte" "c_ida_sens_calc_ic_ya_ydp"
    external c_sens_calc_ic_y :
      ('a, 'k) session ->
      ('a, 'k) nvector option -> ('a, 'k) nvector option -> float -> unit
      = "c_ida_sens_calc_ic_y"
    val calc_ic_ya_yd' :
      ('a, 'b) session ->
      ?y:('a, 'b) nvector ->
      ?y':('a, 'b) nvector ->
      ?ys:('a, 'b) nvector array ->
      ?y's:('a, 'b) nvector array -> ('a, 'b) nvector -> float -> unit
    val calc_ic_y :
      ('a, 'b) session ->
      ?y:('a, 'b) nvector -> ?ys:('a, 'b) nvector -> float -> unit
    external get_num_nonlin_solv_iters : ('a, 'k) session -> int
      = "c_idas_sens_get_num_nonlin_solv_iters"
    external get_num_nonlin_solv_conv_fails : ('a, 'k) session -> int
      = "c_idas_sens_get_num_nonlin_solv_conv_fails"
    external get_nonlin_solv_stats : ('a, 'k) session -> int * int
      = "c_idas_sens_get_nonlin_solv_stats"
    module Quadrature :
      sig
        exception QuadSensNotInitialized
        exception QuadSensRhsFuncFailure
        exception FirstQuadSensRhsFuncErr
        exception RepeatedQuadSensRhsFuncErr
        type 'a quadsensrhsfn =
            float ->
            'a ->
            'a ->
            'a array -> 'a array -> 'a -> 'a array -> 'a -> 'a -> 'a -> unit
        external c_quadsens_init :
          ('a, 'k) session -> bool -> ('a, 'k) nvector array -> unit
          = "c_idas_quadsens_init"
        val init :
          ('a, 'b) session ->
          'a quadsensrhsfn option -> ('a, 'b) nvector array -> unit
        external c_reinit :
          ('a, 'k) session -> ('a, 'k) nvector array -> unit
          = "c_idas_quadsens_reinit"
        val reinit : ('a, 'b) session -> ('a, 'b) nvector array -> unit
        type ('a, 'k) tolerance =
            NoStepSizeControl
          | SStolerances of float * Sundials.RealArray.t
          | SVtolerances of float * ('a, 'k) nvector array
          | EEtolerances
        external set_err_con : ('a, 'k) session -> bool -> unit
          = "c_idas_quadsens_set_err_con"
        external ss_tolerances :
          ('a, 'k) session -> float -> Sundials.RealArray.t -> unit
          = "c_idas_quadsens_ss_tolerances"
        external sv_tolerances :
          ('a, 'k) session -> float -> ('a, 'k) nvector array -> unit
          = "c_idas_quadsens_sv_tolerances"
        external ee_tolerances : ('a, 'k) session -> unit
          = "c_idas_quadsens_ee_tolerances"
        val set_tolerances : ('a, 'b) session -> ('a, 'b) tolerance -> unit
        external c_get : ('a, 'k) session -> ('a, 'k) nvector array -> float
          = "c_idas_quadsens_get"
        val get : ('a, 'b) session -> ('a, 'b) nvector array -> float
        external get1 : ('a, 'k) session -> int -> ('a, 'k) nvector -> float
          = "c_idas_quadsens_get1"
        external c_get_dky :
          ('a, 'k) session -> float -> int -> ('a, 'k) nvector array -> unit
          = "c_idas_quadsens_get_dky"
        val get_dky :
          ('a, 'b) session -> float -> int -> ('a, 'b) nvector array -> unit
        external get_dky1 :
          ('a, 'k) session -> float -> int -> int -> ('a, 'k) nvector -> unit
          = "c_idas_quadsens_get_dky1"
        external get_num_rhs_evals : ('a, 'k) session -> int
          = "c_idas_quadsens_get_num_rhs_evals"
        external get_num_err_test_fails : ('a, 'k) session -> int
          = "c_idas_quadsens_get_num_err_test_fails"
        external c_get_err_weights :
          ('a, 'k) session -> ('a, 'k) nvector array -> unit
          = "c_idas_quadsens_get_err_weights"
        val get_err_weights :
          ('a, 'b) session -> ('a, 'b) nvector array -> unit
        external get_stats : ('a, 'k) session -> int * int
          = "c_idas_quadsens_get_stats"
      end
  end
module Adjoint :
  sig
    exception AdjointNotInitialized
    exception NoForwardCall
    exception ForwardReinitializationFailed
    exception ForwardFailed
    exception NoBackwardProblem
    exception BadFinalTime
    exception BadOutputTime
    type ('a, 'k) bsession = ('a, 'k) Ida_impl.bsession
    type serial_bsession = (real_array, Nvector_serial.kind) bsession
    type interpolation = IPolynomial | IHermite
    val init : ('a, 'b) session -> int -> interpolation -> unit
    val fwdsensext : ('a, 'b) session -> ('a, 'b) fsensext
    val set_var_types : ('a, 'b) bsession -> ('a, 'b) nvector -> unit
    val set_suppress_alg : ('a, 'b) bsession -> bool -> unit
    val calc_ic :
      ('a, 'b) bsession ->
      ?yb:('a, 'b) nvector ->
      ?y'b:('a, 'b) nvector ->
      float -> ('a, 'b) nvector -> ('a, 'b) nvector -> unit
    val calc_ic_sens :
      ('a, 'b) bsession ->
      ?yb:('a, 'b) nvector ->
      ?y'b:('a, 'b) nvector ->
      float ->
      ('a, 'b) nvector ->
      ('a, 'b) nvector ->
      ('a, 'b) nvector array -> ('a, 'b) nvector array -> unit
    external forward_normal :
      ('a, 'k) session ->
      float ->
      ('a, 'k) nvector ->
      ('a, 'k) nvector -> float * int * Sundials.solver_result
      = "c_idas_adj_forward_normal"
    external forward_one_step :
      ('a, 'k) session ->
      float ->
      ('a, 'k) nvector ->
      ('a, 'k) nvector -> float * int * Sundials.solver_result
      = "c_idas_adj_forward_one_step"
    type 'a bresfn =
      'a B.bresfn =
        Basic of 'a B.resfnb
      | WithSens of 'a B.resfnbs
    type 'a single_tmp = 'a
    type 'a triple_tmp = 'a * 'a * 'a
    type ('t, 'a) jacobian_arg =
      ('t, 'a) B.jacobian_arg = {
      jac_t : float;
      jac_y : 'a;
      jac_y' : 'a;
      jac_yb : 'a;
      jac_y'b : 'a;
      jac_resb : 'a;
      jac_coef : float;
      jac_tmp : 't;
    }
    type bandrange = Ida_impl.bandrange = { mupper : int; mlower : int; }
    type ('data, 'kind) linear_solver =
        ('data, 'kind) Ida_impl.blinear_solver
    type serial_linear_solver =
        (real_array, Nvector_serial.kind) linear_solver
    type ('data, 'kind) iter =
        Newton of ('data, 'kind) linear_solver
      | Functional
    type ('a, 'k) tolerance =
        SStolerances of float * float
      | SVtolerances of float * ('a, 'k) nvector
    external ss_tolerances :
      ('a, 'k) session -> int -> float -> float -> unit
      = "c_idas_adj_ss_tolerances"
    external sv_tolerances :
      ('a, 'k) session -> int -> float -> ('a, 'k) nvector -> unit
      = "c_idas_adj_sv_tolerances"
    val set_tolerances : ('a, 'b) bsession -> ('a, 'b) tolerance -> unit
    val set_linear_solver :
      ('a, 'b) bsession ->
      (('a, 'b) bsession -> 'c -> 'd -> 'e) -> 'c -> 'd -> 'e
    external bsession_finalize : ('a, 'k) session -> unit
      = "c_idas_adj_bsession_finalize"
    external c_init_backward :
      ('a, 'k) session ->
      ('a, 'k) session Weak.t ->
      float ->
      ('a, 'k) nvector ->
      ('a, 'k) nvector -> bool -> ida_mem * int * c_weak_ref * ida_file
      = "c_idas_adj_init_backward_byte" "c_idas_adj_init_backward"
    val init_backward :
      ('a, 'b) session ->
      (('a, 'b) bsession -> ('a, 'b) nvector -> ('a, 'b) nvector -> 'c) ->
      ('a, 'b) tolerance ->
      'a bresfn ->
      float -> ('a, 'b) nvector -> ('a, 'b) nvector -> ('a, 'b) bsession
    external c_reinit :
      ('a, 'k) session ->
      int -> float -> ('a, 'k) nvector -> ('a, 'k) nvector -> unit
      = "c_idas_adj_reinit"
    val reinit :
      ('a, 'b) bsession ->
      float -> ('a, 'b) nvector -> ('a, 'b) nvector -> unit
    external backward_normal : ('a, 'k) session -> float -> unit
      = "c_idas_adj_backward_normal"
    external backward_one_step : ('a, 'k) session -> float -> unit
      = "c_idas_adj_backward_one_step"
    external c_get :
      ('a, 'k) session ->
      int -> ('a, 'k) nvector -> ('a, 'k) nvector -> float = "c_idas_adj_get"
    val get :
      ('a, 'b) bsession -> ('a, 'b) nvector -> ('a, 'b) nvector -> float
    val get_dky :
      ('a, 'b) bsession -> float -> int -> ('a, 'b) Ida.nvector -> unit
    external c_set_max_ord : ('a, 'k) session -> int -> int -> unit
      = "c_idas_adj_set_max_ord"
    val set_max_ord : ('a, 'b) bsession -> int -> unit
    external c_set_max_num_steps : ('a, 'k) session -> int -> int -> unit
      = "c_idas_adj_set_max_num_steps"
    val set_max_num_steps : ('a, 'b) bsession -> int -> unit
    external c_set_init_step : ('a, 'k) session -> int -> float -> unit
      = "c_idas_adj_set_init_step"
    val set_init_step : ('a, 'b) bsession -> float -> unit
    external c_set_max_step : ('a, 'k) session -> int -> float -> unit
      = "c_idas_adj_set_max_step"
    val set_max_step : ('a, 'b) bsession -> float -> unit
    module Dls :
      sig
        type dense_jac_fn =
            (real_array triple_tmp, real_array) jacobian_arg ->
            Dls.DenseMatrix.t -> unit
        type band_jac_fn =
            bandrange ->
            (real_array triple_tmp, real_array) jacobian_arg ->
            Dls.BandMatrix.t -> unit
        external c_dls_dense : serial_session -> int -> int -> bool -> unit
          = "c_idas_adj_dls_dense"
        external c_dls_band :
          serial_session * int -> int -> int -> int -> bool -> unit
          = "c_idas_adj_dls_band"
        val dense :
          B.dense_jac_fn option ->
          (real_array, Nvector_serial.kind) bsession ->
          (Sundials.RealArray.t, 'a) Sundials.nvector -> 'b -> unit
        type ('data, 'kind) linear_solver =
            ('data, 'kind) bsession -> ('data, 'kind) nvector -> unit
        val band :
          bandrange ->
          B.band_jac_fn option ->
          (real_array, Nvector_serial.kind) bsession ->
          (Sundials.RealArray.t, 'a) Sundials.nvector -> 'b -> unit
      end
    module Spils :
      sig
        type gramschmidt_type =
          Spils.gramschmidt_type =
            ModifiedGS
          | ClassicalGS
        type preconditioning_type =
          Spils.preconditioning_type =
            PrecNone
          | PrecLeft
          | PrecRight
          | PrecBoth
        type 'a callbacks =
          'a B.spils_callbacks = {
          prec_solve_fn :
            (('a single_tmp, 'a) jacobian_arg -> 'a -> 'a -> float -> unit)
            option;
          prec_setup_fn : (('a triple_tmp, 'a) jacobian_arg -> unit) option;
          jac_times_vec_fn :
            (('a single_tmp, 'a) jacobian_arg -> 'a -> 'a -> unit) option;
        }
        val no_precond : 'a callbacks
        external c_spils_set_preconditioner :
          ('a, 'k) session -> int -> bool -> bool -> unit
          = "c_idas_adj_spils_set_preconditioner"
        external c_spils_spgmr : ('a, 'k) session -> int -> int -> unit
          = "c_idas_adj_spils_spgmr"
        external c_spils_spbcg : ('a, 'k) session -> int -> int -> unit
          = "c_idas_adj_spils_spbcg"
        external c_spils_sptfqmr : ('a, 'k) session -> int -> int -> unit
          = "c_idas_adj_spils_sptfqmr"
        val set_precond :
          ('a, 'b) bsession ->
          ('c, 'd) session -> int -> 'a callbacks -> unit
        val spgmr :
          int option -> 'a callbacks -> ('a, 'b) bsession -> 'c -> unit
        val spbcg :
          int option -> 'a callbacks -> ('a, 'b) bsession -> 'c -> unit
        val sptfqmr :
          int option -> 'a callbacks -> ('a, 'b) bsession -> 'c -> unit
        external set_gs_type :
          ('a, 'k) bsession -> Spils.gramschmidt_type -> unit
          = "c_idas_adj_spils_set_gs_type"
        external set_eps_lin : ('a, 'k) bsession -> float -> unit
          = "c_idas_adj_spils_set_eps_lin"
        external c_set_maxl : ('a, 'k) bsession -> int -> unit
          = "c_idas_adj_spils_set_maxl"
        val set_maxl : ('a, 'b) bsession -> int option -> unit
        val get_work_space : ('a, 'b) bsession -> int * int
        val get_num_lin_iters : ('a, 'b) bsession -> int
        val get_num_conv_fails : ('a, 'b) bsession -> int
        val get_num_prec_evals : ('a, 'b) bsession -> int
        val get_num_prec_solves : ('a, 'b) bsession -> int
        val get_num_jtimes_evals : ('a, 'b) bsession -> int
        val get_num_res_evals : ('a, 'b) bsession -> int
      end
    val get_work_space : ('a, 'b) bsession -> int * int
    val get_num_steps : ('a, 'b) bsession -> int
    val get_num_res_evals : ('a, 'b) bsession -> int
    val get_num_lin_solv_setups : ('a, 'b) bsession -> int
    val get_num_err_test_fails : ('a, 'b) bsession -> int
    val get_last_order : ('a, 'b) bsession -> int
    val get_current_order : ('a, 'b) bsession -> int
    val get_last_step : ('a, 'b) bsession -> float
    val get_current_step : ('a, 'b) bsession -> float
    val get_actual_init_step : ('a, 'b) bsession -> float
    val get_current_time : ('a, 'b) bsession -> float
    val get_tol_scale_factor : ('a, 'b) bsession -> float
    val get_err_weights : ('a, 'b) bsession -> ('a, 'b) Ida.nvector -> unit
    val get_est_local_errors :
      ('a, 'b) bsession -> ('a, 'b) Ida.nvector -> unit
    val get_integrator_stats : ('a, 'b) bsession -> Ida.integrator_stats
    val print_integrator_stats : ('a, 'b) bsession -> unit
    val get_num_nonlin_solv_iters : ('a, 'b) bsession -> int
    val get_num_nonlin_solv_conv_fails : ('a, 'b) bsession -> int
    val get_nonlin_solv_stats : ('a, 'b) bsession -> int * int
    module Quadrature :
      sig
        type 'a bquadrhsfn =
          'a B.bquadrhsfn =
            Basic of (float -> 'a -> 'a -> 'a -> 'a -> 'a -> unit)
          | WithSens of
              (float ->
               'a -> 'a -> 'a array -> 'a array -> 'a -> 'a -> 'a -> unit)
        external c_quad_initb :
          ('a, 'k) session -> int -> ('a, 'k) nvector -> unit
          = "c_idas_adjquad_initb"
        external c_quad_initbs :
          ('a, 'k) session -> int -> ('a, 'k) nvector -> unit
          = "c_idas_adjquad_initbs"
        val init :
          ('a, 'b) bsession -> 'a bquadrhsfn -> ('a, 'b) nvector -> unit
        external c_reinit :
          ('a, 'k) session -> int -> ('a, 'k) nvector -> unit
          = "c_idas_adjquad_reinit"
        val reinit : ('a, 'b) bsession -> ('a, 'b) nvector -> unit
        external c_get : ('a, 'k) session -> int -> ('a, 'k) nvector -> float
          = "c_idas_adjquad_get"
        val get : ('a, 'b) bsession -> ('a, 'b) nvector -> float
        type ('a, 'k) tolerance =
            NoStepSizeControl
          | SStolerances of float * float
          | SVtolerances of float * ('a, 'k) nvector
        external set_err_con : ('a, 'k) session -> int -> bool -> unit
          = "c_idas_adjquad_set_err_con"
        external sv_tolerances :
          ('a, 'k) session -> int -> float -> ('a, 'k) nvector -> unit
          = "c_idas_adjquad_sv_tolerances"
        external ss_tolerances :
          ('a, 'k) session -> int -> float -> float -> unit
          = "c_idas_adjquad_ss_tolerances"
        val set_tolerances : ('a, 'b) bsession -> ('a, 'b) tolerance -> unit
        val get_num_rhs_evals : ('a, 'b) bsession -> int
        val get_num_err_test_fails : ('a, 'b) bsession -> int
        val get_err_weights : ('a, 'b) bsession -> ('a, 'b) nvector -> unit
        val get_stats : ('a, 'b) bsession -> int * int
      end
  end
