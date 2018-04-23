(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

include Ida_impl

(* "Simulate" Linear Solvers in Sundials < 3.0.0 *)
let in_compat_mode =
  match Sundials.sundials_version with
  | 2,_,_ -> true
  | _ -> false

external c_alloc_nvector_array : int -> 'a array
    = "c_idas_alloc_nvector_array"

let add_fwdsensext s =
  match s.sensext with
  | FwdSensExt se -> ()
  | BwdSensExt _ -> failwith "Quadrature.add_fwdsensext: internal error"
  | NoSensExt ->
      s.sensext <- FwdSensExt {
        num_sensitivities = 0;
        sensarray1        = c_alloc_nvector_array 0;
        sensarray2        = c_alloc_nvector_array 0;
        sensarray3        = c_alloc_nvector_array 0;
        quadrhsfn         = dummy_quadrhsfn;
        checkquadvec      = (fun _ -> raise Nvector.IncompatibleNvector);
        has_quad          = true;
        senspvals         = None;
        sensresfn         = dummy_sensresfn;
        quadsensrhsfn     = dummy_quadsensrhsfn;
        bsessions         = [];
      }

let num_sensitivities s =
  match s.sensext with
  | FwdSensExt se -> se.num_sensitivities
  | BwdSensExt se -> se.bnum_sensitivities
  | _ -> 0

let ocheck checkfn oy =
  match oy with
  | Some y -> checkfn y
  | None -> ()

module Quadrature = struct (* {{{ *)
  include QuadratureTypes

  exception QuadNotInitialized
  exception QuadRhsFuncFailure
  exception FirstQuadRhsFuncFailure
  exception RepeatedQuadRhsFuncFailure

  let fwdsensext s =
    match s.sensext with
    | FwdSensExt se -> se
    | _ -> raise QuadNotInitialized

  external c_quad_init : ('a, 'k) session -> ('a, 'k) Nvector.t -> unit
      = "c_idas_quad_init"

  let init session f yQ0 =
    add_fwdsensext session;
    let s = fwdsensext session in
    s.quadrhsfn <- f;
    s.checkquadvec <- Nvector.check yQ0;
    c_quad_init session yQ0;
    s.has_quad <- true

  external c_reinit : ('a, 'k) session -> ('a, 'k) Nvector.t -> unit
    = "c_idas_quad_reinit"

  let reinit s v0 =
    let se = fwdsensext s in
    if Sundials_config.safe then se.checkquadvec v0;
    c_reinit s v0

  external set_err_con    : ('a, 'k) session -> bool -> unit
      = "c_idas_quad_set_err_con"

  external sv_tolerances
      : ('a, 'k) session -> float -> ('a, 'k) Nvector.t -> unit
      = "c_idas_quad_sv_tolerances"

  external ss_tolerances  : ('a, 'k) session -> float -> float -> unit
      = "c_idas_quad_ss_tolerances"

  type ('a, 'k) tolerance =
      NoStepSizeControl
    | SStolerances of float * float
    | SVtolerances of float * ('a, 'k) Nvector.t

  let set_tolerances s tol =
    let se = fwdsensext s in
    match tol with
    | NoStepSizeControl -> set_err_con s false
    | SStolerances (rel, abs) -> (ss_tolerances s rel abs;
                                  set_err_con s true)
    | SVtolerances (rel, abs) -> (if Sundials_config.safe then
                                    se.checkquadvec abs;
                                  sv_tolerances s rel abs;
                                  set_err_con s true)

  external c_get : ('a, 'k) session -> ('a, 'k) Nvector.t -> float
      = "c_idas_quad_get"

  let get s v =
    let se = fwdsensext s in
    if Sundials_config.safe then se.checkquadvec v;
    c_get s v

  external c_get_dky
      : ('a, 'k) session -> float -> int -> ('a, 'k) Nvector.t -> unit
      = "c_idas_quad_get_dky"

  let get_dky s dky =
    let se = fwdsensext s in
    if Sundials_config.safe then se.checkquadvec dky;
    fun t k -> c_get_dky s t k dky

  external get_num_rhs_evals       : ('a, 'k) session -> int
      = "c_idas_quad_get_num_rhs_evals"

  external get_num_err_test_fails  : ('a, 'k) session -> int
      = "c_idas_quad_get_num_err_test_fails"

  external c_get_err_weights : ('a, 'k) session -> ('a, 'k) Nvector.t -> unit
      = "c_idas_quad_get_err_weights"

  let get_err_weights s v =
    let se = fwdsensext s in
    if Sundials_config.safe then se.checkquadvec v;
    c_get_err_weights s v

  external get_stats : ('a, 'k) session -> int * int
      = "c_idas_quad_get_stats"

end (* }}} *)

module Sensitivity = struct (* {{{ *)
  include SensitivityTypes

  exception SensNotInitialized
  exception SensResFuncFailure
  exception FirstSensResFuncFailure
  exception RepeatedSensResFuncFailure
  exception BadSensIdentifier

  let fwdsensext s =
    match s.sensext with
    | FwdSensExt se -> se
    | _ -> raise SensNotInitialized

  type ('a, 'k) tolerance =
      SStolerances of float * Sundials.RealArray.t
    | SVtolerances of float * ('a, 'k) Nvector.t array
    | EEtolerances

  external set_err_con : ('a, 'k) session -> bool -> unit
      = "c_idas_sens_set_err_con"

  external ss_tolerances
      : ('a, 'k) session -> float -> Sundials.RealArray.t -> unit
      = "c_idas_sens_ss_tolerances"

  external ee_tolerances  : ('a, 'k) session -> unit
      = "c_idas_sens_ee_tolerances"

  external sv_tolerances
      : ('a, 'k) session -> float -> ('a, 'k) Nvector.t array -> unit
      = "c_idas_sens_sv_tolerances"

  let set_tolerances s tol =
    let ns = num_sensitivities s in
    match tol with
    | SStolerances (rel, abs) -> begin
          if Sundials_config.safe && Bigarray.Array1.dim abs <> ns
          then invalid_arg "set_tolerances: abstol has the wrong length";
          ss_tolerances s rel abs
        end
    | SVtolerances (rel, abs) -> begin
          if Sundials_config.safe then
            (if Array.length abs <> ns
             then invalid_arg "set_tolerances: abstol has the wrong length";
             Array.iter s.checkvec abs);
          sv_tolerances s rel abs
        end
    | EEtolerances -> ee_tolerances s

  type sens_method =
      Simultaneous
    | Staggered

  type sens_params = {
      pvals  : Sundials.RealArray.t option;
      pbar   : Sundials.RealArray.t option;
      plist  : int array option;
    }

  let no_sens_params = { pvals = None; pbar = None; plist = None }

  external c_sens_init : ('a, 'k) session -> sens_method -> bool
                         -> ('a, 'k) Nvector.t array
                         -> ('a, 'k) Nvector.t array -> unit
    = "c_idas_sens_init"

  external c_set_params : ('a, 'k) session -> sens_params -> unit
      = "c_idas_sens_set_params"

  let check_sens_params ns {pvals; pbar; plist} =
    if Sundials_config.safe then
      begin
        let np = match pvals with None -> 0
                                | Some p -> Bigarray.Array1.dim p in
        let check_pi v =
          if v < 0 || v >= np
          then invalid_arg "set_params: plist has an invalid entry"
        in
        if 0 <> np && np < ns then
          invalid_arg "set_params: pvals is too short";
        (match pbar with
         | None -> ()
         | Some p ->
           if Bigarray.Array1.dim p <> ns
           then invalid_arg "set_params: pbar has the wrong length");
        (match plist with
         | None -> ()
         | Some p ->
           if Array.length p <> ns
           then invalid_arg "set_params: plist has the wrong length"
           else Array.iter check_pi p)
      end

  let init s tol fmethod ?(sens_params=no_sens_params) ?fs y0 y'0 =
    if Sundials_config.safe then
      (Array.iter s.checkvec y0;
       Array.iter s.checkvec y'0);
    add_fwdsensext s;
    let se = fwdsensext s in
    let ns = Array.length y0 in
    if Sundials_config.safe then
      (if ns = 0 then
         invalid_arg "init: require at least one sensitivity parameter";
       if ns <> Array.length y'0 then
         invalid_arg "init: y0 and y'0 have inconsistent lengths");
    check_sens_params ns sens_params;
    c_sens_init s fmethod (fs <> None) y0 y'0;
    (match fs with
     | Some f -> se.sensresfn <- f
     | None -> ());
    se.num_sensitivities <- ns;
    c_set_params s sens_params;
    se.senspvals <- sens_params.pvals;
    se.sensarray1 <- c_alloc_nvector_array ns;
    se.sensarray2 <- c_alloc_nvector_array ns;
    se.sensarray3 <- c_alloc_nvector_array ns;
    set_tolerances s tol

  external c_reinit
    : ('a, 'k) session -> sens_method
      -> ('a, 'k) Nvector.t array -> ('a, 'k) Nvector.t array -> unit
    = "c_idas_sens_reinit"

  let reinit s sm s0 s'0 =
    let ns = num_sensitivities s in
    if Sundials_config.safe then
      (if Array.length s0 <> ns || Array.length s'0 <> ns
       then invalid_arg "reinit: wrong number of sensitivity vectors";
       Array.iter s.checkvec s0;
       Array.iter s.checkvec s'0);
    c_reinit s sm s0 s'0

  external toggle_off : ('a, 'k) session -> unit
    = "c_idas_sens_toggle_off"

  external c_get : ('a, 'k) session -> ('a, 'k) Nvector.t array -> float
    = "c_idas_sens_get"

  let get s ys =
    if Sundials_config.safe then
      (if Array.length ys <> num_sensitivities s
       then invalid_arg "get: wrong number of sensitivity vectors";
       Array.iter s.checkvec ys);
    c_get s ys

  external c_get_dky
    : ('a, 'k) session -> float -> int -> ('a, 'k) Nvector.t array -> unit
    = "c_idas_sens_get_dky"

  let get_dky s dkys =
    if Sundials_config.safe then
      (if Array.length dkys <> num_sensitivities s
       then invalid_arg "get_dky: wrong number of sensitivity vectors";
       Array.iter s.checkvec dkys);
    fun t k -> c_get_dky s t k dkys

  external c_get1 : ('a, 'k) session -> int -> ('a, 'k) Nvector.t -> float
    = "c_idas_sens_get1"

  let get1 s ys =
    if Sundials_config.safe then s.checkvec ys;
    fun i -> c_get1 s i ys

  external c_get_dky1
    : ('a, 'k) session -> float -> int -> int -> ('a, 'k) Nvector.t -> unit
    = "c_idas_sens_get_dky1"

  let get_dky1 s dkys =
    if Sundials_config.safe then s.checkvec dkys;
    fun t k i -> c_get_dky1 s t k i dkys

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
    num_sens_evals :int;
    num_res_evals : int;
    num_err_test_fails : int;
    num_lin_solv_setups :int;
  }

  external get_stats : ('a, 'k) session -> sensitivity_stats
    = "c_idas_sens_get_stats"

  external c_get_err_weights
    : ('a, 'k) session -> ('a, 'k) Nvector.t array -> unit
    = "c_idas_sens_get_err_weights"

  let get_err_weights s esweight =
    if Sundials_config.safe then
      (if Array.length esweight <> num_sensitivities s
       then invalid_arg "get_err_weights: wrong number of vectors";
       Array.iter s.checkvec esweight);
    c_get_err_weights s esweight

  external c_sens_calc_ic_ya_yd' :
       ('a,'k) session
    -> ('a,'k) Nvector.t option
    -> ('a,'k) Nvector.t option
    -> ('a,'k) Nvector.t array option
    -> ('a,'k) Nvector.t array option
    -> float
    -> unit
    = "c_ida_sens_calc_ic_ya_ydp_byte"
      "c_ida_sens_calc_ic_ya_ydp"

  external c_sens_calc_ic_y :
    ('a,'k) session
    -> ('a,'k) Nvector.t option
    -> ('a,'k) Nvector.t array option
    -> float -> unit
    = "c_ida_sens_calc_ic_y"

  let calc_ic_ya_yd' session ?y ?y' ?s ?s' ?varid tout1 =
    let num_sens = num_sensitivities session in
    if Sundials_config.safe then
      (ocheck session.checkvec y;
       ocheck session.checkvec y');
    (match s with
     | Some s ->
         if Sundials_config.safe then
           (if Array.length s <> num_sens
            then invalid_arg "calc_ic_ya_yd': wrong number of vectors in ~s";
            Array.iter session.checkvec s)
     | _ -> ());
    (match s' with
     | Some s' ->
         if Sundials_config.safe then
           (if Array.length s' <> num_sens
            then invalid_arg
                 "calc_ic_ya_yd': wrong number of vectors in ~s'";
            Array.iter session.checkvec s')
     | _ -> ());
    (match varid with
     | None -> if not session.id_set then raise Ida.IdNotSet
     | Some x -> Ida.set_id session x);
    c_sens_calc_ic_ya_yd' session y y' s s' tout1

  (* Note: my understanding is that CalcIC with IDA_Y_INIT corrects
     the non-derivatives of the sensitivity variables while holding
     the derivatives constant, so there's no point querying the
     values of the corrected derivatives.  *)
  let calc_ic_y session ?y ?s tout1 =
    let num_sens = num_sensitivities session in
    if Sundials_config.safe then
      (ocheck session.checkvec y;
       match s with
       | Some s when Array.length s <> num_sens ->
         invalid_arg "calc_ic_y: wrong number of vectors in ~s"
       | _ -> ());
    c_sens_calc_ic_y session y s tout1

  external get_num_nonlin_solv_iters : ('a, 'k) session -> int
    = "c_idas_sens_get_num_nonlin_solv_iters"

  external get_num_nonlin_solv_conv_fails : ('a, 'k) session -> int
    = "c_idas_sens_get_num_nonlin_solv_conv_fails"

  external get_nonlin_solv_stats : ('a, 'k) session -> int * int
    = "c_idas_sens_get_nonlin_solv_stats"

  module Quadrature = struct (* {{{ *)
    include QuadratureTypes

    exception QuadSensNotInitialized
    exception QuadSensRhsFuncFailure
    exception FirstQuadSensRhsFuncFailure
    exception RepeatedQuadSensRhsFuncFailure

    external c_quadsens_init
      : ('a, 'k) session -> bool -> ('a, 'k) Nvector.t array -> unit
      = "c_idas_quadsens_init"

    let init s ?fqs v0 =
      let se = fwdsensext s in
      if not se.has_quad then raise Quadrature.QuadNotInitialized;
      if Sundials_config.safe && Array.length v0 <> se.num_sensitivities
      then invalid_arg "init: wrong number of vectors";
      if Sundials_config.safe then Array.iter se.checkquadvec v0;
      match fqs with
      | Some f -> se.quadsensrhsfn <- f;
                  c_quadsens_init s true v0
      | None -> c_quadsens_init s false v0

    external c_reinit : ('a, 'k) session -> ('a, 'k) Nvector.t array -> unit
      = "c_idas_quadsens_reinit"

    let reinit s v =
      let se = fwdsensext s in
      if Sundials_config.safe then
        (if Array.length v <> se.num_sensitivities
         then invalid_arg "reinit: wrong number of vectors";
         Array.iter se.checkquadvec v);
      c_reinit s v

    type ('a, 'k) tolerance =
        NoStepSizeControl
      | SStolerances of float * Sundials.RealArray.t
      | SVtolerances of float * ('a, 'k) Nvector.t array
      | EEtolerances

    external set_err_con : ('a, 'k) session -> bool -> unit
      = "c_idas_quadsens_set_err_con"

    external ss_tolerances
      : ('a, 'k) session -> float -> Sundials.RealArray.t -> unit
      = "c_idas_quadsens_ss_tolerances"

    external sv_tolerances
      : ('a, 'k) session -> float -> ('a, 'k) Nvector.t array -> unit
      = "c_idas_quadsens_sv_tolerances"

    external ee_tolerances  : ('a, 'k) session -> unit
      = "c_idas_quadsens_ee_tolerances"

    let set_tolerances s tol =
      let se = fwdsensext s in
      match tol with
      | NoStepSizeControl -> set_err_con s false
      | SStolerances (rel, abs) -> begin
          if Sundials_config.safe &&
             Bigarray.Array1.dim abs <> se.num_sensitivities
          then invalid_arg "set_tolerances: abstol has the wrong length";
          ss_tolerances s rel abs;
          set_err_con s true
        end
      | SVtolerances (rel, abs) -> begin
          if Sundials_config.safe then
            (if Array.length abs <> se.num_sensitivities
             then invalid_arg "set_tolerances: abstol has the wrong length";
             Array.iter se.checkquadvec abs);
          sv_tolerances s rel abs;
          set_err_con s true
        end
      | EEtolerances -> (ee_tolerances s;
                         set_err_con s true)

    external c_get : ('a, 'k) session -> ('a, 'k) Nvector.t array -> float
      = "c_idas_quadsens_get"

    let get s ys =
      let se = fwdsensext s in
      if Sundials_config.safe then
        (if Array.length ys <> se.num_sensitivities
         then invalid_arg "get: wrong number of vectors";
         Array.iter se.checkquadvec ys);
      c_get s ys

    external c_get1 : ('a, 'k) session -> int -> ('a, 'k) Nvector.t -> float
      = "c_idas_quadsens_get1"

    let get1 s yqs =
      let se = fwdsensext s in
      if Sundials_config.safe then se.checkquadvec yqs;
      fun i -> c_get1 s i yqs

    external c_get_dky
      : ('a, 'k) session -> float -> int -> ('a, 'k) Nvector.t array -> unit
      = "c_idas_quadsens_get_dky"

    let get_dky s ys =
      let se = fwdsensext s in
      if Sundials_config.safe then
        (if Array.length ys <> se.num_sensitivities
         then invalid_arg "get_dky: wrong number of vectors";
         Array.iter se.checkquadvec ys);
      fun t k -> c_get_dky s t k ys

    external c_get_dky1 : ('a, 'k) session -> float -> int -> int
      -> ('a, 'k) Nvector.t -> unit
      = "c_idas_quadsens_get_dky1"

    let get_dky1 s dkyqs =
      let se = fwdsensext s in
      if Sundials_config.safe then se.checkquadvec dkyqs;
      fun t k i -> c_get_dky1 s t k i dkyqs

    external get_num_rhs_evals       : ('a, 'k) session -> int
      = "c_idas_quadsens_get_num_rhs_evals"

    external get_num_err_test_fails  : ('a, 'k) session -> int
      = "c_idas_quadsens_get_num_err_test_fails"

    external c_get_err_weights
      : ('a, 'k) session -> ('a, 'k) Nvector.t array -> unit
      = "c_idas_quadsens_get_err_weights"

    let get_err_weights s esweight =
      let se = fwdsensext s in
      if Sundials_config.safe then
        (if Array.length esweight <> se.num_sensitivities
         then invalid_arg "get_err_weights: wrong number of vectors";
         Array.iter se.checkquadvec esweight);
      c_get_err_weights s esweight

    external get_stats : ('a, 'k) session -> int * int
      = "c_idas_quadsens_get_stats"
  end (* }}} *)
end (* }}} *)

module Adjoint = struct (* {{{ *)
  include AdjointTypes

  exception AdjointNotInitialized
  exception NoForwardCall
  exception ForwardReinitFailure
  exception ForwardFailure
  exception NoBackwardProblem
  exception BadFinalTime
  exception BadOutputTime

  type interpolation = IPolynomial | IHermite

  external c_init : ('a, 'k) session -> int -> interpolation -> unit
      = "c_idas_adj_init"

  let init s nd interptype =
    add_fwdsensext s;
    c_init s nd interptype

  let fwdsensext s =
    match s.sensext with
    | FwdSensExt se -> se
    | _ -> raise AdjointNotInitialized

  external c_set_id
    : ('a,'k) session -> int -> ('a,'k) Nvector.t -> unit
    = "c_idas_adj_set_id"

  let set_id b ids =
    let bs = tosession b in
    if Sundials_config.safe then bs.checkvec ids;
    let parent, which = parent_and_which b in
    c_set_id parent which ids;
    bs.id_set <- true

  external c_set_suppress_alg : ('a,'k) session -> int -> bool -> unit
    = "c_idas_adj_set_suppress_alg"

  let set_suppress_alg b ?varid v =
    (match varid with
     | None -> if v && not (tosession b).id_set then raise Ida.IdNotSet
     | Some x -> set_id b x);
    let parent, which = parent_and_which b in
    c_set_suppress_alg parent which v

  external c_adj_calc_ic :
    ('a,'k) session
    -> int
    -> float
    -> ('a,'k) Nvector.t
    -> ('a,'k) Nvector.t
    -> unit
    = "c_idas_adj_calc_ic"

  external c_adj_calc_ic_sens :
    ('a,'k) session
    -> int
    -> float
    -> ('a,'k) Nvector.t
    -> ('a,'k) Nvector.t
    -> ('a,'k) Nvector.t array
    -> ('a,'k) Nvector.t array
    -> unit
    = "c_idas_adj_calc_ic_sens_byte"
      "c_idas_adj_calc_ic_sens"

  external c_adj_get_consistent_ic :
    ('a,'k) session -> int
    -> ('a,'k) Nvector.t option
    -> ('a,'k) Nvector.t option
    -> unit
    = "c_idas_adj_get_consistent_ic"

  let calc_ic bsession ?yb ?yb' tout1 y0 y0' =
    let checkvec = (tosession bsession).checkvec in
    let parent, which = parent_and_which bsession in
    if Sundials_config.safe then
      (checkvec y0;
       checkvec y0';
       ocheck checkvec yb;
       ocheck checkvec yb');
    c_adj_calc_ic parent which tout1 y0 y0';
    c_adj_get_consistent_ic parent which yb yb'

  let calc_ic_sens bsession ?yb ?yb' ?varid tout1 y0 y0' ys0 ys0' =
    let bs = tosession bsession in
    let num_sens = num_sensitivities bs in
    if Sundials_config.safe then
      (if Array.length ys0 <> num_sens then
         invalid_arg "calc_ic_sens: wrong number of vectors in ys0";
       if Array.length ys0' <> num_sens then
         invalid_arg "calc_ic_sens: wrong number of vectors in y's0";
       bs.checkvec y0;
       bs.checkvec y0';
       Array.iter bs.checkvec ys0;
       Array.iter bs.checkvec ys0';
       ocheck bs.checkvec yb;
       ocheck bs.checkvec yb');
    let parent, which = parent_and_which bsession in
    (match varid with
     | None -> if not bs.id_set then raise Ida.IdNotSet
     | Some x -> set_id bsession x);
    c_adj_calc_ic_sens parent which tout1 y0 y0' ys0 ys0';
    c_adj_get_consistent_ic parent which yb yb'

  external c_forward_normal : ('a, 'k) session -> float
                            -> ('a, 'k) Nvector.t -> ('a, 'k) Nvector.t
                            -> float * int * Ida.solver_result
      = "c_idas_adj_forward_normal"

  let forward_normal s t y y' =
    if Sundials_config.safe then
      (s.checkvec y;
       s.checkvec y');
    c_forward_normal s t y y'

  external c_forward_one_step : ('a, 'k) session -> float
                              -> ('a, 'k) Nvector.t -> ('a, 'k) Nvector.t
                              -> float * int * Ida.solver_result
      = "c_idas_adj_forward_one_step"

  let forward_one_step s t y y' =
    if Sundials_config.safe then
      (s.checkvec y;
       s.checkvec y');
    c_forward_one_step s t y y'

  type 'a triple = 'a * 'a * 'a

  type ('data, 'kind) iter =
    | Newton of ('data, 'kind) linear_solver
    | Functional

  type ('a, 'k) tolerance =
    | SStolerances of float * float
    | SVtolerances of float * ('a, 'k) Nvector.t

  external ss_tolerances
      : ('a, 'k) session -> int -> float -> float -> unit
      = "c_idas_adj_ss_tolerances"

  external sv_tolerances
      : ('a, 'k) session -> int -> float -> ('a, 'k) Nvector.t -> unit
      = "c_idas_adj_sv_tolerances"

  let set_tolerances bs tol =
    let parent, which = parent_and_which bs in
    match tol with
    | SStolerances (rel, abs) -> ss_tolerances parent which rel abs
    | SVtolerances (rel, abs) -> (if Sundials_config.safe then
                                    (tosession bs).checkvec abs;
                                  sv_tolerances parent which rel abs)

  let bwdsensext = function (Bsession bs) ->
    match bs.sensext with
    | BwdSensExt se -> se
    | _ -> raise AdjointNotInitialized

  let set_linear_solver bs solver nv =
    (tosession bs).ls_callbacks <- NoCallbacks;
    (tosession bs).ls_precfns <- NoPrecFns;
    solver bs nv

  external backward_normal : ('a, 'k) session -> float -> unit
      = "c_idas_adj_backward_normal"

  external backward_one_step : ('a, 'k) session -> float -> unit
      = "c_idas_adj_backward_one_step"

  external c_get : ('a, 'k) session -> int
                   -> ('a, 'k) Nvector.t -> ('a, 'k) Nvector.t -> float
      = "c_idas_adj_get"

  let get bs yb ypb =
    if Sundials_config.safe then
      (let checkvec = (tosession bs).checkvec in
       checkvec yb;
       checkvec ypb);
    let parent, which = parent_and_which bs in
    c_get parent which yb ypb

  let get_dky bs = Ida.get_dky (tosession bs)

  external c_get_y : ('d, 'k) session -> float -> ('d, 'k) Nvector.t
                      -> ('d, 'k) Nvector.t -> unit
      = "c_idas_adj_get_y"

  let get_y s y yp =
    if Sundials_config.safe then (s.checkvec y; s.checkvec yp);
    fun t -> c_get_y s t y yp

  external set_no_sensitivity : ('a, 'k) session -> unit
      = "c_idas_adj_set_no_sensi"

  external c_set_max_ord : ('a, 'k) session -> int -> int -> unit
      = "c_idas_adj_set_max_ord"

  let set_max_ord bs maxordb =
    let parent, which = parent_and_which bs in
    c_set_max_ord parent which maxordb

  external c_set_max_num_steps : ('a, 'k) session -> int -> int -> unit
      = "c_idas_adj_set_max_num_steps"

  let set_max_num_steps bs mxstepsb =
    let parent, which = parent_and_which bs in
    c_set_max_num_steps parent which mxstepsb

  external c_set_init_step : ('a, 'k) session -> int -> float -> unit
      = "c_idas_adj_set_init_step"

  let set_init_step bs hinb =
    let parent, which = parent_and_which bs in
    c_set_init_step parent which hinb

  external c_set_max_step : ('a, 'k) session -> int -> float -> unit
      = "c_idas_adj_set_max_step"

  let set_max_step bs hmaxb =
    let parent, which = parent_and_which bs in
    c_set_max_step parent which hmaxb

  module Dls = struct (* {{{ *)
    include DirectTypes

    (* Sundials < 3.0.0 *)
    external c_dls_dense
      : 'k serial_session -> int -> int -> bool -> bool -> unit
      = "c_idas_adj_dls_dense"

    (* Sundials < 3.0.0 *)
    external c_dls_lapack_dense
      : 'k serial_session -> int -> int -> bool -> bool -> unit
      = "c_idas_adj_dls_lapack_dense"

    (* Sundials < 3.0.0 *)
    external c_dls_band : ('k serial_session * int) -> (int * int * int)
                            -> bool -> bool -> unit
      = "c_idas_adj_dls_band"

    (* Sundials < 3.0.0 *)
    external c_dls_lapack_band : ('k serial_session * int) -> (int * int * int)
                                  -> bool -> bool -> unit
      = "c_idas_adj_dls_lapack_band"

    (* Sundials < 3.0.0 *)
    external c_klub
      : 'k serial_session * int -> 's Matrix.Sparse.sformat
        -> int -> int -> bool -> unit
      = "c_idas_klub_init"

    (* Sundials < 3.0.0 *)
    external c_klu_set_ordering
      : 'k serial_session -> Lsolver.Direct.Klu.ordering -> unit
      = "c_ida_klu_set_ordering"

    (* Sundials < 3.0.0 *)
    external c_klu_reinit : 'k serial_session -> int -> int -> unit
      = "c_ida_klu_reinit"

    (* Sundials < 3.0.0 *)
    let klu_set_ordering session ordering =
      match session.ls_callbacks with
      | SlsKluCallback _ | BSlsKluCallback _ | BSlsKluCallbackSens _ ->
          c_klu_set_ordering session ordering
      | _ -> ()

    (* Sundials < 3.0.0 *)
    let klu_reinit session n onnz =
      match session.ls_callbacks with
      | SlsKluCallback _ | BSlsKluCallback _ | BSlsKluCallbackSens _ ->
          c_klu_reinit session n (match onnz with None -> 0 | Some nnz -> nnz)
      | _ -> ()

    (* Sundials < 3.0.0 *)
    external c_superlumtb : ('k serial_session * int)
                            -> int -> int -> int -> bool -> unit
      = "c_idas_superlumtb_init"

    (* Sundials < 3.0.0 *)
    external c_superlumt_set_ordering
      : 'k serial_session -> Lsolver.Direct.Superlumt.ordering -> unit
      = "c_ida_superlumt_set_ordering"

    (* Sundials < 3.0.0 *)
    let superlumt_set_ordering session ordering =
      match session.ls_callbacks with
      | SlsSuperlumtCallback _ | BSlsSuperlumtCallback _
      | BSlsSuperlumtCallbackSens _ ->
          c_superlumt_set_ordering session ordering
      | _ -> ()

    module LSD = Lsolver_impl.Direct

    (* Sundials < 3.0.0 *)
    let make_compat (type s) (type tag) hasjac usesens
        (solver : (s, 'nd, 'nk, tag) LSD.solver)
        (mat : ('k, s, 'nd, 'nk) Matrix.t) bs =
      let parent, which = parent_and_which bs in
      match solver with
      | LSD.Dense ->
          let m, n = Matrix.(Dense.size (unwrap mat)) in
          if m <> n then raise Lsolver.MatrixNotSquare;
          c_dls_dense parent which m hasjac usesens
      | LSD.LapackDense ->
          let m, n = Matrix.(Dense.size (unwrap mat)) in
          if m <> n then raise Lsolver.MatrixNotSquare;
          c_dls_lapack_dense parent which m hasjac usesens

      | LSD.Band ->
          let open Matrix.Band in
          let { n; mu; ml } = dims (Matrix.unwrap mat) in
          c_dls_band (parent, which) (n, mu, ml) hasjac usesens
      | LSD.LapackBand ->
          let open Matrix.Band in
          let { n; mu; ml } = dims (Matrix.unwrap mat) in
          c_dls_lapack_band (parent, which) (n, mu, ml) hasjac usesens

      | LSD.Klu sinfo ->
          if not Sundials_config.klu_enabled
            then raise Sundials.NotImplementedBySundialsVersion;
          let smat = Matrix.unwrap mat in
          let m, n = Matrix.Sparse.size smat in
          let nnz, _ = Matrix.Sparse.dims smat in
          if m <> n then raise Lsolver.MatrixNotSquare;
          let open Lsolver_impl.Klu in
          let session = tosession bs in
          sinfo.set_ordering <- klu_set_ordering session;
          sinfo.reinit <- klu_reinit session;
          c_klub (parent, which) (Matrix.Sparse.sformat smat) m nnz usesens;
          (match sinfo.ordering with None -> ()
                                   | Some o -> c_klu_set_ordering session o)

      | LSD.Superlumt sinfo ->
          if not Sundials_config.superlumt_enabled
            then raise Sundials.NotImplementedBySundialsVersion;
          let smat = Matrix.unwrap mat in
          let m, n = Matrix.Sparse.size smat in
          let nnz, _ = Matrix.Sparse.dims smat in
          if m <> n then raise Lsolver.MatrixNotSquare;
          let open Lsolver_impl.Superlumt in
          let session = tosession bs in
          sinfo.set_ordering <- superlumt_set_ordering session;
          c_superlumtb (parent, which) m nnz sinfo.num_threads usesens;
          (match sinfo.ordering with None -> ()
                                   | Some o -> c_superlumt_set_ordering session o)

    | LSD.Custom _ ->
        assert false

    let check_dqjac jac mat = let open Matrix in
      match get_id mat with
      | Dense | Band -> ()
      | _ -> if jac = None then invalid_arg "A Jacobian function is required"

    let set_ls_callbacks (type m) (type tag)
          ?(jac : m jac_fn option)
          (solver : (m, 'nd, 'nk, tag) LSD.solver)
          (mat : ('mk, m, 'nd, 'nk) Matrix.t) session =
      let none = (None : m option) in
      begin match solver with
      | LSD.Dense ->
          session.ls_callbacks <- (match jac with
            | None ->
                BDlsDenseCallback { jacfn = no_callback; jmat = none }
            | Some (NoSens f) ->
                BDlsDenseCallback { jacfn = f; jmat = none }
            | Some (WithSens f) ->
                BDlsDenseCallbackSens { jacfn_sens = f; jmat = none })
      | LSD.LapackDense ->
          session.ls_callbacks <- (match jac with
            | None ->
                BDlsDenseCallback { jacfn = no_callback; jmat = none }
            | Some (NoSens f) ->
                BDlsDenseCallback { jacfn = f; jmat = none }
            | Some (WithSens f) ->
                BDlsDenseCallbackSens { jacfn_sens = f; jmat = none })
      | LSD.Band ->
          session.ls_callbacks <- (match jac with
            | None ->
                BDlsBandCallback { jacfn = no_callback; jmat = none }
            | Some (NoSens f) ->
                BDlsBandCallback { jacfn = f; jmat = none }
            | Some (WithSens f) ->
                BDlsBandCallbackSens { jacfn_sens = f; jmat = none })
      | LSD.LapackBand ->
          session.ls_callbacks <- (match jac with
            | None ->
                BDlsBandCallback { jacfn = no_callback; jmat = none }
            | Some (NoSens f) ->
                BDlsBandCallback { jacfn = f; jmat = none }
            | Some (WithSens f) ->
                BDlsBandCallbackSens { jacfn_sens = f; jmat = none })
      | LSD.Klu _ ->
          session.ls_callbacks <- (match jac with
            | None -> invalid_arg "Klu requires Jacobian function";
            | Some (NoSens f) ->
                BSlsKluCallback { jacfn = f; jmat = none }
            | Some (WithSens f) ->
                BSlsKluCallbackSens { jacfn_sens = f; jmat = none })
      | LSD.Superlumt _ ->
          session.ls_callbacks <- (match jac with
            | None -> invalid_arg "Superlumt requires Jacobian function";
            | Some (NoSens f) ->
                BSlsSuperlumtCallback { jacfn = f; jmat = none }
            | Some (WithSens f) ->
                BSlsSuperlumtCallbackSens { jacfn_sens = f; jmat = none })
      | LSD.Custom _ ->
          session.ls_callbacks <- (match jac with
            | None ->
                (match Matrix.get_id mat with
                 | Matrix.Dense | Matrix.Band -> ()
                 | _ -> invalid_arg "A Jacobian function is required");
                BDirectCustomCallback { jacfn = no_callback; jmat = none }
            | Some (NoSens f) ->
                BDirectCustomCallback { jacfn = f; jmat = none }
            | Some (WithSens f) ->
                BDirectCustomCallbackSens { jacfn_sens = f; jmat = none })
      end;
      session.ls_precfns <- NoPrecFns

    (* Sundials >= 3.0.0 *)
    external c_dls_set_linear_solver
      : 'k serial_session * int
         -> ('m, Nvector_serial.data, 'k) LSD.cptr
         -> ('mk, 'm, Nvector_serial.data, 'k) Matrix.t
         -> bool
         -> bool
         -> unit
      = "c_idas_adj_dls_set_linear_solver"

    let solver ?jac ((LSD.S { LSD.rawptr; LSD.solver; LSD.matrix }) as ls)
               bs nv =
      let session = tosession bs in
      let parent, which = parent_and_which bs in
      let use_sens = match jac with Some (WithSens _) -> true | _ -> false in
      set_ls_callbacks ?jac solver matrix session;
      if in_compat_mode then make_compat (jac <> None) use_sens solver matrix bs
      else c_dls_set_linear_solver (parent, which) rawptr matrix
                                                      (jac <> None) use_sens;
      LSD.attach ls;
      session.ls_solver <- Lsolver_impl.DirectSolver ls

    (* Sundials < 3.0.0 *)
    let invalidate_callback s =
      if in_compat_mode then
        match s.ls_callbacks with
        | BDlsDenseCallback ({ jmat = Some d } as cb) ->
            Matrix.Dense.invalidate d;
            cb.jmat <- None
        | BDlsDenseCallbackSens ({ jmat = Some d } as cb) ->
            Matrix.Dense.invalidate d;
            cb.jmat <- None
        | BDlsBandCallback  ({ jmat = Some d } as cb) ->
            Matrix.Band.invalidate d;
            cb.jmat <- None
        | BDlsBandCallbackSens  ({ jmat = Some d } as cb) ->
            Matrix.Band.invalidate d;
            cb.jmat <- None
        | BSlsKluCallback ({ jmat = Some d } as cb) ->
            Matrix.Sparse.invalidate d;
            cb.jmat <- None
        | BSlsKluCallbackSens ({ jmat = Some d } as cb) ->
            Matrix.Sparse.invalidate d;
            cb.jmat <- None
        | BSlsSuperlumtCallback ({ jmat = Some d } as cb) ->
            Matrix.Sparse.invalidate d;
            cb.jmat <- None
        | BSlsSuperlumtCallbackSens ({ jmat = Some d } as cb) ->
            Matrix.Sparse.invalidate d;
            cb.jmat <- None
        | _ -> ()

    let get_work_space bs = Ida.Dls.get_work_space (tosession bs)
    let get_num_jac_evals bs = Ida.Dls.get_num_jac_evals (tosession bs)
    let get_num_res_evals bs = Ida.Dls.get_num_res_evals (tosession bs)
  end (* }}} *)

  module Spils = struct (* {{{ *)
    include SpilsTypes

    (* Sundials < 3.0.0 *)
    external c_spgmr
      : ('a, 'k) session -> int -> int -> unit
      = "c_idas_adj_spils_spgmr"

    (* Sundials < 3.0.0 *)
    external c_spbcgs
      : ('a, 'k) session -> int -> int -> unit
      = "c_idas_adj_spils_spbcgs"

    (* Sundials < 3.0.0 *)
    external c_sptfqmr
      : ('a, 'k) session -> int -> int -> unit
      = "c_idas_adj_spils_sptfqmr"

    (* Sundials < 3.0.0 *)
    external c_set_gs_type
        : ('a, 'k) session -> int -> Lsolver.Iterative.gramschmidt_type -> unit
        = "c_idas_adj_spils_set_gs_type"

    (* Sundials < 3.0.0 *)
    external c_set_maxl : ('a, 'k) session -> int -> int -> unit
        = "c_idas_adj_spils_set_maxl"

    let old_set_maxl bs maxl =
      ls_check_spils (tosession bs);
      let parent, which = parent_and_which bs in
      c_set_maxl parent which maxl

    let old_set_gs_type bs t =
      ls_check_spils (tosession bs);
      let parent, which = parent_and_which bs in
      c_set_gs_type parent which t

    external c_set_jac_times
      : ('a, 'k) session -> int -> bool -> bool -> bool -> unit
      = "c_idas_adj_spils_set_jac_times"

    external c_set_preconditioner
      : ('a, 'k) session -> int -> bool -> bool -> unit
      = "c_idas_adj_spils_set_preconditioner"

    let init_preconditioner solve setup bs parent which nv =
      c_set_preconditioner parent which (setup <> None) false;
      (tosession bs).ls_precfns <- BPrecFns { prec_solve_fn = solve;
                                              prec_setup_fn = setup }

    let prec_none = Lsolver_impl.Iterative.(PrecNone,
                      fun bs _ _ _ -> (tosession bs).ls_precfns <- NoPrecFns)
    let prec_left ?setup solve  = Lsolver_impl.Iterative.(PrecLeft,
                                    init_preconditioner solve setup)

    let check_prec_type prec_type =
      let open Lsolver_impl.Iterative in
      match prec_type with
      | PrecNone | PrecLeft -> true
      | PrecRight | PrecBoth -> false

    let init_preconditioner_with_sens solve setup bs parent which nv =
      c_set_preconditioner parent which (setup <> None) true;
      (tosession bs).ls_precfns <- BPrecFnsSens
            { prec_solve_fn_sens = solve; prec_setup_fn_sens = setup }

    let prec_left_with_sens ?setup solve  = Lsolver_impl.Iterative.(PrecLeft,
                                      init_preconditioner_with_sens solve setup)

    type 'd jac_times_vec_fn =
      | NoSens of 'd jac_times_setup_fn_no_sens option
                  * 'd jac_times_vec_fn_no_sens
      | WithSens of 'd jac_times_setup_fn_with_sens option
                    * 'd jac_times_vec_fn_with_sens

    external c_spils_set_linear_solver
      : ('a, 'k) session -> int -> ('a, 'k) Lsolver_impl.Iterative.cptr -> unit
      = "c_idas_adj_spils_set_linear_solver"

    let solver (type s)
          ({ Lsolver_impl.Iterative.rawptr;
             Lsolver_impl.Iterative.solver;
             Lsolver_impl.Iterative.compat =
               ({ Lsolver_impl.Iterative.maxl;
                  Lsolver_impl.Iterative.gs_type } as compat) } as lsolver)
          ?jac_times_vec (prec_type, set_prec) bs nv =
      let session = tosession bs in
      let parent, which = parent_and_which bs in
      if in_compat_mode then begin
        match jac_times_vec with
        | Some (NoSens (Some _, _)) | Some (WithSens (Some _, _)) ->
            raise Sundials.NotImplementedBySundialsVersion;
        | _ -> ();
        let open Lsolver_impl.Iterative in
        lsolver.check_prec_type <- check_prec_type;
        (match (solver : ('nd, 'nk, s) solver) with
         | Spgmr ->
             c_spgmr parent which maxl;
             (match gs_type with None -> () | Some t ->
                 c_set_gs_type parent which t);
             compat.set_gs_type <- old_set_gs_type bs
         | Spbcgs ->
             c_spbcgs parent which maxl;
             compat.set_maxl <- old_set_maxl bs
         | Sptfqmr ->
             c_sptfqmr parent which maxl;
             compat.set_maxl <- old_set_maxl bs
         | _ -> raise Sundials.NotImplementedBySundialsVersion);
        session.ls_solver <- Lsolver_impl.IterativeSolver lsolver;
        set_prec bs parent which nv;
        (match jac_times_vec with
         | Some (NoSens (ojs, jt)) ->
             c_set_jac_times parent which (ojs <> None) true false;
             session.ls_callbacks <- BSpilsCallback (Some jt, ojs)
         | Some (WithSens (ojs, jt)) ->
             c_set_jac_times parent which (ojs <> None) true true;
             session.ls_callbacks <- BSpilsCallbackSens (Some jt, ojs)
         | None ->
             session.ls_callbacks <- BSpilsCallbackSens (None, None))
      end else
        c_spils_set_linear_solver parent which rawptr;
        Lsolver_impl.Iterative.attach lsolver;
        session.ls_solver <- Lsolver_impl.IterativeSolver lsolver;
        Lsolver_impl.Iterative.(c_set_prec_type rawptr solver prec_type false);
        set_prec bs parent which nv;
        let has_setup, has_times, use_sens =
          match jac_times_vec with
          | Some (NoSens (ojs, jt)) ->
              session.ls_callbacks <- BSpilsCallback (Some jt, ojs);
              ojs <> None, true, false
          | Some (WithSens (ojs, jt)) ->
              session.ls_callbacks <- BSpilsCallbackSens (Some jt, ojs);
              ojs <> None, true, true
          | None ->
              session.ls_callbacks <- BSpilsCallbackSens (None, None);
              false, false, false
        in
        c_set_jac_times parent which has_setup has_times use_sens

    let set_preconditioner bs ?setup solve =
      match (tosession bs).ls_callbacks with
      | BSpilsCallback _ | BSpilsCallbackSens _ ->
          let parent, which = parent_and_which bs in
          c_set_preconditioner parent which (setup <> None) false;
          (tosession bs).ls_precfns
            <- BPrecFns { prec_setup_fn = setup; prec_solve_fn = solve }
      | _ -> raise Sundials.InvalidLinearSolver

    let set_preconditioner_with_sens bs ?setup solve =
      match (tosession bs).ls_callbacks with
      | BSpilsCallback _ | BSpilsCallbackSens _ ->
          let parent, which = parent_and_which bs in
          c_set_preconditioner parent which (setup <> None) true;
          (tosession bs).ls_precfns
            <- BPrecFnsSens { prec_setup_fn_sens = setup;
                              prec_solve_fn_sens = solve }
      | _ -> raise Sundials.InvalidLinearSolver

    let set_jac_times bs jtv =
      match (tosession bs).ls_callbacks with
      | BSpilsCallback _ | BSpilsCallbackSens _ ->
          let parent, which = parent_and_which bs in
          (match jtv with
           | NoSens (ojs, jt) ->
             if in_compat_mode && ojs <> None then
               raise Sundials.NotImplementedBySundialsVersion;
             c_set_jac_times parent which (ojs <> None) true false;
             (tosession bs).ls_callbacks <- BSpilsCallback (Some jt, ojs)
           | WithSens (ojs, jt) ->
             if in_compat_mode && ojs <> None then
               raise Sundials.NotImplementedBySundialsVersion;
             c_set_jac_times parent which (ojs <> None) true true;
             (tosession bs).ls_callbacks <- BSpilsCallbackSens (Some jt, ojs))
      | _ -> raise Sundials.InvalidLinearSolver

    let clear_jac_times bs =
      let parent, which = parent_and_which bs in
      match (tosession bs).ls_callbacks with
      | BSpilsCallback _ ->
          c_set_jac_times parent which false false false;
          (tosession bs).ls_callbacks <- BSpilsCallback (None, None)
      | BSpilsCallbackSens _ ->
          c_set_jac_times parent which false false true;
          (tosession bs).ls_callbacks <- BSpilsCallbackSens (None, None)
      | _ -> raise Sundials.InvalidLinearSolver

    external set_eps_lin : ('a, 'k) bsession -> float -> unit
        = "c_idas_adj_spils_set_eps_lin"

    let set_eps_lin bs epsl =
      ls_check_spils (tosession bs);
      set_eps_lin bs epsl

    let get_work_space bs =
      Ida.Spils.get_work_space (tosession bs)

    let get_num_lin_iters bs =
      Ida.Spils.get_num_lin_iters (tosession bs)

    let get_num_conv_fails bs =
      Ida.Spils.get_num_conv_fails (tosession bs)

    let get_num_prec_evals bs =
      Ida.Spils.get_num_prec_evals (tosession bs)

    let get_num_prec_solves bs =
      Ida.Spils.get_num_prec_solves (tosession bs)

    let get_num_jtsetup_evals bs =
      Ida.Spils.get_num_jtsetup_evals (tosession bs)

    let get_num_jtimes_evals bs =
      Ida.Spils.get_num_jtimes_evals (tosession bs)

    let get_num_res_evals bs =
      Ida.Spils.get_num_res_evals (tosession bs)
  end (* }}} *)

  module Alternate = struct (* {{{ *)

    type ('data, 'kind) linit = ('data, 'kind) bsession -> unit

    type ('data, 'kind) lsetup =
      ('data, 'kind) bsession
      -> 'data Ida.Alternate.lsetup_args
      -> unit

    type ('data, 'kind) lsolve =
      ('data, 'kind) bsession
      -> 'data Ida.Alternate.lsolve_args
      -> 'data
      -> unit

    type ('data, 'kind) callbacks =
      {
        linit  : ('data, 'kind) linit option;
        lsetup : ('data, 'kind) lsetup option;
        lsolve : ('data, 'kind) lsolve;
      }

    let wrapo fo =
      match fo with
      | None -> None
      | Some f -> Some (fun s -> f (Bsession s))

    let wrap_callbacks { linit; lsetup; lsolve } =
      Ida.Alternate.({
          linit  = wrapo linit;
          lsetup = wrapo lsetup;
          lsolve = (fun s -> lsolve (Bsession s));
        })

    external c_set_alternate
      : ('data, 'kind) session -> int -> bool -> bool -> unit
      = "c_ida_adj_set_alternate"

    let get_cj bs = Ida.Alternate.get_cj (tosession bs)
    let get_cjratio bs = Ida.Alternate.get_cjratio (tosession bs)

    let solver f bs nv =
      let { linit; lsetup; lsolve } as cb = f bs nv in
      let s = tosession bs in
      let parent, which = parent_and_which bs in
      c_set_alternate parent which (linit <> None) (lsetup <> None);
      s.ls_precfns <- NoPrecFns;
      s.ls_callbacks <- AlternateCallback (wrap_callbacks cb)

  end (* }}} *)
  external c_bsession_finalize : ('a, 'k) session -> unit
      = "c_idas_adj_bsession_finalize"

  let bsession_finalize s =
    Dls.invalidate_callback s;
    c_bsession_finalize s

  external c_init_backward
      : ('a, 'k) session -> ('a, 'k) session Weak.t
        -> float
        -> ('a, 'k) Nvector.t
        -> ('a, 'k) Nvector.t
        -> bool
        -> (ida_mem * int * c_weak_ref)
      = "c_idas_adj_init_backward_byte"
        "c_idas_adj_init_backward"

  let init_backward s linsolv tol mf ?varid t0 y0 y'0 =
    let { bsessions } as se = fwdsensext s in
    let ns = num_sensitivities s in
    let checkvec = Nvector.check y0 in
    if Sundials_config.safe then checkvec y'0;
    let weakref = Weak.create 1 in
    let ida_mem, which, backref =
      match mf with
      | NoSens _ -> c_init_backward s weakref t0 y0 y'0 false
      | WithSens _ -> c_init_backward s weakref t0 y0 y'0 true
    in
    (* ida_mem and backref have to be immediately captured in a session and
       associated with the finalizer before we do anything else.  *)
    let bs = Bsession {
            ida          = ida_mem;
            backref      = backref;
            nroots       = 0;
            checkvec     = checkvec;

            exn_temp     = None;
            id_set       = false;

            resfn        = dummy_resfn;
            rootsfn      = dummy_rootsfn;
            errh         = dummy_errh;
            errw         = dummy_errw;
            ls_solver    = Lsolver_impl.NoSolver;
            ls_callbacks = NoCallbacks;
            ls_precfns   = NoPrecFns;

            sensext    = BwdSensExt {
              parent   = s;
              which    = which;

              bnum_sensitivities = ns;
              bsensarray1 = c_alloc_nvector_array ns;
              bsensarray2 = c_alloc_nvector_array ns;

              bresfn      = (match mf with
                             | NoSens f -> f
                             | _ -> dummy_bresfn_no_sens);

              bresfn_sens = (match mf with
                             | WithSens f -> f
                             | _ -> dummy_bresfn_with_sens);

              bquadrhsfn = dummy_bquadrhsfn_no_sens;
              bquadrhsfn_sens = dummy_bquadrhsfn_with_sens;
              checkbquadvec = (fun _ -> raise Nvector.IncompatibleNvector);
            };
          } in
    Gc.finalise bsession_finalize (tosession bs);
    Weak.set weakref 0 (Some (tosession bs));
    (* Now the session is safe to use.  If any of the following fails and
       raises an exception, the GC will take care of freeing ida_mem and
       backref. *)
    (match varid with
       None -> ()
     | Some x -> set_id bs x);
    set_linear_solver bs linsolv y0;
    set_tolerances bs tol;
    se.bsessions <- (tosession bs) :: bsessions;
    bs

  external c_reinit
      : ('a, 'k) session -> int -> float
        -> ('a, 'k) Nvector.t
        -> ('a, 'k) Nvector.t
        -> unit
      = "c_idas_adj_reinit"

  let reinit bs ?linsolv tb0 yb0 y'b0 =
    if Sundials_config.safe then
      (let checkvec = (tosession bs).checkvec in
       checkvec yb0;
       checkvec y'b0);
    let parent, which = parent_and_which bs in
    c_reinit parent which tb0 yb0 y'b0;
    (match linsolv with
     | Some linsolv -> set_linear_solver bs linsolv yb0
     | None -> ())

  let get_work_space bs = Ida.get_work_space (tosession bs)

  let get_num_steps bs = Ida.get_num_steps (tosession bs)

  let get_num_res_evals bs = Ida.get_num_res_evals (tosession bs)

  let get_num_lin_solv_setups bs =
    Ida.get_num_lin_solv_setups (tosession bs)

  let get_num_err_test_fails bs =
    Ida.get_num_err_test_fails (tosession bs)

  let get_last_order bs = Ida.get_last_order (tosession bs)

  let get_current_order bs = Ida.get_current_order (tosession bs)

  let get_last_step bs = Ida.get_last_step (tosession bs)

  let get_current_step bs = Ida.get_current_step (tosession bs)

  let get_actual_init_step bs =
    Ida.get_actual_init_step (tosession bs)

  let get_current_time bs = Ida.get_current_time (tosession bs)

  let get_tol_scale_factor bs =
    Ida.get_tol_scale_factor (tosession bs)

  let get_err_weights bs = Ida.get_err_weights (tosession bs)

  let get_est_local_errors bs =
    Ida.get_est_local_errors (tosession bs)

  let get_integrator_stats bs =
    Ida.get_integrator_stats (tosession bs)

  let print_integrator_stats bs oc =
    Ida.print_integrator_stats (tosession bs) oc

  let get_num_nonlin_solv_iters bs =
    Ida.get_num_nonlin_solv_iters (tosession bs)

  let get_num_nonlin_solv_conv_fails bs =
    Ida.get_num_nonlin_solv_conv_fails (tosession bs)

  let get_nonlin_solv_stats bs =
    Ida.get_nonlin_solv_stats (tosession bs)

  module Quadrature = struct (* {{{ *)
    include QuadratureTypes

    external c_quad_initb
        : ('a, 'k) session -> int -> ('a, 'k) Nvector.t -> unit
        = "c_idas_adjquad_initb"
    external c_quad_initbs
        : ('a, 'k) session -> int -> ('a, 'k) Nvector.t -> unit
        = "c_idas_adjquad_initbs"

    let init bs mf y0 =
      let parent, which = parent_and_which bs in
      let se = bwdsensext bs in
      se.checkbquadvec <- Nvector.check y0;
      match mf with
       | NoSens f -> (se.bquadrhsfn <- f;
                      c_quad_initb parent which y0)
       | WithSens f -> (se.bquadrhsfn_sens <- f;
                        c_quad_initbs parent which y0)

    external c_reinit
        : ('a, 'k) session -> int -> ('a, 'k) Nvector.t -> unit
        = "c_idas_adjquad_reinit"

    let reinit bs yqb0 =
      let parent, which = parent_and_which bs in
      let se = bwdsensext bs in
      if Sundials_config.safe then se.checkbquadvec yqb0;
      c_reinit parent which yqb0

    external c_get : ('a, 'k) session -> int -> ('a, 'k) Nvector.t -> float
        = "c_idas_adjquad_get"

    let get bs yqb =
      let parent, which = parent_and_which bs in
      let se = bwdsensext bs in
      if Sundials_config.safe then se.checkbquadvec yqb;
      c_get parent which yqb

    type ('a, 'k) tolerance =
        NoStepSizeControl
      | SStolerances of float * float
      | SVtolerances of float * ('a, 'k) Nvector.t

    external set_err_con : ('a, 'k) session -> int -> bool -> unit
        = "c_idas_adjquad_set_err_con"

    external sv_tolerances
        : ('a, 'k) session -> int -> float -> ('a, 'k) Nvector.t -> unit
        = "c_idas_adjquad_sv_tolerances"

    external ss_tolerances
        : ('a, 'k) session -> int -> float -> float -> unit
        = "c_idas_adjquad_ss_tolerances"

    let set_tolerances bs tol =
      let parent, which = parent_and_which bs in
      match tol with
      | NoStepSizeControl -> set_err_con parent which false
      | SStolerances (rel, abs) -> (ss_tolerances parent which rel abs;
                                    set_err_con parent which true)
      | SVtolerances (rel, abs) -> (let se = bwdsensext bs in
                                    if Sundials_config.safe then
                                      se.checkbquadvec abs;
                                    sv_tolerances parent which rel abs;
                                    set_err_con parent which true)

    let get_num_rhs_evals bs =
      Quadrature.get_num_rhs_evals (tosession bs)

    let get_num_err_test_fails bs =
      Quadrature.get_num_err_test_fails (tosession bs)

    let get_err_weights bs =
      Quadrature.get_err_weights (tosession bs)

    let get_stats bs = Quadrature.get_stats (tosession bs)
  end (* }}} *)
end (* }}} *)

(* Let C code know about some of the values in this module.  *)
external c_init_module : exn array -> unit =
  "c_idas_init_module"

let _ =
  c_init_module

    (* Exceptions must be listed in the same order as
       idas_exn_index.  *)
    [|Quadrature.QuadNotInitialized;
      Quadrature.QuadRhsFuncFailure;
      Quadrature.FirstQuadRhsFuncFailure;
      Quadrature.RepeatedQuadRhsFuncFailure;

      Sensitivity.SensNotInitialized;
      Sensitivity.SensResFuncFailure;
      Sensitivity.FirstSensResFuncFailure;
      Sensitivity.RepeatedSensResFuncFailure;
      Sensitivity.BadSensIdentifier;

      Sensitivity.Quadrature.QuadSensNotInitialized;
      Sensitivity.Quadrature.QuadSensRhsFuncFailure;
      Sensitivity.Quadrature.FirstQuadSensRhsFuncFailure;
      Sensitivity.Quadrature.RepeatedQuadSensRhsFuncFailure;

      Adjoint.AdjointNotInitialized;
      Adjoint.NoForwardCall;
      Adjoint.ForwardReinitFailure;
      Adjoint.ForwardFailure;
      Adjoint.NoBackwardProblem;
      Adjoint.BadFinalTime;
      Adjoint.BadOutputTime;
    |]

