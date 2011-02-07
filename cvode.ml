(* Aug 2010, Timothy Bourke (INRIA) *)

(*
 * NB: The order of variant constructors and record fields is important!
 *     If these types are changed or augmented, the corresponding declarations
 *     in cvode_serial.h (and code in cvode_serial.c) must also be updated.
 *)

module type GENERIC =
  sig

    val extra_time_precision : bool ref
    val print_time : string * string -> float -> unit

    val big_real : float
    type real_array =
      (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
    val new_real_array : int -> real_array

    type real_array2 =
      (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t
    val new_real_array2 : int -> int -> real_array2

    module Carray :
      sig
        type t = real_array

        val kind : (float, Bigarray.float64_elt) Bigarray.kind
        val layout : Bigarray.c_layout Bigarray.layout

        val empty : t
        val create : int -> t
        val of_array : float array -> t
        val fill : t -> float -> unit
        val length : t -> int

        val print_with_time : float -> t -> unit
        val print_with_time' : float -> t -> unit

        val app : (float -> unit) -> t -> unit
        val appi : (int -> float -> unit) -> t -> unit

        val map : (float -> float) -> t -> unit
        val mapi : (int -> float -> float) -> t -> unit

        val clamp : float -> t -> unit

        val vmax_norm : t -> float (* N_VMaxNorm *)
      end

    type rootval_array = Carray.t
    type int_array = (int32, Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array1.t
    val new_int_array  : int -> int_array

    module Roots :
      sig
        type t = int_array
        val empty : t
        val create : int -> t
        val print : t -> unit
        val print' : t -> unit
        val get : t -> int -> bool
        val get' : t -> int -> int
        val set : t -> int -> bool -> unit
        val length : t -> int
        val reset : t -> unit
        val exists : t -> bool

        val app : (bool -> unit) -> t -> unit
        val appi : (int -> bool -> unit) -> t -> unit
      end

    type lmm =
      | Adams
      | BDF

    type preconditioning_type =
      | PrecNone
      | PrecLeft
      | PrecRight
      | PrecBoth

    type bandrange = { mupper : int; mlower : int }

    val sprange_default_maxl : int
    type sprange = { pretype : preconditioning_type; maxl: int }

    type linear_solver =
      | Dense
      | LapackDense
      | Band of bandrange
      | LapackBand of bandrange
      | Diag
      | Spgmr of sprange
      | Spbcg of sprange
      | Sptfqmr of sprange
      | BandedSpgmr of sprange * bandrange
      | BandedSpbcg of sprange * bandrange
      | BandedSptfqmr of sprange * bandrange

    type iter =
      | Newton of linear_solver
      | Functional

    type solver_result =
      | Continue
      | RootsFound
      | StopTimeReached

    type root_direction =
      | Increasing
      | Decreasing
      | IncreasingOrDecreasing

    type error_details = {
        error_code : int;
        module_name : string;
        function_name : string;
        error_message : string;
      }

    (* Solver exceptions *)
    exception IllInput
    exception TooClose
    exception TooMuchWork
    exception TooMuchAccuracy
    exception ErrFailure
    exception ConvergenceFailure
    exception LinearInitFailure
    exception LinearSetupFailure
    exception LinearSolveFailure
    exception RhsFuncErr
    exception FirstRhsFuncFailure
    exception RepeatedRhsFuncErr
    exception UnrecoverableRhsFuncErr
    exception RootFuncFailure

    (* get_dky exceptions *)
    exception BadK
    exception BadT
    exception BadDky

    val no_roots : (int * ('a -> 'b -> 'c -> unit))

    (* Throw inside the f callback if the derivatives cannot be calculated at
       the given time. *)
    exception RecoverableFailure

    type integrator_stats = {
        num_steps : int;
        num_rhs_evals : int;
        num_lin_solv_setups : int;
        num_err_test_fails : int;
        last_order : int;
        current_order : int;
        actual_init_step : float;
        last_step : float;
        current_step : float;
        current_time : float
      }

    (* direct linear solvers functions *)

    (* Thrown by GETRF routines for a zero diagonal element at the given
       column index. *)
    exception ZeroDiagonalElement of int

    module Densematrix :
      sig
        type t

        val new_dense_mat  : int * int -> t
        val print_mat      : t -> unit

        val set_to_zero    : t -> unit
        val add_identity   : t -> unit
        val dense_copy     : t -> t -> unit
        val dense_scale    : float -> t -> unit
        val dense_getrf    : t -> int_array -> unit
        val dense_getrs    : t -> int_array -> real_array -> unit
        val dense_potrf    : t -> unit
        val dense_potrs    : t -> real_array -> unit
        val dense_geqrf    : t -> real_array -> real_array -> unit

        type ormqr = {
              beta : real_array;
              vn   : real_array;
              vm   : real_array;
              work : real_array;
            }

        val dense_ormqr    : t -> ormqr -> unit

        val get : t -> (int * int) -> float
        val set : t -> (int * int) -> float -> unit

        module Direct :
          sig
            type t

            val new_dense_mat  : int * int -> t

            val get : t -> (int * int) -> float
            val set : t -> (int * int) -> float -> unit

            val dense_copy  : t -> t -> int * int -> unit
            val dense_scale : float -> t -> int * int -> unit
            val dense_add_identity : t -> int -> unit
            val dense_getrf : t -> int * int -> int_array -> unit
            val dense_getrs : t -> int -> int_array -> real_array -> unit
            val dense_potrf : t -> int -> unit
            val dense_potrs : t -> int -> real_array -> unit
            val dense_geqrf : t -> int * int -> real_array -> real_array -> unit
            val dense_ormqr : t -> int * int -> ormqr -> unit
          end

      end

    module Bandmatrix :
      sig
        type t

        val new_band_mat : int * int * int * int -> t (* n, mu, ml, smu *)
        val print_mat : t -> unit

        val set_to_zero    : t -> unit
        val add_identity   : t -> unit

        val band_copy : t -> t -> int -> int -> unit
        val band_scale : float -> t -> unit
        val band_gbtrf : t -> int_array -> unit
        val band_gbtrs : t -> int_array -> real_array -> unit

        val get : t -> (int * int) -> float
        val set : t -> (int * int) -> float -> unit

        module Col :
          sig
            type c

            val get_col : t -> int -> c

            val get : c -> int -> int -> float
            val set : c -> int -> int -> float -> unit
          end

        module Direct :
          sig
            type t

            val new_band_mat : int * int * int -> t (* n smu ml *)

            val get : t -> (int * int) -> float
            val set : t -> (int * int) -> float -> unit

            val band_copy : t -> t -> int -> int -> int -> int -> int -> unit
                        (*  a    b    n     a_smu  b_smu  copymu  copyml *)

            val band_scale : float -> t -> int -> int -> int -> int -> unit
                        (*  c         a    n      mu     ml     smu *)

            val band_add_identity : t -> int -> int -> unit
                        (*          a    n      smu *)

            val band_gbtrf : t -> int -> int -> int -> int -> int_array -> unit
                        (*   a    n      mu     ml     smu    p *)

            val band_gbtrs
                : t -> int -> int -> int -> int_array -> real_array -> unit
                (*a    n      smu    ml     p            b *)
          end
      end
  end

module Generic =
  struct

    let extra_time_precision = ref false

    let print_time (s1, s2) t =
      if !extra_time_precision
      then Printf.printf "%s%.15e%s" s1 t s2
      else Printf.printf "%s%e%s" s1 t s2

    external get_big_real : unit -> float
        = "ml_cvode_big_real"
    let big_real = get_big_real ()
    type real_array =
      (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t
    let new_real_array =
      Bigarray.Array1.create Bigarray.float64 Bigarray.c_layout

    type real_array2 =
      (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t
    let new_real_array2 =
      Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout

    module Carray =
      struct
        type t = real_array

        let kind = Bigarray.float64
        let layout = Bigarray.c_layout
        let empty = Bigarray.Array1.create kind layout 0

        let create = Bigarray.Array1.create kind layout
        let of_array = Bigarray.Array1.of_array kind layout

        let fill = Bigarray.Array1.fill

        let length = Bigarray.Array1.dim

        let app f v =
          for i = 0 to (length v - 1) do
            f v.{i}
          done

        let map f v =
          for i = 0 to (length v - 1) do
            v.{i} <- f v.{i}
          done

        let appi f v =
          for i = 0 to (length v - 1) do
            f i v.{i}
          done

        let mapi f v =
          for i = 0 to (length v - 1) do
            v.{i} <- f i v.{i}
          done

        let print_with_time' t v =
          print_time ("", "") t;
          app (Printf.printf "\t% .8f") v;
          print_newline ()

        let print_with_time t v =
          print_time ("", "") t;
          app (Printf.printf "\t% e") v;
          print_newline ()

        let clamp thres =
          let cf v = if abs_float v <= thres then 0.0 else v
          in
          if thres = 0.0 then (fun x -> ())
          else map cf

        external vmax_norm : t -> float
          = "c_ba_vmax_norm"
      end

    (* root arrays *)

    type int_array = (int32, Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array1.t
    let create_int_array = Bigarray.Array1.create Bigarray.int32 Carray.layout
    let new_int_array  = create_int_array

    type rootval_array = Carray.t

    module Roots =
      struct
        type t = int_array

        let reset v = Bigarray.Array1.fill v 0l

        let create n =
          let a = create_int_array n in
          reset a;
          a

        let empty = create 0

        let length = Bigarray.Array1.dim

        let get roots i = roots.{i} <> 0l
        let get' roots i = Int32.to_int roots.{i}

        let set a i v = Bigarray.Array1.set a i (if v then 1l else 0l)

        let appi f v =
          for i = 0 to (length v - 1) do
            f i (v.{i} <> 0l)
          done

        let app f v =
          for i = 0 to (length v - 1) do
            f (v.{i} <> 0l)
          done

        let print vs =
          app (fun v -> print_string (if v then "\t1" else "\t0")) vs;
          print_newline ()

        let print' vs =
          Carray.appi (fun i v -> Printf.printf "\t% ld" v) vs;
          print_newline ()

        let fold_left f a vs =
          let rec check (i, a) =
            if i < 0 then a
            else check (i - 1, f a (Int32.to_int vs.{i}))
          in
          check (Bigarray.Array1.dim vs - 1, a)

        let exists = fold_left (fun a x -> a || x <> 0) false

      end

    type lmm =
      | Adams
      | BDF

    type preconditioning_type =
      | PrecNone
      | PrecLeft
      | PrecRight
      | PrecBoth

    type bandrange = { mupper : int; mlower : int }

    let sprange_default_maxl = 0
    type sprange = { pretype : preconditioning_type; maxl: int }

    type linear_solver =
      | Dense
      | LapackDense
      | Band of bandrange
      | LapackBand of bandrange
      | Diag
      | Spgmr of sprange
      | Spbcg of sprange
      | Sptfqmr of sprange
      | BandedSpgmr of sprange * bandrange
      | BandedSpbcg of sprange * bandrange
      | BandedSptfqmr of sprange * bandrange

    type iter =
      | Newton of linear_solver
      | Functional

    type solver_result =
      | Continue
      | RootsFound
      | StopTimeReached

    type root_direction =
      | Increasing
      | Decreasing
      | IncreasingOrDecreasing

    let int_of_root_direction x =
      match x with
      | Increasing -> 1l
      | Decreasing -> -1l
      | IncreasingOrDecreasing -> 0l

    type error_details = {
        error_code : int;
        module_name : string;
        function_name : string;
        error_message : string;
      }

    (* Solver exceptions *)
    exception IllInput
    exception TooClose
    exception TooMuchWork
    exception TooMuchAccuracy
    exception ErrFailure
    exception ConvergenceFailure
    exception LinearInitFailure
    exception LinearSetupFailure
    exception LinearSolveFailure
    exception RhsFuncErr
    exception FirstRhsFuncFailure
    exception RepeatedRhsFuncErr
    exception UnrecoverableRhsFuncErr
    exception RootFuncFailure

    (* get_dky exceptions *)
    exception BadK
    exception BadT
    exception BadDky

    let no_roots = (0, (fun _ _ _ -> ()))

    (* Throw inside the f callback if the derivatives cannot be calculated at
       the given time. *)
    exception RecoverableFailure

    type integrator_stats = {
        num_steps : int;
        num_rhs_evals : int;
        num_lin_solv_setups : int;
        num_err_test_fails : int;
        last_order : int;
        current_order : int;
        actual_init_step : float;
        last_step : float;
        current_step : float;
        current_time : float
      }

    exception StopTimeReached

    exception ZeroDiagonalElement of int

    let _ =
      List.iter (fun (nm, ex) -> Callback.register_exception nm ex)
      [
        ("cvode_RecoverableFailure",      RecoverableFailure);

        ("cvode_StopTimeReached",         StopTimeReached);
        ("cvode_IllInput",                IllInput);
        ("cvode_TooClose",                TooClose);
        ("cvode_TooMuchWork",             TooMuchWork);
        ("cvode_TooMuchAccuracy",         TooMuchAccuracy);
        ("cvode_ErrFailure",              ErrFailure);
        ("cvode_ConvergenceFailure",      ConvergenceFailure);
        ("cvode_LinearInitFailure",       LinearInitFailure);
        ("cvode_LinearSetupFailure",      LinearSetupFailure);
        ("cvode_LinearSolveFailure",      LinearSolveFailure);
        ("cvode_RhsFuncErr",              RhsFuncErr);
        ("cvode_FirstRhsFuncFailure",     FirstRhsFuncFailure);
        ("cvode_RepeatedRhsFuncErr",      RepeatedRhsFuncErr);
        ("cvode_UnrecoverableRhsFuncErr", UnrecoverableRhsFuncErr);
        ("cvode_RootFuncFailure",         RootFuncFailure);

        ("cvode_BadK",                    BadK);
        ("cvode_BadT",                    BadT);
        ("cvode_BadDky",                  BadDky);

        ("cvode_ZeroDiagonalElement",     ZeroDiagonalElement 0);
      ]

    (* passing callbacks to c *)

    type handler =
      | RhsFn
      | RootsFn
      | ErrorHandler
      | ErrorWeight
      | JacFn
      | BandJacFn
      | PreSetupFn
      | PreSolveFn
      | JacTimesFn

    let handler_name h = match h with
      | RhsFn        -> "cvode_serial_callback_rhsfn"
      | RootsFn      -> "cvode_serial_callback_rootsfn"
      | ErrorHandler -> "cvode_serial_callback_errorhandler"
      | ErrorWeight  -> "cvode_serial_callback_errorweight"
      | JacFn        -> "cvode_serial_callback_jacfn"
      | BandJacFn    -> "cvode_serial_callback_bandjacfn"
      | PreSetupFn   -> "cvode_serial_callback_presetupfn"
      | PreSolveFn   -> "cvode_serial_callback_presolvefn"
      | JacTimesFn   -> "cvode_serial_callback_jactimesfn"

    (* direct linear solvers functions *)

    (* note: uses DENSE_ELEM rather than the more efficient DENSE_COL. *)
    module Densematrix =
      struct
        type t

        external new_dense_mat  : int * int -> t
            = "c_densematrix_new_dense_mat"

        external print_mat      : t -> unit
            = "c_densematrix_print_mat"

        external set_to_zero    : t -> unit
            = "c_densematrix_set_to_zero"

        external add_identity   : t -> unit
            = "c_densematrix_add_identity"

        external dense_copy     : t -> t -> unit
            = "c_densematrix_dense_copy"

        external dense_scale    : float -> t -> unit
            = "c_densematrix_dense_scale"

        external dense_getrf    : t -> int_array -> unit
            = "c_densematrix_getrf"

        external dense_getrs    : t -> int_array -> real_array -> unit
            = "c_densematrix_getrs"

        external dense_potrf    : t -> unit
            = "c_densematrix_potrf"

        external dense_potrs    : t -> real_array -> unit
            = "c_densematrix_potrs"

        external dense_geqrf    : t -> real_array -> real_array -> unit
            = "c_densematrix_geqrf"

        type ormqr = {
              beta : real_array;
              vn   : real_array;
              vm   : real_array;
              work : real_array;
            }

        external dense_ormqr : t -> ormqr -> unit
            = "c_densematrix_ormqr"

        external get : t -> (int * int) -> float
            = "c_densematrix_get"

        external set : t -> (int * int) -> float -> unit
            = "c_densematrix_set"

        module Direct =
          struct
            type t

            external new_dense_mat  : int * int -> t
                = "c_densematrix_direct_new_dense_mat"

            external get : t -> (int * int) -> float
                = "c_densematrix_direct_get"

            external set : t -> (int * int) -> float -> unit
                = "c_densematrix_direct_set"

            external dense_copy  : t -> t -> int * int -> unit
                = "c_densematrix_direct_copy"

            external dense_scale : float -> t -> int * int -> unit
                = "c_densematrix_direct_scale"

            external dense_add_identity : t -> int -> unit
                = "c_densematrix_direct_add_identity"

            external dense_getrf : t -> int * int -> int_array -> unit
                = "c_densematrix_direct_getrf"

            external dense_getrs : t -> int -> int_array -> real_array -> unit
                = "c_densematrix_direct_getrs"

            external dense_potrf : t -> int -> unit
                = "c_densematrix_direct_potrf"

            external dense_potrs : t -> int -> real_array -> unit
                = "c_densematrix_direct_potrs"

            external dense_geqrf : t -> int * int -> real_array -> real_array -> unit
                = "c_densematrix_direct_geqrf"

            external dense_ormqr : t -> int * int -> ormqr -> unit
                = "c_densematrix_direct_ormqr"
          end
      end

    (* note: uses BAND_ELEM rather than the more efficient BAND_COL/BAND_COL_ELEM *)
    module Bandmatrix =
      struct
        type t

        external new_band_mat : int * int * int * int -> t
            = "c_bandmatrix_new_band_mat"

        external print_mat : t -> unit
            = "c_densematrix_print_mat"
              (* NB: same as densematrix *)

        external set_to_zero    : t -> unit
            = "c_densematrix_set_to_zero"
              (* NB: same as densematrix *)

        external add_identity : t -> unit
            = "c_densematrix_add_identity"
              (* NB: same as densematrix *)

        external band_copy : t -> t -> int -> int -> unit
            = "c_bandmatrix_copy"

        external band_scale : float -> t -> unit
            = "c_bandmatrix_scale"

        external band_gbtrf : t -> int_array -> unit
            = "c_bandmatrix_gbtrf"

        external band_gbtrs : t -> int_array -> real_array -> unit
            = "c_bandmatrix_gbtrs"

        external get : t -> (int * int) -> float
            = "c_bandmatrix_get"

        external set : t -> (int * int) -> float -> unit
            = "c_bandmatrix_set"

        module Col =
          struct
            type c

            external get_col : t -> int -> c
                = "c_bandmatrix_col_get_col"

            external get : c -> int -> int -> float
                = "c_bandmatrix_col_get"

            external set : c -> int -> int -> float -> unit
                = "c_bandmatrix_col_set"
          end

        module Direct =
          struct
            type t

            external new_band_mat : int * int * int -> t
                = "c_bandmatrix_direct_new_band_mat"

            external get : t -> (int * int) -> float
                = "c_densematrix_direct_get"
                (* NB: same as densematrix_direct *)

            external set : t -> (int * int) -> float -> unit
                = "c_densematrix_direct_set"
                (* NB: same as densematrix_direct *)

            external band_copy' : t -> t -> int * int * int * int * int -> unit
                = "c_bandmatrix_direct_copy"

            let band_copy a b n a_smu b_smu copymu copyml
                = band_copy' a b (n, a_smu, b_smu, copymu, copyml)

            external band_scale' : float -> t -> int * int * int * int -> unit
                = "c_bandmatrix_direct_scale"

            let band_scale c a n mu ml smu = band_scale' c a (n, mu, ml, smu)

            external band_add_identity : t -> int -> int -> unit
                = "c_bandmatrix_direct_add_identity"

            external band_gbtrf' : t -> int * int * int * int -> int_array -> unit
                = "c_bandmatrix_direct_gbtrf"

            let band_gbtrf a n mu ml smu p = band_gbtrf' a (n, mu, ml, smu) p

            external band_gbtrs'
                : t -> int * int * int -> int_array -> real_array -> unit
                = "c_bandmatrix_direct_gbtrs"

            let band_gbtrs a n smu ml p b = band_gbtrs' a (n, smu, ml) p b
          end
      end
  end

module Serial =
  struct
    include Generic

    type nvec = Carray.t
    type val_array = Carray.t
    type der_array = Carray.t

    type session

    (* interface *)

    external external_register_handler : session -> handler -> unit
        = "c_register_handler"

    let register_handler s h f =
      Callback.register (handler_name h) f;
      external_register_handler s h

    external external_init
        : lmm -> iter -> val_array -> int -> float -> session
        = "c_ba_init"

    external nroots         : session -> int
        = "c_nroots"

    external neqs           : session -> int
        = "c_neqs"

    external reinit
        : session -> float -> val_array -> unit
        = "c_ba_reinit"

    external sv_tolerances  : session -> float -> nvec -> unit
        = "c_ba_sv_tolerances"
    external ss_tolerances  : session -> float -> float -> unit
        = "c_ss_tolerances"
    external wf_tolerances  : session -> (val_array -> nvec -> unit) -> unit
        = "c_ba_wf_tolerances"

    let wf_tolerances s efun =
      register_handler s ErrorWeight efun;
      wf_tolerances s efun

    external get_root_info  : session -> Roots.t -> unit
        = "c_get_root_info"

    external free           : session -> unit
        = "c_free"

    external normal
        : session -> float -> val_array -> float * solver_result
        = "c_ba_normal"

    external one_step
        : session -> float -> val_array -> float * solver_result
        = "c_ba_one_step"

    external get_dky
        : session -> float -> int -> nvec -> unit
        = "c_ba_get_dky"

    let init' lmm iter f (num_roots, roots) y0 t0 =
      let s = external_init lmm iter y0 num_roots t0 in
      register_handler s RhsFn f;
      register_handler s RootsFn roots;
      s

    let init lmm iter f roots y0 = init' lmm iter f roots y0 0.0

    external get_integrator_stats   : session -> integrator_stats
        = "c_get_integrator_stats"

    external last_step_size         : session -> float
        = "c_last_step_size"

    external next_step_size         : session -> float
        = "c_next_step_size"

    external get_work_space         : session -> int * int
        = "c_get_work_space"

    external get_num_steps          : session -> int
        = "c_get_num_steps"

    external get_num_rhs_evals      : session -> int
        = "c_get_num_rhs_evals"

    external get_num_lin_solv_setups : session -> int
        = "c_get_num_lin_solv_setups"

    external get_num_err_test_fails : session -> int
        = "c_get_num_err_test_fails"

    external get_last_order         : session -> int
        = "c_get_last_order"

    external get_current_order      : session -> int
        = "c_get_current_order"

    external get_actual_init_step   : session -> float
        = "c_get_actual_init_step"

    external get_last_step          : session -> float
        = "c_get_last_step"

    external get_current_step       : session -> float
        = "c_get_current_step"

    external get_current_time       : session -> float
        = "c_get_current_time"

    let print_integrator_stats s =
      let stats = get_integrator_stats s
      in
        Printf.printf "num_steps = %d\n"           stats.num_steps;
        Printf.printf "num_rhs_evals = %d\n"       stats.num_rhs_evals;
        Printf.printf "num_lin_solv_setups = %d\n" stats.num_lin_solv_setups;
        Printf.printf "num_err_test_fails = %d\n"  stats.num_err_test_fails;
        Printf.printf "last_order = %d\n"          stats.last_order;
        Printf.printf "current_order = %d\n"       stats.current_order;
        Printf.printf "actual_init_step = %e\n"    stats.actual_init_step;
        Printf.printf "last_step = %e\n"           stats.last_step;
        Printf.printf "current_step = %e\n"        stats.current_step;
        Printf.printf "current_time = %e\n"        stats.current_time;

    external set_error_file         : session -> string -> bool -> unit
        = "c_set_error_file"

    external enable_err_handler_fn  : session -> unit
        = "c_enable_err_handler_fn"

    let set_err_handler_fn s errh =
      register_handler s ErrorHandler errh;
      enable_err_handler_fn s

    external set_max_ord            : session -> int -> unit
        = "c_set_max_ord"
    external set_max_num_steps      : session -> int -> unit
        = "c_set_max_num_steps"
    external set_max_hnil_warns     : session -> int -> unit
        = "c_set_max_hnil_warns"
    external set_stab_lim_det       : session -> bool -> unit
        = "c_set_stab_lim_det"
    external set_init_step          : session -> float -> unit
        = "c_set_init_step"
    external set_min_step           : session -> float -> unit
        = "c_set_min_step"
    external set_max_step           : session -> float -> unit
        = "c_set_max_step"
    external set_stop_time          : session -> float -> unit
        = "c_set_stop_time"
    external set_max_err_test_fails : session -> int -> unit
        = "c_set_max_err_test_fails"
    external set_max_nonlin_iters   : session -> int -> unit
        = "c_set_max_nonlin_iters"
    external set_max_conv_fails     : session -> int -> unit
        = "c_set_max_conv_fails"
    external set_nonlin_conv_coef   : session -> float -> unit
        = "c_set_nonlin_conv_coef"
    external set_iter_type          : session -> iter -> unit
        = "c_set_iter_type"

    external set_root_direction'    : session -> int_array -> unit
        = "c_set_root_direction"

    let set_root_direction s rda =
      let n = nroots s in
      let rdirs = create_int_array n in
      if (n > Array.length rda)
        then Bigarray.Array1.fill rdirs
                (int_of_root_direction IncreasingOrDecreasing);
      Array.iteri (fun i v -> rdirs.{i} <- int_of_root_direction v) rda;
      set_root_direction' s rdirs

    let set_all_root_directions s rd =
      let n = nroots s in
      if (n > 0) then begin
        let rdirs = create_int_array n in
        Bigarray.Array1.fill rdirs (int_of_root_direction rd);
        set_root_direction' s rdirs
      end; ()

    external set_no_inactive_root_warn      : session -> unit
        = "c_set_no_inactive_root_warn"

    external get_num_stab_lim_order_reds    : session -> int
        = "c_get_num_stab_lim_order_reds"

    external get_tol_scale_factor           : session -> float
        = "c_get_tol_scale_factor"

    external get_err_weights                : session -> nvec -> unit
        = "c_ba_get_err_weights"

    external get_est_local_errors           : session -> nvec -> unit
        = "c_ba_get_est_local_errors"

    external get_num_nonlin_solv_iters      : session -> int
        = "c_get_num_nonlin_solv_iters"

    external get_num_nonlin_solv_conv_fails : session -> int
        = "c_get_num_nonlin_solv_conv_fails"

    external get_num_g_evals                : session -> int
        = "c_get_num_g_evals"

    type 't jacobian_arg =
      {
        jac_t   : float;
        jac_y   : val_array;
        jac_fy  : val_array;
        jac_tmp : 't
      }

    type triple_tmp = val_array * val_array * val_array

    module Dls =
      struct
        external enable_dense_jac_fn  : session -> unit
            = "c_ba_dls_enable_dense_jac_fn"

        let set_dense_jac_fn s f =
            register_handler s JacFn f;
            enable_dense_jac_fn s

        external enable_band_jac_fn   : session -> unit
            = "c_ba_dls_enable_band_jac_fn"

        let set_band_jac_fn s f =
            register_handler s BandJacFn f;
            enable_band_jac_fn s

        external get_num_jac_evals    : session -> int
            = "c_dls_get_num_jac_evals"

        external get_num_rhs_evals    : session -> int
            = "c_dls_get_num_rhs_evals"
      end

    module Diag =
      struct
        external get_num_rhs_evals    : session -> int
            = "c_diag_get_num_rhs_evals"
      end

    module BandPrec =
      struct
        external get_num_rhs_evals    : session -> int
            = "c_bandprec_get_num_rhs_evals"
      end

    module Spils =
      struct
        type solve_arg =
          {
            rhs   : val_array;
            gamma : float;
            delta : float;
            left  : bool;
          }

        type single_tmp = val_array

        type gramschmidt_type =
          | ModifiedGS
          | ClassicalGS

        external enable_preconditioner  : session -> unit
            = "c_ba_enable_preconditioner"

        let set_preconditioner s fsetup fsolve =
            register_handler s PreSetupFn fsetup;
            register_handler s PreSolveFn fsolve;
            enable_preconditioner s

        external enable_jac_times_vec_fn : session -> unit
            = "c_ba_enable_jac_times_vec_fn"

        let set_jac_times_vec_fn s f =
            register_handler s JacTimesFn f;
            enable_jac_times_vec_fn s

        external set_prec_type
            : session -> preconditioning_type -> unit
            = "c_set_prec_type"

        external set_gs_type : session -> gramschmidt_type -> unit
            = "c_set_gs_type"

        external set_eps_lin : session -> float -> unit
            = "c_set_eps_lin"

        external set_maxl   : session -> int -> unit
            = "c_set_maxl"

        external get_num_lin_iters      : session -> int
            = "c_spils_get_num_lin_iters"

        external get_num_conv_fails     : session -> int
            = "c_spils_get_num_conv_fails"

        external get_work_space         : session -> int * int
            = "c_spils_get_work_space"

        external get_num_prec_evals     : session -> int
            = "c_spils_get_num_prec_evals"

        external get_num_prec_solves    : session -> int
            = "c_spils_get_num_prec_solves"

        external get_num_jtimes_evals   : session -> int
            = "c_spils_get_num_jtimes_evals"

        external get_num_rhs_evals      : session -> int
            = "c_spils_get_num_rhs_evals"

      end
  end

module Nvector =
  struct
    include Generic

    type 'a nvector = 'a Nvector.nvector
    type 'a session

    (* interface *)

    external external_register_handler : 'a session -> handler -> unit
        = "c_register_handler"

    let register_handler s h f =
      Callback.register (handler_name h) f;
      external_register_handler s h

    external external_init
        : lmm -> iter -> 'a nvector -> int -> float -> 'a session
        = "c_nvec_init"

    external nroots : 'a session -> int
        = "c_nroots"

    external neqs   : 'a session -> int
        = "c_neqs"

    external reinit
        : 'a session -> float -> 'a nvector -> unit
        = "c_nvec_reinit"

    external sv_tolerances  : 'a session -> float -> 'a nvector -> unit
        = "c_nvec_sv_tolerances"
    external ss_tolerances  : 'a session -> float -> float -> unit
        = "c_ss_tolerances"
    external wf_tolerances  : 'a session -> ('a -> 'a -> unit)
                             -> unit
        = "c_nvec_wf_tolerances"

    let wf_tolerances s efun =
      register_handler s ErrorWeight efun;
      wf_tolerances s efun

    external get_root_info  : 'a session -> Roots.t -> unit
        = "c_get_root_info"

    external free           : 'a session -> unit
        = "c_free"

    external normal
        : 'a session -> float -> 'a nvector -> float * solver_result
        = "c_nvec_normal"

    external one_step
        : 'a session -> float -> 'a nvector -> float * solver_result
        = "c_nvec_one_step"

    external get_dky
        : 'a session -> float -> int -> 'a nvector -> unit
        = "c_nvec_get_dky"

    let init' lmm iter f (num_roots, roots) y0 t0 =
      let s = external_init lmm iter y0 num_roots t0 in
      register_handler s RhsFn f;
      register_handler s RootsFn roots;
      s

    let init lmm iter f roots y0 = init' lmm iter f roots y0 0.0

    external get_integrator_stats : 'a session -> integrator_stats
        = "c_get_integrator_stats"

    external last_step_size         : 'a session -> float
        = "c_last_step_size"

    external next_step_size         : 'a session -> float
        = "c_next_step_size"
 
    external get_work_space         : 'a session -> int * int
        = "c_get_work_space"

    external get_num_steps          : 'a session -> int
        = "c_get_num_steps"

    external get_num_rhs_evals      : 'a session -> int
        = "c_get_num_rhs_evals"

    external get_num_lin_solv_setups : 'a session -> int
        = "c_get_num_lin_solv_setups"

    external get_num_err_test_fails : 'a session -> int
        = "c_get_num_err_test_fails"

    external get_last_order         : 'a session -> int
        = "c_get_last_order"

    external get_current_order      : 'a session -> int
        = "c_get_current_order"

    external get_actual_init_step   : 'a session -> float
        = "c_get_actual_init_step"

    external get_last_step          : 'a session -> float
        = "c_get_last_step"

    external get_current_step       : 'a session -> float
        = "c_get_current_step"

    external get_current_time       : 'a session -> float
        = "c_get_current_time"

    let print_integrator_stats s =
      let stats = get_integrator_stats s
      in
        Printf.printf "num_steps = %d\n"           stats.num_steps;
        Printf.printf "num_rhs_evals = %d\n"       stats.num_rhs_evals;
        Printf.printf "num_lin_solv_setups = %d\n" stats.num_lin_solv_setups;
        Printf.printf "num_err_test_fails = %d\n"  stats.num_err_test_fails;
        Printf.printf "last_order = %d\n"          stats.last_order;
        Printf.printf "current_order = %d\n"       stats.current_order;
        Printf.printf "actual_init_step = %e\n"    stats.actual_init_step;
        Printf.printf "last_step = %e\n"           stats.last_step;
        Printf.printf "current_step = %e\n"        stats.current_step;
        Printf.printf "current_time = %e\n"        stats.current_time;

    external set_error_file : 'a session -> string -> bool -> unit
        = "c_set_error_file"

    external enable_err_handler_fn : 'a session -> unit
        = "c_enable_err_handler_fn"

    let set_err_handler_fn s errh =
      register_handler s ErrorHandler errh;
      enable_err_handler_fn s

    external set_max_ord            : 'a session -> int -> unit
        = "c_set_max_ord"
    external set_max_num_steps      : 'a session -> int -> unit
        = "c_set_max_num_steps"
    external set_max_hnil_warns     : 'a session -> int -> unit
        = "c_set_max_hnil_warns"
    external set_stab_lim_det       : 'a session -> bool -> unit
        = "c_set_stab_lim_det"
    external set_init_step          : 'a session -> float -> unit
        = "c_set_init_step"
    external set_min_step           : 'a session -> float -> unit
        = "c_set_min_step"
    external set_max_step           : 'a session -> float -> unit
        = "c_set_max_step"
    external set_stop_time          : 'a session -> float -> unit
        = "c_set_stop_time"
    external set_max_err_test_fails : 'a session -> int -> unit
        = "c_set_max_err_test_fails"
    external set_max_nonlin_iters   : 'a session -> int -> unit
        = "c_set_max_nonlin_iters"
    external set_max_conv_fails     : 'a session -> int -> unit
        = "c_set_max_conv_fails"
    external set_nonlin_conv_coef   : 'a session -> float -> unit
        = "c_set_nonlin_conv_coef"
    external set_iter_type          : 'a session -> iter -> unit
        = "c_set_iter_type"

    external set_root_direction' : 'a session -> int_array -> unit
        = "c_set_root_direction"

    let set_root_direction s rda =
      let n = nroots s in
      let rdirs = create_int_array n in
      if (n > Array.length rda)
        then Bigarray.Array1.fill rdirs
                (int_of_root_direction IncreasingOrDecreasing);
      Array.iteri (fun i v -> rdirs.{i} <- int_of_root_direction v) rda;
      set_root_direction' s rdirs

    let set_all_root_directions s rd =
      let rdirs = create_int_array (nroots s) in
      Bigarray.Array1.fill rdirs (int_of_root_direction rd);
      set_root_direction' s rdirs

    external set_no_inactive_root_warn      : 'a session -> unit
        = "c_set_no_inactive_root_warn"

    external get_num_stab_lim_order_reds    : 'a session -> int
        = "c_get_num_stab_lim_order_reds"

    external get_tol_scale_factor           : 'a session -> float
        = "c_get_tol_scale_factor"

    external get_err_weights                : 'a session -> 'a nvector -> unit
        = "c_nvec_get_err_weights"

    external get_est_local_errors           : 'a session -> 'a nvector -> unit
        = "c_nvec_get_est_local_errors"

    external get_num_nonlin_solv_iters      : 'a session -> int
        = "c_get_num_nonlin_solv_iters"

    external get_num_nonlin_solv_conv_fails : 'a session -> int
        = "c_get_num_nonlin_solv_conv_fails"

    external get_num_g_evals                : 'a session -> int
        = "c_get_num_g_evals"

    type ('t, 'a) jacobian_arg =
      {
        jac_t   : float;
        jac_y   : 'a;
        jac_fy  : 'a;
        jac_tmp : 't
      }

    type 'a triple_tmp = 'a * 'a * 'a

    module Dls =
      struct
        external enable_dense_jac_fn    : 'a session -> unit
            = "c_nvec_dls_enable_dense_jac_fn"

        let set_dense_jac_fn s f =
            register_handler s JacFn f;
            enable_dense_jac_fn s

        external enable_band_jac_fn     : 'a session -> unit
            = "c_nvec_dls_enable_band_jac_fn"

        let set_band_jac_fn s f =
            register_handler s BandJacFn f;
            enable_band_jac_fn s

        external get_num_jac_evals      : 'a session -> int
            = "c_dls_get_num_jac_evals"

        external get_num_rhs_evals      : 'a session -> int
            = "c_dls_get_num_rhs_evals"
      end

    module Diag =
      struct
        external get_num_rhs_evals      : 'a session -> int
            = "c_diag_get_num_rhs_evals"
      end

    module BandPrec =
      struct
        external get_num_rhs_evals      : 'a session -> int
            = "c_bandprec_get_num_rhs_evals"
      end

    module Spils =
      struct
        type 'a solve_arg =
          {
            rhs   : 'a;
            gamma : float;
            delta : float;
            left  : bool;
          }

        type 'a single_tmp = 'a nvector

        type gramschmidt_type =
          | ModifiedGS
          | ClassicalGS

        external enable_preconditioner      : 'a session -> unit
            = "c_nvec_enable_preconditioner"

        let set_preconditioner s fsetup fsolve =
            register_handler s PreSetupFn fsetup;
            register_handler s PreSolveFn fsolve;
            enable_preconditioner s

        external enable_jac_times_vec_fn    : 'a session -> unit
            = "c_nvec_enable_jac_times_vec_fn"

        let set_jac_times_vec_fn s f =
            register_handler s JacTimesFn f;
            enable_jac_times_vec_fn s

        external set_prec_type : 'a session -> preconditioning_type -> unit
            = "c_set_prec_type"

        external set_gs_type : 'a session -> gramschmidt_type -> unit
            = "c_set_gs_type"

        external set_eps_lin            : 'a session -> float -> unit
            = "c_set_eps_lin"

        external set_maxl               : 'a session -> int -> unit
            = "c_set_maxl"

        external get_num_lin_iters      : 'a session -> int
            = "c_spils_get_num_lin_iters"

        external get_num_conv_fails     : 'a session -> int
            = "c_spils_get_num_conv_fails"

        external get_work_space         : 'a session -> int * int
            = "c_spils_get_work_space"

        external get_num_prec_evals     : 'a session -> int
            = "c_spils_get_num_prec_evals"

        external get_num_prec_solves    : 'a session -> int
            = "c_spils_get_num_prec_solves"

        external get_num_jtimes_evals   : 'a session -> int
            = "c_spils_get_num_jtimes_evals"

        external get_num_rhs_evals      : 'a session -> int
            = "c_spils_get_num_rhs_evals"

      end
  end

