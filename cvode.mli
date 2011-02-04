(* Aug 2010, Timothy Bourke (INRIA) *)

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
    type int_array =
      (int32, Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array1.t
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
            type t = real_array2

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
            type t = real_array2

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

module Serial :
  sig
    include GENERIC

    type nvec = Carray.t
    type val_array = Carray.t
    type der_array = Carray.t

    type session

    val init :
        lmm
        -> iter
        -> (float -> val_array -> der_array -> unit)
        -> (int * (float -> val_array -> rootval_array -> unit))
        -> val_array
        -> session

    val init' :
        lmm
        -> iter
        -> (float -> val_array -> der_array -> unit)
        -> (int * (float -> val_array -> rootval_array -> unit))
        -> val_array
        -> float (* start time *)
        -> session

    val nroots : session -> int
    val neqs : session -> int

    val reinit : session -> float -> val_array -> unit

    val sv_tolerances : session -> float -> nvec -> unit
    val ss_tolerances : session -> float -> float -> unit
    val wf_tolerances : session -> (val_array -> nvec -> unit) -> unit

    val get_root_info : session -> Roots.t -> unit

    val normal : session -> float -> val_array -> float * solver_result
    val one_step : session -> float -> val_array -> float * solver_result

    val free : session -> unit
    val get_dky : session -> float -> int -> nvec -> unit

    val get_work_space          : session -> int * int
    val get_num_steps           : session -> int
    val get_num_rhs_evals       : session -> int
    val get_num_lin_solv_setups : session -> int
    val get_num_err_test_fails  : session -> int
    val get_last_order          : session -> int
    val get_current_order       : session -> int

    val get_actual_init_step    : session -> float
    val get_last_step           : session -> float
    val get_current_step        : session -> float
    val get_current_time        : session -> float

    val get_integrator_stats    : session -> integrator_stats
    val print_integrator_stats  : session -> unit
    val last_step_size          : session -> float
    val next_step_size          : session -> float

    (* optional input functions *)

    (* path to an error file, whether to truncate (true) or append (false) *)
    val set_error_file : session -> string -> bool -> unit

    (* Error handler function *)
    val set_err_handler_fn : session -> (error_details -> unit) -> unit

    (* Maximum order for BDF or Adams method *)
    val set_max_ord : session -> int -> unit

    (* Maximum no. of internal steps before tout *)
    val set_max_num_steps : session -> int -> unit

    (* Maximum no. of warnings for tn + h = tn *)
    val set_max_hnil_warns : session -> int -> unit

    (* Flag to activate stability limit detection *)
    val set_stab_lim_det : session -> bool -> unit

    (* Initial step size *)
    val set_init_step : session -> float -> unit

    (* Minimum absolute step size *)
    val set_min_step : session -> float -> unit

    (* Maximum absolute step size *)
    val set_max_step : session -> float -> unit

    (* Value of tstop *)
    val set_stop_time : session -> float -> unit

    (* Maximum no. of error test failures *)
    val set_max_err_test_fails : session -> int -> unit

    (* Maximum no. of nonlinear iterations *)
    val set_max_nonlin_iters : session -> int -> unit

    (* Maximum no. of convergence failures *)
    val set_max_conv_fails : session -> int -> unit

    (* Coefficient in the nonlinear convergence test *)
    val set_nonlin_conv_coef : session -> float -> unit

    (* Nonlinear iteration type *)
    val set_iter_type : session -> iter -> unit

    (* Direction of zero-crossing *)
    val set_root_direction : session -> root_direction array -> unit
    val set_all_root_directions : session -> root_direction -> unit

    (* Disable rootfinding warnings *)
    val set_no_inactive_root_warn : session -> unit

    (* Optional output functions *)

    (* No. of order reductions due to stability limit detection *)
    val get_num_stab_lim_order_reds : session -> int

    (* Suggested factor for tolerance scaling *)
    val get_tol_scale_factor : session -> float

    (* Error weight vector for state variables *)
    val get_err_weights : session -> nvec -> unit

    (* Estimated local error vector *)
    val get_est_local_errors : session -> nvec -> unit

    (* No. of nonlinear solver iterations  *)
    val get_num_nonlin_solv_iters : session -> int

    (* No. of nonlinear convergence failures  *)
    val get_num_nonlin_solv_conv_fails : session -> int

    (* No. of calls to user root function  *)
    val get_num_g_evals : session -> int

    (* direct linear solvers functions *)

    type 't jacobian_arg =
      {
        jac_t   : float;
        jac_y   : val_array;
        jac_fy  : val_array;
        jac_tmp : 't
      }

    type triple_tmp = val_array * val_array * val_array

    module Dls :
      sig
        val set_dense_jac_fn :
             session
          -> (triple_tmp jacobian_arg -> Densematrix.t -> unit)
          -> unit

        val set_band_jac_fn :
             session
          -> (triple_tmp jacobian_arg -> int -> int -> Bandmatrix.t -> unit)
          -> unit

        (* No. of Jacobian evaluations *)
        val get_num_jac_evals : session -> int

        (* No. of r.h.s. calls for finite diff. Jacobian evals. *)
        val get_num_rhs_evals : session -> int
      end

    module Diag :
      sig
        (* No. of r.h.s. calls for finite diff. Jacobian evals. *)
        val get_num_rhs_evals : session -> int
      end

    module BandPrec :
      sig
        (* No. of r.h.s. calls for finite diff. banded Jacobian evals. *)
        val get_num_rhs_evals : session -> int
      end

    (* iterative linear solvers *)
    module Spils :
      sig
        type solve_arg =
          {
            rhs   : val_array;
            gamma : float;
            delta : float;
            left  : bool; (* true: left, false: right *)
          }

        type single_tmp = val_array

        type gramschmidt_type =
          | ModifiedGS
          | ClassicalGS

        val set_preconditioner :
          session
          -> (triple_tmp jacobian_arg -> bool -> float -> bool)
          -> (single_tmp jacobian_arg -> solve_arg -> val_array -> unit)
          -> unit

        val set_jac_times_vec_fn :
          session
          -> (single_tmp jacobian_arg
              -> val_array (* v *)
              -> val_array (* Jv *)
              -> unit)
          -> unit

        val set_prec_type : session -> preconditioning_type -> unit

        val set_gs_type :
          session -> gramschmidt_type -> unit

        val set_eps_lin : session -> float -> unit

        val set_maxl : session -> int -> unit

        val get_work_space       : session -> int * int
        val get_num_prec_evals   : session -> int
        val get_num_prec_solves  : session -> int
        val get_num_lin_iters    : session -> int
        val get_num_conv_fails   : session -> int
        val get_num_jtimes_evals : session -> int
        val get_num_rhs_evals    : session -> int
      end
  end

module Nvector :
  sig
    include GENERIC

    type 'a nvector = 'a Nvector.nvector
    type 'a session

    val init :
        lmm
        -> iter
        -> (float -> 'a -> 'a -> unit)
        -> (int * (float -> 'a -> rootval_array -> unit))
        -> 'a nvector
        -> 'a session

    val init' :
        lmm
        -> iter
        -> (float -> 'a -> 'a -> unit)
        -> (int * (float -> 'a -> rootval_array -> unit))
        -> 'a nvector
        -> float (* start time *)
        -> 'a session

    val nroots : 'a session -> int
    val neqs : 'a session -> int

    val reinit : 'a session -> float -> 'a nvector -> unit

    val sv_tolerances : 'a session -> float -> 'a nvector -> unit
    val ss_tolerances : 'a session -> float -> float -> unit
    val wf_tolerances : 'a session -> ('a -> 'a -> unit) -> unit

    val get_root_info : 'a session -> Roots.t -> unit

    val normal : 'a session -> float -> 'a nvector -> float * solver_result
    val one_step : 'a session -> float -> 'a nvector -> float * solver_result

    val free : 'a session -> unit
    val get_dky : 'a session -> float -> int -> 'a nvector -> unit

    val get_work_space          : 'a session -> int * int
    val get_num_steps           : 'a session -> int
    val get_num_rhs_evals       : 'a session -> int
    val get_num_lin_solv_setups : 'a session -> int
    val get_num_err_test_fails  : 'a session -> int
    val get_last_order          : 'a session -> int
    val get_current_order       : 'a session -> int

    val get_actual_init_step    : 'a session -> float
    val get_last_step           : 'a session -> float
    val get_current_step        : 'a session -> float
    val get_current_time        : 'a session -> float

    val get_integrator_stats    : 'a session -> integrator_stats
    val print_integrator_stats  : 'a session -> unit
    val last_step_size          : 'a session -> float
    val next_step_size          : 'a session -> float

    (* optional input functions *)

    (* path to an error file, whether to truncate (true) or append (false) *)
    val set_error_file : 'a session -> string -> bool -> unit

    (* Error handler function *)
    val set_err_handler_fn : 'a session -> (error_details -> unit) -> unit

    (* Maximum order for BDF or Adams method *)
    val set_max_ord : 'a session -> int -> unit

    (* Maximum no. of internal steps before tout *)
    val set_max_num_steps : 'a session -> int -> unit

    (* Maximum no. of warnings for tn + h = tn *)
    val set_max_hnil_warns : 'a session -> int -> unit

    (* Flag to activate stability limit detection *)
    val set_stab_lim_det : 'a session -> bool -> unit

    (* Initial step size *)
    val set_init_step : 'a session -> float -> unit

    (* Minimum absolute step size *)
    val set_min_step : 'a session -> float -> unit

    (* Maximum absolute step size *)
    val set_max_step : 'a session -> float -> unit

    (* Value of tstop *)
    val set_stop_time : 'a session -> float -> unit

    (* Maximum no. of error test failures *)
    val set_max_err_test_fails : 'a session -> int -> unit

    (* Maximum no. of nonlinear iterations *)
    val set_max_nonlin_iters : 'a session -> int -> unit

    (* Maximum no. of convergence failures *)
    val set_max_conv_fails : 'a session -> int -> unit

    (* Coefficient in the nonlinear convergence test *)
    val set_nonlin_conv_coef : 'a session -> float -> unit

    (* Nonlinear iteration type *)
    val set_iter_type : 'a session -> iter -> unit

    (* Direction of zero-crossing *)
    val set_root_direction : 'a session -> root_direction array -> unit
    val set_all_root_directions : 'a session -> root_direction -> unit

    (* Disable rootfinding warnings *)
    val set_no_inactive_root_warn : 'a session -> unit

    (* Optional output functions *)

    (* No. of order reductions due to stability limit detection *)
    val get_num_stab_lim_order_reds : 'a session -> int

    (* Suggested factor for tolerance scaling *)
    val get_tol_scale_factor : 'a session -> float

    (* Error weight vector for state variables *)
    val get_err_weights : 'a session -> 'a nvector -> unit

    (* Estimated local error vector *)
    val get_est_local_errors : 'a session -> 'a nvector -> unit

    (* No. of nonlinear solver iterations  *)
    val get_num_nonlin_solv_iters : 'a session -> int

    (* No. of nonlinear convergence failures  *)
    val get_num_nonlin_solv_conv_fails : 'a session -> int

    (* No. of calls to user root function  *)
    val get_num_g_evals : 'a session -> int

    (* direct linear solvers functions *)

    type ('t, 'a) jacobian_arg =
      {
        jac_t   : float;
        jac_y   : 'a;
        jac_fy  : 'a;
        jac_tmp : 't
      }

    type 'a triple_tmp = 'a * 'a * 'a

    module Dls :
      sig
        val set_dense_jac_fn :
             'a session
          -> (('a triple_tmp, 'a) jacobian_arg -> Densematrix.t -> unit)
          -> unit

        val set_band_jac_fn :
             'a session
          -> (('a triple_tmp, 'a) jacobian_arg -> int -> int -> Bandmatrix.t
              -> unit)
          -> unit

        (* No. of Jacobian evaluations *)
        val get_num_jac_evals : 'a session -> int

        (* No. of r.h.s. calls for finite diff. Jacobian evals. *)
        val get_num_rhs_evals : 'a session -> int
      end

    module Diag :
      sig
        (* No. of r.h.s. calls for finite diff. Jacobian evals. *)
        val get_num_rhs_evals : 'a session -> int
      end

    module BandPrec :
      sig
        (* No. of r.h.s. calls for finite diff. banded Jacobian evals. *)
        val get_num_rhs_evals : 'a session -> int
      end

    (* iterative linear solvers *)
    module Spils :
      sig
        type 'a solve_arg =
          {
            rhs   : 'a;
            gamma : float;
            delta : float;
            left  : bool; (* true: left, false: right *)
          }

        type 'a single_tmp = 'a nvector

        type gramschmidt_type =
          | ModifiedGS
          | ClassicalGS

        val set_preconditioner :
          'a session
          -> (('a triple_tmp, 'a) jacobian_arg -> bool -> float -> bool)
          -> (('a single_tmp, 'a) jacobian_arg -> 'a solve_arg -> 'a nvector
              -> unit)
          -> unit

        val set_jac_times_vec_fn :
          'a session
          -> (('a single_tmp, 'a) jacobian_arg
              -> 'a (* v *)
              -> 'a (* Jv *)
              -> unit)
          -> unit

        val set_prec_type : 'a session -> preconditioning_type -> unit

        val set_gs_type :
          'a session -> gramschmidt_type -> unit

        val set_eps_lin : 'a session -> float -> unit

        val set_maxl : 'a session -> int -> unit

        val get_work_space       : 'a session -> int * int
        val get_num_prec_evals   : 'a session -> int
        val get_num_prec_solves  : 'a session -> int
        val get_num_lin_iters    : 'a session -> int
        val get_num_conv_fails   : 'a session -> int
        val get_num_jtimes_evals : 'a session -> int
        val get_num_rhs_evals    : 'a session -> int
      end
  end

