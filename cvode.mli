(***********************************************************************)
(*                                                                     *)
(*              Ocaml interface to Sundials CVODE solver               *)
(*                                                                     *)
(*       Timothy Bourke (INRIA Rennes) and Marc Pouzet (LIENS)         *)
(*                                                                     *)
(*  Copyright 2011 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under the terms of the GNU Library General Public License, with    *)
(*  the special exception on linking described in file LICENSE.        *)
(*                                                                     *)
(***********************************************************************)

(** Vector-independent types and values for the CVODE solver.

 @version VERSION()
 @author Timothy Bourke (INRIA)
 @author Marc Pouzet (LIENS)
 *)

include module type of Sundials

(**
 Specify a linear multistep method in {!Cvode_serial.init} and
 {!Cvode_nvector.init}.

 @cvode <node3#ss:ivp_sol> IVP Solution
 @cvode <node5#sss:cvodemalloc> CVodeCreate
 *)
type lmm =
  | Adams   (** Non-stiff systems; Adams-Moulton formulas *)
  | BDF     (** Stiff systems;     Backward Differentiation Formulas *)

(**
 Specify a solution method in {!Cvode_serial.init} and {!Cvode_nvector.init}.

 @see <CVODE_DOC_ROOT(node3#ss:ivp_sol)> IVP Solution
 @see <CVODE_DOC_ROOT(node5#sss:cvodemalloc)> CVodeCreate
 *)
type iter =
  | Newton of linear_solver (** Newton iteration with a given linear solver *)
  | Functional              (** Functional iteration (non-stiff systems only) *)

(**
 Specify a linear solver in {!Cvode_serial.init}, {!Cvode_serial.set_iter_type},
 {!Cvode_nvector.init}, and {!Cvode_nvector.set_iter_type}.

 The Lapack solvers require that both Sundials and the Ocaml interface were
 built to link with a LAPACK library.

 The Banded Krylov solvers imply an additional call to
 {{:CVODE_DOC_ROOT(node5#sss:cvbandpre)} CVBandPrecInit}.

 @see <CVODE_DOC_ROOT(node5#sss:lin_solv_init)> Linear Solver Specification
                                                 Functions
 *)
and linear_solver =
  | Dense                                   (** Direct with dense matrix *)
  | LapackDense                             (** Direct with dense matrix,
                                                with Lapack *)

  | Band of bandrange                       (** Direct with banded matrix *)
  | LapackBand of bandrange                 (** Direct with banded matrix
                                                with Lapack *)

  | Diag                                    (** Diagonal approximation
                                                of the Jacobian *)
  | Spgmr of sprange                        (** Krylov Spils solver: SPGMR *)
  | Spbcg of sprange                        (** Krylov Spils solver: SPBCG *)
  | Sptfqmr of sprange                      (** Krylov Spils solver: SPFQMR *)

  | BandedSpgmr of sprange * bandrange      (** Krylov Spils solver
                                                with banded matrix: SPGMR *)
  | BandedSpbcg of sprange * bandrange      (** Krylov Spils solver
                                                with banded matrix: SPBCG *)
  | BandedSptfqmr of sprange * bandrange    (** Krylov Spils solver
                                                with banded matrix: SPTFQMR *)

and bandrange = {
    mupper : int;
    mlower : int
  }

and sprange = {
    pretype : preconditioning_type;
    maxl: int
  }

and preconditioning_type =
  | PrecNone
  | PrecLeft
  | PrecRight
  | PrecBoth

(**
 Possible values returned when a solver step function ({!Cvode_serial.normal},
 {!Cvode_serial.one_step}, {!Cvode_nvector.normal}, {!Cvode_nvector.one_step})
 succeeds.
 Failures are indicated by exceptions.

 @see <CVODE_DOC_ROOT(node5#sss:cvode)> CVode
 *)
type solver_result =
  | Continue            (** CV_SUCCESS *)
  | RootsFound          (** CV_ROOT_RETURN *)
  | StopTimeReached     (** CV_TSTOP_RETURN *)


(**
 Exceptions that may be thrown when solver step function
 ({!Cvode_serial.normal}, {!Cvode_serial.one_step}, {!Cvode_nvector.normal},
  {!Cvode_nvector.one_step}) fails.

 @see <CVODE_DOC_ROOT(node5#sss:cvode)> CVode
 *)
exception IllInput                  (** CV_ILL_INPUT *)
exception TooClose                  (** CV_TOO_CLOSE *)
exception TooMuchWork               (** CV_TOO_MUCH_WORK *)
exception TooMuchAccuracy           (** CV_TOO_MUCH_ACC *)
exception ErrFailure                (** CV_ERR_FAIL *)
exception ConvergenceFailure        (** CV_CONV_FAIL *)
exception LinearInitFailure         (** CV_LINIT_FAIL *)
exception LinearSetupFailure        (** CV_LSETUP_FAIL *)
exception LinearSolveFailure        (** CV_LSOLVE_FAIL *)
exception RhsFuncErr                (** CV_RHSFUNC_FAIL *)
exception FirstRhsFuncFailure       (** CV_FIRST_RHSFUNC_FAIL *)
exception RepeatedRhsFuncErr        (** CV_REPTD_RHSFUNC_ERR *)
exception UnrecoverableRhsFuncErr   (** CV_UNREC_RHSFUNC_ERR *)
exception RootFuncFailure           (** CV_RTFUNC_FAIL *)

(**
 Type of values passed to a registered error handler function
 ({!Cvode_serial.set_err_handler_fn}, {!Cvode_nvector.set_err_handler_fn}).

 @see <CVODE_DOC_ROOT(node5#sss:optin_main)> CVodeSetErrHandlerFn
 *)
type error_details = {
    error_code : int;
    module_name : string;
    function_name : string;
    error_message : string;
  }


(**
 Values for root directions. Passed to {!Cvode_serial.set_root_direction},
 {!Cvode_serial.set_all_root_directions}, {!Cvode_nvector.set_root_direction},
 and {!Cvode_nvector.set_all_root_directions}.

 @see <CVODE_DOC_ROOT(node5#sss:optin_root)> CVodeSetRootDirection
 *)
type root_direction =
  | Increasing              (** +1 *)
  | Decreasing              (** -1 *)
  | IncreasingOrDecreasing  (**  0 *)

(**
 These exceptions may be thrown by {!Cvode_serial.get_dky} and
 {!Cvode_nvector.get_dky}.

 @see <CVODE_DOC_ROOT(node5#ss:optional_dky)> CVodeGetDky
 *)
exception BadK      (** k is not in the range 0, 1, ..., q{_u} (CV_BAD_K) *)
exception BadT      (** t is not in the interval
                        [t{_n} - h{_u}, t{_n}] (CV_BAD_T) *)
exception BadDky    (** invalid dky argument (CV_BAD_DKY) *)

(**
 This is a convenience value for signalling to {!Cvode_serial.init} and
 {!Cvode_nvector.init} that there are no roots (zero-crossings) to
 monitor.
 *)
val no_roots : (int * ('a -> 'b -> 'c -> unit))

(**
 This exception may be thrown inside the RHS callback function (f)
 if one or more derivatives cannot be calculated at the given time. *)
exception RecoverableFailure

(**
 This record holds the results returned by
 {!Cvode_serial.get_integrator_stats} and
 {!Cvode_nvector.get_integrator_stats}.
 @see <CVODE_DOC_ROOT(node5#sss:optout_main)> CVodeGetIntegratorStats
 *)
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

(**
  {2 Data structures for Direct Linear Solvers}
  @see <CVODE_DOC_ROOT(node9#s:dls)>  The DLS Modules
 *)

(**
 Thrown by the getrf functions if a zero diagonal element is encountered during
 factorization. The argument indicates the column index (from 1).

 @see <CVODE_DOC_ROOT(node9#ss:dense)> DenseGETRF/denseGETRF 
 *)
exception ZeroDiagonalElement of int

(** {3 Dense linear solver}
 @see <CVODE_DOC_ROOT(node9#ss:dense)> The DENSE Module
 *)

module Densematrix :
  sig
    (**
    This type is essentially a [DlsMat] returned from a call to
    {!new_dense_mat}.

     @see <CVODE_DOC_ROOT(node9#s:dls)>  Type DlsMat
     @see <CVODE_DOC_ROOT(node9#ss:dense)> NewDenseMat 
     *)
    type t

    (**
     @see <CVODE_DOC_ROOT(node9#ss:dense)> NewDenseMat
     *)
    val new_dense_mat  : int * int -> t

    (**
     @see <CVODE_DOC_ROOT(node9#ss:dense)> PrintMat
     *)
    val print_mat      : t -> unit

    (**
     @see <CVODE_DOC_ROOT(node9#ss:dense)> SetToZero
     *)
    val set_to_zero    : t -> unit

    (**
     @see <CVODE_DOC_ROOT(node9#ss:dense)> AddIdentity
     *)
    val add_identity   : t -> unit

    (**
     @see <CVODE_DOC_ROOT(node9#ss:dense)> DenseCopy
     *)
    val copy     : t -> t -> unit

    (**
     @see <CVODE_DOC_ROOT(node9#ss:dense)> DenseScale
     *)
    val scale    : float -> t -> unit

    (**
     @raise ZeroDiagonalElement Zero found in matrix diagonal
     @see <CVODE_DOC_ROOT(node9#ss:dense)> DenseGETRF
     *)
    val getrf    : t -> int_array -> unit

    (**
     @see <CVODE_DOC_ROOT(node9#ss:dense)> DenseGETRS
     *)
    val getrs    : t -> int_array -> real_array -> unit

    (**
     @see <CVODE_DOC_ROOT(node9#ss:dense)> DensePOTRF
     *)
    val potrf    : t -> unit

    (**
     @see <CVODE_DOC_ROOT(node9#ss:dense)> DensePOTRS
     *)
    val potrs    : t -> real_array -> unit

    (**
     @see <CVODE_DOC_ROOT(node9#ss:dense)> DenseGEQRF
     *)
    val geqrf    : t -> real_array -> real_array -> unit

    (**
     @see <CVODE_DOC_ROOT(node9#ss:dense)> DenseORMQR
     *)
    type ormqr = {
          beta : real_array;
          vn   : real_array;
          vm   : real_array;
          work : real_array;
        }

    (**
     @see <CVODE_DOC_ROOT(node9#ss:dense)> DenseORMQR
     *)
    val ormqr : t -> ormqr -> unit

    (**
     @see <CVODE_DOC_ROOT(node9#s:dls)> DENSE_ELEM
     *)
    val get : t -> (int * int) -> float

    (**
     @see <CVODE_DOC_ROOT(node9#s:dls)> DENSE_ELEM
     *)
    val set : t -> (int * int) -> float -> unit

    module Direct :
      sig
        (**
        This type is essentially a [realtype **] returned from a call to
        {!new_dense_mat}.

         @see <CVODE_DOC_ROOT(node9#ss:dense)> Small dense matrices
         @see <CVODE_DOC_ROOT(node9#ss:dense)> newDenseMat 
         *)
        type t

        (**
         @see <CVODE_DOC_ROOT(node9#ss:dense)> newDenseMat
         *)
        val new_dense_mat  : int * int -> t

        (**
         *)
        val get : t -> (int * int) -> float

        (**
         *)
        val set : t -> (int * int) -> float -> unit

        (**
         @see <CVODE_DOC_ROOT(node9#ss:dense)> denseCopy
         *)
        val copy  : t -> t -> int * int -> unit

        (*
         @see <CVODE_DOC_ROOT(node9#ss:dense)> denseScale
         *)
        val scale : float -> t -> int * int -> unit

        (**
         @see <CVODE_DOC_ROOT(node9#ss:dense)> denseAddIdentity
         *)
        val add_identity : t -> int -> unit

        (**
         @raise ZeroDiagonalElement Zero found in matrix diagonal
         @see <CVODE_DOC_ROOT(node9#ss:dense)> denseGETRF
         *)
        val getrf : t -> int * int -> int_array -> unit

        (**
         @see <CVODE_DOC_ROOT(node9#ss:dense)> denseGETRS
         *)
        val getrs : t -> int -> int_array -> real_array -> unit

        (**
         @see <CVODE_DOC_ROOT(node9#ss:dense)> densePOTRF
         *)
        val potrf : t -> int -> unit

        (**
         @see <CVODE_DOC_ROOT(node9#ss:dense)> densePOTRS
         *)
        val potrs : t -> int -> real_array -> unit

        (**
         @see <CVODE_DOC_ROOT(node9#ss:dense)> denseGEQRF
         *)
        val geqrf : t -> int * int -> real_array -> real_array -> unit

        (**
         @see <CVODE_DOC_ROOT(node9#ss:dense)> denseORMQR
         *)
        val ormqr : t -> int * int -> ormqr -> unit
      end
  end

(** {3 Banded linear solver}
 @see <CVODE_DOC_ROOT(node9#ss:band)> The BAND Module
 *)

module Bandmatrix :
  sig
    (**
    This type is essentially a [DlsMat] returned from a call to
    {!new_band_mat}.

     @see <CVODE_DOC_ROOT(node9#s:dls)>  Type DlsMat
     @see <CVODE_DOC_ROOT(node9#ss:band)> NewBandMat 
     *)
    type t

    (**
     TODO: n, mu, ml, smu
     @see <CVODE_DOC_ROOT(node9#ss:band)> NewBandMat
     *)
    val new_band_mat : int * int * int * int -> t

    (**
     @see <CVODE_DOC_ROOT(node9#ss:band)> PrintMat
     *)
    val print_mat : t -> unit


    (**
     @see <CVODE_DOC_ROOT(node9#ss:band)> SetToZero
     *)
    val set_to_zero    : t -> unit

    (**
     @see <CVODE_DOC_ROOT(node9#ss:band)> AddIdentity
     *)
    val add_identity   : t -> unit

    (**
     @see <CVODE_DOC_ROOT(node9#ss:band)> BandCopy
     *)
    val copy : t -> t -> int -> int -> unit

    (**
     @see <CVODE_DOC_ROOT(node9#ss:band)> BandScale
     *)
    val scale : float -> t -> unit

    (**
     @see <CVODE_DOC_ROOT(node9#ss:band)> BandGBTRF
     *)
    val gbtrf : t -> int_array -> unit

    (**
     @see <CVODE_DOC_ROOT(node9#ss:band)> BandGBTRS
     *)
    val gbtrs : t -> int_array -> real_array -> unit

    (**
     @see <CVODE_DOC_ROOT(node9#s:dls)> BAND_ELEM
     *)
    val get : t -> (int * int) -> float

    (**
     @see <CVODE_DOC_ROOT(node9#s:dls)> BAND_ELEM
     *)
    val set : t -> (int * int) -> float -> unit

    module Col :
      sig
        (**
         TODO: explain this type
         @see <CVODE_DOC_ROOT(node9#s:dls)> BAND_COL
         *)
        type c

        (**
         @see <CVODE_DOC_ROOT(node9#s:dls)> BAND_COL
         *)
        val get_col : t -> int -> c

        (**
         @see <CVODE_DOC_ROOT(node9#s:dls)> BAND_COL_ELEM
         *)
        val get : c -> int -> int -> float

        (**
         @see <CVODE_DOC_ROOT(node9#s:dls)> BAND_COL_ELEM
         *)
        val set : c -> int -> int -> float -> unit
      end

    module Direct :
      sig
        (**
        This type is essentially a [realtype **] returned from a call to
        {!new_band_mat}.

         @see <CVODE_DOC_ROOT(node9#ss:band)> Small band matrices
         @see <CVODE_DOC_ROOT(node9#ss:band)> NewBandMat 
         *)
        type t

        (**
         TODO: n smu ml
         @see <CVODE_DOC_ROOT(node9#ss:band)> newBandMat
         *)
        val new_band_mat : int * int * int -> t

        (**
         *)
        val get : t -> (int * int) -> float

        (**
         *)
        val set : t -> (int * int) -> float -> unit

        (**
         TODO:  a    b    n     a_smu  b_smu  copymu  copyml
         @see <CVODE_DOC_ROOT(node9#ss:band)> bandCopy
         *)
        val copy : t -> t -> int -> int -> int -> int -> int -> unit

        (**
         TODO:  c         a    n      mu     ml     smu
         @see <CVODE_DOC_ROOT(node9#ss:band)> bandScale
         *)
        val scale : float -> t -> int -> int -> int -> int -> unit

        (**
         TODO: a    n      smu
         @see <CVODE_DOC_ROOT(node9#ss:band)> bandAddIdentity
         *)
        val add_identity : t -> int -> int -> unit

        (**
         TODO: a    n      mu     ml     smu    p
         @see <CVODE_DOC_ROOT(node9#ss:band)> bandGBTRF
         *)
        val gbtrf : t -> int -> int -> int -> int -> int_array -> unit

        (**
         TODO: a    n      smu    ml     p            b
         @see <CVODE_DOC_ROOT(node9#ss:band)> bandGBTRS
         *)
        val gbtrs
            : t -> int -> int -> int -> int_array -> real_array -> unit
      end
  end

