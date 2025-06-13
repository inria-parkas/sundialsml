(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

(** Generic definitions, arrays, matrices, linear solvers, nonlinear solvers,
    and utility functions.

 @version VERSION()
 @author Timothy Bourke (Inria/ENS)
 @author Jun Inoue (Inria/ENS)
 @author Marc Pouzet (UPMC/ENS/Inria)
 *)

(** Installation specific constants. *)
module Config = Sundials_Config

(** Index values for sparse matrices. *)
module Index = Sundials_Index

(** A rudimentary interface to C streams for logging in Sundials.

    Files for error and diagnostic information. File values are passed
    to functions like {!Cvode.set_error_file} and {!Kinsol.set_info_file} to
    log solver errors and diagnostics. *)
module Logfile : sig (* {{{ *)

  (** An open log file. *)
  type t = Sundials_impl.Logfile.t

  (** The stderr file.
      This log file has its own buffer, distinct from that of
      [Stdlib.stderr]; take care to flush before and after. *)
  val stderr : t

  (** The stdout file.
      This log file has its own buffer, distinct from that of
      [Stdlib.stdout]; take care to flush before and after. *)
  val stdout : t

  (** Opens the named file. When [trunc] is false, the default, writes are
      appended to the file. When [trunc] is true, the opened file is
      truncated to zero length. Files are closed on garbage collection. *)
  val openfile : ?trunc:bool -> string -> t

  (** Writes the given string to an open log file. *)
  val output_string : t -> string -> unit

  (** Writes the given byte sequence to an open log file. *)
  val output_bytes : t -> bytes -> unit

  (** Flushes the given file. *)
  val flush : t -> unit

  (** Closes the given file. *)
  val close : t -> unit

end (* }}} *)

(** Performance profiling

    The underlying Sundials library must be built with
    [SUNDIALS_BUILD_WITH_PROFILING] set to on. Profiling is
    light-weight but can still reduce performance.

    The [SUNPROFILER_PRINT] environment variable determines whether
    profiler information is printed when a context is freed (by the garbage
    collector).

    The profiling functions (silently) do nothing when profiling is not
    available.

    @profiler SUNProfiler
    @since 6.0.0 *)
module Profiler : sig (* {{{ *)

  (** A Sundials profiler.

      @profiler SUNProfiler *)
  type t = Sundials_impl.Profiler.t

  (** Indicates whether the underlying library was built with profiling
      enabled. *)
  val enabled : bool

  (** Creates a new profiler with the given name.

      @profiler SUNProfiler_Create *)
  val make : string -> t

  (** Starts timing the region indicated by the given name.

      @profiler SUNProfiler_Begin *)
  external start : t -> string -> unit = "sunml_profiler_begin" [@@noalloc]

  (** Ends timing the region indicated by the given name.

      @profiler SUNProfiler_End *)
  external finish : t -> string -> unit = "sunml_profiler_end" [@@noalloc]

  (** Prints out a profiling summary.

      @profiler SUNProfiler_Print *)
  val print : t -> Logfile.t -> unit

  (** Reset the region timings and counters to zero.

      @profiler SUNProfiler_Reset
      @since 6.2.0 *)
  val reset : t -> unit

end (* }}} *)

(** Status logging

    The underlying Sundials library must be built with
    [SUNDIALS_LOGGING_LEVEL] greater than zero.

    The logging functions (silently) do nothing when logging is not
    available.

    @logger SUNLogger
    @since 6.2.0 *)
module Logger : sig (* {{{ *)

  (** A Sundials logger.

      @logger SUNLogger *)
  type t = Sundials_impl.Logger.t

  (** Identifies the logging level to use for a given message
      ({cconst SUNLogLevel}) or that was set when the underlying library was
      compiled.

     @logger SUNLogLevel
     @logger <#enabling-logging> SUNDIALS_LOGGING_LEVEL *)
  type level =
    | Error     (** Error-level logging messages
        {cconst SUN_LOGLEVEL_ERROR}/{cconst SUNDIALS_LOGGING_ERROR} *)
    | Warning   (** Warning-level logging messages
        {cconst SUN_LOGLEVEL_WARNING}/{cconst SUNDIALS_LOGGING_WARNING} *)
    | Info      (** Info-level logging messages
        {cconst SUN_LOGLEVEL_INFO}/{cconst SUNDIALS_LOGGING_INFO} *)
    | Debug     (** Debug-level logging messages
        {cconst SUN_LOGLEVEL_DEBUG}/{cconst SUNDIALS_LOGGING_DEBUG} *)

  (** Returns true if the first level implies the second. For instance,
      [level_implies Warning Error = true]. *)
  val level_implies : level -> level -> bool

  (** Indicates what level of debugging was specified when the underlying
      library was compiled. Returns [None] if logging is not enabled.

     @logger <#enabling-logging> SUNDIALS_LOGGING_LEVEL *)
  val logging_level : level option

  (** Creates a new logger.

      @logger SUNLogger_Create
      @logger SUNLogger_SetErrorFilename
      @logger SUNLogger_SetWarningFilename
      @logger SUNLogger_SetInfoFilename
      @logger SUNLogger_SetDebugFilename *)
  val make :
       ?error_filename:string
    -> ?warning_filename:string
    -> ?info_filename:string
    -> ?debug_filename:string
    -> unit
    -> t

  (** Creates a new logger and opens the output streams/files from
      environment variables. The environment variables are:
      {cconst SUNLOGGER_ERROR_FILENAME}
      {cconst SUNLOGGER_WARNING_FILENAME}
      {cconst SUNLOGGER_INFO_FILENAME}
      {cconst SUNLOGGER_DEBUG_FILENAME}.

      @logger SUNLogger_CreateFromEnv *)
  val make_from_env : unit -> t

  (** Sets the filename for error output.

      @logger SUNLogger_SetErrorFilename *)
  val set_error_filename : t -> string -> unit

  (** Sets the filename for warning output.

      @logger SUNLogger_SetWarningFilename *)
  val set_warning_filename : t -> string -> unit

  (** Sets the filename for info output.

      @logger SUNLogger_SetInfoFilename *)
  val set_info_filename : t -> string -> unit

  (** Sets the filename for debug output.

      @logger SUNLogger_SetDebugFilename *)
  val set_debug_filename : t -> string -> unit

  (** Queues a message to the output log level.

      @logger SUNLogger_QueueMsg *)
  val queue_msg : t -> level -> scope:string -> label:string -> string -> unit

  (** Flush the message queue(s). If [level] is not explicitly given, then all
      queues are flushed ({cconst SUN_LOGLEVEL_ALL}).

      @logger SUNLogger_Flush *)
  val flush : ?level:level -> t -> unit

end (* }}} *)

(** Contexts for creating Sundials values

    Every function that creates a Sundials value (integrator, nvector,
    linear solver, nonlinear solver, etcetera) does so within a context.
    All such functions have an optional [?context] argument. When this
    argument is not given explicitly, it defaults to the context returned
    by {!default}.

    @context SUNContext
    @since 6.0.0 *)
module Context : sig (* {{{ *)

  (** A context required to create Sundials values.

      @context SUNContext *)
  type t = Sundials_impl.Context.t

  (** The default context when creating values.

      @context SUNContext_Create *)
  val default : unit -> t

  (** Create a new context, optionally specifying the profiler to use.

      @context SUNContext_Create
      @context SUNContext_SetProfiler
      @context SUNContext_SetLogger *)
  val make : ?profiler:Profiler.t -> ?logger:Logger.t -> unit -> t

  (** Indicates that an external library (i.e., caliper) is being use for
      profiling. *)
  exception ExternalProfilerInUse

  (** Return the profiler associated with a context.

      @context SUNContext_GetProfiler
      @raise ExternalProfilerInUse If an external library is used for profiling *)
  val get_profiler : t -> Profiler.t

  (** Sets the profiler associated with a context.

      @context SUNContext_SetProfiler *)
  val set_profiler : t -> Profiler.t -> unit

  (** Return the logger associated with a context.

      @context SUNContext_GetLogger *)
  val get_logger : t -> Logger.t

  (** Sets the logger associated with a context.

      @context SUNContext_SetLogger *)
  val set_logger : t -> Logger.t -> unit

end (* }}} *)

(** Accuracy-based adaptivity controllers that estimate the step sizes, and
    other parameters, of integrators.

    @adaptcontroller <index.html> Time Step Adaptivity Controllers
    @since 6.7.0 *)
module AdaptController = Sundials_AdaptController

(** {2:exceptions Exceptions} *)

(** Indicates a recoverable failure within a callback function.
    Any other exception normally indicates an unrecoverable failure. *)
exception RecoverableFailure

(** Raised by error-weight functions on non-positive error weights. See
    {{!Cvode.tolerance}[Cvode.WFtolerances]} or
    {{!Ida.tolerance}[Ida.WFtolerances]}. *)
exception NonPositiveEwt

(** {2:values Generic values} *)

(** Specifies the output format for statistics. See
    {!Cvode.print_all_stats}, {!Ida.print_all_stats},
    {!Kinsol.print_all_stats}, {!Arkode.ARKStep.print_all_stats},
    {!Arkode.MRIStep.print_all_stats}, or {!Arkode.ERKStep.print_all_stats}. *)
type output_format =
  | OutputTable     (* A table of values ({cconst SUN_OUTPUTFORMAT_TABLE}) *)
  | OutputCSV       (* Comma-separated list of key and value pairs
                       ({cconst SUN_OUTPUTFORMAT_CSV}) *)

(** A callback function into the underlying library. *)
type 'f cfun = 'f Sundials_impl.Callback.cfun

(** Use a callback function provided by the underlying library. *)
val invoke : 'f cfun -> 'f

(** {2:arrays Arrays} *)

(** Vectors of floats (one-dimensional bigarrays). *)
module RealArray = Sundials_RealArray

(** Matrices of floats (two-dimensional bigarrays plus
    extra information for Sundials. *)
module RealArray2 = Sundials_RealArray2

(** Vectors of integers (one-dimensional bigarrays). *)
module LintArray = Sundials_LintArray

(** Read-only polymorphic arrays. *)
module ROArray = Sundials_ROArray

(** {2:roots Arrays of roots (zero-crossings)} *)

(** Vectors of root (zero-crossing) statuses. *)
module Roots : sig (* {{{ *)
  (** Arrays that communicate the occurrence of zero-crossings. The
      underlying representation is hidden to isolate compatability
      issues related to integers. *)
  type t

  (** Values indicating the status of root functions.
      @cvode CVodeGetRootInfo
      @ida IdaGetRootInfo *)
  type r =
    | NoRoot      (** No root was found on the corresponding function (0). *)
    | Rising      (** The corresponding root function is increasing (1). *)
    | Falling     (** The corresponding root function is decreasing (-1). *)

  (** [create n] returns an array with [n] elements each set to
      {{!r}NoRoot}. *)
  val create : int -> t

  (** [make n x] returns an array with [n] elements each set to [x]. *)
  val make : int -> r -> t

  (** [init n f] returns an array with [n] elements, with element [i] set
      to [f i]. *)
  val init : int -> (int -> r) -> t

  (** Returns the length of an array. *)
  val length : t -> int

  (** Pretty-print a root array using the
      {{:OCAML_DOC_ROOT(Format.html)} Format} module. *)
  val pp : Format.formatter -> t -> unit

  (** Pretty-print a root array using the
      {{:OCAML_DOC_ROOT(Format.html)} Format} module.
      The defaults are: [start="\["], [stop="\]"], [sep="; "], and
      [item] prints '_' for {{!r}NoRoot}, 'R' for {{!r}Rising}, and 'F' for
      {{!r}Falling}. *)
  val ppi : ?start:string -> ?stop:string -> ?sep:string
            -> ?item:(Format.formatter -> int -> r -> unit)
            -> unit
            -> Format.formatter -> t -> unit

  (** Returns [true] only if the specified element is either {{!r}Rising} or
      {{!r}Falling}. *)
  val detected : t -> int -> bool

  (** Returns [true] only if the specified element is {{!r}Rising}. *)
  val rising : t -> int -> bool

  (** Returns [true] only if the specified element is {{!r}Falling}. *)
  val falling : t -> int -> bool

  (** [get r i] returns the [i]th element of [r]. *)
  val get : t -> int -> r

  (** [set r i v] sets the [i]th element of [r] to [v]. *)
  val set : t -> int -> r -> unit

  (** [set_noroot r i] sets the [i]th element of [r] to {{!r}NoRoot}. *)
  val set_noroot : t -> int -> unit

  (** [set_rising r i] sets the [i]th element of [r] to {{!r}Rising}. *)
  val set_rising : t -> int -> unit

  (** [set_falling r i] sets the [i]th element of [r] to {{!r}Falling}. *)
  val set_falling : t -> int -> unit

  (** [fill a x] sets all elements in [a] to [x]. *)
  val fill : t -> r -> unit

  (** Creates a new array with the same contents as an existing one. *)
  val copy : t -> t

  (** Returns [0] for {{!r}NoRoot}, [1] for {{!r}Rising}, and [-1] for
      {{!r}Falling}. *)
  val int_of_root : r -> int

  (** Resets all elements to {{!r}NoRoot}. *)
  val reset : t -> unit

  (** [true] if any elements are equal to {{!r}Rising} or {{!r}Falling}. *)
  val exists : t -> bool

  (** [iter f r] successively applies [f] to each element in [r]. *)
  val iter : (r -> unit) -> t -> unit

  (** [iteri f r] successively applies [f] to the indexes and
      elements of [r]. *)
  val iteri : (int -> r -> unit) -> t -> unit

  (** Creates an array by copying the contents of a [r list]. *)
  val of_list : r list -> t

  (** Copies into a list. *)
  val to_list : t -> r list

  (** Creates a new value from the contents of an
      {{:OCAML_DOC_ROOT(Array.html)} array}. *)
  val of_array : r array -> t

  (** Creates a new array from the contents of a given value. *)
  val to_array : t -> r array
end (* }}} *)

(** Vectors of root (zero-crossing) directions. *)
module RootDirs : sig (* {{{ *)
  type t
  (** Arrays that communicate which zero-crossings are sought. The
      underlying representation is hidden to isolate compatability
      issues related to integers. *)

  (** Values indicating which types of roots are sought.

      @cvode CVodeSetRootDirection
      @ida IdaSetRootDirection *)
  type d =
    | Increasing              (** Only look for rising zero-crossings. *)
    | Decreasing              (** Only look for falling zero-crossings. *)
    | IncreasingOrDecreasing  (** Look for any zero-crossing. *)

  (** [make n x] returns an array with [n] elements each set to [x]. *)
  val make : int -> d -> t

  (** [create n] returns an array with [n] elements each set to
      {{!d}IncreasingOrDecreasing}. *)
  val create : int -> t

  (** [init n f] returns an array with [n] elements, with element [i] set
      to [f i]. *)
  val init : int -> (int -> d) -> t

  (** Pretty-print a root direction array using the
      {{:OCAML_DOC_ROOT(Format.html)} Format} module. *)
  val pp : Format.formatter -> t -> unit

  (** Pretty-print a root direction array using the
      {{:OCAML_DOC_ROOT(Format.html)} Format} module.
      The defaults are: [start="\["], [stop="\]"], [sep="; "], and
      [item] prints 'R' for {{!d}Increasing}, 'F' for {{!d}Decreasing}, and
      'E' (either) for {{!d}IncreasingOrDecreasing}. *)
  val ppi : ?start:string -> ?stop:string -> ?sep:string
            -> ?item:(Format.formatter -> int -> d -> unit)
            -> unit
            -> Format.formatter -> t -> unit

  (** [copy n a] returns an array with [n] elements, initialized from
      the contents of a. If [n > Array.length a] then the extra space is
      initialized to {{!d}IncreasingOrDecreasing}. *)
  val copy : int -> d array -> t

  (** Returns the length of an array *)
  val length : t -> int

  (** [get r i] returns the [i]th element of [r]. *)
  val get : t -> int -> d

  (** [set r i v] sets the [i]th element of [r] to [v]. *)
  val set : t -> int -> d -> unit

  (** [fill_all a x] sets the values of [a] to [x] everywhere. *)
  val fill : t -> d -> unit

  (** [blitn ~src ?spos ~dst ?dpos len] copies [len] elements of [src] at
      offset [spos] to [dst] at offset [dpos].
      The [spos] and [dpos] arguments are optional and default to zero.

      @raise Invalid_argument "RootDirs.blitn" if [spos], [dpos], and
      [len] do not specify valid subarrays of [src] and [dst]. *)
  val blitn : src:t -> ?spos:int -> dst:t -> ?dpos:int -> int -> unit

  (** Copy the first array into the second one.
      See {{:OCAML_DOC_ROOT(Bigarray.Genarray.html#VALblit)}
      [Bigarray.Genarray.blit]} for more details. *)
  val blit : src:t -> dst:t -> unit

  (** Creates an array by copying the contents of a [d list]. *)
  val of_list : d list -> t

  (** Copies into a list. *)
  val to_list : t -> d list

  (** Creates a new value from the contents of an
      {{:OCAML_DOC_ROOT(Array.html)} array}. *)
  val of_array : d array -> t

  (** Creates a new array from the contents of a given value. *)
  val to_array : t -> d array
end (* }}} *)

(** {2:constraints Constraints} *)

(** Symbolic names for variable constraints. These names describe
    the constants passed to {!Ida.set_constraints} and
    {!Kinsol.set_constraints}.

    @ida IDASetConstraints
    @kinsol KINSetConstraints *)
module Constraint : sig (* {{{ *)
  (** The constant [0.0]. *)
  val unconstrained : float

  (** The constant [1.0]. *)
  val geq_zero : float

  (** The constant [-1.0]. *)
  val leq_zero : float

  (** The constant [2.0]. *)
  val gt_zero : float

  (** The constant [-2.0]. *)
  val lt_zero : float

  (** For pattern-matching on constraints. See {!of_float}. *)
  type t =
  | Unconstrained   (** [true] *)
  | GeqZero         (** [>= 0] *)
  | LeqZero         (** [<= 0] *)
  | GtZero          (** [> 0]  *)
  | LtZero          (** [< 0]  *)

  (** Map constraint values to floating-point constants. *)
  val to_float : t -> float

  (** Map floating-point constants to constraint values.

      @raise Invalid_argument The given value is not a legal constraint. *)
  val of_float : float -> t
end (* }}} *)

(** {2:matlin Matrices, Linear Solvers, and Nonlinear Solvers} *)

(** Generic matrices.

    @matrix <index.html#matrix-data-structures> Matrix Data Structures
    @since 3.0.0 *)
module Matrix = Sundials_Matrix

(** Generic linear solvers.

    Sundials provides a set of functions for instantiating linear solvers from
    two families: {!module:Sundials_LinearSolver.Direct} and
    {!module:Sundials_LinearSolver.Iterative}. Any instance may
    be associated with at most one solver session.

    @linsol <SUNLinSol_API_link.html#the-sunlinearsolver-api> The SUNLinearSolver API
    @linsol <index.html#linear-algebraic-solvers> Linear Algebraic Solvers
    @since 3.0.0 *)
module LinearSolver = Sundials_LinearSolver

(** Generic nonlinear solvers.

    Sundials provides generic nonlinear solvers of two main types:
    {!module:Sundials_NonlinearSolver.Newton} and
    {!module:Sundials_NonlinearSolver.FixedPoint}. An instance of a nonlinear
    solver may only be associated with at most one integrator session at
    a time.

    This module supports calling both Sundials and custom OCaml nonlinear
    solvers from both Sundials integrators and OCaml applications.

    @nonlinsol <SUNNonlinSol_API_link.html#the-sunnonlinearsolver-api> The SUNNonlinearSolver API
    @since 4.0.0 *)
module NonlinearSolver = Sundials_NonlinearSolver

(** Shared definitions and miscellaneous utility functions. *)
module Util : sig (* {{{ *)

  (** Information passed to registered error handler functions.
      See {!Cvode.set_err_handler_fn}, {!Ida.set_err_handler_fn}, and
      {!Kinsol.set_err_handler_fn}.

      @cvode CVodeErrHandlerFn
      @ida IDAErrHandlerFn
      @kinsol KINErrHandlerFn *)
  type error_details = {
      error_code : int;
      module_name : string;        (** IDA, CVODE, CVSPGMR, etc. *)
      function_name : string;
      error_message : string;
    }

  (** {2:misc Miscellaneous utility functions} *)

  (** [format_float fmt f] formats [f] according to the format string [fmt].
      It uses the low-level [caml_format_float] function. *)
  val format_float : string -> float -> string

  (** Returns the bit-level representation of a float in hexadecimal as a string.
      Equivalent to [format_float "%a"]. *)
  val floata : float -> string

  (** Returns true if the relative difference of the two arguments is less
      than or equal to the tolerance.
      The tolerance defaults to 10 times {!Sundials_Config.unit_roundoff}.
      This function handles the cases where the arguments are near zero and
      where either is {{:OCAML_DOC_ROOT(Float.html)} Float.infinity}
      or {{:OCAML_DOC_ROOT(Float.html)} Float.nan}.

      @nodoc SUNRCompare
      @nodoc SUNRCompareTol
      @since 5.8.0 *)
  val compare_float : ?tol:float -> float -> float -> bool

  (** The (inclusive) upper bound on values returned from {!rand}. *)
  val rand_max : int

  (** Returns a random integer between 0 and {!rand_max} inclusive using
      the standard C library. *)
  external rand : unit -> int = "sunml_rand"

  (** Resets the seed used to generate the sequence of numbers returned by
      {!rand}. *)
  external srand : int -> unit = "sunml_srand"

  (* Returns [None] if POSIX timers are not available, otherwise returns the
     timer resolution [(n, f)], meaning that [n] nanoseconds are represented
     by the floating-point value [f]. *)
  external get_time_precision : unit -> (int * float) option
    = "sunml_get_timing_precision"

  (* Returns the current value of a monotonic timer in the precision given by
     {!get_time_precision}, or 0 if POSIX timers are not available. Allows
     measuring an interval by taking the difference of two monotonic times. *)
  external get_monotonic_time : unit -> (float [@unboxed])
    = "sunml_get_time_byte" "sunml_get_time"

end (* }}} *)

