external crash : string -> 'a = "sunml_crash"
module Vptr :
  sig type 'a vptr val make : 'a -> 'a vptr val unwrap : 'a vptr -> 'a end
module Callback :
  sig
    type 'f cfunptr
    type 'f cfun = { cptr : 'f cfunptr; call : 'f; }
    val invoke : 'a cfun -> 'a
  end
module Version :
  sig
    val in_compat_mode2 : bool
    val in_compat_mode2_3 : bool
    val lt400 : bool
    val lt500 : bool
    val lt530 : bool
    val lt540 : bool
    val lt580 : bool
    val lt600 : bool
    val lt620 : bool
    val lt640 : bool
    val lt660 : bool
    val has_nvector_get_id : bool
  end
module Logfile :
  sig
    type t
    external c_stderr : unit -> t = "sunml_sundials_stderr"
    external c_stdout : unit -> t = "sunml_sundials_stdout"
    external fopen : string -> bool -> t = "sunml_sundials_fopen"
    val stderr : t
    val stdout : t
    val openfile : ?trunc:bool -> string -> t
    external output_string : t -> string -> unit = "sunml_sundials_write"
    external output_bytes : t -> bytes -> unit = "sunml_sundials_write"
    external flush : t -> unit = "sunml_sundials_fflush"
    external close : t -> unit = "sunml_sundials_close"
  end
module Profiler :
  sig type t external make : string -> t = "sunml_profiler_make" end
module Logger : sig type t end
module Context :
  sig
    type cptr
    type t = {
      cptr : cptr;
      mutable profiler : Profiler.t option;
      mutable logger : Logger.t;
    }
    exception ExternalProfilerInUse
    external c_make : unit -> cptr * Logger.t = "sunml_context_make"
    external c_set_profiler : cptr -> Profiler.t -> unit
      = "sunml_context_set_profiler"
    val set_profiler : t -> Profiler.t -> unit
    external c_set_logger : cptr -> Logger.t -> unit
      = "sunml_context_set_logger"
    val set_logger : t -> Logger.t -> unit
    val get_logger : t -> Logger.t
    val make : ?profiler:Profiler.t -> ?logger:Logger.t -> unit -> t
    val default_context : t Weak.t
    val default : unit -> t
    val get : t option -> t
    val get_profiler : t -> Profiler.t
  end
