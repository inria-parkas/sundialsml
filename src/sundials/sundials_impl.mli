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
    val has_nvector_get_id : bool
  end
module Context :
  sig
    type cptr
    type t = { cptr : cptr; }
    external c_make : unit -> cptr = "sunml_context_make"
    val make : unit -> t
    val default_context : t Weak.t
    val default : unit -> t
    val get : t option -> t
  end
