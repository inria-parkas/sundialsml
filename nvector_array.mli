
module type ARRAY_NVECTOR =
  sig
    type t

    val array_nvec_ops  : t Nvector.Mutable.nvector_ops
    val make            : int -> float -> t Nvector.nvector
    val wrap            : t -> t Nvector.nvector
    val data            : t Nvector.nvector -> t
  end

module Bigarray :
  ARRAY_NVECTOR
  with
    type t = (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t

include ARRAY_NVECTOR with type t = float array

(*
val array_nvec_ops  : float array Nvector.Mutable.nvector_ops
val make            : int -> float -> float array Nvector.nvector
val wrap            : float array -> float array Nvector.nvector
val data            : 'a Nvector.nvector -> 'a
*)

