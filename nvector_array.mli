
module type ARRAY_NVECTOR =
  sig
    type t

    val array_nvec_ops  : t Nvector.Mutable.nvector_ops
    val make            : int -> float -> t Nvector.nvector
    val wrap            : t -> t Nvector.nvector
    val data            : t Nvector.nvector -> t
  end

include ARRAY_NVECTOR with type t = float array

module Bigarray :
  ARRAY_NVECTOR
  with
    type t = (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t

