
type kind
type t = (Sundials.RealArray.t, kind) Sundials.nvector

external wrap : Sundials.RealArray.t -> t
  = "ml_nvec_wrap_serial"

let make n iv = wrap (Sundials.RealArray.init n iv)

