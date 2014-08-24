
type kind
type t = (Sundials.RealArray.t, kind) Sundials.nvector

external wrap : Sundials.RealArray.t -> t
  = "ml_nvec_wrap_serial"

let make n iv = wrap (Sundials.RealArray.make n iv)

module Ops = struct
  type t = (Sundials.RealArray.t, kind) Sundials.nvector

  let n_vclone nv =
    let data = Sundials.unvec nv in
    wrap (Sundials.RealArray.clone data)

  external n_vlinearsum    : float -> t -> float -> t -> t -> unit
    = "ml_nvec_ser_n_vlinearsum"

  external n_vconst        : float -> t -> unit
    = "ml_nvec_ser_n_vconst"

  external n_vprod         : t -> t -> t -> unit
    = "ml_nvec_ser_n_vprod"

  external n_vdiv          : t -> t -> t -> unit
    = "ml_nvec_ser_n_vdiv"

  external n_vscale        : float -> t -> t -> unit
    = "ml_nvec_ser_n_vscale"

  external n_vabs          : t -> t -> unit
    = "ml_nvec_ser_n_vabs"

  external n_vinv          : t -> t -> unit
    = "ml_nvec_ser_n_vinv"

  external n_vaddconst     : t -> float -> t -> unit
    = "ml_nvec_ser_n_vaddconst"

  external n_vdotprod      : t -> t -> float
    = "ml_nvec_ser_n_vdotprod"

  external n_vmaxnorm      : t -> float
    = "ml_nvec_ser_n_vmaxnorm"

  external n_vwrmsnorm     : t -> t -> float
    = "ml_nvec_ser_n_vwrmsnorm"

  external n_vwrmsnormmask : t -> t -> t -> float
    = "ml_nvec_ser_n_vwrmsnormmask"

  external n_vmin          : t -> float
    = "ml_nvec_ser_n_vmin"

  external n_vwl2norm      : t -> t -> float
    = "ml_nvec_ser_n_vwl2norm"

  external n_vl1norm       : t -> float
    = "ml_nvec_ser_n_vl1norm"

  external n_vcompare      : float -> t -> t -> unit
    = "ml_nvec_ser_n_vcompare"

  external n_vinvtest      : t -> t -> bool
    = "ml_nvec_ser_n_vinvtest"

  external n_vconstrmask   : t -> t -> t -> bool
    = "ml_nvec_ser_n_vconstrmask"

  external n_vminquotient  : t -> t -> float
    = "ml_nvec_ser_n_vminquotient"

end

module ArrayOps = Nvector_array.MakeOps (struct
    type data = Sundials.RealArray.t

    let get       = Bigarray.Array1.get
    let set       = Bigarray.Array1.set
    let fill      = Bigarray.Array1.fill

    let make      = Sundials.RealArray.make
    let length    = Sundials.RealArray.length
    let clone     = Sundials.RealArray.clone
    let fold_left = Sundials.RealArray.fold_left
  end)
module DataOps = ArrayOps.DataOps

