
type kind
type data = Sundials.RealArray.t * int * Mpi.communicator
type t = (data, kind) Sundials.nvector

exception IncorrectGlobalSize

let _ = Callback.register_exception
          "nvector_parallel_IncorrectGlobalSize" IncorrectGlobalSize

external wrap : (Sundials.RealArray.t * int * Mpi.communicator) -> t
  = "ml_nvec_wrap_parallel"

let make nl ng comm iv = wrap (Sundials.RealArray.init nl iv, ng, comm)

let unwrap nv =
  let data, _, _ = Sundials.unvec nv in
  data

let global_length nv =
  let _, gl, _ = Sundials.unvec nv in
  gl

let communicator nv =
  let _, _, comm = Sundials.unvec nv in
  comm

module Ops = struct
  let n_vclone nv =
    let loc, glen, comm = Sundials.unvec nv in
    wrap (Sundials.RealArray.clone loc, glen, comm)

  external n_vlinearsum    : float -> t -> float -> t -> t -> unit
    = "ml_nvec_par_n_vlinearsum"

  external n_vconst        : float -> t -> unit
    = "ml_nvec_par_n_vconst"

  external n_vprod         : t -> t -> t -> unit
    = "ml_nvec_par_n_vprod"

  external n_vdiv          : t -> t -> t -> unit
    = "ml_nvec_par_n_vdiv"

  external n_vscale        : float -> t -> t -> unit
    = "ml_nvec_par_n_vscale"

  external n_vabs          : t -> t -> unit
    = "ml_nvec_par_n_vabs"

  external n_vinv          : t -> t -> unit
    = "ml_nvec_par_n_vinv"

  external n_vaddconst     : t -> float -> t -> unit
    = "ml_nvec_par_n_vaddconst"

  external n_vdotprod      : t -> t -> float
    = "ml_nvec_par_n_vdotprod"

  external n_vmaxnorm      : t -> float
    = "ml_nvec_par_n_vmaxnorm"

  external n_vwrmsnorm     : t -> t -> float
    = "ml_nvec_par_n_vwrmsnorm"

  external n_vwrmsnormmask : t -> t -> t -> float
    = "ml_nvec_par_n_vwrmsnormmask"

  external n_vmin          : t -> float
    = "ml_nvec_par_n_vmin"

  external n_vwl2norm      : t -> t -> float
    = "ml_nvec_par_n_vwl2norm"

  external n_vl1norm       : t -> float
    = "ml_nvec_par_n_vl1norm"

  external n_vcompare      : float -> t -> t -> unit
    = "ml_nvec_par_n_vcompare"

  external n_vinvtest      : t -> t -> bool
    = "ml_nvec_par_n_vinvtest"

  external n_vconstrmask   : t -> t -> t -> bool
    = "ml_nvec_par_n_vconstrmask"

  external n_vminquotient  : t -> t -> float
    = "ml_nvec_par_n_vminquotient"

end

