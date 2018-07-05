
type data = RealArray.t
type kind = [`OpenMP|Nvector_serial.kind]
type t = (data, kind) Nvector.t

external c_wrap : int -> RealArray.t
                    -> (RealArray.t -> bool) -> t
  = "ml_nvec_wrap_openmp"

let wrap nthreads v =
  let len = RealArray.length v in
  c_wrap nthreads v (fun v' -> len = RealArray.length v')

let unwrap = Nvector.unwrap

let pp fmt v = RealArray.pp fmt (unwrap v)

let make nthreads n iv = wrap nthreads (RealArray.make n iv)

external num_threads : t -> int
  = "ml_nvec_openmp_num_threads"

module Ops = struct
  type t = (RealArray.t, kind) Nvector.t

  let n_vclone nv =
    let data = Nvector.unwrap nv in
    wrap (num_threads nv) (RealArray.copy data)

  external n_vlinearsum    : float -> t -> float -> t -> t -> unit
    = "ml_nvec_openmp_n_vlinearsum"

  external n_vconst        : float -> t -> unit
    = "ml_nvec_openmp_n_vconst"

  external n_vprod         : t -> t -> t -> unit
    = "ml_nvec_openmp_n_vprod"

  external n_vdiv          : t -> t -> t -> unit
    = "ml_nvec_openmp_n_vdiv"

  external n_vscale        : float -> t -> t -> unit
    = "ml_nvec_openmp_n_vscale"

  external n_vabs          : t -> t -> unit
    = "ml_nvec_openmp_n_vabs"

  external n_vinv          : t -> t -> unit
    = "ml_nvec_openmp_n_vinv"

  external n_vaddconst     : t -> float -> t -> unit
    = "ml_nvec_openmp_n_vaddconst"

  external n_vdotprod      : t -> t -> float
    = "ml_nvec_openmp_n_vdotprod"

  external n_vmaxnorm      : t -> float
    = "ml_nvec_openmp_n_vmaxnorm"

  external n_vwrmsnorm     : t -> t -> float
    = "ml_nvec_openmp_n_vwrmsnorm"

  external n_vwrmsnormmask : t -> t -> t -> float
    = "ml_nvec_openmp_n_vwrmsnormmask"

  external n_vmin          : t -> float
    = "ml_nvec_openmp_n_vmin"

  external n_vwl2norm      : t -> t -> float
    = "ml_nvec_openmp_n_vwl2norm"

  external n_vl1norm       : t -> float
    = "ml_nvec_openmp_n_vl1norm"

  external n_vcompare      : float -> t -> t -> unit
    = "ml_nvec_openmp_n_vcompare"

  external n_vinvtest      : t -> t -> bool
    = "ml_nvec_openmp_n_vinvtest"

  external n_vconstrmask   : t -> t -> t -> bool
    = "ml_nvec_openmp_n_vconstrmask"

  external n_vminquotient  : t -> t -> float
    = "ml_nvec_openmp_n_vminquotient"

end

