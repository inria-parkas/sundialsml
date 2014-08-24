
module type NVECTOR_OPS =
  sig
    type t

    val n_vclone        : t -> t
    val n_vlinearsum    : float -> t -> float -> t -> t -> unit
    val n_vconst        : float -> t -> unit
    val n_vprod         : t -> t -> t -> unit
    val n_vdiv          : t -> t -> t -> unit
    val n_vscale        : float -> t -> t -> unit
    val n_vabs          : t -> t -> unit
    val n_vinv          : t -> t -> unit
    val n_vaddconst     : t -> float -> t -> unit
    val n_vdotprod      : t -> t -> float
    val n_vmaxnorm      : t -> float
    val n_vwrmsnorm     : t -> t -> float
    val n_vwrmsnormmask : t -> t -> t -> float
    val n_vmin          : t -> float
    val n_vwl2norm      : t -> t -> float
    val n_vl1norm       : t -> float
    val n_vcompare      : float -> t -> t -> unit
    val n_vinvtest      : t -> t -> bool
    val n_vconstrmask   : t -> t -> t -> bool
    val n_vminquotient  : t -> t -> float
  end

