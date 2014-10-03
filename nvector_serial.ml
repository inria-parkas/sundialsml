
type data = Sundials.RealArray.t
type kind
type t = (data, kind) Sundials.nvector

external wrap : Sundials.RealArray.t -> t
  = "ml_nvec_wrap_serial"

let make n iv = wrap (Sundials.RealArray.make n iv)

module Ops = struct
  type t = (Sundials.RealArray.t, kind) Sundials.nvector

  let n_vclone nv =
    let data = Nvector.unwrap nv in
    wrap (Sundials.RealArray.copy data)

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

(* (* Too slow! *)
module ArrayOps = Nvector_array.MakeOps (struct
    type data = Sundials.RealArray.t

    let get       = Bigarray.Array1.get
    let set       = Bigarray.Array1.set
    let fill      = Bigarray.Array1.fill

    let make      = Sundials.RealArray.make
    let length    = Sundials.RealArray.length
    let clone     = Sundials.RealArray.clone
  end)
module DataOps = ArrayOps.DataOps
*)

module DataOps =
  struct

    module A = Bigarray.Array1
    type t = Sundials.RealArray.t

    let n_vclone     = Sundials.RealArray.copy

    let arr_vaxpy a (x : t) (y : t) =
      if a = 1.0 then
        for i = 0 to A.dim x - 1 do
          A.set y i (A.get y i +. A.get x i)
        done
      else if a = -1.0 then
        for i = 0 to A.dim x - 1 do
          A.set y i (A.get y i -. A.get x i)
        done
      else
        for i = 0 to A.dim x - 1 do
          A.set y i (A.get y i +. a *. A.get x i)
        done

    let n_vlinearsum a (x : t) b (y : t) (z : t) =
      if b = 1.0 && z == y then
        arr_vaxpy a x y
      else if a = 1.0 && z == x then
        arr_vaxpy b y x
      else if a = 1.0 && b = 1.0 then
        for i = 0 to A.dim x - 1 do
          A.set z i (A.get x i +. A.get y i)
        done
      else if (a = 1.0 && b = -1.0) || (a = -1.0 && b == 1.0) then
        let v1, v2 = if (a = 1.0 && b = -1.0) then y, x else x, y in
        for i = 0 to A.dim v1 - 1 do
          A.set z i (A.get v1 i -. A.get v2 i)
        done
      else if a = 1.0 || b = 1.0 then
        let c, v1, v2 = if a = 1.0 then b, y, x else a, x, y in
        for i = 0 to A.dim v1 - 1 do
          A.set z i (c *. A.get v1 i +. A.get v2 i)
        done
      else if a = -1.0 || b = -1.0 then
        let c, v1, v2 = if a = -1.0 then b, y, x else a, x, y in
        for i = 0 to A.dim v1 - 1 do
          A.set z i (c *. A.get v1 i -. A.get v2 i)
        done
      else if a = b then
        for i = 0 to A.dim x - 1 do
          A.set z i (a *. (A.get x i +. A.get y i))
        done
      else if a = -.b then
        for i = 0 to A.dim x - 1 do
          A.set z i (a *. (A.get x i -. A.get y i))
        done
      else
        for i = 0 to A.dim x - 1 do
          A.set z i (a *. A.get x i +. b *. A.get y i)
        done

    let n_vconst c a = A.fill a c

    let n_vscale c (x : t) (z : t) =
      if c = 1.0 then
        for i = 0 to A.dim x - 1 do
          A.set z i (A.get x i)
        done
      else if c = -1.0 then
        for i = 0 to A.dim x - 1 do
          A.set z i (-. A.get x i)
        done
      else
        for i = 0 to A.dim x - 1 do
          A.set z i (c *. A.get x i)
        done

    let n_vaddconst (x : t) b (z : t) =
      for i = 0 to A.dim x - 1 do
        A.set z i (A.get x i +. b)
      done

    let n_vmaxnorm (x : t) =
      let max = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        let ax = abs_float (A.get x i) in
        if ax > !max then max := ax
      done;
      !max

    let n_vwrmsnorm (x : t) (w : t) =
      let a = ref 0.0 in
      let lx = A.dim x in
      for i = 0 to lx - 1 do
        a := !a +. ((A.get x i) *. (A.get w i) *. (A.get x i) *. (A.get w i))
      done;
      sqrt (!a /. float lx)

    let n_vwrmsnormmask (x : t) (w : t) (id : t) =
      let a = ref 0.0 in
      let lx = A.dim x in
      for i = 0 to lx - 1 do
        if A.get id i > 0.0 then
          a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
      done;
      sqrt (!a /. float lx)

    let n_vmin (x : t) =
      let min = ref max_float in
      for i = 0 to A.dim x - 1 do
        let xv = A.get x i in
        if xv < !min then min := xv
      done;
      !min

    let n_vdotprod (x : t) (y : t) =
      let a = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        a := !a +. (A.get x i *. A.get y i)
      done;
      !a

    let n_vcompare c (x : t) (z : t) =
      for i = 0 to A.dim x - 1 do
        A.set z i (if abs_float (A.get x i) >= c then 1.0 else 0.0)
      done

    let n_vinvtest (x : t) (z : t) =
      let r = ref true in
      for i = 0 to A.dim x do
        if A.get x i = 0.0 then r := false
        else A.set z i (1.0 /. (A.get x i))
      done;
      !r

    let n_vwl2norm (x : t) (w : t) =
      let a = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
      done;
      sqrt !a

    let n_vl1norm (x : t) =
      let a = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        a := !a +. abs_float (A.get x i)
      done;
      !a

    let n_vconstrmask c (x : t) (m : t) =
      let test = ref true in
      let check b = if b then 0.0 else (test := false; 1.0) in
      let f c x =
        match c with
        |  2.0 -> check (x >  0.0)
        |  1.0 -> check (x >= 0.0)
        | -1.0 -> check (x <= 0.0)
        | -2.0 -> check (x <  0.0)
        |  0.0 -> 0.0
        |    _ -> assert false
      in
      for i = 0 to A.dim c - 1 do
        A.set m i (f (A.get c i) (A.get x i))
      done;
      !test

    let n_vminquotient (n : t) (d : t) =
      let m = ref Sundials.big_real in
      for i = 0 to A.dim n - 1 do
        if (A.get d i) <> 0.0 then
          m := min !m (A.get n i /. A.get d i)
      done;
      !m

    let n_vprod (x : t) (y : t) (z : t) =
      for i = 0 to A.dim x - 1 do
        A.set z i (A.get x i *. A.get y i)
      done

    let n_vdiv (x : t) (y : t) (z : t) =
      for i = 0 to A.dim x - 1 do
        A.set z i (A.get x i /. A.get y i)
      done

    let n_vabs (x : t) (z : t) =
      for i = 0 to A.dim x - 1 do
        A.set z i (abs_float (A.get x i))
      done

    let n_vinv (x : t) (z : t) =
      for i = 0 to A.dim x - 1 do
        A.set z i (1.0 /. (A.get x i))
      done

  end

