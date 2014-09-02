
type kind
type data = Sundials.RealArray.t * int * Mpi.communicator
type t = (data, kind) Sundials.nvector

exception IncorrectGlobalSize

let _ = Callback.register_exception
          "nvector_parallel_IncorrectGlobalSize" IncorrectGlobalSize

external wrap : (Sundials.RealArray.t * int * Mpi.communicator) -> t
  = "ml_nvec_wrap_parallel"

let make nl ng comm iv = wrap (Sundials.RealArray.make nl iv, ng, comm)

let clone nv =
  let loc, glen, comm = Sundials.unvec nv in
  wrap (Sundials.RealArray.clone loc, glen, comm)

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
  type t = (data, kind) Sundials.nvector

  let n_vclone = clone

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

module MakeOps =
  functor (A : sig
      type local_data
      val get       : local_data -> int -> float
      val set       : local_data -> int -> float -> unit
      val fill      : local_data -> float -> unit
      val make      : int -> float -> local_data
      val clone     : local_data -> local_data
      val length    : local_data -> int
    end) ->
  struct
    type t = A.local_data * int * Mpi.communicator

    (* let n_vclone (d, gl, comm) = (A.clone d, gl, comm) *)
    let n_vclone (d, gl, comm) = (d, gl, comm)

    let arr_vaxpy a x y =
      if a = 1.0 then
        for i = 0 to A.length x - 1 do
          A.set y i (A.get y i +. A.get x i)
        done
      else if a = -1.0 then
        for i = 0 to A.length x - 1 do
          A.set y i (A.get y i -. A.get x i)
        done
      else
        for i = 0 to A.length x - 1 do
          A.set y i (A.get y i +. a *. A.get x i)
        done

    let n_vlinearsum a (x, _, _) b (y, _, _) (z, _, _) =
      if b = 1.0 && z == y then
        arr_vaxpy a x y
      else if a = 1.0 && z == x then
        arr_vaxpy b y x
      else if a = 1.0 && b = 1.0 then
        for i = 0 to A.length x - 1 do
          A.set z i (A.get x i +. A.get y i)
        done
      else if (a = 1.0 && b = -1.0) || (a = -1.0 && b == 1.0) then
        let v1, v2 = if (a = 1.0 && b = -1.0) then y, x else x, y in
        for i = 0 to A.length v1 - 1 do
          A.set z i (A.get v1 i -. A.get v2 i)
        done
      else if a = 1.0 || b = 1.0 then
        let c, v1, v2 = if a = 1.0 then b, y, x else a, x, y in
        for i = 0 to A.length v1 - 1 do
          A.set z i (c *. A.get v1 i +. A.get v2 i)
        done
      else if a = -1.0 || b = -1.0 then
        let c, v1, v2 = if a = -1.0 then b, y, x else a, x, y in
        for i = 0 to A.length v1 - 1 do
          A.set z i (c *. A.get v1 i -. A.get v2 i)
        done
      else if a = b then
        for i = 0 to A.length x - 1 do
          A.set z i (a *. (A.get x i +. A.get y i))
        done
      else if a = -.b then
        for i = 0 to A.length x - 1 do
          A.set z i (a *. (A.get x i -. A.get y i))
        done
      else
        for i = 0 to A.length x - 1 do
          A.set z i (a *. A.get x i +. b *. A.get y i)
        done

    let n_vconst c (a, _, _) = A.fill a c

    let n_vscale c (x, _, _) (z, _, _) =
      if c = 1.0 then
        for i = 0 to A.length x - 1 do
          A.set z i (A.get x i)
        done
      else if c = -1.0 then
        for i = 0 to A.length x - 1 do
          A.set z i (-. A.get x i)
        done
      else
        for i = 0 to A.length x - 1 do
          A.set z i (c *. A.get x i)
        done

    let n_vaddconst (x, _, _) b (z, _, _) =
      for i = 0 to A.length x - 1 do
        A.set z i (A.get x i +. b)
      done

    let n_vmaxnorm (x, _, comm) =
      let lmax = ref 0.0 in
      for i = 0 to A.length x - 1 do
        let ax = abs_float (A.get x i) in
        if ax > !lmax then lmax := ax
      done;
      Mpi.allreduce_float !lmax Mpi.Float_max comm

    let n_vwrmsnorm (x, n_global, comm) (w, _, _) =
      let lsum = ref 0.0 in
      for i = 0 to A.length x - 1 do
        lsum := !lsum
                  +. ((A.get x i) *. (A.get w i) *. (A.get x i) *. (A.get w i))
      done;
      let gsum = Mpi.allreduce_float !lsum Mpi.Float_sum comm in
      sqrt (gsum /. float n_global)

    let n_vwrmsnormmask (x, n_global, comm) (w, _, _) (id, _, _) =
      let lsum = ref 0.0 in
      for i = 0 to A.length x - 1 do
        if A.get id i > 0.0 then
          lsum := !lsum +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
      done;
      let gsum = Mpi.allreduce_float !lsum Mpi.Float_sum comm in
      sqrt (gsum /. float n_global)

    let n_vmin (x, _, comm) =
      let lmin = ref max_float in
      for i = 0 to A.length x - 1 do
        let xv = A.get x i in
        if xv < !lmin then lmin := xv
      done;
      Mpi.allreduce_float !lmin Mpi.Float_min comm

    let n_vdotprod (x, _, comm) (y, _, _) =
      let lsum = ref 0.0 in
      for i = 0 to A.length x - 1 do
        lsum := !lsum +. (A.get x i *. A.get y i)
      done;
      Mpi.allreduce_float !lsum Mpi.Float_sum comm

    let n_vcompare c (x, _, _) (z, _, _) =
      for i = 0 to A.length x - 1 do
        A.set z i (if abs_float (A.get x i) >= c then 1.0 else 0.0)
      done

    let n_vinvtest (x, _, comm) (z, _, _) =
      let lval = ref 1.0 in
      for i = 0 to A.length x do
        if A.get x i = 0.0 then lval := 0.0
        else A.set z i (1.0 /. (A.get x i))
      done;
      (Mpi.allreduce_float !lval Mpi.Float_min comm = 0.0)

    let n_vwl2norm (x, _, comm) (w, _, _) =
      let lsum = ref 0.0 in
      for i = 0 to A.length x - 1 do
        lsum := !lsum +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
      done;
      let gsum = Mpi.allreduce_float !lsum Mpi.Float_sum comm in
      sqrt gsum

    let n_vl1norm (x, _, comm) =
      let lsum = ref 0.0 in
      for i = 0 to A.length x - 1 do
        lsum := !lsum +. abs_float (A.get x i)
      done;
      Mpi.allreduce_float !lsum Mpi.Float_sum comm

    let n_vconstrmask (c, _, _) (x, _, comm) (m, _, _) =
      let test = ref 1.0 in
      let check b = if b then 0.0 else (test := 0.0; 1.0) in
      let f c x =
        match c with
        |  2.0 -> check (x >  0.0)
        |  1.0 -> check (x >= 0.0)
        | -1.0 -> check (x <= 0.0)
        | -2.0 -> check (x <  0.0)
        |  0.0 -> 0.0
        |    _ -> assert false
      in
      for i = 0 to A.length c - 1 do
        A.set m i (f (A.get c i) (A.get x i))
      done;
      (Mpi.allreduce_float !test Mpi.Float_min comm = 1.0)

    let n_vminquotient (num, _, comm) (denom, _, _) =
      let lmin = ref Sundials.big_real in
      for i = 0 to A.length num - 1 do
        if (A.get denom i) <> 0.0 then
          lmin := min !lmin (A.get num i /. A.get denom i)
      done;
      Mpi.allreduce_float !lmin Mpi.Float_min comm

    let n_vprod (x, _, _) (y, _, _) (z, _, _) =
      for i = 0 to A.length x - 1 do
        A.set z i (A.get x i *. A.get y i)
      done

    let n_vdiv (x, _, _) (y, _, _) (z, _, _) =
      for i = 0 to A.length x - 1 do
        A.set z i (A.get x i /. A.get y i)
      done

    let n_vabs (x, _, _) (z, _, _) =
      for i = 0 to A.length x - 1 do
        A.set z i (abs_float (A.get x i))
      done

    let n_vinv (x, _, _) (z, _, _) =
      for i = 0 to A.length x - 1 do
        A.set z i (1.0 /. (A.get x i))
      done
  end

(* (* Too slow *)
module SlowerDataOps = MakeOps (struct
    type local_data = Sundials.RealArray.t

    let get       = Bigarray.Array1.get
    let set       = Bigarray.Array1.set
    let fill      = Bigarray.Array1.fill

    let make      = Sundials.RealArray.make
    let length    = Sundials.RealArray.length
    let clone     = Sundials.RealArray.clone
  end)
*)

module DataOps =
  struct
    module A = Bigarray.Array1

    let make      = Sundials.RealArray.make
    let clone     = Sundials.RealArray.clone

    type t = Sundials.RealArray.t * int * Mpi.communicator
    type d = Sundials.RealArray.t

    (* let n_vclone (d, gl, comm) = (A.clone d, gl, comm) *)
    let n_vclone (d, gl, comm) = (d, gl, comm)

    let arr_vaxpy a (x : d) (y : d) =
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

    let n_vlinearsum a ((x : d), _, _) b ((y : d), _, _) ((z : d), _, _) =
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

    let n_vconst c ((a : d), _, _) = A.fill a c

    let n_vscale c ((x : d), _, _) ((z : d), _, _) =
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

    let n_vaddconst ((x : d), _, _) b ((z : d), _, _) =
      for i = 0 to A.dim x - 1 do
        A.set z i (A.get x i +. b)
      done

    let n_vmaxnorm ((x : d), _, comm) =
      let lmax = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        let ax = abs_float (A.get x i) in
        if ax > !lmax then lmax := ax
      done;
      Mpi.allreduce_float !lmax Mpi.Float_max comm

    let n_vwrmsnorm ((x : d), n_global, comm) ((w : d), _, _) =
      let lsum = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        lsum := !lsum
                  +. ((A.get x i) *. (A.get w i) *. (A.get x i) *. (A.get w i))
      done;
      let gsum = Mpi.allreduce_float !lsum Mpi.Float_sum comm in
      sqrt (gsum /. float n_global)

    let n_vwrmsnormmask ((x : d), n_global, comm) ((w : d), _, _) ((id : d), _, _) =
      let lsum = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        if A.get id i > 0.0 then
          lsum := !lsum +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
      done;
      let gsum = Mpi.allreduce_float !lsum Mpi.Float_sum comm in
      sqrt (gsum /. float n_global)

    let n_vmin ((x : d), _, comm) =
      let lmin = ref max_float in
      for i = 0 to A.dim x - 1 do
        let xv = A.get x i in
        if xv < !lmin then lmin := xv
      done;
      Mpi.allreduce_float !lmin Mpi.Float_min comm

    let n_vdotprod ((x : d), _, comm) ((y : d), _, _) =
      let lsum = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        lsum := !lsum +. (A.get x i *. A.get y i)
      done;
      Mpi.allreduce_float !lsum Mpi.Float_sum comm

    let n_vcompare c ((x : d), _, _) ((z : d), _, _) =
      for i = 0 to A.dim x - 1 do
        A.set z i (if abs_float (A.get x i) >= c then 1.0 else 0.0)
      done

    let n_vinvtest ((x : d), _, comm) ((z : d), _, _) =
      let lval = ref 1.0 in
      for i = 0 to A.dim x do
        if A.get x i = 0.0 then lval := 0.0
        else A.set z i (1.0 /. (A.get x i))
      done;
      (Mpi.allreduce_float !lval Mpi.Float_min comm = 0.0)

    let n_vwl2norm ((x : d), _, comm) ((w : d), _, _) =
      let lsum = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        lsum := !lsum +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
      done;
      let gsum = Mpi.allreduce_float !lsum Mpi.Float_sum comm in
      sqrt gsum

    let n_vl1norm ((x : d), _, comm) =
      let lsum = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        lsum := !lsum +. abs_float (A.get x i)
      done;
      Mpi.allreduce_float !lsum Mpi.Float_sum comm

    let n_vconstrmask ((c : d), _, _) ((x : d), _, comm) ((m : d), _, _) =
      let test = ref 1.0 in
      let check b = if b then 0.0 else (test := 0.0; 1.0) in
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
      (Mpi.allreduce_float !test Mpi.Float_min comm = 1.0)

    let n_vminquotient ((num : d), _, comm) ((denom : d), _, _) =
      let lmin = ref Sundials.big_real in
      for i = 0 to A.dim num - 1 do
        if (A.get denom i) <> 0.0 then
          lmin := min !lmin (A.get num i /. A.get denom i)
      done;
      Mpi.allreduce_float !lmin Mpi.Float_min comm

    let n_vprod ((x : d), _, _) ((y : d), _, _) ((z : d), _, _) =
      for i = 0 to A.dim x - 1 do
        A.set z i (A.get x i *. A.get y i)
      done

    let n_vdiv ((x : d), _, _) ((y : d), _, _) ((z : d), _, _) =
      for i = 0 to A.dim x - 1 do
        A.set z i (A.get x i /. A.get y i)
      done

    let n_vabs ((x : d), _, _) ((z : d), _, _) =
      for i = 0 to A.dim x - 1 do
        A.set z i (abs_float (A.get x i))
      done

    let n_vinv ((x : d), _, _) ((z : d), _, _) =
      for i = 0 to A.dim x - 1 do
        A.set z i (1.0 /. (A.get x i))
      done
  end

