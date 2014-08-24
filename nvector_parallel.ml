
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
      val fold_left : ('a -> float -> 'a) -> 'a -> local_data -> 'a
    end) ->
  struct
    type t = A.local_data * int * Mpi.communicator

    (* let n_vclone (d, gl, comm) = (A.clone d, gl, comm) *)
    let n_vclone (d, gl, comm) = (d, gl, comm)

    let lift_map f x y =
      for i = 0 to A.length x - 1 do
        A.set x i (f (A.get x i) (A.get y i))
      done

    let lift_bop f x y z =
      for i = 0 to A.length x - 1 do
        A.set z i (f (A.get x i) (A.get y i))
      done

    let lift_op f x z =
      for i = 0 to A.length x - 1 do
        A.set z i (f (A.get x i))
      done

    let zip_fold_left f iv x y =
      let a = ref iv in
      for i = 0 to A.length x - 1 do
        a := f !a (A.get x i) (A.get y i)
      done;
      !a

    let triple_zip_fold_left f iv x y z =
      let a = ref iv in
      for i = 0 to A.length x - 1 do
        a := f !a (A.get x i) (A.get y i) (A.get z i)
      done;
      !a

    let arr_vaxpy a x y =
      if a = 1.0 then
        lift_map (fun y x -> y +. x) y x
      else if a = -1.0 then
        lift_map (fun y x -> y -. x) y x
      else
        lift_map (fun y x -> y +. a *. x) y x

    let n_vlinearsum a (x, _, _) b (y, _, _) (z, _, _) =
      if b = 1.0 && z == y then
        arr_vaxpy a x y
      else if a = 1.0 && z == x then
        arr_vaxpy b y x
      else if a = 1.0 && b = 1.0 then
        lift_bop (fun x y -> x +. y) x y z
      else if (a = 1.0 && b = -1.0) || (a = -1.0 && b == 1.0) then
        let v1, v2 = if (a = 1.0 && b = -1.0) then y, x else x, y in
        lift_bop (fun x y -> x -. y) v1 v2 z
      else if a = 1.0 || b = 1.0 then
        let c, v1, v2 = if a = 1.0 then b, y, x else a, x, y in
        lift_bop (fun x y -> c *. x +. y) v1 v2 z
      else if a = -1.0 || b = -1.0 then
        let c, v1, v2 = if a = -1.0 then b, y, x else a, x, y in
        lift_bop (fun x y -> a *. x -. y) v1 v2 z
      else if a = b then
        lift_bop (fun x y -> a *. (x +. y)) x y z
      else if a = -.b then
        lift_bop (fun x y -> a *. (x -. y)) x y z
      else
        lift_bop (fun x y -> a *. x +. b *. y) x y z

    let n_vconst c (a, _, _) = A.fill a c

    let n_vscale c (x, _, _) (z, _, _) =
      if c = 1.0 then
        lift_op (fun x -> x) x z
      else if c = -1.0 then
        lift_op (fun x -> -. x) x z
      else
        lift_op (fun x -> c *. x) x z

    let n_vaddconst (x, _, _) b (z, _, _) = lift_op (fun x -> x +. b) x z

    let n_vmaxnorm (x, _, comm) =
      let f max x =
        let ax = abs_float x in
        if ax > max then ax else max
      in
      let lmax = A.fold_left f 0.0 x in
      Mpi.allreduce_float lmax Mpi.Float_max comm

    let n_vwrmsnorm (x, n_global, comm) (w, _, _) =
      let f a x w = a +. ((x *. w) ** 2.0) in
      let lsum = zip_fold_left f 0.0 x w in
      let gsum = Mpi.allreduce_float lsum Mpi.Float_sum comm in
      sqrt (gsum /. float n_global)

    let n_vwrmsnormmask (x, n_global, comm) (w, _, _) (id, _, _) =
      let f a id x w = if id > 0.0 then a +. ((x *. w) ** 2.0) else a in
      let lsum = triple_zip_fold_left f 0.0 id x w in
      let gsum = Mpi.allreduce_float lsum Mpi.Float_sum comm in
      sqrt (gsum /. float n_global)

    let n_vmin (x, _, comm) =
      let f min x = if x < min then x else min in
      let lmin = A.fold_left f max_float x in
      Mpi.allreduce_float lmin Mpi.Float_min comm

    let n_vdotprod (x, _, comm) (y, _, _) =
      let f a x y = a +. x *. y in
      let lsum = zip_fold_left f 0.0 x y in
      Mpi.allreduce_float lsum Mpi.Float_sum comm

    let n_vcompare c (x, _, _) (z, _, _) =
      lift_op (fun x -> if abs_float x >= c then 1.0 else 0.0) x z

    let n_vinvtest (x, _, comm) (z, _, _) =
      let l = A.length x in
      let rec f r i =
        if i = l then r
        else if (A.get x i) = 0.0 then f 0.0 (i + 1)
        else (A.set z i (1.0 /. (A.get x i)); f r (i + 1))
      in
      let lval = f 1.0 0 in
      (Mpi.allreduce_float lval Mpi.Float_min comm = 0.0)

    let n_vwl2norm (x, _, comm) (w, _, _) =
      let f a x w = a +. ((x *. w) ** 2.0) in
      let lsum = zip_fold_left f 0.0 x w in
      let gsum = Mpi.allreduce_float lsum Mpi.Float_sum comm in
      sqrt gsum

    let n_vl1norm (x, _, comm) =
      let f a x = a +. abs_float x in
      let lsum = A.fold_left f 0.0 x in
      Mpi.allreduce_float lsum Mpi.Float_sum comm

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
      lift_bop f c x m;
      (Mpi.allreduce_float !test Mpi.Float_min comm = 1.0)

    let n_vminquotient (num, _, comm) (denom, _, _) =
      let f m n d =
        if d = 0.0 then m
        else min m (n /. d)
      in
      let lmin = zip_fold_left f Sundials.big_real num denom in
      Mpi.allreduce_float lmin Mpi.Float_min comm

    let n_vprod (x, _, _) (y, _, _) (z, _, _) = lift_bop ( *. ) x y z
    let n_vdiv (x, _, _) (y, _, _) (z, _, _) = lift_bop ( /. ) x y z
    let n_vabs (x, _, _) (z, _, _) = lift_op abs_float x z
    let n_vinv (x, _, _) (z, _, _) = lift_op (fun x -> 1.0 /. x) x z
  end

module DataOps = MakeOps (struct
    type local_data = Sundials.RealArray.t

    let get       = Bigarray.Array1.get
    let set       = Bigarray.Array1.set
    let fill      = Bigarray.Array1.fill

    let make      = Sundials.RealArray.make
    let length    = Sundials.RealArray.length
    let clone     = Sundials.RealArray.clone
    let fold_left = Sundials.RealArray.fold_left
  end)

