open Sundials

type kind
type data = RealArray.t * int * Mpi.communicator
type t = (data, kind) Nvector.t

exception IncorrectGlobalSize

external c_wrap : (RealArray.t * int * Mpi.communicator)
                  -> (RealArray.t * int * Mpi.communicator -> bool)
                  -> t
  = "sunml_nvec_wrap_parallel"

(* Selectively enable and disable fused and array operations *)
external c_enablefusedops_parallel                          : t -> bool -> unit
  = "sunml_nvec_par_enablefusedops"
external c_enablelinearcombination_parallel                 : t -> bool -> unit
  = "sunml_nvec_par_enablelinearcombination"
external c_enablescaleaddmulti_parallel                     : t -> bool -> unit
  = "sunml_nvec_par_enablescaleaddmulti"
external c_enabledotprodmulti_parallel                      : t -> bool -> unit
  = "sunml_nvec_par_enabledotprodmulti"
external c_enablelinearsumvectorarray_parallel              : t -> bool -> unit
  = "sunml_nvec_par_enablelinearsumvectorarray"
external c_enablescalevectorarray_parallel                  : t -> bool -> unit
  = "sunml_nvec_par_enablescalevectorarray"
external c_enableconstvectorarray_parallel                  : t -> bool -> unit
  = "sunml_nvec_par_enableconstvectorarray"
external c_enablewrmsnormvectorarray_parallel               : t -> bool -> unit
  = "sunml_nvec_par_enablewrmsnormvectorarray"
external c_enablewrmsnormmaskvectorarray_parallel           : t -> bool -> unit
  = "sunml_nvec_par_enablewrmsnormmaskvectorarray"
external c_enablescaleaddmultivectorarray_parallel          : t -> bool -> unit
  = "sunml_nvec_par_enablescaleaddmultivectorarray"
external c_enablelinearcombinationvectorarray_parallel      : t -> bool -> unit
  = "sunml_nvec_par_enablelinearcombinationvectorarray"

let wrap ?(with_fused_ops=false) ((nl, ng, comm) as v) =
  let nl_len = RealArray.length nl in
  let check (nl', ng', comm') =
    (nl_len = RealArray.length nl') && (ng <= ng') && (comm == comm')
  in
  let nv = c_wrap v check in
  if with_fused_ops then c_enablefusedops_parallel nv true;
  nv

let make ?with_fused_ops nl ng comm iv =
  wrap ?with_fused_ops (RealArray.make nl iv, ng, comm)

let clone nv =
  let loc, glen, comm = Nvector.unwrap nv in
  wrap (RealArray.copy loc, glen, comm)

let unwrap = Nvector.unwrap

let pp fmt nv =
  let data, _, _ = Nvector.unwrap nv in
  RealArray.pp fmt data

let local_array nv =
  let data, _, _ = Nvector.unwrap nv in
  data

let local_length nv = RealArray.length (local_array nv)

let global_length nv =
  let _, gl, _ = Nvector.unwrap nv in
  gl

let communicator nv =
  let _, _, comm = Nvector.unwrap nv in
  comm

external c_get_communicator : ('d, 'k) Nvector.t -> Mpi.communicator option
  = "sunml_nvec_par_n_vgetcommunicator"

let get_communicator nv =
  match c_get_communicator nv with
  | None -> raise Nvector.IncompatibleNvector
  | Some c -> c

let do_enable f nv v =
  match v with
  | None -> ()
  | Some v -> f nv v

let enable
   ?with_fused_ops
   ?with_linear_combination
   ?with_scale_add_multi
   ?with_dot_prod_multi
   ?with_linear_sum_vector_array
   ?with_scale_vector_array
   ?with_const_vector_array
   ?with_wrms_norm_vector_array
   ?with_wrms_norm_mask_vector_array
   ?with_scale_add_multi_vector_array
   ?with_linear_combination_vector_array
   nv
  = do_enable c_enablefusedops_parallel nv
              with_fused_ops;
    do_enable c_enablelinearcombination_parallel nv
              with_linear_combination;
    do_enable c_enablescaleaddmulti_parallel nv
              with_scale_add_multi;
    do_enable c_enabledotprodmulti_parallel nv
              with_dot_prod_multi;
    do_enable c_enablelinearsumvectorarray_parallel nv
              with_linear_sum_vector_array;
    do_enable c_enablescalevectorarray_parallel nv
              with_scale_vector_array;
    do_enable c_enableconstvectorarray_parallel nv
              with_const_vector_array;
    do_enable c_enablewrmsnormvectorarray_parallel nv
              with_wrms_norm_vector_array;
    do_enable c_enablewrmsnormmaskvectorarray_parallel nv
              with_wrms_norm_mask_vector_array;
    do_enable c_enablescaleaddmultivectorarray_parallel nv
              with_scale_add_multi_vector_array;
    do_enable c_enablelinearcombinationvectorarray_parallel nv
              with_linear_combination_vector_array

external hide_communicator
  : Mpi.communicator -> Nvector_custom.communicator
  = "%identity"

module Ops = struct (* {{{ *)
  type t = (data, kind) Nvector.t

  let n_vclone = clone

  external n_vlinearsum    : float -> t -> float -> t -> t -> unit
    = "sunml_nvec_par_n_vlinearsum"

  external n_vconst        : float -> t -> unit
    = "sunml_nvec_par_n_vconst"

  external n_vprod         : t -> t -> t -> unit
    = "sunml_nvec_par_n_vprod"

  external n_vdiv          : t -> t -> t -> unit
    = "sunml_nvec_par_n_vdiv"

  external n_vscale        : float -> t -> t -> unit
    = "sunml_nvec_par_n_vscale"

  external n_vabs          : t -> t -> unit
    = "sunml_nvec_par_n_vabs"

  external n_vinv          : t -> t -> unit
    = "sunml_nvec_par_n_vinv"

  external n_vaddconst     : t -> float -> t -> unit
    = "sunml_nvec_par_n_vaddconst"

  external n_vdotprod      : t -> t -> float
    = "sunml_nvec_par_n_vdotprod"

  external n_vmaxnorm      : t -> float
    = "sunml_nvec_par_n_vmaxnorm"

  external n_vwrmsnorm     : t -> t -> float
    = "sunml_nvec_par_n_vwrmsnorm"

  external n_vwrmsnormmask : t -> t -> t -> float
    = "sunml_nvec_par_n_vwrmsnormmask"

  external n_vmin          : t -> float
    = "sunml_nvec_par_n_vmin"

  external n_vwl2norm      : t -> t -> float
    = "sunml_nvec_par_n_vwl2norm"

  external n_vl1norm       : t -> float
    = "sunml_nvec_par_n_vl1norm"

  external n_vcompare      : float -> t -> t -> unit
    = "sunml_nvec_par_n_vcompare"

  external n_vinvtest      : t -> t -> bool
    = "sunml_nvec_par_n_vinvtest"

  external n_vconstrmask   : t -> t -> t -> bool
    = "sunml_nvec_par_n_vconstrmask"

  external n_vminquotient  : t -> t -> float
    = "sunml_nvec_par_n_vminquotient"

  external n_vspace  : t -> int * int
    = "sunml_nvec_par_n_vspace"

  external n_vgetlength  : t -> int
    = "sunml_nvec_par_n_vgetlength"

  external n_vlinearcombination : RealArray.t -> t array -> t -> unit
    = "sunml_nvec_par_n_vlinearcombination"

  external n_vscaleaddmulti : RealArray.t -> t -> t array -> t array -> unit
    = "sunml_nvec_par_n_vscaleaddmulti"

  external n_vdotprodmulti : t -> t array -> RealArray.t -> unit
    = "sunml_nvec_par_n_vdotprodmulti"

  external n_vlinearsumvectorarray
    : float -> t array -> float -> t array -> t array -> unit
    = "sunml_nvec_par_n_vlinearsumvectorarray"

  external n_vscalevectorarray
    : RealArray.t -> t array -> t array -> unit
    = "sunml_nvec_par_n_vscalevectorarray"

  external n_vconstvectorarray
    : float -> t array -> unit
    = "sunml_nvec_par_n_vconstvectorarray"

  external n_vwrmsnormvectorarray
    : t array -> t array -> RealArray.t -> unit
    = "sunml_nvec_par_n_vwrmsnormvectorarray"

  external n_vwrmsnormmaskvectorarray
    : t array -> t array -> t -> RealArray.t -> unit
    = "sunml_nvec_par_n_vwrmsnormmaskvectorarray"

  external n_vscaleaddmultivectorarray
    : RealArray.t -> t array -> t array array -> t array array -> unit
    = "sunml_nvec_par_n_vscaleaddmultivectorarray"

  external n_vlinearcombinationvectorarray
    : RealArray.t -> t array array -> t array -> unit
    = "sunml_nvec_par_n_vlinearcombinationvectorarray"

  module Local = struct

    external n_vdotprod
      : t -> t -> float
      = "sunml_nvec_par_n_vdotprodlocal"

    external n_vmaxnorm
      : t -> float
      = "sunml_nvec_par_n_vmaxnormlocal"

    external n_vmin
      : t -> float
      = "sunml_nvec_par_n_vminlocal"

    external n_vl1norm
      : t -> float
      = "sunml_nvec_par_n_vl1normlocal"

    external n_vinvtest
      : t -> t -> bool
      = "sunml_nvec_par_n_vinvtestlocal"

    external n_vconstrmask
      : t -> t -> t -> bool
      = "sunml_nvec_par_n_vconstrmasklocal"

    external n_vminquotient
      : t -> t -> float
      = "sunml_nvec_par_n_vminquotientlocal"

    external n_vwsqrsum
      : t -> t -> float
      = "sunml_nvec_par_n_vwsqrsumlocallocal"

    external n_vwsqrsummask
      : t -> t -> t -> float
      = "sunml_nvec_par_n_vwsqrsummasklocal"
  end
end (* }}} *)

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
  struct (* {{{ *)
    type t = A.local_data * int * Mpi.communicator

    (* let n_vclone (d, gl, comm) = (A.clone d, gl, comm) *)
    let n_vclone (d, gl, comm) = (A.clone d, gl, comm)

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
          A.set z i (A.get v2 i -. A.get v1 i)
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
      Mpi.allreduce_float !lmax Mpi.Max comm

    let n_vwrmsnorm (x, n_global, comm) (w, _, _) =
      let lsum = ref 0.0 in
      for i = 0 to A.length x - 1 do
        lsum := !lsum
                  +. ((A.get x i) *. (A.get w i) *. (A.get x i) *. (A.get w i))
      done;
      let gsum = Mpi.allreduce_float !lsum Mpi.Sum comm in
      sqrt (gsum /. float n_global)

    let n_vwrmsnormmask (x, n_global, comm) (w, _, _) (id, _, _) =
      let lsum = ref 0.0 in
      for i = 0 to A.length x - 1 do
        if A.get id i > 0.0 then
          lsum := !lsum +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
      done;
      let gsum = Mpi.allreduce_float !lsum Mpi.Sum comm in
      sqrt (gsum /. float n_global)

    let n_vmin (x, _, comm) =
      let lmin = ref max_float in
      for i = 0 to A.length x - 1 do
        let xv = A.get x i in
        if xv < !lmin then lmin := xv
      done;
      Mpi.allreduce_float !lmin Mpi.Min comm

    let n_vdotprod (x, _, comm) (y, _, _) =
      let lsum = ref 0.0 in
      for i = 0 to A.length x - 1 do
        lsum := !lsum +. (A.get x i *. A.get y i)
      done;
      Mpi.allreduce_float !lsum Mpi.Sum comm

    let n_vcompare c (x, _, _) (z, _, _) =
      for i = 0 to A.length x - 1 do
        A.set z i (if abs_float (A.get x i) >= c then 1.0 else 0.0)
      done

    let n_vinvtest (x, _, comm) (z, _, _) =
      let lval = ref 1 in
      for i = 0 to A.length x - 1 do
        if A.get x i = 0.0 then lval := 0
        else A.set z i (1.0 /. (A.get x i))
      done;
      (Mpi.allreduce_int !lval Mpi.Min comm <> 0)

    let n_vwl2norm (x, _, comm) (w, _, _) =
      let lsum = ref 0.0 in
      for i = 0 to A.length x - 1 do
        lsum := !lsum +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
      done;
      let gsum = Mpi.allreduce_float !lsum Mpi.Sum comm in
      sqrt gsum

    let n_vl1norm (x, _, comm) =
      let lsum = ref 0.0 in
      for i = 0 to A.length x - 1 do
        lsum := !lsum +. abs_float (A.get x i)
      done;
      Mpi.allreduce_float !lsum Mpi.Sum comm

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
      (Mpi.allreduce_float !test Mpi.Min comm = 1.0)

    let n_vminquotient (num, _, comm) (denom, _, _) =
      let lmin = ref Config.big_real in
      for i = 0 to A.length num - 1 do
        if (A.get denom i) <> 0.0 then
          lmin := min !lmin (A.get num i /. A.get denom i)
      done;
      Mpi.allreduce_float !lmin Mpi.Min comm

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

    let n_vspace (_, ng, comm) = (ng, 2 * Mpi.comm_size comm)

    let n_vgetlength (_, ng, _) = ng

    (* fused and array operations *)

    let n_vlinearcombination (ca : RealArray.t) (xa : t array) (z : t) =
      let nvec = Array.length xa in
      if nvec = 1 then n_vscale ca.{0} xa.(0) z
      else if nvec = 2 then n_vlinearsum ca.{0} xa.(0) ca.{1} xa.(1) z
      else
        let zd, _, _ = z in
        let n = A.length zd in
        if xa.(0) == z then begin
          let c0 = ca.{0} in
          if c0 <> 1.0 then
            for j = 0 to n - 1 do
              A.set zd j (A.get zd j *. c0)
            done;
          for i = 1 to nvec - 1 do
            let ci, (xd, _, _) = ca.{i}, xa.(i) in
            for j = 0 to n - 1 do
              A.set zd j (A.get zd j +. ci *. A.get xd j)
            done
          done
        end
        else begin
          let c0, (xd, _, _) = ca.{0}, xa.(0) in
          for j = 0 to n - 1 do
            A.set zd j (c0 *. A.get xd j)
          done;
          for i = 1 to nvec - 1 do
            let ci, (xd, _, _) = ca.{i}, xa.(i) in
            for j = 0 to n - 1 do
              A.set zd j (A.get zd j +. ci *. A.get xd j)
            done
          done
        end

    let n_vscaleaddmulti (aa : RealArray.t) (x : t) (ya : t array) (za : t array) =
      let nvec = Array.length ya in
      if nvec = 1 then n_vlinearsum aa.{0} x 1.0 ya.(0) za.(0)
      else
        let xd, _, _ = x in
        let n = A.length xd in
        if ya == za then
          for i = 0 to nvec - 1 do
            let a, (yd, _, _) = aa.{i}, ya.(i) in
            for j = 0 to n - 1 do
              A.set yd j (A.get yd j +. a *. A.get xd j)
            done
          done
        else
          for i = 0 to nvec - 1 do
            let ai, (yd, _, _), (zd, _, _) = aa.{i}, ya.(i), za.(i) in
            for j = 0 to n - 1 do
              A.set zd j (ai *. A.get xd j +. A.get yd j)
            done
          done

    let n_vdotprodmulti (x : t) (ya : t array) (dp : RealArray.t) =
      let nvec = Array.length ya in
      if nvec = 1 then dp.{0} <- n_vdotprod x ya.(0)
      else
        let xd, _, comm = x in
        let n = A.length xd in
        for i = 0 to nvec - 1 do
          let yd, _, _ = ya.(i) in
          dp.{i} <- 0.0;
          for j = 0 to n - 1 do
            dp.{i} <- dp.{i} +. A.get xd j *. A.get yd j
          done
        done;
        (* Note: ocamlmpi does not provide MPI_IN_PLACE *)
        Mpi.(allreduce_bigarray1 dp dp Sum comm)

    let arr_vaxpy_array a (xa : t array) (ya : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.length xad in
      if a = 1.0 then
        for i = 0 to nvec - 1 do
          let xd, _, _ = xa.(i) in
          let yd, _, _ = ya.(i) in
          for j = 0 to n - 1 do
            A.set yd j (A.get yd j +. A.get xd j)
          done
        done
      else if a = -1.0 then
        for i = 0 to nvec - 1 do
          let xd, _, _ = xa.(i) in
          let yd, _, _ = ya.(i) in
          for j = 0 to n - 1 do
            A.set yd j (A.get yd j -. A.get xd j)
          done
        done
      else
        for i = 0 to nvec - 1 do
          let xd, _, _ = xa.(i) in
          let yd, _, _ = ya.(i) in
          for j = 0 to n - 1 do
            A.set yd j (A.get yd j +. a *. A.get xd j)
          done
        done

    let v_sumvectorarray (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.length xad in
      for i = 0 to nvec - 1 do
        let xd, _, _ = xa.(i) in
        let yd, _, _ = ya.(i) in
        let zd, _, _ = za.(i) in
        for j = 0 to n - 1 do
          A.set zd j (A.get xd j +. A.get yd j)
        done
      done

    let v_diffvectorarray (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.length xad in
      for i = 0 to nvec - 1 do
        let xd, _, _ = xa.(i) in
        let yd, _, _ = ya.(i) in
        let zd, _, _ = za.(i) in
        for j = 0 to n - 1 do
          A.set zd j (A.get xd j -. A.get yd j)
        done
      done

    let v_lin1vectorarray a (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.length xad in
      for i = 0 to nvec - 1 do
        let xd, _, _ = xa.(i) in
        let yd, _, _ = ya.(i) in
        let zd, _, _ = za.(i) in
        for j = 0 to n - 1 do
          A.set zd j (a *. A.get xd j +. A.get yd j)
        done
      done

    let v_lin2vectorarray a (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.length xad in
      for i = 0 to nvec - 1 do
        let xd, _, _ = xa.(i) in
        let yd, _, _ = ya.(i) in
        let zd, _, _ = za.(i) in
        for j = 0 to n - 1 do
          A.set zd j (a *. A.get xd j -. A.get yd j)
        done
      done

    let v_scalesumvectorarray c (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.length xad in
      for i = 0 to nvec - 1 do
        let xd, _, _ = xa.(i) in
        let yd, _, _ = ya.(i) in
        let zd, _, _ = za.(i) in
        for j = 0 to n - 1 do
          A.set zd j (c *. (A.get xd j +. A.get yd j))
        done
      done

    let v_scalediffvectorarray c (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.length xad in
      for i = 0 to nvec - 1 do
        let xd, _, _ = xa.(i) in
        let yd, _, _ = ya.(i) in
        let zd, _, _ = za.(i) in
        for j = 0 to n - 1 do
          A.set zd j (c *. (A.get xd j -. A.get yd j))
        done
      done

    let n_vlinearsumvectorarray a (xa : t array) b (ya : t array) (za : t array) =
      let nvec = Array.length ya in
      if nvec = 1 then n_vlinearsum a xa.(0) b ya.(0) za.(0)
      else if b =  1.0 && (za == ya) then arr_vaxpy_array a xa ya
      else if a =  1.0 && (za == xa) then arr_vaxpy_array a ya xa
      else if a =  1.0 && b =  1.0 then v_sumvectorarray xa ya za
      else if a =  1.0 && b = -1.0 then v_diffvectorarray xa ya za
      else if a = -1.0 && b =  1.0 then v_diffvectorarray ya xa za
      else if a =  1.0 then v_lin1vectorarray b ya xa za
      else if b =  1.0 then v_lin1vectorarray a xa ya za
      else if a = -1.0 then v_lin2vectorarray b ya xa za
      else if b = -1.0 then v_lin2vectorarray a xa ya za
      else if a = b then v_scalesumvectorarray a xa ya za
      else if a = -. b then v_scalediffvectorarray a xa ya za
      else
        let xad, _, _ = xa.(0) in
        let n = A.length xad in
        for i = 0 to nvec - 1 do
          let xd, _, _ = xa.(i) in
          let yd, _, _ = ya.(i) in
          let zd, _, _ = za.(i) in
          for j = 0 to n - 1 do
            A.set zd j (a *. A.get xd j +. b *. A.get yd j)
          done
        done

    let n_vscalevectorarray (c : RealArray.t) (xa : t array) (za : t array) =
      let nvec = Array.length xa in
      if nvec = 1 then n_vscale c.{0} xa.(0) za.(0)
      else
        let xad, _, _ = xa.(0) in
        let n = A.length xad in
        if xa == za then
          for i = 0 to nvec - 1 do
            let (xd, _, _), c = xa.(i), c.{i} in
            for j = 0 to n - 1 do
              A.set xd j (c *. A.get xd j)
            done
          done
        else
          for i = 0 to nvec - 1 do
            let c = c.{i} in
            let xd, _, _ = xa.(i) in
            let zd, _, _ = za.(i) in
            for j = 0 to n - 1 do
              A.set zd j (c *. A.get xd j)
            done
          done

    let n_vconstvectorarray c (za : t array) =
      let nvec = Array.length za in
      if nvec = 1 then n_vconst c za.(0)
      else
        let zad, _, _ = za.(0) in
        let n = A.length zad in
        for i = 0 to nvec - 1 do
          let zd, _, _ = za.(i) in
          for j = 0 to n - 1 do
            A.set zd j c
          done
        done

    let n_vwrmsnormvectorarray (xa : t array) (wa : t array) (nrm : RealArray.t) =
      let nvec = Array.length xa in
      if nvec = 1 then nrm.{0} <- n_vwrmsnorm xa.(0) wa.(0)
      else
        let xad, ng, comm = xa.(0) in
        let nl = A.length xad in
        let nf = float ng in
        for i = 0 to nvec - 1 do
          let xd, _, _ = xa.(i) in
          let wd, _, _ = wa.(i) in
          nrm.{i} <- 0.0;
          for j = 0 to nl - 1 do
            let s = A.get xd j *. A.get wd j in
            nrm.{i} <- nrm.{i} +. s *. s
          done
        done;
        (* Note: ocamlmpi does not provide MPI_IN_PLACE *)
        Mpi.(allreduce_bigarray1 nrm nrm Sum comm);
        for i = 0 to nvec - 1 do
          nrm.{i} <- sqrt (nrm.{i}/.nf)
        done

    let n_vwrmsnormmaskvectorarray (xa : t array) (wa : t array) (id : t)
                                   (nrm : RealArray.t) =
      let nvec = Array.length xa in
      if nvec = 1 then nrm.{0} <- n_vwrmsnormmask xa.(0) wa.(0) id
      else
        let xad, ng, comm = xa.(0) in
        let idd, _, _ = id in
        let nl = A.length xad in
        let nf = float ng in
        for i = 0 to nvec - 1 do
          let xd, _, _ = xa.(i) in
          let wd, _, _ = wa.(i) in
          nrm.{i} <- 0.0;
          for j = 0 to nl - 1 do
            if A.get idd j > 0.0 then begin
              let s = A.get xd j *. A.get wd j in
              nrm.{i} <- nrm.{i} +. s *. s
            end
          done;
        done;
        (* Note: ocamlmpi does not provide MPI_IN_PLACE *)
        Mpi.(allreduce_bigarray1 nrm nrm Sum comm);
        for i = 0 to nvec - 1 do
          nrm.{i} <- sqrt (nrm.{i}/.nf)
        done

    let n_vscaleaddmultivectorarray (ra : RealArray.t) (xa : t array)
                                    (yaa : t array array) (zaa : t array array) =
      let nsum = Array.length yaa in
      let nvec = Array.length yaa.(0) in
      if nvec = 1 then begin
        if nsum = 1 then n_vlinearsum ra.{0} xa.(0) 1.0 yaa.(0).(0) zaa.(0).(0)
        else
          let yya = Array.init nsum (fun j -> yaa.(j).(0)) in
          let zza = Array.init nsum (fun j -> zaa.(j).(0)) in
          n_vscaleaddmulti ra xa.(0) yya zza
      end
      else if nsum = 1 then n_vlinearsumvectorarray ra.{0} xa 1.0 yaa.(0) zaa.(0)
      else
        let xad, _, _ = xa.(0) in
        let n = A.length xad in
        if (yaa == zaa) then
          for i = 0 to nvec - 1 do
            let xd, _, _ = xa.(i) in
            for j = 0 to nsum - 1 do
              let a, (yd, _, _) = ra.{j}, yaa.(j).(i) in
              for k = 0 to n - 1 do
                A.set yd k (A.get yd k +. a *. A.get xd k)
              done
            done
          done
        else
          for i = 0 to nvec - 1 do
            let xd, _, _ = xa.(i) in
            for j = 0 to nsum - 1 do
              let a = ra.{j} in
              let yd, _, _ = yaa.(j).(i) in
              let zd, _, _ = zaa.(j).(i) in
              for k = 0 to n - 1 do
                A.set zd k (a *. A.get xd k +. A.get yd k)
              done
            done
          done

    let n_vlinearcombinationvectorarray (ca : RealArray.t) (xaa : t array array)
                                        (za : t array) =
      let nsum = Array.length xaa in
      let nvec = Array.length xaa.(0) in
      if nvec = 1 then begin
        if nsum = 1 then n_vscale ca.{0} xaa.(0).(0) za.(0)
        else if nsum = 2
             then n_vlinearsum ca.{0} xaa.(0).(0) ca.{1} xaa.(1).(0) za.(0)
        else
          let ya = Array.init nsum (fun i -> xaa.(i).(0)) in
          n_vlinearcombination ca ya za.(0)
      end
      else
        if nsum = 1 then
          let ctmp = RealArray.make nvec ca.{0} in
          n_vscalevectorarray ctmp xaa.(0) za
        else if nsum = 2 then
          n_vlinearsumvectorarray ca.{0} xaa.(0) ca.{1} xaa.(1) za
        else
          let zad, _, _ = za.(0) in
          let n = A.length zad in
          if xaa.(0) == za then begin
            if ca.{0} = 1.0 then
              for j = 0 to nvec - 1 do
                let zd, _, _ = za.(j) in
                for i = 1 to nsum - 1 do
                  let c, (xd, _, _) = ca.{i}, xaa.(i).(j) in
                  for k = 0 to n - 1 do
                    A.set zd k (A.get zd k +. c *. A.get xd k)
                  done
                done
              done
            else
              let c0 = ca.{0} in
              for j = 0 to nvec - 1 do
                let zd, _, _ = za.(j) in
                for k = 0 to n - 1 do
                  A.set zd k (A.get zd k *. c0)
                done;
                for i = 1 to nsum - 1 do
                  let c, (xd, _, _) = ca.{i}, xaa.(i).(j) in
                  for k = 0 to n - 1 do
                    A.set zd k (A.get zd k +. c *. A.get xd k)
                  done
                done
              done
          end
          else
            let c0 = ca.{0} in
            for j = 0 to nvec - 1 do
              let xd, _, _ = xaa.(0).(j) in
              let zd, _, _ = za.(j) in
              for k = 0 to n - 1 do
                A.set zd k (c0 *. A.get xd k)
              done;
              for i = 1 to nsum - 1 do
                let c, (xd, _, _) = ca.{i}, xaa.(i).(j) in
                for k = 0 to n - 1 do
                  A.set zd k (A.get zd k +. c *. A.get xd k)
                done
              done
            done

    module Local = struct
      let n_vdotprod     = n_vdotprod
      let n_vmaxnorm     = n_vmaxnorm
      let n_vmin         = n_vmin
      let n_vl1norm      = n_vl1norm
      let n_vinvtest     = n_vinvtest
      let n_vconstrmask  = n_vconstrmask
      let n_vminquotient = n_vminquotient

      let n_vwsqrsum (x, _, _) (w, _, _) =
        let a = ref 0.0 in
        let lx = A.length x in
        for i = 0 to lx - 1 do
          a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
        done;
        !a

      let n_vwsqrsummask (x, _, _) (w, _, _) (id, _, _) =
        let a = ref 0.0 in
        let lx = A.length x in
        for i = 0 to lx - 1 do
          if A.get id i > 0.0 then
            a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
        done;
        !a
    end
  end (* }}} *)

(* (* Too slow *)
module SlowerDataOps = MakeOps (struct
    type local_data = RealArray.t

    let get       = Bigarray.Array1.get
    let set       = Bigarray.Array1.set
    let fill      = Bigarray.Array1.fill

    let make      = RealArray.make
    let length    = RealArray.length
    let clone     = RealArray.clone
  end)
*)

module DataOps =
  struct (* {{{ *)
    module A = Bigarray.Array1

    let make      = RealArray.make
    let clone     = RealArray.copy

    type t = RealArray.t * int * Mpi.communicator
    type d = RealArray.t

    let n_vclone (d, gl, comm) = (clone d, gl, comm)

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
          A.set z i (A.get v2 i -. A.get v1 i)
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
      Mpi.allreduce_float !lmax Mpi.Max comm

    let n_vwrmsnorm ((x : d), n_global, comm) ((w : d), _, _) =
      let lsum = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        lsum := !lsum
                  +. ((A.get x i) *. (A.get w i) *. (A.get x i) *. (A.get w i))
      done;
      let gsum = Mpi.allreduce_float !lsum Mpi.Sum comm in
      sqrt (gsum /. float n_global)

    let n_vwrmsnormmask ((x : d), n_global, comm) ((w : d), _, _) ((id : d), _, _) =
      let lsum = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        if A.get id i > 0.0 then
          lsum := !lsum +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
      done;
      let gsum = Mpi.allreduce_float !lsum Mpi.Sum comm in
      sqrt (gsum /. float n_global)

    let n_vmin ((x : d), _, comm) =
      let lmin = ref max_float in
      for i = 0 to A.dim x - 1 do
        let xv = A.get x i in
        if xv < !lmin then lmin := xv
      done;
      Mpi.allreduce_float !lmin Mpi.Min comm

    let n_vdotprod ((x : d), _, comm) ((y : d), _, _) =
      let lsum = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        lsum := !lsum +. (A.get x i *. A.get y i)
      done;
      Mpi.allreduce_float !lsum Mpi.Sum comm

    let n_vcompare c ((x : d), _, _) ((z : d), _, _) =
      for i = 0 to A.dim x - 1 do
        A.set z i (if abs_float (A.get x i) >= c then 1.0 else 0.0)
      done

    let n_vinvtest ((x : d), _, comm) ((z : d), _, _) =
      let lval = ref 1 in
      for i = 0 to A.dim x - 1 do
        if A.get x i = 0.0 then lval := 0
        else A.set z i (1.0 /. (A.get x i))
      done;
      (Mpi.allreduce_int !lval Mpi.Min comm <> 0)

    let n_vwl2norm ((x : d), _, comm) ((w : d), _, _) =
      let lsum = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        lsum := !lsum +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
      done;
      let gsum = Mpi.allreduce_float !lsum Mpi.Sum comm in
      sqrt gsum

    let n_vl1norm ((x : d), _, comm) =
      let lsum = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        lsum := !lsum +. abs_float (A.get x i)
      done;
      Mpi.allreduce_float !lsum Mpi.Sum comm

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
      (Mpi.allreduce_float !test Mpi.Min comm = 1.0)

    let n_vminquotient ((num : d), _, comm) ((denom : d), _, _) =
      let lmin = ref Config.big_real in
      for i = 0 to A.dim num - 1 do
        if (A.get denom i) <> 0.0 then
          lmin := min !lmin (A.get num i /. A.get denom i)
      done;
      Mpi.allreduce_float !lmin Mpi.Min comm

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

    let n_vspace (_, ng, comm) = (ng, 2 * Mpi.comm_size comm)

    let n_vgetlength (_, ng, _) = ng

    (* fused and array operations *)

    let n_vlinearcombination (ca : RealArray.t) (xa : t array) (z : t) =
      let nvec = Array.length xa in
      if nvec = 1 then n_vscale ca.{0} xa.(0) z
      else if nvec = 2 then n_vlinearsum ca.{0} xa.(0) ca.{1} xa.(1) z
      else
        let zd, _, _ = z in
        let n = RealArray.length zd in
        if xa.(0) == z then begin
          let c0 = ca.{0} in
          if c0 <> 1.0 then
            for j = 0 to n - 1 do
              A.set zd j (A.get zd j *. c0)
            done;
          for i = 1 to nvec - 1 do
            let ci, (xd, _, _) = ca.{i}, xa.(i) in
            for j = 0 to n - 1 do
              A.set zd j (A.get zd j +. ci *. A.get xd j)
            done
          done
        end
        else begin
          let c0, (xd, _, _) = ca.{0}, xa.(0) in
          for j = 0 to n - 1 do
            A.set zd j (c0 *. A.get xd j)
          done;
          for i = 1 to nvec - 1 do
            let ci, (xd, _, _) = ca.{i}, xa.(i) in
            for j = 0 to n - 1 do
              A.set zd j (A.get zd j +. ci *. A.get xd j)
            done
          done
        end

    let n_vscaleaddmulti (aa : RealArray.t) (x : t) (ya : t array) (za : t array) =
      let nvec = Array.length ya in
      if nvec = 1 then n_vlinearsum aa.{0} x 1.0 ya.(0) za.(0)
      else
        let xd, _, _ = x in
        let n = RealArray.length xd in
        if ya == za then
          for i = 0 to nvec - 1 do
            let a, (yd, _, _) = aa.{i}, ya.(i) in
            for j = 0 to n - 1 do
              A.set yd j (A.get yd j +. a *. A.get xd j)
            done
          done
        else
          for i = 0 to nvec - 1 do
            let ai, (yd, _, _), (zd, _, _) = aa.{i}, ya.(i), za.(i) in
            for j = 0 to n - 1 do
              A.set zd j (ai *. A.get xd j +. A.get yd j)
            done
          done

    let n_vdotprodmulti (x : t) (ya : t array) (dp : RealArray.t) =
      let nvec = Array.length ya in
      if nvec = 1 then dp.{0} <- n_vdotprod x ya.(0)
      else
        let xd, _, comm = x in
        let n = RealArray.length xd in
        for i = 0 to nvec - 1 do
          let yd, _, _ = ya.(i) in
          dp.{i} <- 0.0;
          for j = 0 to n - 1 do
            dp.{i} <- dp.{i} +. A.get xd j *. A.get yd j
          done
        done;
        (* Note: ocamlmpi does not provide MPI_IN_PLACE *)
        Mpi.(allreduce_bigarray1 dp dp Sum comm)

    let arr_vaxpy_array a (xa : t array) (ya : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.dim xad in
      if a = 1.0 then
        for i = 0 to nvec - 1 do
          let xd, _, _ = xa.(i) in
          let yd, _, _ = ya.(i) in
          for j = 0 to n - 1 do
            A.set yd j (A.get yd j +. A.get xd j)
          done
        done
      else if a = -1.0 then
        for i = 0 to nvec - 1 do
          let xd, _, _ = xa.(i) in
          let yd, _, _ = ya.(i) in
          for j = 0 to n - 1 do
            A.set yd j (A.get yd j -. A.get xd j)
          done
        done
      else
        for i = 0 to nvec - 1 do
          let xd, _, _ = xa.(i) in
          let yd, _, _ = ya.(i) in
          for j = 0 to n - 1 do
            A.set yd j (A.get yd j +. a *. A.get xd j)
          done
        done

    let v_sumvectorarray (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.dim xad in
      for i = 0 to nvec - 1 do
        let xd, _, _ = xa.(i) in
        let yd, _, _ = ya.(i) in
        let zd, _, _ = za.(i) in
        for j = 0 to n - 1 do
          A.set zd j (A.get xd j +. A.get yd j)
        done
      done

    let v_diffvectorarray (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.dim xad in
      for i = 0 to nvec - 1 do
        let xd, _, _ = xa.(i) in
        let yd, _, _ = ya.(i) in
        let zd, _, _ = za.(i) in
        for j = 0 to n - 1 do
          A.set zd j (A.get xd j -. A.get yd j)
        done
      done

    let v_lin1vectorarray a (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.dim xad in
      for i = 0 to nvec - 1 do
        let xd, _, _ = xa.(i) in
        let yd, _, _ = ya.(i) in
        let zd, _, _ = za.(i) in
        for j = 0 to n - 1 do
          A.set zd j (a *. A.get xd j +. A.get yd j)
        done
      done

    let v_lin2vectorarray a (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.dim xad in
      for i = 0 to nvec - 1 do
        let xd, _, _ = xa.(i) in
        let yd, _, _ = ya.(i) in
        let zd, _, _ = za.(i) in
        for j = 0 to n - 1 do
          A.set zd j (a *. A.get xd j -. A.get yd j)
        done
      done

    let v_scalesumvectorarray c (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.dim xad in
      for i = 0 to nvec - 1 do
        let xd, _, _ = xa.(i) in
        let yd, _, _ = ya.(i) in
        let zd, _, _ = za.(i) in
        for j = 0 to n - 1 do
          A.set zd j (c *. (A.get xd j +. A.get yd j))
        done
      done

    let v_scalediffvectorarray c (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.dim xad in
      for i = 0 to nvec - 1 do
        let xd, _, _ = xa.(i) in
        let yd, _, _ = ya.(i) in
        let zd, _, _ = za.(i) in
        for j = 0 to n - 1 do
          A.set zd j (c *. (A.get xd j -. A.get yd j))
        done
      done

    let n_vlinearsumvectorarray a (xa : t array) b (ya : t array) (za : t array) =
      let nvec = Array.length ya in
      if nvec = 1 then n_vlinearsum a xa.(0) b ya.(0) za.(0)
      else if b =  1.0 && (za == ya) then arr_vaxpy_array a xa ya
      else if a =  1.0 && (za == xa) then arr_vaxpy_array a ya xa
      else if a =  1.0 && b =  1.0 then v_sumvectorarray xa ya za
      else if a =  1.0 && b = -1.0 then v_diffvectorarray xa ya za
      else if a = -1.0 && b =  1.0 then v_diffvectorarray ya xa za
      else if a =  1.0 then v_lin1vectorarray b ya xa za
      else if b =  1.0 then v_lin1vectorarray a xa ya za
      else if a = -1.0 then v_lin2vectorarray b ya xa za
      else if b = -1.0 then v_lin2vectorarray a xa ya za
      else if a = b then v_scalesumvectorarray a xa ya za
      else if a = -. b then v_scalediffvectorarray a xa ya za
      else
        let xad, _, _ = xa.(0) in
        let n = A.dim xad in
        for i = 0 to nvec - 1 do
          let xd, _, _ = xa.(i) in
          let yd, _, _ = ya.(i) in
          let zd, _, _ = za.(i) in
          for j = 0 to n - 1 do
            A.set zd j (a *. A.get xd j +. b *. A.get yd j)
          done
        done

    let n_vscalevectorarray (c : RealArray.t) (xa : t array) (za : t array) =
      let nvec = Array.length xa in
      if nvec = 1 then n_vscale c.{0} xa.(0) za.(0)
      else
        let xad, _, _ = xa.(0) in
        let n = A.dim xad in
        if xa == za then
          for i = 0 to nvec - 1 do
            let (xd, _, _), c = xa.(i), c.{i} in
            for j = 0 to n - 1 do
              A.set xd j (c *. A.get xd j)
            done
          done
        else
          for i = 0 to nvec - 1 do
            let c = c.{i} in
            let xd, _, _ = xa.(i) in
            let zd, _, _ = za.(i) in
            for j = 0 to n - 1 do
              A.set zd j (c *. A.get xd j)
            done
          done

    let n_vconstvectorarray c (za : t array) =
      let nvec = Array.length za in
      if nvec = 1 then n_vconst c za.(0)
      else
        let zad, _, _ = za.(0) in
        let n = A.dim zad in
        for i = 0 to nvec - 1 do
          let zd, _, _ = za.(i) in
          for j = 0 to n - 1 do
            A.set zd j c
          done
        done

    let n_vwrmsnormvectorarray (xa : t array) (wa : t array) (nrm : RealArray.t) =
      let nvec = Array.length xa in
      if nvec = 1 then nrm.{0} <- n_vwrmsnorm xa.(0) wa.(0)
      else
        let xad, ng, comm = xa.(0) in
        let nl = A.dim xad in
        let nf = float ng in
        for i = 0 to nvec - 1 do
          let xd, _, _ = xa.(i) in
          let wd, _, _ = wa.(i) in
          nrm.{i} <- 0.0;
          for j = 0 to nl - 1 do
            let s = A.get xd j *. A.get wd j in
            nrm.{i} <- nrm.{i} +. s *. s
          done
        done;
        (* Note: ocamlmpi does not provide MPI_IN_PLACE *)
        Mpi.(allreduce_bigarray1 nrm nrm Sum comm);
        for i = 0 to nvec - 1 do
          nrm.{i} <- sqrt (nrm.{i}/.nf)
        done

    let n_vwrmsnormmaskvectorarray (xa : t array) (wa : t array) (id : t)
                                   (nrm : RealArray.t) =
      let nvec = Array.length xa in
      if nvec = 1 then nrm.{0} <- n_vwrmsnormmask xa.(0) wa.(0) id
      else
        let xad, ng, comm = xa.(0) in
        let idd, _, _ = id in
        let nl = A.dim xad in
        let nf = float ng in
        for i = 0 to nvec - 1 do
          let xd, _, _ = xa.(i) in
          let wd, _, _ = wa.(i) in
          nrm.{i} <- 0.0;
          for j = 0 to nl - 1 do
            if A.get idd j > 0.0 then begin
              let s = A.get xd j *. A.get wd j in
              nrm.{i} <- nrm.{i} +. s *. s
            end
          done;
        done;
        (* Note: ocamlmpi does not provide MPI_IN_PLACE *)
        Mpi.(allreduce_bigarray1 nrm nrm Sum comm);
        for i = 0 to nvec - 1 do
          nrm.{i} <- sqrt (nrm.{i}/.nf)
        done

    let n_vscaleaddmultivectorarray (ra : RealArray.t) (xa : t array)
                                    (yaa : t array array) (zaa : t array array) =
      let nsum = Array.length yaa in
      let nvec = Array.length yaa.(0) in
      if nvec = 1 then begin
        if nsum = 1 then n_vlinearsum ra.{0} xa.(0) 1.0 yaa.(0).(0) zaa.(0).(0)
        else
          let yya = Array.init nsum (fun j -> yaa.(j).(0)) in
          let zza = Array.init nsum (fun j -> zaa.(j).(0)) in
          n_vscaleaddmulti ra xa.(0) yya zza
      end
      else if nsum = 1 then n_vlinearsumvectorarray ra.{0} xa 1.0 yaa.(0) zaa.(0)
      else
        let xad, _, _ = xa.(0) in
        let n = A.dim xad in
        if (yaa == zaa) then
          for i = 0 to nvec - 1 do
            let xd, _, _ = xa.(i) in
            for j = 0 to nsum - 1 do
              let a, (yd, _, _) = ra.{j}, yaa.(j).(i) in
              for k = 0 to n - 1 do
                A.set yd k (A.get yd k +. a *. A.get xd k)
              done
            done
          done
        else
          for i = 0 to nvec - 1 do
            let xd, _, _ = xa.(i) in
            for j = 0 to nsum - 1 do
              let a = ra.{j} in
              let yd, _, _ = yaa.(j).(i) in
              let zd, _, _ = zaa.(j).(i) in
              for k = 0 to n - 1 do
                A.set zd k (a *. A.get xd k +. A.get yd k)
              done
            done
          done

    let n_vlinearcombinationvectorarray (ca : RealArray.t) (xaa : t array array)
                                        (za : t array) =
      let nsum = Array.length xaa in
      let nvec = Array.length xaa.(0) in
      if nvec = 1 then begin
        if nsum = 1 then n_vscale ca.{0} xaa.(0).(0) za.(0)
        else if nsum = 2
             then n_vlinearsum ca.{0} xaa.(0).(0) ca.{1} xaa.(1).(0) za.(0)
        else
          let ya = Array.init nsum (fun i -> xaa.(i).(0)) in
          n_vlinearcombination ca ya za.(0)
      end
      else
        if nsum = 1 then
          let ctmp = RealArray.make nvec ca.{0} in
          n_vscalevectorarray ctmp xaa.(0) za
        else if nsum = 2 then
          n_vlinearsumvectorarray ca.{0} xaa.(0) ca.{1} xaa.(1) za
        else
          let zad, _, _ = za.(0) in
          let n = A.dim zad in
          if xaa.(0) == za then begin
            if ca.{0} = 1.0 then
              for j = 0 to nvec - 1 do
                let zd, _, _ = za.(j) in
                for i = 1 to nsum - 1 do
                  let c, (xd, _, _) = ca.{i}, xaa.(i).(j) in
                  for k = 0 to n - 1 do
                    A.set zd k (A.get zd k +. c *. A.get xd k)
                  done
                done
              done
            else
              let c0 = ca.{0} in
              for j = 0 to nvec - 1 do
                let zd, _, _ = za.(j) in
                for k = 0 to n - 1 do
                  A.set zd k (A.get zd k *. c0)
                done;
                for i = 1 to nsum - 1 do
                  let c, (xd, _, _) = ca.{i}, xaa.(i).(j) in
                  for k = 0 to n - 1 do
                    A.set zd k (A.get zd k +. c *. A.get xd k)
                  done
                done
              done
          end
          else
            let c0 = ca.{0} in
            for j = 0 to nvec - 1 do
              let xd, _, _ = xaa.(0).(j) in
              let zd, _, _ = za.(j) in
              for k = 0 to n - 1 do
                A.set zd k (c0 *. A.get xd k)
              done;
              for i = 1 to nsum - 1 do
                let c, (xd, _, _) = ca.{i}, xaa.(i).(j) in
                for k = 0 to n - 1 do
                  A.set zd k (A.get zd k +. c *. A.get xd k)
                done
              done
            done

    module Local = struct
      let n_vdotprod     = n_vdotprod
      let n_vmaxnorm     = n_vmaxnorm
      let n_vmin         = n_vmin
      let n_vl1norm      = n_vl1norm
      let n_vinvtest     = n_vinvtest
      let n_vconstrmask  = n_vconstrmask
      let n_vminquotient = n_vminquotient

      let n_vwsqrsum ((x : d), _, _) ((w : d), _, _) =
        let a = ref 0.0 in
        let lx = A.dim x in
        for i = 0 to lx - 1 do
          a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
        done;
        !a

      let n_vwsqrsummask ((x : d), _, _) ((w : d), _, _) ((id : d), _, _) =
        let a = ref 0.0 in
        let lx = A.dim x in
        for i = 0 to lx - 1 do
          if A.get id i > 0.0 then
            a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
        done;
        !a
    end
  end (* }}} *)


(* Let C code know about some of the values in this module.  *)
external c_init_module : exn array -> unit =
  "sunml_nvector_parallel_init_module"

let _ =
  c_init_module
    (* Exceptions must be listed in the same order as
       nvector_parallel_exn_index.  *)
    [|IncorrectGlobalSize|]
