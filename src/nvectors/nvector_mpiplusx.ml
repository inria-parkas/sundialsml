(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2021 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

open Sundials

(* Accessed directly from nvector_many_ml.c. *)
type data = Nvector.any * Mpi.communicator

type kind

type t = (data, kind) Nvector.t

type Nvector.gdata += MpiPlusX of data

external c_wrap
  : Nvector.any * Mpi.communicator * int
    -> (t -> bool)
    -> (t -> t)
    -> t
  = "sunml_nvec_wrap_mpiplusx"

let unwrap = Nvector.unwrap

let subvector_mpi_rank nv =
  match Nvector_parallel.get_communicator nv with
  | None -> 0
  | Some comm -> Mpi.comm_rank comm

let sumlens comm nv =
  let local_length =
    if subvector_mpi_rank nv = 0 then Nvector.Ops.n_vgetlength nv else 0
  in
  Mpi.(allreduce_int local_length Sum comm)

let rec wrap comm nv =
  if Sundials_impl.Versions.sundials_lt500
    then raise Config.NotImplementedBySundialsVersion;
  let check mnv' =
    let nv', _ = unwrap mnv' in
    try Nvector.check nv nv'; true
    with Nvector.IncompatibleNvector -> false
  in
  c_wrap (nv, comm, sumlens comm nv) check clone

and clone mnv =
  let nv, comm = unwrap mnv in
  wrap comm (Nvector.clone nv)

external c_ident_or_congruent : Mpi.communicator -> Mpi.communicator -> bool
  = "sunml_nvector_parallel_compare_comms" [@@noalloc]

let communicator nv = snd (Nvector.unwrap nv)

module Ops : Nvector.NVECTOR_OPS with type t = t =
struct (* {{{ *)
  type t = (data, kind) Nvector.t
  let check = Nvector.check

  let n_vclone = Nvector.clone

  (* The underlying representation is a ManyVector with one element, so we
     just use its routines directly. *)

  external c_n_vlinearsum    : float -> t -> float -> t -> t -> unit
    = "sunml_nvec_mpimany_n_vlinearsum" [@@noalloc]

  let n_vlinearsum a (x : t) b (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_n_vlinearsum a x b y z

  external n_vconst          : float -> t -> unit
    = "sunml_nvec_mpimany_n_vconst" [@@noalloc]

  external c_n_vprod         : t -> t -> t -> unit
    = "sunml_nvec_mpimany_n_vprod" [@@noalloc]

  let n_vprod (x : t) (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_n_vprod x y z

  external c_n_vdiv          : t -> t -> t -> unit
    = "sunml_nvec_mpimany_n_vdiv" [@@noalloc]

  let n_vdiv (x : t) (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_n_vdiv x y z

  external c_n_vscale        : float -> t -> t -> unit
    = "sunml_nvec_mpimany_n_vscale" [@@noalloc]

  let n_vscale c (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_n_vscale c x z

  external c_n_vabs          : t -> t -> unit
    = "sunml_nvec_mpimany_n_vabs" [@@noalloc]

  let n_vabs (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_n_vabs x z

  external c_n_vinv          : t -> t -> unit
    = "sunml_nvec_mpimany_n_vinv" [@@noalloc]

  external c_n_vaddconst     : t -> float -> t -> unit
    = "sunml_nvec_mpimany_n_vaddconst" [@@noalloc]

  external c_n_vdotprod      : t -> t -> float
    = "sunml_nvec_mpimany_n_vdotprod"

  external n_vmaxnorm        : t -> float
    = "sunml_nvec_mpimany_n_vmaxnorm"

  external c_n_vwrmsnorm     : t -> t -> float
    = "sunml_nvec_mpimany_n_vwrmsnorm"

  external c_n_vwrmsnormmask : t -> t -> t -> float
    = "sunml_nvec_mpimany_n_vwrmsnormmask"

  external n_vmin            : t -> float
    = "sunml_nvec_mpimany_n_vmin"

  external c_n_vwl2norm      : t -> t -> float
    = "sunml_nvec_mpimany_n_vwl2norm"

  external n_vl1norm         : t -> float
    = "sunml_nvec_mpimany_n_vl1norm"

  external c_n_vcompare      : float -> t -> t -> unit
    = "sunml_nvec_mpimany_n_vcompare" [@@noalloc]

  external c_n_vinvtest      : t -> t -> bool
    = "sunml_nvec_mpimany_n_vinvtest" [@@noalloc]

  let n_vinv (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_n_vinv x z

  external c_n_vconstrmask   : t -> t -> t -> bool
    = "sunml_nvec_mpimany_n_vconstrmask" [@@noalloc]

  let n_vconstrmask (c : t) (x : t) (m : t) =
    if Sundials_configuration.safe then (check c x; check c m);
    c_n_vconstrmask c x m

  external c_n_vminquotient  : t -> t -> float
    = "sunml_nvec_mpimany_n_vminquotient"

  let n_vminquotient (n : t) (d : t) =
    if Sundials_configuration.safe then check n d;
    c_n_vminquotient n d

  external n_vspace          : t -> int * int
    = "sunml_nvec_mpimany_n_vspace" [@@noalloc]

  external n_vgetlength      : t -> int
    = "sunml_nvec_mpimany_n_vgetlength" [@@noalloc]

  external c_n_vlinearcombination : RealArray.t -> t array -> t -> unit
    = "sunml_nvec_mpimany_n_vlinearcombination"

  let n_vlinearcombination ca (xa : t array) (z : t) =
    if Sundials_configuration.safe then Array.iter (check z) xa;
    c_n_vlinearcombination ca xa z

  external c_n_vscaleaddmulti : RealArray.t -> t -> t array -> t array -> unit
    = "sunml_nvec_mpimany_n_vscaleaddmulti"

  let n_vscaleaddmulti aa (x : t) (ya : t array) (za : t array) =
    if Sundials_configuration.safe then
      (Array.iter (check x) ya; Array.iter (check x) za);
    c_n_vscaleaddmulti aa x ya za

  external c_n_vdotprodmulti  : t -> t array -> RealArray.t -> unit
    = "sunml_nvec_mpimany_n_vdotprodmulti"

  let n_vdotprodmulti (x : t) (ya : t array) dp =
    if Sundials_configuration.safe then Array.iter (check x) ya;
    c_n_vdotprodmulti x ya dp

  external c_n_vlinearsumvectorarray
    : float -> t array -> float -> t array -> t array -> unit
    = "sunml_nvec_mpimany_n_vlinearsumvectorarray"

  let n_vlinearsumvectorarray a (xa : t array) b (ya : t array) (za : t array) =
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (check x) ya;
          Array.iter (check x) za);
    c_n_vlinearsumvectorarray a xa b ya za

  external c_n_vscalevectorarray
    : RealArray.t -> t array -> t array -> unit
    = "sunml_nvec_mpimany_n_vscalevectorarray"

  let n_vscalevectorarray c (xa : t array) (za : t array) =
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (check x) za);
    c_n_vscalevectorarray c xa za

  external c_n_vconstvectorarray
    : float -> t array -> unit
    = "sunml_nvec_mpimany_n_vconstvectorarray"

  let n_vconstvectorarray c (za : t array) =
    if Sundials_configuration.safe
    then (let z = Array.get za 0 in
          Array.iter (check z) za);
    c_n_vconstvectorarray c za

  external c_n_vwrmsnormvectorarray
    : t array -> t array -> RealArray.t -> unit
    = "sunml_nvec_mpimany_n_vwrmsnormvectorarray"

  let n_vwrmsnormvectorarray (xa : t array) (wa : t array) nrm =
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (check x) wa);
    c_n_vwrmsnormvectorarray xa wa nrm

  external c_n_vwrmsnormmaskvectorarray
    : t array -> t array -> t -> RealArray.t -> unit
    = "sunml_nvec_mpimany_n_vwrmsnormmaskvectorarray"

  let n_vwrmsnormmaskvectorarray (xa : t array) (wa : t array) (id : t) nrm =
    if Sundials_configuration.safe
    then (Array.iter (check id) xa;
          Array.iter (check id) wa);
    c_n_vwrmsnormmaskvectorarray xa wa id nrm

  let n_vscaleaddmultivectorarray ra (xa : t array) (yaa : t array array)
                                     (zaa : t array array) =
    raise Nvector.OperationNotProvided

  let n_vlinearcombinationvectorarray ca (xaa : t array array) (za : t array) =
    raise Nvector.OperationNotProvided

  let n_vaddconst (x : t) b (z : t) =
    if Sundials_configuration.safe then check x z;
    c_n_vaddconst x b z

  let n_vwrmsnorm (x : t) (w : t) =
    if Sundials_configuration.safe then check x w;
    c_n_vwrmsnorm x w

  let n_vwrmsnormmask (x : t) (w : t) (id : t) =
    if Sundials_configuration.safe then (check x w; check x id);
    c_n_vwrmsnormmask x w id

  let n_vdotprod (x : t) (y : t) =
    if Sundials_configuration.safe then check x y;
    c_n_vdotprod x y

  let n_vcompare c (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_n_vcompare c x z

  let n_vinvtest (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_n_vinvtest x z

  let n_vwl2norm (x : t) (w : t) =
    if Sundials_configuration.safe then check x w;
    c_n_vwl2norm x w

  module Local = struct
    external c_n_vdotprod      : t -> t -> float
      = "sunml_nvec_mpimany_n_vdotprodlocal"

    let n_vdotprod (x : t) (y : t) =
      if Sundials_configuration.safe then check x y;
      c_n_vdotprod x y

    external n_vmaxnorm        : t -> float
      = "sunml_nvec_mpimany_n_vmaxnormlocal"

    external n_vmin            : t -> float
      = "sunml_nvec_mpimany_n_vminlocal"

    external n_vl1norm         : t -> float
      = "sunml_nvec_mpimany_n_vl1normlocal"

    external c_n_vinvtest      : t -> t -> bool
      = "sunml_nvec_mpimany_n_vinvtestlocal" [@@noalloc]

    let n_vinvtest (x : t) (z : t) =
      if Sundials_configuration.safe then check x z;
      c_n_vinvtest x z

    external c_n_vconstrmask   : t -> t -> t -> bool
      = "sunml_nvec_mpimany_n_vconstrmasklocal" [@@noalloc]

    let n_vconstrmask (c : t) (x : t) (m : t) =
      if Sundials_configuration.safe then (check c x; check c m);
      c_n_vconstrmask c x m

    external c_n_vminquotient  : t -> t -> float
      = "sunml_nvec_mpimany_n_vminquotientlocal"

    let n_vminquotient (n : t) (d : t) =
      if Sundials_configuration.safe then check n d;
      c_n_vminquotient n d

    external c_n_vwsqrsum
      : t -> t -> float
      = "sunml_nvec_mpimany_n_vwsqrsumlocal"

    let n_vwsqrsum (x : t) (w : t) =
      if Sundials_configuration.safe then check x w;
      c_n_vwsqrsum x w

    external c_n_vwsqrsummask
      : t -> t -> t -> float
      = "sunml_nvec_mpimany_n_vwsqrsummasklocal"

    let n_vwsqrsummask (x : t) (w : t) (id : t) =
      if Sundials_configuration.safe then (check x w; check x id);
      c_n_vwsqrsummask x w id
  end
end (* }}} *)

module DataOps : Nvector.NVECTOR_OPS with type t = data =
struct (* {{{ *)

  type t = data

  let n_vclone (nv, comm) = (Nvector.clone nv, comm)

  module Local = struct

    let n_vdotprod ((x, _) : t) ((y, _) : t) =
      if Nvector.Local.has_n_vdotprod x
      then Nvector.Ops.Local.n_vdotprod x y
      else
        let contrib = Nvector.Ops.n_vdotprod x y in
        let rank = subvector_mpi_rank x in
        if rank = 0 then contrib else 0.

    let n_vmaxnorm ((x, _) : t) =
      if Nvector.Any.Local.has_n_vmaxnorm x
      then Nvector.Ops.Local.n_vmaxnorm x
      else Nvector.Ops.n_vmaxnorm x

    let n_vmin ((x, _) : t) =
      if Nvector.Any.Local.has_n_vmin x
      then Nvector.Ops.Local.n_vmin x
      else Nvector.Ops.n_vmin x

    let n_vl1norm ((x, _) : t) =
      if Nvector.Local.has_n_vl1norm x
      then Nvector.Ops.Local.n_vl1norm x
      else
        let contrib = Nvector.Ops.n_vl1norm x in
        let rank = subvector_mpi_rank x in
        if rank = 0 then contrib else 0.

    let n_vinvtest ((x, _) : t) ((z, _) : t) =
      if Nvector.Any.Local.has_n_vinvtest x
      then Nvector.Ops.Local.n_vinvtest x z
      else Nvector.Ops.n_vinvtest x z

    let n_vconstrmask ((c, _) : t) ((x, _) : t) ((m, _) : t) =
      if Nvector.Any.Local.has_n_vconstrmask c
      then Nvector.Ops.Local.n_vconstrmask c x m
      else Nvector.Ops.n_vconstrmask c x m

    let n_vminquotient ((n, _) : t) ((d, _) : t) =
      if Nvector.Any.Local.has_n_vminquotient n
      then Nvector.Ops.Local.n_vminquotient n d
      else Nvector.Ops.n_vminquotient n d

    let n_vwsqrsum ((x, _) : t) ((w, _) : t) =
      if Nvector.Local.has_n_vwsqrsum x
      then Nvector.Ops.Local.n_vwsqrsum x w
      else
        let contrib = Nvector.Ops.n_vwrmsnorm x w in
        let rank = subvector_mpi_rank x in
        if rank = 0
        then contrib *. contrib *. float (Nvector.Ops.n_vgetlength x)
        else 0.

    let n_vwsqrsummask ((x, _) : t) ((w, _) : t) ((id, _) : t) =
      if Nvector.Local.has_n_vwsqrsummask x
      then Nvector.Ops.Local.n_vwsqrsummask x w id
      else
        let contrib = Nvector.Ops.n_vwrmsnormmask x w id in
        let rank = subvector_mpi_rank x in
        if rank = 0
        then contrib *. contrib *. float (Nvector.Ops.n_vgetlength x)
        else 0.
  end

  let n_vlinearsum a ((x, _) : t) b ((y, _) : t) ((z, _) : t) =
    Nvector.Ops.n_vlinearsum a x b y z

  let n_vconst c ((z, _) : t) = Nvector.Ops.n_vconst c z

  let n_vscale c ((x, _) : t) ((z, _) : t) = Nvector.Ops.n_vscale c x z

  let n_vaddconst ((x, _) : t) b ((z, _) : t) = Nvector.Ops.n_vaddconst x b z

  let n_vmaxnorm ((x, comm) as xv : t) =
    let lmax = Local.n_vmaxnorm xv in
    Mpi.(allreduce_float lmax Max comm)

  let n_vgetlength ((x, comm) : t) = sumlens comm x

  let n_vwrmsnorm (x : t) (w : t) =
    let lx = n_vgetlength x in
    if lx <> n_vgetlength w then raise Nvector.IncompatibleNvector;
    let gsum = Local.n_vwsqrsum x w in
    sqrt (gsum /. float lx)

  let n_vwrmsnormmask (x : t) (w : t) (id : t) =
    let lx = n_vgetlength x in
    if lx <> n_vgetlength w || lx <> n_vgetlength id
    then raise Nvector.IncompatibleNvector;
    let gsum = Local.n_vwsqrsummask x w id in
    sqrt (gsum /. float lx)

  let n_vmin ((x, comm) as xv : t) =
    let lmin = Local.n_vmin xv in
    Mpi.(allreduce_float lmin Min comm)

  let n_vdotprod ((x, comm) as xv : t) (yv : t) =
    let lsum = Local.n_vdotprod xv yv in
    Mpi.(allreduce_float lsum Sum comm)

  let n_vcompare c ((x, _) : t) ((z, _) : t) = Nvector.Ops.n_vcompare c x z

  let n_vinvtest ((x, comm) as xv : t) (zv : t) =
    let v = if Local.n_vinvtest xv zv then 1. else 0. in
    Mpi.(allreduce_float v Min comm) <> 0.

  let n_vwl2norm (xv : t) (wv : t) =
    let gsum = Local.n_vwsqrsum xv wv in
    sqrt (gsum)

  let n_vl1norm ((x, comm) as xv : t) =
    let lsum = Local.n_vl1norm xv in
    Mpi.(allreduce_float lsum Sum comm)

  let n_vconstrmask (cv : t) ((x, comm) as xv : t) (mv : t) =
    let v = if Local.n_vconstrmask cv xv mv then 1. else 0. in
    Mpi.(allreduce_float v Min comm) <> 0.

  let n_vminquotient ((n, comm) as nv : t) (dv : t) =
    let lmin = Local.n_vminquotient nv dv in
    Mpi.(allreduce_float lmin Min comm)

  let n_vprod ((x, _) : t) ((y, _) : t) ((z, _) : t) = Nvector.Ops.n_vprod x y z

  let n_vdiv ((x, _) : t) ((y, _) : t) ((z, _) : t) = Nvector.Ops.n_vdiv x y z

  let n_vabs ((x, _) : t) ((z, _) : t) = Nvector.Ops.n_vabs x z

  let n_vinv ((x, _) : t) ((z, _) : t) = Nvector.Ops.n_vinv x z

  let n_vspace ((x, _) : t) = Nvector.Ops.n_vspace x

  (* fused and array operations *)

  let pnvs (nv, _) = nv

  let n_vlinearcombination (ca : RealArray.t) (xa : t array) ((z, _) : t) =
    let xdata = Array.map pnvs xa in
    Nvector.Ops.n_vlinearcombination ca xdata z

  let n_vscaleaddmulti (aa : RealArray.t) ((x, _) : t)
                                          (ya : t array) (za : t array) =
    Nvector.Ops.n_vscaleaddmulti aa x (Array.map pnvs ya) (Array.map pnvs za)

  let n_vdotprodmulti (x : t) (ya : t array) (dp : RealArray.t) =
    let f i yi = dp.{i} <- n_vdotprod x yi in
    Array.iteri f ya

  let n_vlinearsumvectorarray a (xa : t array) b (ya : t array) (za : t array) =
    Nvector.Ops.n_vlinearsumvectorarray a (Array.map pnvs xa)
                                        b (Array.map pnvs ya)
                                        (Array.map pnvs za)

  let n_vscalevectorarray (c : RealArray.t) (xa : t array) (za : t array) =
    Nvector.Ops.n_vscalevectorarray c (Array.map pnvs xa) (Array.map pnvs za)

  let n_vconstvectorarray c (za : t array) =
    Nvector.Ops.n_vconstvectorarray c (Array.map pnvs za)

  let n_vwrmsnormvectorarray (xa : t array) (wa : t array) (nrm : RealArray.t) =
    let xlen = Array.length xa in
    let wlen = Array.length wa in
    if xlen <> wlen || xlen <> RealArray.length nrm
      then invalid_arg "n_vwrmsnormvectorarray: arrays must have the same length";
    for i = 0 to xlen - 1 do
      let len = n_vgetlength xa.(i) in
      nrm.{i} <- sqrt (Local.n_vwsqrsum xa.(i) wa.(i) /. float len)
    done

  let n_vwrmsnormmaskvectorarray (xa : t array) (wa : t array) (id : t)
                                 (nrm : RealArray.t) =
    let xlen = Array.length xa in
    let wlen = Array.length wa in
    if xlen <> wlen || xlen <> RealArray.length nrm
      then invalid_arg "n_vwrmsnormvectorarray: arrays must have the same length";
    for i = 0 to xlen - 1 do
      let len = n_vgetlength xa.(i) in
      nrm.{i} <- sqrt (Local.n_vwsqrsummask xa.(i) wa.(i) id /. float len)
    done

  let n_vscaleaddmultivectorarray (ra : RealArray.t) (xa : t array)
                                  (yaa : t array array) (zaa : t array array) =
    raise Config.NotImplementedBySundialsVersion

  let n_vlinearcombinationvectorarray (ca : RealArray.t) (xaa : t array array)
                                      (za : t array) =
    raise Config.NotImplementedBySundialsVersion

end (* }}} *)

module Any = struct (* {{{ *)

  external c_any_wrap
    : extension_constructor
      -> Nvector.any * Mpi.communicator * int
      -> (Nvector.any -> bool)
      -> (Nvector.any -> Nvector.any)
      -> Nvector.any
    = "sunml_nvec_anywrap_mpimany"

  let rec wrap comm nv =
    if Sundials_impl.Versions.sundials_lt500
      then raise Config.NotImplementedBySundialsVersion;
    let check mnv' =
      match unwrap mnv' with
      | MpiPlusX (nv', _) ->
          (try
             Nvector.check nv nv';
             Nvector.get_id nv' = Nvector.MpiPlusX
           with Nvector.IncompatibleNvector -> false)
      | _ -> false
    in
    c_any_wrap [%extension_constructor MpiPlusX]
               (nv, comm, sumlens comm nv)
               check clone

  and clone mnv =
    let nv, comm = match unwrap mnv with
                   | MpiPlusX v -> v
                   | _ -> assert false
    in
    wrap comm (Nvector.clone nv)

end (* }}} *)

