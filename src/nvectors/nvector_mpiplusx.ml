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

(* Selectively enable and disable fused and array operations *)
(* NOTE: intentionally linked to the *_manyvector primitives! *)
external c_enablefusedops_manyvector                : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_mpimany_enablefusedops"
external c_enablelinearcombination_manyvector       : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_mpimany_enablelinearcombination"
external c_enablescaleaddmulti_manyvector           : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_mpimany_enablescaleaddmulti"
external c_enabledotprodmulti_manyvector            : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_mpimany_enabledotprodmulti"
external c_enablelinearsumvectorarray_manyvector    : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_mpimany_enablelinearsumvectorarray"
external c_enablescalevectorarray_manyvector        : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_mpimany_enablescalevectorarray"
external c_enableconstvectorarray_manyvector        : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_mpimany_enableconstvectorarray"
external c_enablewrmsnormvectorarray_manyvector     : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_mpimany_enablewrmsnormvectorarray"
external c_enablewrmsnormmaskvectorarray_manyvector : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_mpimany_enablewrmsnormmaskvectorarray"

let subvector_mpi_rank nv =
  match Nvector_parallel.get_communicator nv with
  | None -> 0
  | Some comm -> Mpi.comm_rank comm

let sumlens comm nv =
  let local_length =
    if subvector_mpi_rank nv = 0 then Nvector.Ops.getlength nv else 0
  in
  Mpi.(allreduce_int local_length Sum comm)

let rec wrap ?(with_fused_ops=false) comm nv =
  if Sundials_impl.Version.lt500
    then raise Config.NotImplementedBySundialsVersion;
  let check mnv' =
    let nv', _ = unwrap mnv' in
    try Nvector.check nv nv'; true
    with Nvector.IncompatibleNvector -> false
  in
  let nv = c_wrap (nv, comm, sumlens comm nv) check clone in
  if with_fused_ops then c_enablefusedops_manyvector nv true;
  nv

and clone mnv =
  let nv, comm = unwrap mnv in
  wrap comm (Nvector.clone nv)

external c_ident_or_congruent : Mpi.communicator -> Mpi.communicator -> bool
  = "sunml_nvector_parallel_compare_comms" [@@noalloc]

let communicator nv = snd (Nvector.unwrap nv)

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
   nv
  = do_enable c_enablefusedops_manyvector nv
              with_fused_ops;
    do_enable c_enablelinearcombination_manyvector nv
              with_linear_combination;
    do_enable c_enablescaleaddmulti_manyvector nv
              with_scale_add_multi;
    do_enable c_enabledotprodmulti_manyvector nv
              with_dot_prod_multi;
    do_enable c_enablelinearsumvectorarray_manyvector nv
              with_linear_sum_vector_array;
    do_enable c_enablescalevectorarray_manyvector nv
              with_scale_vector_array;
    do_enable c_enableconstvectorarray_manyvector nv
              with_const_vector_array;
    do_enable c_enablewrmsnormvectorarray_manyvector nv
              with_wrms_norm_vector_array;
    do_enable c_enablewrmsnormmaskvectorarray_manyvector nv
              with_wrms_norm_mask_vector_array

module Ops : Nvector.NVECTOR_OPS with type t = t =
struct (* {{{ *)
  type t = (data, kind) Nvector.t
  let check = Nvector.check

  let clone = Nvector.clone

  (* The underlying representation is a ManyVector with one element, so we
     just use its routines directly. *)

  external c_linearsum    : float -> t -> float -> t -> t -> unit
    = "sunml_nvec_mpimany_linearsum" [@@noalloc]

  let linearsum a (x : t) b (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_linearsum a x b y z

  external const          : float -> t -> unit
    = "sunml_nvec_mpimany_const" [@@noalloc]

  external c_prod         : t -> t -> t -> unit
    = "sunml_nvec_mpimany_prod" [@@noalloc]

  let prod (x : t) (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_prod x y z

  external c_div          : t -> t -> t -> unit
    = "sunml_nvec_mpimany_div" [@@noalloc]

  let div (x : t) (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_div x y z

  external c_scale        : float -> t -> t -> unit
    = "sunml_nvec_mpimany_scale" [@@noalloc]

  let scale c (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_scale c x z

  external c_abs          : t -> t -> unit
    = "sunml_nvec_mpimany_abs" [@@noalloc]

  let abs (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_abs x z

  external c_inv          : t -> t -> unit
    = "sunml_nvec_mpimany_inv" [@@noalloc]

  external c_addconst     : t -> float -> t -> unit
    = "sunml_nvec_mpimany_addconst" [@@noalloc]

  external c_dotprod      : t -> t -> float
    = "sunml_nvec_mpimany_dotprod"

  external maxnorm        : t -> float
    = "sunml_nvec_mpimany_maxnorm"

  external c_wrmsnorm     : t -> t -> float
    = "sunml_nvec_mpimany_wrmsnorm"

  external c_wrmsnormmask : t -> t -> t -> float
    = "sunml_nvec_mpimany_wrmsnormmask"

  external min            : t -> float
    = "sunml_nvec_mpimany_min"

  external c_wl2norm      : t -> t -> float
    = "sunml_nvec_mpimany_wl2norm"

  external l1norm         : t -> float
    = "sunml_nvec_mpimany_l1norm"

  external c_compare      : float -> t -> t -> unit
    = "sunml_nvec_mpimany_compare" [@@noalloc]

  external c_invtest      : t -> t -> bool
    = "sunml_nvec_mpimany_invtest" [@@noalloc]

  let inv (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_inv x z

  external c_constrmask   : t -> t -> t -> bool
    = "sunml_nvec_mpimany_constrmask" [@@noalloc]

  let constrmask (c : t) (x : t) (m : t) =
    if Sundials_configuration.safe then (check c x; check c m);
    c_constrmask c x m

  external c_minquotient  : t -> t -> float
    = "sunml_nvec_mpimany_minquotient"

  let minquotient (n : t) (d : t) =
    if Sundials_configuration.safe then check n d;
    c_minquotient n d

  external space          : t -> int * int
    = "sunml_nvec_mpimany_space" [@@noalloc]

  external getlength      : t -> int
    = "sunml_nvec_mpimany_getlength" [@@noalloc]

  external c_print_file : t -> Sundials.Logfile.t option -> unit
    = "sunml_nvec_mpimany_print_file"

  let print ?logfile (x : t) = c_print_file x logfile

  external c_linearcombination : RealArray.t -> t array -> t -> unit
    = "sunml_nvec_mpimany_linearcombination"

  let linearcombination ca (xa : t array) (z : t) =
    if Sundials_configuration.safe then Array.iter (check z) xa;
    c_linearcombination ca xa z

  let same_len' n ya =
    if n <> Array.length ya then invalid_arg "arrays of unequal length"
  let same_len xa ya = same_len' (Array.length xa) ya

  external c_scaleaddmulti : RealArray.t -> t -> t array -> t array -> unit
    = "sunml_nvec_mpimany_scaleaddmulti"

  let scaleaddmulti aa (x : t) (ya : t array) (za : t array) =
    if Sundials_configuration.safe then
      (Array.iter (check x) ya; Array.iter (check x) za;
       let nv = RealArray.length aa in
       same_len' nv ya; same_len' nv za);
    c_scaleaddmulti aa x ya za

  external c_dotprodmulti  : t -> t array -> RealArray.t -> unit
    = "sunml_nvec_mpimany_dotprodmulti"

  let dotprodmulti (x : t) (ya : t array) dp =
    if Sundials_configuration.safe then
      (let nv = RealArray.length dp in
       same_len' nv ya;
       Array.iter (check x) ya);
    c_dotprodmulti x ya dp

  external c_linearsumvectorarray
    : float -> t array -> float -> t array -> t array -> unit
    = "sunml_nvec_mpimany_linearsumvectorarray"

  let linearsumvectorarray a (xa : t array) b (ya : t array) (za : t array) =
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (check x) ya;
          Array.iter (check x) za;
          same_len xa ya; same_len xa za);
    c_linearsumvectorarray a xa b ya za

  external c_scalevectorarray
    : RealArray.t -> t array -> t array -> unit
    = "sunml_nvec_mpimany_scalevectorarray"

  let scalevectorarray c (xa : t array) (za : t array) =
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (check x) za;
          same_len xa za);
    c_scalevectorarray c xa za

  external c_constvectorarray
    : float -> t array -> unit
    = "sunml_nvec_mpimany_constvectorarray"

  let constvectorarray c (za : t array) =
    if Sundials_configuration.safe
    then (let z = Array.get za 0 in
          Array.iter (check z) za);
    c_constvectorarray c za

  external c_wrmsnormvectorarray
    : t array -> t array -> RealArray.t -> unit
    = "sunml_nvec_mpimany_wrmsnormvectorarray"

  let wrmsnormvectorarray (xa : t array) (wa : t array) nrm =
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (check x) wa;
         same_len xa wa);
    c_wrmsnormvectorarray xa wa nrm

  external c_wrmsnormmaskvectorarray
    : t array -> t array -> t -> RealArray.t -> unit
    = "sunml_nvec_mpimany_wrmsnormmaskvectorarray"

  let wrmsnormmaskvectorarray (xa : t array) (wa : t array) (id : t) nrm =
    if Sundials_configuration.safe
    then (Array.iter (check id) xa;
          Array.iter (check id) wa;
          same_len xa wa);
    c_wrmsnormmaskvectorarray xa wa id nrm

  (* The generic nvector routine compensates for a missing
     nvscaleaddmultivectorarray operation. *)
  external c_scaleaddmultivectorarray
    : RealArray.t -> t array -> t array array -> t array array -> unit
    = "sunml_nvec_scaleaddmultivectorarray"

  let scaleaddmultivectorarray ra (xa : t array) (yaa : t array array)
                                     (zaa : t array array) =
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          let ns = RealArray.length ra in
          let nv = Array.length xa in
          same_len' ns yaa;
          same_len' ns zaa;
          Array.iter (check x) xa;
          Array.iter (fun ya -> same_len' nv ya; Array.iter (check x) ya) yaa;
          Array.iter (fun za -> same_len' nv za; Array.iter (check x) za) zaa;
          same_len yaa zaa);
    (* MpiPlusX Nvectors do not provide this operation. *)
    c_scaleaddmultivectorarray ra xa yaa zaa

  (* The generic nvector routine compensates for a missing
     nvscaleaddmultivectorarray operation. *)
  external c_linearcombinationvectorarray 
    : RealArray.t -> t array array -> t array -> unit
    = "sunml_nvec_linearcombinationvectorarray"

  let linearcombinationvectorarray ca (xaa : t array array) (za : t array) =
    if Sundials_configuration.safe
    then (let z = Array.get za 0 in
          let ns = RealArray.length ca in
          let nv = Array.length za in
          same_len' ns xaa;
          Array.iter (check z) za;
          Array.iter (fun xa -> same_len' nv xa; Array.iter (check z) xa) xaa);
    (* MpiPlusX Nvectors do not provide this operation. *)
    c_linearcombinationvectorarray ca xaa za

  let addconst (x : t) b (z : t) =
    if Sundials_configuration.safe then check x z;
    c_addconst x b z

  let wrmsnorm (x : t) (w : t) =
    if Sundials_configuration.safe then check x w;
    c_wrmsnorm x w

  let wrmsnormmask (x : t) (w : t) (id : t) =
    if Sundials_configuration.safe then (check x w; check x id);
    c_wrmsnormmask x w id

  let dotprod (x : t) (y : t) =
    if Sundials_configuration.safe then check x y;
    c_dotprod x y

  let compare c (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_compare c x z

  let invtest (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_invtest x z

  let wl2norm (x : t) (w : t) =
    if Sundials_configuration.safe then check x w;
    c_wl2norm x w

  module Local = struct
    external c_dotprod      : t -> t -> float
      = "sunml_nvec_mpimany_dotprodlocal"

    let dotprod (x : t) (y : t) =
      if Sundials_configuration.safe then check x y;
      c_dotprod x y

    external maxnorm        : t -> float
      = "sunml_nvec_mpimany_maxnormlocal"

    external min            : t -> float
      = "sunml_nvec_mpimany_minlocal"

    external l1norm         : t -> float
      = "sunml_nvec_mpimany_l1normlocal"

    external c_invtest      : t -> t -> bool
      = "sunml_nvec_mpimany_invtestlocal" [@@noalloc]

    let invtest (x : t) (z : t) =
      if Sundials_configuration.safe then check x z;
      c_invtest x z

    external c_constrmask   : t -> t -> t -> bool
      = "sunml_nvec_mpimany_constrmasklocal" [@@noalloc]

    let constrmask (c : t) (x : t) (m : t) =
      if Sundials_configuration.safe then (check c x; check c m);
      c_constrmask c x m

    external c_minquotient  : t -> t -> float
      = "sunml_nvec_mpimany_minquotientlocal"

    let minquotient (n : t) (d : t) =
      if Sundials_configuration.safe then check n d;
      c_minquotient n d

    external c_wsqrsum
      : t -> t -> float
      = "sunml_nvec_mpimany_wsqrsumlocal"

    let wsqrsum (x : t) (w : t) =
      if Sundials_configuration.safe then check x w;
      c_wsqrsum x w

    external c_wsqrsummask
      : t -> t -> t -> float
      = "sunml_nvec_mpimany_wsqrsummasklocal"

    let wsqrsummask (x : t) (w : t) (id : t) =
      if Sundials_configuration.safe then (check x w; check x id);
      c_wsqrsummask x w id
  end
end (* }}} *)

module DataOps : Nvector.NVECTOR_OPS with type t = data =
struct (* {{{ *)

  type t = data

  let clone (nv, comm) = (Nvector.clone nv, comm)

  module Local = struct

    let dotprod ((x, _) : t) ((y, _) : t) =
      if Nvector.Ops.Local.has_dotprod x
      then Nvector.Ops.Local.dotprod x y
      else
        let contrib = Nvector.Ops.dotprod x y in
        let rank = subvector_mpi_rank x in
        if rank = 0 then contrib else 0.

    let maxnorm ((x, _) : t) =
      if Nvector.Ops.Local.has_maxnorm x
      then Nvector.Ops.Local.maxnorm x
      else Nvector.Ops.maxnorm x

    let min ((x, _) : t) =
      if Nvector.Ops.Local.has_min x
      then Nvector.Ops.Local.min x
      else Nvector.Ops.min x

    let l1norm ((x, _) : t) =
      if Nvector.Ops.Local.has_l1norm x
      then Nvector.Ops.Local.l1norm x
      else
        let contrib = Nvector.Ops.l1norm x in
        let rank = subvector_mpi_rank x in
        if rank = 0 then contrib else 0.

    let invtest ((x, _) : t) ((z, _) : t) =
      if Nvector.Ops.Local.has_invtest x
      then Nvector.Ops.Local.invtest x z
      else Nvector.Ops.invtest x z

    let constrmask ((c, _) : t) ((x, _) : t) ((m, _) : t) =
      if Nvector.Ops.Local.has_constrmask c
      then Nvector.Ops.Local.constrmask c x m
      else Nvector.Ops.constrmask c x m

    let minquotient ((n, _) : t) ((d, _) : t) =
      if Nvector.Ops.Local.has_minquotient n
      then Nvector.Ops.Local.minquotient n d
      else Nvector.Ops.minquotient n d

    let wsqrsum ((x, _) : t) ((w, _) : t) =
      if Nvector.Ops.Local.has_wsqrsum x
      then Nvector.Ops.Local.wsqrsum x w
      else
        let contrib = Nvector.Ops.wrmsnorm x w in
        let rank = subvector_mpi_rank x in
        if rank = 0
        then contrib *. contrib *. float (Nvector.Ops.getlength x)
        else 0.

    let wsqrsummask ((x, _) : t) ((w, _) : t) ((id, _) : t) =
      if Nvector.Ops.Local.has_wsqrsummask x
      then Nvector.Ops.Local.wsqrsummask x w id
      else
        let contrib = Nvector.Ops.wrmsnormmask x w id in
        let rank = subvector_mpi_rank x in
        if rank = 0
        then contrib *. contrib *. float (Nvector.Ops.getlength x)
        else 0.
  end

  let linearsum a ((x, _) : t) b ((y, _) : t) ((z, _) : t) =
    Nvector.Ops.linearsum a x b y z

  let const c ((z, _) : t) = Nvector.Ops.const c z

  let scale c ((x, _) : t) ((z, _) : t) = Nvector.Ops.scale c x z

  let addconst ((x, _) : t) b ((z, _) : t) = Nvector.Ops.addconst x b z

  let maxnorm ((x, comm) as xv : t) =
    let lmax = Local.maxnorm xv in
    Mpi.(allreduce_float lmax Max comm)

  let getlength ((x, comm) : t) = sumlens comm x

  let print ?logfile ((x, _) : t) = Nvector.Ops.print ?logfile x

  let wrmsnorm (x : t) (w : t) =
    let lx = getlength x in
    if lx <> getlength w then raise Nvector.IncompatibleNvector;
    let gsum = Local.wsqrsum x w in
    sqrt (gsum /. float lx)

  let wrmsnormmask (x : t) (w : t) (id : t) =
    let lx = getlength x in
    if lx <> getlength w || lx <> getlength id
    then raise Nvector.IncompatibleNvector;
    let gsum = Local.wsqrsummask x w id in
    sqrt (gsum /. float lx)

  let min ((x, comm) as xv : t) =
    let lmin = Local.min xv in
    Mpi.(allreduce_float lmin Min comm)

  let dotprod ((x, comm) as xv : t) (yv : t) =
    let lsum = Local.dotprod xv yv in
    Mpi.(allreduce_float lsum Sum comm)

  let compare c ((x, _) : t) ((z, _) : t) = Nvector.Ops.compare c x z

  let invtest ((x, comm) as xv : t) (zv : t) =
    let v = if Local.invtest xv zv then 1. else 0. in
    Mpi.(allreduce_float v Min comm) <> 0.

  let wl2norm (xv : t) (wv : t) =
    let gsum = Local.wsqrsum xv wv in
    sqrt (gsum)

  let l1norm ((x, comm) as xv : t) =
    let lsum = Local.l1norm xv in
    Mpi.(allreduce_float lsum Sum comm)

  let constrmask (cv : t) ((x, comm) as xv : t) (mv : t) =
    let v = if Local.constrmask cv xv mv then 1. else 0. in
    Mpi.(allreduce_float v Min comm) <> 0.

  let minquotient ((n, comm) as nv : t) (dv : t) =
    let lmin = Local.minquotient nv dv in
    Mpi.(allreduce_float lmin Min comm)

  let prod ((x, _) : t) ((y, _) : t) ((z, _) : t) = Nvector.Ops.prod x y z

  let div ((x, _) : t) ((y, _) : t) ((z, _) : t) = Nvector.Ops.div x y z

  let abs ((x, _) : t) ((z, _) : t) = Nvector.Ops.abs x z

  let inv ((x, _) : t) ((z, _) : t) = Nvector.Ops.inv x z

  let space ((x, _) : t) = Nvector.Ops.space x

  (* fused and array operations *)

  let pnvs (nv, _) = nv

  let linearcombination (ca : RealArray.t) (xa : t array) ((z, _) : t) =
    let xdata = Array.map pnvs xa in
    Nvector.Ops.linearcombination ca xdata z

  let scaleaddmulti (aa : RealArray.t) ((x, _) : t)
                                          (ya : t array) (za : t array) =
    Nvector.Ops.scaleaddmulti aa x (Array.map pnvs ya) (Array.map pnvs za)

  let dotprodmulti (x : t) (ya : t array) (dp : RealArray.t) =
    let f i yi = dp.{i} <- dotprod x yi in
    Array.iteri f ya

  let linearsumvectorarray a (xa : t array) b (ya : t array) (za : t array) =
    Nvector.Ops.linearsumvectorarray a (Array.map pnvs xa)
                                        b (Array.map pnvs ya)
                                        (Array.map pnvs za)

  let scalevectorarray (c : RealArray.t) (xa : t array) (za : t array) =
    Nvector.Ops.scalevectorarray c (Array.map pnvs xa) (Array.map pnvs za)

  let constvectorarray c (za : t array) =
    Nvector.Ops.constvectorarray c (Array.map pnvs za)

  let wrmsnormvectorarray (xa : t array) (wa : t array) (nrm : RealArray.t) =
    let xlen = Array.length xa in
    let wlen = Array.length wa in
    if xlen <> wlen || xlen <> RealArray.length nrm
      then invalid_arg "wrmsnormvectorarray: arrays must have the same length";
    for i = 0 to xlen - 1 do
      let len = getlength xa.(i) in
      nrm.{i} <- sqrt (Local.wsqrsum xa.(i) wa.(i) /. float len)
    done

  let wrmsnormmaskvectorarray (xa : t array) (wa : t array) (id : t)
                                 (nrm : RealArray.t) =
    let xlen = Array.length xa in
    let wlen = Array.length wa in
    if xlen <> wlen || xlen <> RealArray.length nrm
      then invalid_arg "wrmsnormvectorarray: arrays must have the same length";
    for i = 0 to xlen - 1 do
      let len = getlength xa.(i) in
      nrm.{i} <- sqrt (Local.wsqrsummask xa.(i) wa.(i) id /. float len)
    done

  let scaleaddmultivectorarray (ra : RealArray.t) (xa : t array)
                                  (yaa : t array array) (zaa : t array array) =
    raise Config.NotImplementedBySundialsVersion

  let linearcombinationvectorarray (ca : RealArray.t) (xaa : t array array)
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

  let rec wrap ?(with_fused_ops=false) comm nv =
    if Sundials_impl.Version.lt500
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
    let nv = c_any_wrap [%extension_constructor MpiPlusX]
                        (nv, comm, sumlens comm nv)
                        check clone
    in
    if with_fused_ops then c_enablefusedops_manyvector nv true;
    nv

  and clone mnv =
    let nv, comm = match unwrap mnv with
                   | MpiPlusX v -> v
                   | _ -> assert false
    in
    wrap comm (Nvector.clone nv)

  let unwrap nv =
    match Nvector.unwrap nv with
    | MpiPlusX a -> a
    | _ -> raise Nvector.BadGenericType

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
     nv
    = if Sundials_impl.Version.lt400
        then raise Config.NotImplementedBySundialsVersion;
      if Nvector.get_id nv <> Nvector.MpiManyVector
        then raise Nvector.BadGenericType;
      do_enable c_enablefusedops_manyvector nv
                with_fused_ops;
      do_enable c_enablelinearcombination_manyvector nv
                with_linear_combination;
      do_enable c_enablescaleaddmulti_manyvector nv
                with_scale_add_multi;
      do_enable c_enabledotprodmulti_manyvector nv
                with_dot_prod_multi;
      do_enable c_enablelinearsumvectorarray_manyvector nv
                with_linear_sum_vector_array;
      do_enable c_enablescalevectorarray_manyvector nv
                with_scale_vector_array;
      do_enable c_enableconstvectorarray_manyvector nv
                with_const_vector_array;
      do_enable c_enablewrmsnormvectorarray_manyvector nv
                with_wrms_norm_vector_array;
      do_enable c_enablewrmsnormmaskvectorarray_manyvector nv
                with_wrms_norm_mask_vector_array

end (* }}} *)

