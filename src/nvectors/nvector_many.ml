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
type data = Nvector.any ROArray.t * int

type kind

type t = (data, kind) Nvector.t

type Nvector.gdata += Many of data

external c_wrap : data -> (t -> bool) -> (t -> t) -> t
  = "sunml_nvec_wrap_many"

let unwrap = Nvector.unwrap

let rec wrap_withlen ((nvs, _) as payload) =
  let check nv' =
    let nvs', _ = unwrap nv' in
    try ROArray.iter2 Nvector.check nvs nvs'; true
    with Nvector.IncompatibleNvector -> false
  in
  c_wrap payload check clone

and clone nv =
  let nvs, gl = unwrap nv in
  wrap_withlen (ROArray.map Nvector.clone nvs, gl)

let sumlens =
  let f sum nv = sum + Nvector.Ops.n_vgetlength nv in
  ROArray.fold_left f 0

let wrap nvs =
  if Sundials_impl.Versions.sundials_lt500
    then raise Config.NotImplementedBySundialsVersion;
  wrap_withlen (nvs, sumlens nvs)

let length nv = snd (unwrap nv)

let num_subvectors nv = ROArray.length (fst (unwrap nv))

module Ops : Nvector.NVECTOR_OPS with type t = t =
struct (* {{{ *)
  type t = (data, kind) Nvector.t
  let check = Nvector.check

  let n_vclone = Nvector.clone

  external c_n_vlinearsum    : float -> t -> float -> t -> t -> unit
    = "sunml_nvec_many_n_vlinearsum" [@@noalloc]

  let n_vlinearsum a (x : t) b (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_n_vlinearsum a x b y z

  external n_vconst          : float -> t -> unit
    = "sunml_nvec_many_n_vconst" [@@noalloc]

  external c_n_vprod         : t -> t -> t -> unit
    = "sunml_nvec_many_n_vprod" [@@noalloc]

  let n_vprod (x : t) (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_n_vprod x y z

  external c_n_vdiv          : t -> t -> t -> unit
    = "sunml_nvec_many_n_vdiv" [@@noalloc]

  let n_vdiv (x : t) (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_n_vdiv x y z

  external c_n_vscale        : float -> t -> t -> unit
    = "sunml_nvec_many_n_vscale" [@@noalloc]

  let n_vscale c (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_n_vscale c x z

  external c_n_vabs          : t -> t -> unit
    = "sunml_nvec_many_n_vabs" [@@noalloc]

  let n_vabs (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_n_vabs x z

  external c_n_vinv          : t -> t -> unit
    = "sunml_nvec_many_n_vinv" [@@noalloc]

  let n_vinv (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_n_vinv x z

  external c_n_vaddconst     : t -> float -> t -> unit
    = "sunml_nvec_many_n_vaddconst" [@@noalloc]

  let n_vaddconst (x : t) b (z : t) =
    if Sundials_configuration.safe then check x z;
    c_n_vaddconst x b z

  external c_n_vwrmsnorm     : t -> t -> float
    = "sunml_nvec_many_n_vwrmsnorm"

  let n_vwrmsnorm (x : t) (w : t) =
    if Sundials_configuration.safe then check x w;
    c_n_vwrmsnorm x w

  external c_n_vwrmsnormmask : t -> t -> t -> float
    = "sunml_nvec_many_n_vwrmsnormmask"

  let n_vwrmsnormmask (x : t) (w : t) (id : t) =
    if Sundials_configuration.safe then (check x w; check x id);
    c_n_vwrmsnormmask x w id

  external c_n_vwl2norm      : t -> t -> float
    = "sunml_nvec_many_n_vwl2norm"

  let n_vwl2norm (x : t) (w : t) =
    if Sundials_configuration.safe then check x w;
    c_n_vwl2norm x w

  external c_n_vcompare      : float -> t -> t -> unit
    = "sunml_nvec_many_n_vcompare" [@@noalloc]

  let n_vcompare c (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_n_vcompare c x z

  external c_n_vinvtest      : t -> t -> bool
    = "sunml_nvec_many_n_vinvtestlocal" [@@noalloc]

  let n_vinvtest (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_n_vinvtest x z

  external n_vspace          : t -> int * int
    = "sunml_nvec_many_n_vspace" [@@noalloc]

  external n_vgetlength      : t -> int
    = "sunml_nvec_many_n_vgetlength" [@@noalloc]

  external c_n_vlinearcombination : RealArray.t -> t array -> t -> unit
    = "sunml_nvec_many_n_vlinearcombination"

  let n_vlinearcombination ca (xa : t array) (z : t) =
    if Sundials_configuration.safe then Array.iter (check z) xa;
    c_n_vlinearcombination ca xa z

  external c_n_vscaleaddmulti : RealArray.t -> t -> t array -> t array -> unit
    = "sunml_nvec_many_n_vscaleaddmulti"

  let n_vscaleaddmulti aa (x : t) (ya : t array) (za : t array) =
    if Sundials_configuration.safe then
      (Array.iter (check x) ya; Array.iter (check x) za);
    c_n_vscaleaddmulti aa x ya za

  external c_n_vdotprodmulti  : t -> t array -> RealArray.t -> unit
    = "sunml_nvec_many_n_vdotprodmulti"

  let n_vdotprodmulti (x : t) (ya : t array) dp =
    if Sundials_configuration.safe then Array.iter (check x) ya;
    c_n_vdotprodmulti x ya dp

  external c_n_vlinearsumvectorarray
    : float -> t array -> float -> t array -> t array -> unit
    = "sunml_nvec_many_n_vlinearsumvectorarray"

  let n_vlinearsumvectorarray a (xa : t array) b (ya : t array) (za : t array) =
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (check x) ya;
          Array.iter (check x) za);
    c_n_vlinearsumvectorarray a xa b ya za

  external c_n_vscalevectorarray
    : RealArray.t -> t array -> t array -> unit
    = "sunml_nvec_many_n_vscalevectorarray"

  let n_vscalevectorarray c (xa : t array) (za : t array) =
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (check x) za);
    c_n_vscalevectorarray c xa za

  external c_n_vconstvectorarray
    : float -> t array -> unit
    = "sunml_nvec_many_n_vconstvectorarray"

  let n_vconstvectorarray c (za : t array) =
    if Sundials_configuration.safe
    then (let z = Array.get za 0 in
          Array.iter (check z) za);
    c_n_vconstvectorarray c za

  external c_n_vwrmsnormvectorarray
    : t array -> t array -> RealArray.t -> unit
    = "sunml_nvec_many_n_vwrmsnormvectorarray"

  let n_vwrmsnormvectorarray (xa : t array) (wa : t array) nrm =
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (check x) wa);
    c_n_vwrmsnormvectorarray xa wa nrm

  external c_n_vwrmsnormmaskvectorarray
    : t array -> t array -> t -> RealArray.t -> unit
    = "sunml_nvec_many_n_vwrmsnormmaskvectorarray"

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

  module Local = struct
    external c_n_vdotprod      : t -> t -> float
      = "sunml_nvec_many_n_vdotprodlocal"

    let n_vdotprod (x : t) (y : t) =
      if Sundials_configuration.safe then check x y;
      c_n_vdotprod x y

    external n_vmaxnorm        : t -> float
      = "sunml_nvec_many_n_vmaxnormlocal"

    external n_vmin            : t -> float
      = "sunml_nvec_many_n_vminlocal"

    external n_vl1norm         : t -> float
      = "sunml_nvec_many_n_vl1normlocal"

    external c_n_vinvtest      : t -> t -> bool
      = "sunml_nvec_many_n_vinvtestlocal" [@@noalloc]

    let n_vinvtest (x : t) (z : t) =
      if Sundials_configuration.safe then check x z;
      c_n_vinvtest x z

    external c_n_vconstrmask   : t -> t -> t -> bool
      = "sunml_nvec_many_n_vconstrmasklocal" [@@noalloc]

    let n_vconstrmask (c : t) (x : t) (m : t) =
      if Sundials_configuration.safe then (check c x; check c m);
      c_n_vconstrmask c x m

    external c_n_vminquotient  : t -> t -> float
      = "sunml_nvec_many_n_vminquotientlocal"

    let n_vminquotient (n : t) (d : t) =
      if Sundials_configuration.safe then check n d;
      c_n_vminquotient n d

    external c_n_vwsqrsum
      : t -> t -> float
      = "sunml_nvec_many_n_vwsqrsumlocal"

    let n_vwsqrsum (x : t) (w : t) =
      if Sundials_configuration.safe then check x w;
      c_n_vwsqrsum x w

    external c_n_vwsqrsummask
      : t -> t -> t -> float
      = "sunml_nvec_many_n_vwsqrsummasklocal"

    let n_vwsqrsummask (x : t) (w : t) (id : t) =
      if Sundials_configuration.safe then (check x w; check x id);
      c_n_vwsqrsummask x w id
  end

  let n_vdotprod = Local.n_vdotprod
  let n_vmaxnorm = Local.n_vmaxnorm
  let n_vmin = Local.n_vmin
  let n_vl1norm = Local.n_vl1norm
  let n_vconstrmask = Local.n_vconstrmask
  let n_vminquotient = Local.n_vminquotient

end (* }}} *)

module DataOps : Nvector.NVECTOR_OPS with type t = data =
struct (* {{{ *)

  type t = data

  let n_vclone (nvecs, gl) = (ROArray.map Nvector.clone nvecs, gl)

  module Local = struct
    let n_vdotprod ((x, _) : t) ((y, _) : t) =
      let f sum xi yi = sum +. Nvector.Ops.n_vdotprod xi yi in
      ROArray.fold_left2 f 0. x y

    let n_vmaxnorm ((x, _) : t) =
      let f max xi =
        let maxl =
          if Nvector.Any.Local.has_n_vmaxnorm xi
          then Nvector.Ops.Local.n_vmaxnorm xi
          else Nvector.Ops.n_vmaxnorm xi
        in
        Float.max max maxl
      in
      ROArray.fold_left f 0. x

    let n_vmin ((x, _) : t) =
      let f min xi =
        let minl =
          if Nvector.Any.Local.has_n_vmin xi
          then Nvector.Ops.Local.n_vmin xi
          else Nvector.Ops.n_vmin xi
        in
        Float.min min minl
      in
      ROArray.fold_left f max_float x

    let n_vl1norm ((x, _) : t) =
      let f sum xi = sum +. Nvector.Ops.n_vl1norm xi in
      ROArray.fold_left f 0. x

    let n_vinvtest ((x, _) : t) ((z, _) : t) =
      let f v xi zi =
        v && (if Nvector.Any.Local.has_n_vinvtest xi
              then Nvector.Ops.Local.n_vinvtest xi zi
              else Nvector.Ops.n_vinvtest xi zi)
      in
      ROArray.fold_left2 f true x z

    let n_vconstrmask ((c, _) : t) ((x, _) : t) ((m, _) : t) =
      let f v ci xi mi =
        v && (if Nvector.Any.Local.has_n_vconstrmask ci
              then Nvector.Ops.Local.n_vconstrmask ci xi mi
              else Nvector.Ops.n_vconstrmask ci xi mi)
      in
      ROArray.fold_left3 f true c x m

    let n_vminquotient ((n, _) : t) ((d, _) : t) =
      let f min ni di =
        let minl =
          if Nvector.Any.Local.has_n_vminquotient ni
          then Nvector.Ops.Local.n_vminquotient ni di
          else Nvector.Ops.n_vminquotient ni di
        in
        Float.min min minl
      in
      ROArray.fold_left2 f max_float n d

    let n_vwsqrsum ((x, _) : t) ((w, _) : t) =
      let f sum xi wi =
        let contrib = Nvector.Ops.n_vwrmsnorm xi wi in
        let n = Nvector.Ops.n_vgetlength xi in
        sum +. contrib *. contrib *. float n
      in
      ROArray.fold_left2 f 0. x w

    let n_vwsqrsummask ((x, _) : t) ((w, _) : t) ((id, _) : t) =
      let f sum xi wi idi =
        let contrib = Nvector.Ops.n_vwrmsnormmask xi wi idi in
        let n = Nvector.Ops.n_vgetlength xi in
        sum +. contrib *. contrib *. float n
      in
      ROArray.fold_left3 f 0. x w id
  end

  let n_vlinearsum a ((x, _) : t) b ((y, _) : t) ((z, _) : t) =
    ROArray.iter3 (fun xi yi zi -> Nvector.Ops.n_vlinearsum a xi b yi zi) x y z

  let n_vconst c ((z, _) : t) =
    ROArray.iter (fun zi -> Nvector.Ops.n_vconst c zi) z

  let n_vscale c ((x, _) : t) ((z, _) : t) =
    ROArray.iter2 (fun xi zi -> Nvector.Ops.n_vscale c xi zi) x z

  let n_vaddconst ((x, _) : t) b ((z, _) : t) =
    ROArray.iter2 (fun xi zi -> Nvector.Ops.n_vaddconst xi b zi) x z

  let n_vwrmsnorm ((_, gx) as xv : t) ((_, gw) as wv : t) =
    if gx <> gw then raise Nvector.IncompatibleNvector;
    let gsum = Local.n_vwsqrsum xv wv in
    sqrt (gsum /. float gx)

  let n_vwrmsnormmask ((_, gx) as xv : t) ((_, gw) as wv : t)
                      ((_, gid) as idv : t) =
    if gx <> gw || gx <> gid then raise Nvector.IncompatibleNvector;
    let gsum = Local.n_vwsqrsummask xv wv idv in
    sqrt (gsum /. float gx)

  let n_vdotprod = Local.n_vdotprod
  let n_vmaxnorm = Local.n_vmaxnorm
  let n_vmin = Local.n_vmin
  let n_vl1norm = Local.n_vl1norm
  let n_vconstrmask = Local.n_vconstrmask
  let n_vminquotient = Local.n_vminquotient

  let n_vcompare c ((x, _) : t) ((z, _) : t) =
    ROArray.iter2 (fun xi zi -> Nvector.Ops.n_vcompare c xi zi) x z

  let n_vinvtest = Local.n_vinvtest

  let n_vwl2norm (xv : t) (wv : t) =
    let gsum = Local.n_vwsqrsum xv wv in
    sqrt (gsum)

  let n_vprod ((x, _) : t) ((y, _) : t) ((z, _) : t) =
    ROArray.iter3 (fun xi yi zi -> Nvector.Ops.n_vprod xi yi zi) x y z

  let n_vdiv ((x, _) : t) ((y, _) : t) ((z, _) : t) =
    ROArray.iter3 (fun xi yi zi -> Nvector.Ops.n_vdiv xi yi zi) x y z

  let n_vabs ((x, _) : t) ((z, _) : t) =
    ROArray.iter2 (fun xi zi -> Nvector.Ops.n_vabs xi zi) x z

  let n_vinv ((x, _) : t) ((z, _) : t) =
    ROArray.iter2 (fun xi zi -> Nvector.Ops.n_vinv xi zi) x z

  let n_vspace ((x, _) : t) =
    let f (lrw, liw) xi =
      let lrwi, liwi = Nvector.Ops.n_vspace xi in
      (lrw + lrwi, liw + liwi)
    in
    ROArray.fold_left f (0, 0) x

  let n_vgetlength ((_, gx) : t) = gx

  (* fused and array operations *)

  let n_vlinearcombination (ca : RealArray.t) (xa : t array) ((z, _) : t) =
    let xdata = Array.map fst xa in
    let xsub = Array.map (fun xd -> ROArray.get xd 0) xdata in
    let f i zi =
      Array.iteri (fun j xd -> xsub.(j) <- ROArray.get xd i) xdata;
      Nvector.Ops.n_vlinearcombination ca xsub zi
    in
    ROArray.iteri f z

  let n_vscaleaddmulti (aa : RealArray.t) ((x, _) : t)
                                          (ya : t array) (za : t array) =
    let ydata = Array.map fst ya in
    let ylen = Array.length ydata in
    let zdata = Array.map fst za in
    let zlen = Array.length zdata in
    if ylen <> zlen
      then invalid_arg "n_vscaleaddmulti: arrays must have the same length";
    let ysub = Array.map (fun yd -> ROArray.get yd 0) ydata in
    let zsub = Array.map (fun zd -> ROArray.get zd 0) zdata in
    let f i xi =
      for j = 0 to ylen - 1 do
        ysub.(j) <- ROArray.get ydata.(j) i;
        zsub.(j) <- ROArray.get zdata.(j) i
      done;
      Nvector.Ops.n_vscaleaddmulti aa xi ysub zsub
    in
    ROArray.iteri f x

  let n_vdotprodmulti (x : t) (ya : t array) (dp : RealArray.t) =
    let f i yi = dp.{i} <- n_vdotprod x yi in
    Array.iteri f ya

  let n_vlinearsumvectorarray a (xa : t array) b (ya : t array) (za : t array) =
    let xdata = Array.map fst xa in
    let xlen = Array.length xdata in
    let ydata = Array.map fst ya in
    let ylen = Array.length ydata in
    let zdata = Array.map fst za in
    let zlen = Array.length zdata in
    if xlen <> ylen || xlen <> zlen
      then invalid_arg "n_vlinearsumvectorarray: arrays must have the same length";
    let xsub = Array.map (fun xd -> ROArray.get xd 0) xdata in
    let ysub = Array.map (fun yd -> ROArray.get yd 0) ydata in
    let zsub = Array.map (fun zd -> ROArray.get zd 0) zdata in
    for i = 0 to ROArray.length xdata.(0) - 1 do
      for j = 0 to xlen - 1 do
        xsub.(j) <- ROArray.get xdata.(j) i;
        ysub.(j) <- ROArray.get ydata.(j) i;
        zsub.(j) <- ROArray.get zdata.(j) i
      done;
      Nvector.Ops.n_vlinearsumvectorarray a xsub b ysub zsub
    done

  let n_vscalevectorarray (c : RealArray.t) (xa : t array) (za : t array) =
    let xdata = Array.map fst xa in
    let xlen = Array.length xdata in
    let zdata = Array.map fst za in
    let zlen = Array.length zdata in
    if xlen <> zlen
      then invalid_arg "n_vscalevectorarray: arrays must have the same length";
    let xsub = Array.map (fun xd -> ROArray.get xd 0) xdata in
    let zsub = Array.map (fun zd -> ROArray.get zd 0) zdata in
    for i = 0 to ROArray.length xdata.(0) - 1 do
      for j = 0 to xlen - 1 do
        xsub.(j) <- ROArray.get xdata.(j) i;
        zsub.(j) <- ROArray.get zdata.(j) i
      done;
      Nvector.Ops.n_vscalevectorarray c xsub zsub
    done

  let n_vconstvectorarray c (za : t array) =
    let zdata = Array.map fst za in
    let zlen = Array.length zdata in
    let zsub = Array.map (fun zd -> ROArray.get zd 0) zdata in
    for i = 0 to ROArray.length zdata.(0) - 1 do
      for j = 0 to zlen - 1 do
        zsub.(j) <- ROArray.get zdata.(j) i
      done;
      Nvector.Ops.n_vconstvectorarray c zsub
    done

  let n_vwrmsnormvectorarray (xa : t array) (wa : t array) (nrm : RealArray.t) =
    let xlen = Array.length xa in
    let wlen = Array.length wa in
    if xlen <> wlen || xlen <> RealArray.length nrm
      then invalid_arg "n_vwrmsnormvectorarray: arrays must have the same length";
    for i = 0 to xlen - 1 do
      let _, glen = xa.(i) in
      nrm.{i} <- sqrt (Local.n_vwsqrsum xa.(i) wa.(i) /. float glen)
    done

  let n_vwrmsnormmaskvectorarray (xa : t array) (wa : t array) (id : t)
                                 (nrm : RealArray.t) =
    let xlen = Array.length xa in
    let wlen = Array.length wa in
    if xlen <> wlen || xlen <> RealArray.length nrm
      then invalid_arg "n_vwrmsnormvectorarray: arrays must have the same length";
    for i = 0 to xlen - 1 do
      let _, glen = xa.(i) in
      nrm.{i} <- sqrt (Local.n_vwsqrsummask xa.(i) wa.(i) id /. float glen)
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
      -> data
      -> (Nvector.any -> bool)
      -> (Nvector.any -> Nvector.any)
      -> Nvector.any
    = "sunml_nvec_anywrap_many"

  let rec wrap_with_len ((nvs, _) as payload) =
    if Sundials_impl.Versions.sundials_lt500
      then raise Config.NotImplementedBySundialsVersion;
    let check nv' =
      match unwrap nv' with
      | Many (nvs', _) ->
          (try
             ROArray.iter2 Nvector.check nvs nvs';
             Nvector.get_id nv' = Nvector.ManyVector
           with Nvector.IncompatibleNvector -> false)
      | _ -> false
    in
    c_any_wrap [%extension_constructor Many] payload check clone

  and clone nv =
    let nvs, gl = match unwrap nv with
                  | Many v -> v
                  | _ -> assert false
    in
    wrap_with_len (ROArray.map Nvector.clone nvs, gl)

  let wrap nvs =
    if Sundials_impl.Versions.sundials_lt500
      then raise Config.NotImplementedBySundialsVersion;
    wrap_with_len (nvs, sumlens nvs)

end (* }}} *)

