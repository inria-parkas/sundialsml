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

external c_wrap
  : data
    -> (t -> bool)
    -> (t -> t)
    -> t
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
  let f sum nv = sum + Nvector.Ops.getlength nv in
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

  let clone = Nvector.clone

  external c_linearsum    : float -> t -> float -> t -> t -> unit
    = "sunml_nvec_many_linearsum" [@@noalloc]

  let linearsum a (x : t) b (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_linearsum a x b y z

  external const          : float -> t -> unit
    = "sunml_nvec_many_const" [@@noalloc]

  external c_prod         : t -> t -> t -> unit
    = "sunml_nvec_many_prod" [@@noalloc]

  let prod (x : t) (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_prod x y z

  external c_div          : t -> t -> t -> unit
    = "sunml_nvec_many_div" [@@noalloc]

  let div (x : t) (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_div x y z

  external c_scale        : float -> t -> t -> unit
    = "sunml_nvec_many_scale" [@@noalloc]

  let scale c (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_scale c x z

  external c_abs          : t -> t -> unit
    = "sunml_nvec_many_abs" [@@noalloc]

  let abs (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_abs x z

  external c_inv          : t -> t -> unit
    = "sunml_nvec_many_inv" [@@noalloc]

  let inv (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_inv x z

  external c_addconst     : t -> float -> t -> unit
    = "sunml_nvec_many_addconst" [@@noalloc]

  let addconst (x : t) b (z : t) =
    if Sundials_configuration.safe then check x z;
    c_addconst x b z

  external c_wrmsnorm     : t -> t -> float
    = "sunml_nvec_many_wrmsnorm"

  let wrmsnorm (x : t) (w : t) =
    if Sundials_configuration.safe then check x w;
    c_wrmsnorm x w

  external c_wrmsnormmask : t -> t -> t -> float
    = "sunml_nvec_many_wrmsnormmask"

  let wrmsnormmask (x : t) (w : t) (id : t) =
    if Sundials_configuration.safe then (check x w; check x id);
    c_wrmsnormmask x w id

  external c_wl2norm      : t -> t -> float
    = "sunml_nvec_many_wl2norm"

  let wl2norm (x : t) (w : t) =
    if Sundials_configuration.safe then check x w;
    c_wl2norm x w

  external c_compare      : float -> t -> t -> unit
    = "sunml_nvec_many_compare" [@@noalloc]

  let compare c (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_compare c x z

  external c_invtest      : t -> t -> bool
    = "sunml_nvec_many_invtestlocal" [@@noalloc]

  let invtest (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_invtest x z

  external space          : t -> int * int
    = "sunml_nvec_many_space" [@@noalloc]

  external getlength      : t -> int
    = "sunml_nvec_many_getlength" [@@noalloc]

  external c_linearcombination : RealArray.t -> t array -> t -> unit
    = "sunml_nvec_many_linearcombination"

  let linearcombination ca (xa : t array) (z : t) =
    if Sundials_configuration.safe then Array.iter (check z) xa;
    c_linearcombination ca xa z

  external c_scaleaddmulti : RealArray.t -> t -> t array -> t array -> unit
    = "sunml_nvec_many_scaleaddmulti"

  let scaleaddmulti aa (x : t) (ya : t array) (za : t array) =
    if Sundials_configuration.safe then
      (Array.iter (check x) ya; Array.iter (check x) za);
    c_scaleaddmulti aa x ya za

  external c_dotprodmulti  : t -> t array -> RealArray.t -> unit
    = "sunml_nvec_many_dotprodmulti"

  let dotprodmulti (x : t) (ya : t array) dp =
    if Sundials_configuration.safe then Array.iter (check x) ya;
    c_dotprodmulti x ya dp

  external c_linearsumvectorarray
    : float -> t array -> float -> t array -> t array -> unit
    = "sunml_nvec_many_linearsumvectorarray"

  let linearsumvectorarray a (xa : t array) b (ya : t array) (za : t array) =
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (check x) ya;
          Array.iter (check x) za);
    c_linearsumvectorarray a xa b ya za

  external c_scalevectorarray
    : RealArray.t -> t array -> t array -> unit
    = "sunml_nvec_many_scalevectorarray"

  let scalevectorarray c (xa : t array) (za : t array) =
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (check x) za);
    c_scalevectorarray c xa za

  external c_constvectorarray
    : float -> t array -> unit
    = "sunml_nvec_many_constvectorarray"

  let constvectorarray c (za : t array) =
    if Sundials_configuration.safe
    then (let z = Array.get za 0 in
          Array.iter (check z) za);
    c_constvectorarray c za

  external c_wrmsnormvectorarray
    : t array -> t array -> RealArray.t -> unit
    = "sunml_nvec_many_wrmsnormvectorarray"

  let wrmsnormvectorarray (xa : t array) (wa : t array) nrm =
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (check x) wa);
    c_wrmsnormvectorarray xa wa nrm

  external c_wrmsnormmaskvectorarray
    : t array -> t array -> t -> RealArray.t -> unit
    = "sunml_nvec_many_wrmsnormmaskvectorarray"

  let wrmsnormmaskvectorarray (xa : t array) (wa : t array) (id : t) nrm =
    if Sundials_configuration.safe
    then (Array.iter (check id) xa;
          Array.iter (check id) wa);
    c_wrmsnormmaskvectorarray xa wa id nrm

  let scaleaddmultivectorarray ra (xa : t array) (yaa : t array array)
                                     (zaa : t array array) =
    raise Nvector.OperationNotProvided

  let linearcombinationvectorarray ca (xaa : t array array) (za : t array) =
    raise Nvector.OperationNotProvided

  module Local = struct
    external c_dotprod      : t -> t -> float
      = "sunml_nvec_many_dotprodlocal"

    let dotprod (x : t) (y : t) =
      if Sundials_configuration.safe then check x y;
      c_dotprod x y

    external maxnorm        : t -> float
      = "sunml_nvec_many_maxnormlocal"

    external min            : t -> float
      = "sunml_nvec_many_minlocal"

    external l1norm         : t -> float
      = "sunml_nvec_many_l1normlocal"

    external c_invtest      : t -> t -> bool
      = "sunml_nvec_many_invtestlocal" [@@noalloc]

    let invtest (x : t) (z : t) =
      if Sundials_configuration.safe then check x z;
      c_invtest x z

    external c_constrmask   : t -> t -> t -> bool
      = "sunml_nvec_many_constrmasklocal" [@@noalloc]

    let constrmask (c : t) (x : t) (m : t) =
      if Sundials_configuration.safe then (check c x; check c m);
      c_constrmask c x m

    external c_minquotient  : t -> t -> float
      = "sunml_nvec_many_minquotientlocal"

    let minquotient (n : t) (d : t) =
      if Sundials_configuration.safe then check n d;
      c_minquotient n d

    external c_wsqrsum
      : t -> t -> float
      = "sunml_nvec_many_wsqrsumlocal"

    let wsqrsum (x : t) (w : t) =
      if Sundials_configuration.safe then check x w;
      c_wsqrsum x w

    external c_wsqrsummask
      : t -> t -> t -> float
      = "sunml_nvec_many_wsqrsummasklocal"

    let wsqrsummask (x : t) (w : t) (id : t) =
      if Sundials_configuration.safe then (check x w; check x id);
      c_wsqrsummask x w id
  end

  let dotprod = Local.dotprod
  let maxnorm = Local.maxnorm
  let min = Local.min
  let l1norm = Local.l1norm
  let constrmask = Local.constrmask
  let minquotient = Local.minquotient

end (* }}} *)

module DataOps : Nvector.NVECTOR_OPS with type t = data =
struct (* {{{ *)

  type t = data

  let clone (nvecs, gl) = (ROArray.map Nvector.clone nvecs, gl)
  let floatmin (x : float) (y : float) = min x y
  let floatmax (x : float) (y : float) = max x y

  module Local = struct
    let dotprod ((x, _) : t) ((y, _) : t) =
      let f sum xi yi = sum +. Nvector.Ops.dotprod xi yi in
      ROArray.fold_left2 f 0. x y

    let maxnorm ((x, _) : t) =
      let f max xi =
        let maxl =
          if Nvector.Any.Local.has_maxnorm xi
          then Nvector.Ops.Local.maxnorm xi
          else Nvector.Ops.maxnorm xi
        in
        floatmax max maxl
      in
      ROArray.fold_left f 0. x

    let min ((x, _) : t) =
      let f min xi =
        let minl =
          if Nvector.Any.Local.has_min xi
          then Nvector.Ops.Local.min xi
          else Nvector.Ops.min xi
        in
        floatmin min minl
      in
      ROArray.fold_left f max_float x

    let l1norm ((x, _) : t) =
      let f sum xi = sum +. Nvector.Ops.l1norm xi in
      ROArray.fold_left f 0. x

    let invtest ((x, _) : t) ((z, _) : t) =
      let f v xi zi =
        v && (if Nvector.Any.Local.has_invtest xi
              then Nvector.Ops.Local.invtest xi zi
              else Nvector.Ops.invtest xi zi)
      in
      ROArray.fold_left2 f true x z

    let constrmask ((c, _) : t) ((x, _) : t) ((m, _) : t) =
      let f v ci xi mi =
        v && (if Nvector.Any.Local.has_constrmask ci
              then Nvector.Ops.Local.constrmask ci xi mi
              else Nvector.Ops.constrmask ci xi mi)
      in
      ROArray.fold_left3 f true c x m

    let minquotient ((n, _) : t) ((d, _) : t) =
      let f min ni di =
        let minl =
          if Nvector.Any.Local.has_minquotient ni
          then Nvector.Ops.Local.minquotient ni di
          else Nvector.Ops.minquotient ni di
        in
        floatmin min minl
      in
      ROArray.fold_left2 f max_float n d

    let wsqrsum ((x, _) : t) ((w, _) : t) =
      let f sum xi wi =
        let contrib = Nvector.Ops.wrmsnorm xi wi in
        let n = Nvector.Ops.getlength xi in
        sum +. contrib *. contrib *. float n
      in
      ROArray.fold_left2 f 0. x w

    let wsqrsummask ((x, _) : t) ((w, _) : t) ((id, _) : t) =
      let f sum xi wi idi =
        let contrib = Nvector.Ops.wrmsnormmask xi wi idi in
        let n = Nvector.Ops.getlength xi in
        sum +. contrib *. contrib *. float n
      in
      ROArray.fold_left3 f 0. x w id
  end

  let linearsum a ((x, _) : t) b ((y, _) : t) ((z, _) : t) =
    ROArray.iter3 (fun xi yi zi -> Nvector.Ops.linearsum a xi b yi zi) x y z

  let const c ((z, _) : t) =
    ROArray.iter (fun zi -> Nvector.Ops.const c zi) z

  let scale c ((x, _) : t) ((z, _) : t) =
    ROArray.iter2 (fun xi zi -> Nvector.Ops.scale c xi zi) x z

  let addconst ((x, _) : t) b ((z, _) : t) =
    ROArray.iter2 (fun xi zi -> Nvector.Ops.addconst xi b zi) x z

  let wrmsnorm ((_, gx) as xv : t) ((_, gw) as wv : t) =
    if gx <> gw then raise Nvector.IncompatibleNvector;
    let gsum = Local.wsqrsum xv wv in
    sqrt (gsum /. float gx)

  let wrmsnormmask ((_, gx) as xv : t) ((_, gw) as wv : t)
                      ((_, gid) as idv : t) =
    if gx <> gw || gx <> gid then raise Nvector.IncompatibleNvector;
    let gsum = Local.wsqrsummask xv wv idv in
    sqrt (gsum /. float gx)

  let dotprod = Local.dotprod
  let maxnorm = Local.maxnorm
  let min = Local.min
  let l1norm = Local.l1norm
  let constrmask = Local.constrmask
  let minquotient = Local.minquotient

  let compare c ((x, _) : t) ((z, _) : t) =
    ROArray.iter2 (fun xi zi -> Nvector.Ops.compare c xi zi) x z

  let invtest = Local.invtest

  let wl2norm (xv : t) (wv : t) =
    let gsum = Local.wsqrsum xv wv in
    sqrt (gsum)

  let prod ((x, _) : t) ((y, _) : t) ((z, _) : t) =
    ROArray.iter3 (fun xi yi zi -> Nvector.Ops.prod xi yi zi) x y z

  let div ((x, _) : t) ((y, _) : t) ((z, _) : t) =
    ROArray.iter3 (fun xi yi zi -> Nvector.Ops.div xi yi zi) x y z

  let abs ((x, _) : t) ((z, _) : t) =
    ROArray.iter2 (fun xi zi -> Nvector.Ops.abs xi zi) x z

  let inv ((x, _) : t) ((z, _) : t) =
    ROArray.iter2 (fun xi zi -> Nvector.Ops.inv xi zi) x z

  let space ((x, _) : t) =
    let f (lrw, liw) xi =
      let lrwi, liwi = Nvector.Ops.space xi in
      (lrw + lrwi, liw + liwi)
    in
    ROArray.fold_left f (0, 0) x

  let getlength ((_, gx) : t) = gx

  (* fused and array operations *)

  let linearcombination (ca : RealArray.t) (xa : t array) ((z, _) : t) =
    let xdata = Array.map fst xa in
    let xsub = Array.map (fun xd -> ROArray.get xd 0) xdata in
    let f i zi =
      Array.iteri (fun j xd -> xsub.(j) <- ROArray.get xd i) xdata;
      Nvector.Ops.linearcombination ca xsub zi
    in
    ROArray.iteri f z

  let scaleaddmulti (aa : RealArray.t) ((x, _) : t)
                                          (ya : t array) (za : t array) =
    let ydata = Array.map fst ya in
    let ylen = Array.length ydata in
    let zdata = Array.map fst za in
    let zlen = Array.length zdata in
    if ylen <> zlen
      then invalid_arg "scaleaddmulti: arrays must have the same length";
    let ysub = Array.map (fun yd -> ROArray.get yd 0) ydata in
    let zsub = Array.map (fun zd -> ROArray.get zd 0) zdata in
    let f i xi =
      for j = 0 to ylen - 1 do
        ysub.(j) <- ROArray.get ydata.(j) i;
        zsub.(j) <- ROArray.get zdata.(j) i
      done;
      Nvector.Ops.scaleaddmulti aa xi ysub zsub
    in
    ROArray.iteri f x

  let dotprodmulti (x : t) (ya : t array) (dp : RealArray.t) =
    let f i yi = dp.{i} <- dotprod x yi in
    Array.iteri f ya

  let linearsumvectorarray a (xa : t array) b (ya : t array) (za : t array) =
    let xdata = Array.map fst xa in
    let xlen = Array.length xdata in
    let ydata = Array.map fst ya in
    let ylen = Array.length ydata in
    let zdata = Array.map fst za in
    let zlen = Array.length zdata in
    if xlen <> ylen || xlen <> zlen
      then invalid_arg "linearsumvectorarray: arrays must have the same length";
    let xsub = Array.map (fun xd -> ROArray.get xd 0) xdata in
    let ysub = Array.map (fun yd -> ROArray.get yd 0) ydata in
    let zsub = Array.map (fun zd -> ROArray.get zd 0) zdata in
    for i = 0 to ROArray.length xdata.(0) - 1 do
      for j = 0 to xlen - 1 do
        xsub.(j) <- ROArray.get xdata.(j) i;
        ysub.(j) <- ROArray.get ydata.(j) i;
        zsub.(j) <- ROArray.get zdata.(j) i
      done;
      Nvector.Ops.linearsumvectorarray a xsub b ysub zsub
    done

  let scalevectorarray (c : RealArray.t) (xa : t array) (za : t array) =
    let xdata = Array.map fst xa in
    let xlen = Array.length xdata in
    let zdata = Array.map fst za in
    let zlen = Array.length zdata in
    if xlen <> zlen
      then invalid_arg "scalevectorarray: arrays must have the same length";
    let xsub = Array.map (fun xd -> ROArray.get xd 0) xdata in
    let zsub = Array.map (fun zd -> ROArray.get zd 0) zdata in
    for i = 0 to ROArray.length xdata.(0) - 1 do
      for j = 0 to xlen - 1 do
        xsub.(j) <- ROArray.get xdata.(j) i;
        zsub.(j) <- ROArray.get zdata.(j) i
      done;
      Nvector.Ops.scalevectorarray c xsub zsub
    done

  let constvectorarray c (za : t array) =
    let zdata = Array.map fst za in
    let zlen = Array.length zdata in
    let zsub = Array.map (fun zd -> ROArray.get zd 0) zdata in
    for i = 0 to ROArray.length zdata.(0) - 1 do
      for j = 0 to zlen - 1 do
        zsub.(j) <- ROArray.get zdata.(j) i
      done;
      Nvector.Ops.constvectorarray c zsub
    done

  let wrmsnormvectorarray (xa : t array) (wa : t array) (nrm : RealArray.t) =
    let xlen = Array.length xa in
    let wlen = Array.length wa in
    if xlen <> wlen || xlen <> RealArray.length nrm
      then invalid_arg "wrmsnormvectorarray: arrays must have the same length";
    for i = 0 to xlen - 1 do
      let _, glen = xa.(i) in
      nrm.{i} <- sqrt (Local.wsqrsum xa.(i) wa.(i) /. float glen)
    done

  let wrmsnormmaskvectorarray (xa : t array) (wa : t array) (id : t)
                                 (nrm : RealArray.t) =
    let xlen = Array.length xa in
    let wlen = Array.length wa in
    if xlen <> wlen || xlen <> RealArray.length nrm
      then invalid_arg "wrmsnormvectorarray: arrays must have the same length";
    for i = 0 to xlen - 1 do
      let _, glen = xa.(i) in
      nrm.{i} <- sqrt (Local.wsqrsummask xa.(i) wa.(i) id /. float glen)
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

  let unwrap nv =
    match Nvector.unwrap nv with
    | Many a -> a
    | _ -> raise Nvector.BadGenericType

end (* }}} *)

