(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

type kind
type 'a t = ('a, kind) Nvector.t

type communicator

(* Must match with nvector_ml.h:nvector_ops_tag *)
type 'a nvector_ops = { (* {{{ *)
  check           : 'a -> 'a -> bool;
  clone           : 'a -> 'a;
  space           : ('a -> (int * int)) option;
  getlength       : 'a -> int;
  print           : ('a -> Sundials.Logfile.t option-> unit) option;
  linearsum       : float -> 'a -> float -> 'a -> 'a -> unit;
  const           : float -> 'a -> unit;
  prod            : 'a -> 'a -> 'a -> unit;
  div             : 'a -> 'a -> 'a -> unit;
  scale           : float -> 'a -> 'a -> unit;
  abs             : 'a -> 'a -> unit;
  inv             : 'a -> 'a -> unit;
  addconst        : 'a -> float -> 'a -> unit;
  maxnorm         : 'a -> float;
  wrmsnorm        : 'a -> 'a -> float;
  min             : 'a -> float;

  dotprod         : 'a -> 'a -> float;
  compare         : float -> 'a -> 'a -> unit;
  invtest         : 'a -> 'a -> bool;

  wl2norm         : ('a -> 'a -> float) option;
  l1norm          : ('a -> float) option;
  wrmsnormmask    : ('a -> 'a -> 'a -> float) option;
  constrmask      : ('a -> 'a -> 'a -> bool) option;
  minquotient     : ('a -> 'a -> float) option;

  getcommunicator : ('a -> communicator) option;

  linearcombination :
    (Sundials.RealArray.t -> 'a array -> 'a -> unit) option;
  scaleaddmulti :
    (Sundials.RealArray.t -> 'a -> 'a array -> 'a array -> unit) option;
  dotprodmulti : ('a -> 'a array -> Sundials.RealArray.t -> unit) option;

  linearsumvectorarray :
    (float -> 'a array -> float -> 'a array -> 'a array -> unit) option;
  scalevectorarray :
    (Sundials.RealArray.t -> 'a array -> 'a array -> unit) option;
  constvectorarray : (float -> 'a array -> unit) option;
  wrmsnormvectorarray :
    ('a array -> 'a array -> Sundials.RealArray.t -> unit) option;
  wrmsnormmaskvectorarray :
    ('a array -> 'a array -> 'a -> Sundials.RealArray.t -> unit) option;
  scaleaddmultivectorarray :
    (Sundials.RealArray.t -> 'a array -> 'a array array -> 'a array array ->
      unit) option;
  linearcombinationvectorarray :
    (Sundials.RealArray.t -> 'a array array -> 'a array -> unit) option;

  dotprod_local      : ('a -> 'a -> float) option;
  maxnorm_local      : ('a -> float) option;
  min_local          : ('a -> float) option;
  l1norm_local       : ('a -> float) option;
  invtest_local      : ('a -> 'a -> bool) option;
  constrmask_local   : ('a -> 'a -> 'a -> bool) option;
  minquotient_local  : ('a -> 'a -> float) option;
  wsqrsum_local      : ('a -> 'a -> float) option;
  wsqrsummask_local  : ('a -> 'a -> 'a -> float) option;

} (* }}} *)

(* Selectively enable and disable fused and array operations *)
external c_enablefusedops_custom                     : ('d, 'k) Nvector.t -> bool -> bool
  = "sunml_nvec_custom_enablefusedops"
external c_enablelinearcombination_custom            : ('d, 'k) Nvector.t -> bool -> bool
  = "sunml_nvec_custom_enablelinearcombination"
external c_enablescaleaddmulti_custom                : ('d, 'k) Nvector.t -> bool -> bool
  = "sunml_nvec_custom_enablescaleaddmulti"
external c_enabledotprodmulti_custom                 : ('d, 'k) Nvector.t -> bool -> bool
  = "sunml_nvec_custom_enabledotprodmulti"
external c_enablelinearsumvectorarray_custom         : ('d, 'k) Nvector.t -> bool -> bool
  = "sunml_nvec_custom_enablelinearsumvectorarray"
external c_enablescalevectorarray_custom             : ('d, 'k) Nvector.t -> bool -> bool
  = "sunml_nvec_custom_enablescalevectorarray"
external c_enableconstvectorarray_custom             : ('d, 'k) Nvector.t -> bool -> bool
  = "sunml_nvec_custom_enableconstvectorarray"
external c_enablewrmsnormvectorarray_custom          : ('d, 'k) Nvector.t -> bool -> bool
  = "sunml_nvec_custom_enablewrmsnormvectorarray"
external c_enablewrmsnormmaskvectorarray_custom      : ('d, 'k) Nvector.t -> bool -> bool
  = "sunml_nvec_custom_enablewrmsnormmaskvectorarray"
external c_enablescaleaddmultivectorarray_custom     : ('d, 'k) Nvector.t -> bool -> bool
  = "sunml_nvec_custom_enablescaleaddmultivectorarray"
external c_enablelinearcombinationvectorarray_custom : ('d, 'k) Nvector.t -> bool -> bool
  = "sunml_nvec_custom_enablelinearcombinationvectorarray"

external c_make_wrap
    : 'a nvector_ops -> 'a -> ('a t -> bool) -> ('t -> 't) -> 'a t
    = "sunml_nvec_wrap_custom"

let do_enable f nv v =
  match v with
  | None -> ()
  | Some v -> if not (f nv v) then raise Nvector.OperationNotProvided

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
  = do_enable c_enablefusedops_custom nv
              with_fused_ops;
    do_enable c_enablelinearcombination_custom nv
              with_linear_combination;
    do_enable c_enablescaleaddmulti_custom nv
              with_scale_add_multi;
    do_enable c_enabledotprodmulti_custom nv
              with_dot_prod_multi;
    do_enable c_enablelinearsumvectorarray_custom nv
              with_linear_sum_vector_array;
    do_enable c_enablescalevectorarray_custom nv
              with_scale_vector_array;
    do_enable c_enableconstvectorarray_custom nv
              with_const_vector_array;
    do_enable c_enablewrmsnormvectorarray_custom nv
              with_wrms_norm_vector_array;
    do_enable c_enablewrmsnormmaskvectorarray_custom nv
              with_wrms_norm_mask_vector_array;
    do_enable c_enablescaleaddmultivectorarray_custom nv
              with_scale_add_multi_vector_array;
    do_enable c_enablelinearcombinationvectorarray_custom nv
              with_linear_combination_vector_array

let uv = Nvector.unwrap

let rec make_wrap ops ?(with_fused_ops=false) v =
  let check nv' = ops.check v (uv nv') in
  let nv = c_make_wrap ops v check (clone ops) in
  if with_fused_ops && not (c_enablefusedops_custom nv true)
    then raise Nvector.OperationNotProvided;
  nv

and clone ops nv =
  let nv' = make_wrap ops (ops.clone (uv nv)) in
  ignore (c_enablelinearcombination_custom nv'
            (Nvector.Ops.has_linearcombination nv));
  ignore (c_enablescaleaddmulti_custom nv'
            (Nvector.Ops.has_scaleaddmulti nv));
  ignore (c_enabledotprodmulti_custom nv'
            (Nvector.Ops.has_dotprodmulti nv));
  ignore (c_enablelinearsumvectorarray_custom nv'
            (Nvector.Ops.has_linearsumvectorarray nv));
  ignore (c_enablescalevectorarray_custom nv'
            (Nvector.Ops.has_scalevectorarray nv));
  ignore (c_enableconstvectorarray_custom nv'
            (Nvector.Ops.has_constvectorarray nv));
  ignore (c_enablewrmsnormvectorarray_custom nv'
            (Nvector.Ops.has_wrmsnormvectorarray nv));
  ignore (c_enablewrmsnormmaskvectorarray_custom nv'
            (Nvector.Ops.has_wrmsnormmaskvectorarray nv));
  ignore (c_enablescaleaddmultivectorarray_custom nv'
            (Nvector.Ops.has_scaleaddmultivectorarray nv));
  ignore (c_enablelinearcombinationvectorarray_custom nv'
            (Nvector.Ops.has_linearcombinationvectorarray nv));
  nv'

let add_tracing msg ops =
  let pr s = print_string msg; print_endline s in
  let { (* {{{ *)
      check           = check;
      clone           = clone;
      space           = space;
      getlength       = getlength;
      print           = print;
      linearsum       = linearsum;
      const           = const;
      prod            = prod;
      div             = div;
      scale           = scale;
      abs             = abs;
      inv             = inv;
      addconst        = addconst;
      maxnorm         = maxnorm;
      wrmsnorm        = wrmsnorm;
      min             = min;

      dotprod         = dotprod;
      compare         = compare;
      invtest         = invtest;

      wl2norm         = wl2norm;
      l1norm          = l1norm;
      wrmsnormmask    = wrmsnormmask;
      constrmask      = constrmask;
      minquotient     = minquotient;

      getcommunicator = getcommunicator;

      linearcombination            = linearcombination;
      scaleaddmulti                = scaleaddmulti;
      dotprodmulti                 = dotprodmulti;

      linearsumvectorarray         = linearsumvectorarray;
      scalevectorarray             = scalevectorarray;
      constvectorarray             = constvectorarray;
      wrmsnormvectorarray          = wrmsnormvectorarray;
      wrmsnormmaskvectorarray      = wrmsnormmaskvectorarray;
      scaleaddmultivectorarray     = scaleaddmultivectorarray;
      linearcombinationvectorarray = linearcombinationvectorarray;

      dotprod_local     = dotprod_local;
      maxnorm_local     = maxnorm_local;
      min_local         = min_local;
      l1norm_local      = l1norm_local;
      invtest_local     = invtest_local;
      constrmask_local  = constrmask_local;
      minquotient_local = minquotient_local;
      wsqrsum_local     = wsqrsum_local;
      wsqrsummask_local = wsqrsummask_local;
    } = ops (* }}} *)
  in
  let fo f f' = match f with None -> None | Some f -> Some (f' f) in

  let tr_nvcheck x =
    pr "check-create";
    let check = check x in
    function y ->
      pr "check-check";
      check y
  and tr_nvclone a = pr "clone"; clone a
  (* ... {{{ *)
  and tr_nvspace = fo space (fun f -> fun a -> (pr "space"; f a))
  and tr_nvgetlength a = pr "getlength"; getlength a
  and tr_nvprint = match print with None -> None
                   | Some f -> Some (fun a lf -> pr "print"; f a lf);
  and tr_nvlinearsum a x b y z = pr "linearsum"; linearsum a x b y z
  and tr_nvconst c z = pr "const"; const c z
  and tr_nvprod x y z = pr "prod"; prod x y z
  and tr_nvdiv x y z = pr "div"; div x y z
  and tr_nvscale c x z = pr "scale"; scale c x z
  and tr_nvabs x z = pr "abs"; abs x z
  and tr_nvinv x z = pr "inv"; inv x z
  and tr_nvaddconst x b z = pr "addconst"; addconst x b z
  and tr_nvmaxnorm x = pr "maxnorm"; maxnorm x
  and tr_nvwrmsnorm x w = pr "wrmsnorm"; wrmsnorm x w
  and tr_nvmin x = pr "min"; min x
  and tr_nvdotprod x y = pr "dotprod"; dotprod x y
  and tr_nvcompare c x z = pr "compare"; compare c x z
  and tr_nvinvtest x z = pr "invtest"; invtest x z

  and tr_nvwl2norm = fo wl2norm (fun f -> fun x w -> pr "wl2norm"; f x w)
  and tr_nvl1norm = fo l1norm (fun f -> fun x -> pr "l1norm"; f x)
  and tr_nvwrmsnormmask =
    fo wrmsnormmask (fun f -> fun x w id -> pr "wrmsnormmask"; f x w id)
  and tr_nvconstrmask =
    fo constrmask (fun f -> fun c x m -> pr "constrmask"; f c x m)
  and tr_nvminquotient =
    fo minquotient (fun f -> fun n d -> pr "minquotient"; f n d)

  and tr_nvgetcommunicator =
    fo getcommunicator (fun f -> fun x -> pr "getcommunicator"; f x)

  and tr_nvlinearcombination =
    fo linearcombination
    (fun f -> fun c x z -> pr "linearcombination"; f c x z)
  and tr_nvscaleaddmulti =
    fo scaleaddmulti
    (fun f -> fun c x y z -> pr "scaleaddmulti"; f c x y z)
  and tr_nvdotprodmulti =
    fo dotprodmulti
    (fun f -> fun x y d -> pr "dotprodmulti"; f x y d)

  and tr_nvlinearsumvectorarray =
    fo linearsumvectorarray
    (fun f -> fun a x b y z -> pr "linearsumvectorarray"; f a x b y z)
  and tr_nvscalevectorarray =
    fo scalevectorarray
    (fun f -> fun c x z -> pr "scalevectorarray"; f c x z)
  and tr_nvconstvectorarray =
    fo constvectorarray
    (fun f -> fun c x -> pr "constvectorarray"; f c x)
  and tr_nvwrmsnormvectorarray =
    fo wrmsnormvectorarray
    (fun f -> fun x w m -> pr "wrmsnormvectorarray"; f x w m)
  and tr_nvwrmsnormmaskvectorarray =
    fo wrmsnormmaskvectorarray
    (fun f -> fun x w m -> pr "wrmsnormmaskvectorarray"; f x w m)
  and tr_nvscaleaddmultivectorarray =
    fo scaleaddmultivectorarray
    (fun f -> fun a x yy zz -> pr "scaleaddmultivectorarray"; f a x yy zz)
  and tr_nvlinearcombinationvectorarray =
    fo linearcombinationvectorarray
    (fun f -> fun c xx z -> pr "linearcombinationvectorarray"; f c xx z)

  and tr_nvdotprod_local     =
    fo dotprod_local (fun f -> fun x t -> pr "dotprod_local"; f x t)
  and tr_nvmaxnorm_local     =
    fo maxnorm_local (fun f -> fun x -> pr "maxnorm_local"; f x)
  and tr_nvmin_local         =
    fo min_local (fun f -> fun x -> pr "min_local"; f x)
  and tr_nvl1norm_local      =
    fo l1norm_local (fun f -> fun x -> pr "l1norm_local"; f x)
  and tr_nvinvtest_local     =
    fo invtest_local (fun f -> fun x z -> pr "invtest_local"; f x z)
  and tr_nvconstrmask_local  =
    fo constrmask_local (fun f -> fun c x m -> pr "constrmask_local"; f c x m)
  and tr_nvminquotient_local =
    fo minquotient_local (fun f -> fun n d -> pr "minquotient_local"; f n d)
  and tr_nvwsqrsum_local     =
    fo wsqrsum_local (fun f -> fun x w -> pr "wsqrsum_local"; f x w)
  and tr_nvwsqrsummask_local =
    fo wsqrsummask_local (fun f -> fun x w id -> pr "wsqrsummask_local"; f x w id)
  (* }}} *)
  in
  {
      (* {{{ *)
      check           = tr_nvcheck;
      clone           = tr_nvclone;
      space           = tr_nvspace;
      getlength       = tr_nvgetlength;
      print           = tr_nvprint;
      linearsum       = tr_nvlinearsum;
      const           = tr_nvconst;
      prod            = tr_nvprod;
      div             = tr_nvdiv;
      scale           = tr_nvscale;
      abs             = tr_nvabs;
      inv             = tr_nvinv;
      addconst        = tr_nvaddconst;
      maxnorm         = tr_nvmaxnorm;
      wrmsnorm        = tr_nvwrmsnorm;
      min             = tr_nvmin;

      dotprod         = tr_nvdotprod;
      compare         = tr_nvcompare;
      invtest         = tr_nvinvtest;

      wl2norm         = tr_nvwl2norm;
      l1norm          = tr_nvl1norm;
      wrmsnormmask    = tr_nvwrmsnormmask;
      constrmask      = tr_nvconstrmask;
      minquotient     = tr_nvminquotient;

      getcommunicator = tr_nvgetcommunicator;

      linearcombination            = tr_nvlinearcombination;
      scaleaddmulti                = tr_nvscaleaddmulti;
      dotprodmulti                 = tr_nvdotprodmulti;

      linearsumvectorarray         = tr_nvlinearsumvectorarray;
      scalevectorarray             = tr_nvscalevectorarray;
      constvectorarray             = tr_nvconstvectorarray;
      wrmsnormvectorarray          = tr_nvwrmsnormvectorarray;
      wrmsnormmaskvectorarray      = tr_nvwrmsnormmaskvectorarray;
      scaleaddmultivectorarray     = tr_nvscaleaddmultivectorarray;
      linearcombinationvectorarray = tr_nvlinearcombinationvectorarray;

      dotprod_local     = tr_nvdotprod_local;
      maxnorm_local     = tr_nvmaxnorm_local;
      min_local         = tr_nvmin_local;
      l1norm_local      = tr_nvl1norm_local;
      invtest_local     = tr_nvinvtest_local;
      constrmask_local  = tr_nvconstrmask_local;
      minquotient_local = tr_nvminquotient_local;
      wsqrsum_local     = tr_nvwsqrsum_local;
      wsqrsummask_local = tr_nvwsqrsummask_local;
      (* }}} *)
  }

module MakeOps = functor (A : sig
    type data
    val ops : data nvector_ops
  end) -> struct
    type _kind = kind
    type kind = _kind
    type data = A.data
    type t = (data, kind) Nvector.t

    let wrap = make_wrap A.ops
    let enable = enable

    module Ops = struct (* {{{ *)
      type t = (data, kind) Nvector.t

      let clone n = wrap (A.ops.clone (uv n))
      let linearsum a x b y z = A.ops.linearsum a (uv x) b (uv y) (uv z)
      let const c z = A.ops.const c (uv z)
      let prod x y z = A.ops.prod (uv x) (uv y) (uv z)
      let div x y z = A.ops.div (uv x) (uv y) (uv z)
      let scale c x z = A.ops.scale c (uv x) (uv z)
      let abs x z = A.ops.abs (uv x) (uv z)
      let inv x z = A.ops.inv (uv x) (uv z)
      let addconst x b z = A.ops.addconst (uv x) b (uv z)
      let dotprod x y = A.ops.dotprod (uv x) (uv y)
      let maxnorm x = A.ops.maxnorm (uv x)
      let wrmsnorm x w = A.ops.wrmsnorm (uv x) (uv w)
      let min x = A.ops.min (uv x)
      let compare c x z = A.ops.compare c (uv x) (uv z)
      let invtest x z = A.ops.invtest (uv x) (uv z)

      let wl2norm =
        match A.ops.wl2norm with
        | None -> (fun x w -> raise Nvector.OperationNotProvided)
        | Some f -> (fun x w -> f (uv x) (uv w))

      let l1norm =
        match A.ops.l1norm with
        | None -> (fun x -> raise Nvector.OperationNotProvided)
        | Some f -> (fun x -> f (uv x))

      let wrmsnormmask =
        match A.ops.wrmsnormmask with
        | None -> (fun x w id -> raise Nvector.OperationNotProvided)
        | Some f -> (fun x w id -> f (uv x) (uv w) (uv id))

      let constrmask =
        match A.ops.constrmask with
        | None -> (fun c x m -> raise Nvector.OperationNotProvided)
        | Some f -> (fun c x m -> f (uv c) (uv x) (uv m))

      let minquotient =
        match A.ops.minquotient with
        | None -> (fun num denom -> raise Nvector.OperationNotProvided)
        | Some f -> (fun num denom -> f (uv num) (uv denom))

      let space =
        match A.ops.space with
        | None -> (fun x -> raise Nvector.OperationNotProvided)
        | Some f -> (fun x -> f (uv x))

      let getlength a = A.ops.getlength (uv a)

      let print =
        match A.ops.print with
        | None -> (fun ?logfile _ -> raise Nvector.OperationNotProvided)
        | Some f -> (fun ?logfile x -> f (uv x) logfile)

      let linearcombination =
        match A.ops.linearcombination with
        | None -> (fun c x z -> raise Nvector.OperationNotProvided)
        | Some f -> (fun c x z -> f c (Array.map uv x) (uv z))

      let scaleaddmulti =
        match A.ops.scaleaddmulti with
        | None -> (fun c x y z -> raise Nvector.OperationNotProvided)
        | Some f ->
            (fun c x y z -> f c (uv x) (Array.map uv y) (Array.map uv z))

      let dotprodmulti =
        match A.ops.dotprodmulti with
        | None -> (fun x y d -> raise Nvector.OperationNotProvided)
        | Some f -> (fun x y d -> f (uv x) (Array.map uv y) d)

      let linearsumvectorarray =
        match A.ops.linearsumvectorarray with
        | None -> (fun a x b y z -> raise Nvector.OperationNotProvided)
        | Some f -> (fun a x b y z -> f a (Array.map uv x) b (Array.map uv y)
                                          (Array.map uv z))

      let scalevectorarray =
        match A.ops.scalevectorarray with
        | None -> (fun c x z -> raise Nvector.OperationNotProvided)
        | Some f -> (fun c x z -> f c (Array.map uv x) (Array.map uv z))

      let constvectorarray =
        match A.ops.constvectorarray with
        | None -> (fun c x -> raise Nvector.OperationNotProvided)
        | Some f -> (fun c x -> f c (Array.map uv x))

      let wrmsnormvectorarray =
        match A.ops.wrmsnormvectorarray with
        | None -> (fun x w m -> raise Nvector.OperationNotProvided)
        | Some f -> (fun x w m -> f (Array.map uv x) (Array.map uv w) m)

      let wrmsnormmaskvectorarray =
        match A.ops.wrmsnormmaskvectorarray with
        | None -> (fun x w id m -> raise Nvector.OperationNotProvided)
        | Some f -> (fun x w id m -> f (Array.map uv x) (Array.map uv w)
                                       (uv id) m)

      let scaleaddmultivectorarray =
        match A.ops.scaleaddmultivectorarray with
        | None -> (fun a x yy zz -> raise Nvector.OperationNotProvided)
        | Some f -> (fun a x yy zz -> f a (Array.map uv x)
                                        (Array.map (Array.map uv) yy)
                                        (Array.map (Array.map uv) zz))

      let linearcombinationvectorarray =
        match A.ops.linearcombinationvectorarray with
        | None -> (fun c xx z -> raise Nvector.OperationNotProvided)
        | Some f -> (fun c xx z -> f c (Array.map (Array.map uv) xx)
                                       (Array.map uv z))

      module Local = struct
        let dotprod =
          match A.ops.dotprod_local with
          | None -> (fun _ _ -> raise Nvector.OperationNotProvided)
          | Some f -> (fun x t -> f (uv x) (uv t))

        let maxnorm =
          match A.ops.maxnorm_local with
          | None -> (fun _ -> raise Nvector.OperationNotProvided)
          | Some f -> (fun x -> f (uv x))

        let min =
          match A.ops.min_local with
          | None -> (fun _ -> raise Nvector.OperationNotProvided)
          | Some f -> (fun x -> f (uv x))

        let l1norm =
          match A.ops.l1norm_local with
          | None -> (fun _ -> raise Nvector.OperationNotProvided)
          | Some f -> (fun x -> f (uv x))

        let invtest =
          match A.ops.invtest_local with
          | None -> (fun _ _ -> raise Nvector.OperationNotProvided)
          | Some f -> (fun x z -> f (uv x) (uv z))

        let constrmask =
          match A.ops.constrmask_local with
          | None -> (fun _ _ _ -> raise Nvector.OperationNotProvided)
          | Some f -> (fun c x m -> f (uv c) (uv x) (uv m))

        let minquotient =
          match A.ops.minquotient_local with
          | None -> (fun _ _ -> raise Nvector.OperationNotProvided)
          | Some f -> (fun n d -> f (uv n) (uv d))

        let wsqrsum =
          match A.ops.wsqrsum_local with
          | None -> (fun _ _ -> raise Nvector.OperationNotProvided)
          | Some f -> (fun x w -> f (uv x) (uv w))

        let wsqrsummask =
          match A.ops.wsqrsummask_local with
          | None -> (fun _ _ _ -> raise Nvector.OperationNotProvided)
          | Some f -> (fun x w id -> f (uv x) (uv w) (uv id))
      end
    end (* }}} *)

    module DataOps = struct (* {{{ *)
      type t = data

      let clone        = A.ops.clone
      let linearsum    = A.ops.linearsum
      let const        = A.ops.const
      let prod         = A.ops.prod
      let div          = A.ops.div
      let scale        = A.ops.scale
      let abs          = A.ops.abs
      let inv          = A.ops.inv
      let addconst     = A.ops.addconst
      let dotprod      = A.ops.dotprod
      let maxnorm      = A.ops.maxnorm
      let wrmsnorm     = A.ops.wrmsnorm
      let min          = A.ops.min
      let compare      = A.ops.compare
      let invtest      = A.ops.invtest

      let wl2norm =
        match A.ops.wl2norm with
        | None -> (fun x w -> raise Nvector.OperationNotProvided)
        | Some f -> f

      let l1norm =
        match A.ops.l1norm with
        | None -> (fun x -> raise Nvector.OperationNotProvided)
        | Some f -> f

      let wrmsnormmask =
        match A.ops.wrmsnormmask with
        | None -> (fun x w id -> raise Nvector.OperationNotProvided)
        | Some f -> f

      let constrmask =
        match A.ops.constrmask with
        | None -> (fun c x m -> raise Nvector.OperationNotProvided)
        | Some f -> f

      let minquotient =
        match A.ops.minquotient with
        | None -> (fun num denom -> raise Nvector.OperationNotProvided)
        | Some f -> f

      let space =
        match A.ops.space with
        | None -> (fun x -> raise Nvector.OperationNotProvided)
        | Some f -> f

      let getlength = A.ops.getlength

      let print =
        match A.ops.print with
        | None -> (fun ?logfile _ -> raise Nvector.OperationNotProvided)
        | Some f -> (fun ?logfile x -> f x logfile)

      let linearcombination =
        match A.ops.linearcombination with
        | None -> (fun c x z -> raise Nvector.OperationNotProvided)
        | Some f -> f

      let scaleaddmulti =
        match A.ops.scaleaddmulti with
        | None -> (fun c x y z -> raise Nvector.OperationNotProvided)
        | Some f -> f

      let dotprodmulti =
        match A.ops.dotprodmulti with
        | None -> (fun x y d -> raise Nvector.OperationNotProvided)
        | Some f -> f

      let linearsumvectorarray =
        match A.ops.linearsumvectorarray with
        | None -> (fun a x b y z -> raise Nvector.OperationNotProvided)
        | Some f -> f

      let scalevectorarray =
        match A.ops.scalevectorarray with
        | None -> (fun c x z -> raise Nvector.OperationNotProvided)
        | Some f -> f

      let constvectorarray =
        match A.ops.constvectorarray with
        | None -> (fun c x -> raise Nvector.OperationNotProvided)
        | Some f -> f

      let wrmsnormvectorarray =
        match A.ops.wrmsnormvectorarray with
        | None -> (fun x w m -> raise Nvector.OperationNotProvided)
        | Some f -> f

      let wrmsnormmaskvectorarray =
        match A.ops.wrmsnormmaskvectorarray with
        | None -> (fun x w id m -> raise Nvector.OperationNotProvided)
        | Some f -> f

      let scaleaddmultivectorarray =
        match A.ops.scaleaddmultivectorarray with
        | None -> (fun a x yy zz -> raise Nvector.OperationNotProvided)
        | Some f -> f

      let linearcombinationvectorarray =
        match A.ops.linearcombinationvectorarray with
        | None -> (fun c xx z -> raise Nvector.OperationNotProvided)
        | Some f -> f

      module Local = struct
        let dotprod =
          match A.ops.dotprod_local with
          | None -> (fun _ _ -> raise Nvector.OperationNotProvided)
          | Some f -> f

        let maxnorm =
          match A.ops.maxnorm_local with
          | None -> (fun _ -> raise Nvector.OperationNotProvided)
          | Some f -> f

        let min =
          match A.ops.min_local with
          | None -> (fun _ -> raise Nvector.OperationNotProvided)
          | Some f -> f

        let l1norm =
          match A.ops.l1norm_local with
          | None -> (fun _ -> raise Nvector.OperationNotProvided)
          | Some f -> f

        let invtest =
          match A.ops.invtest_local with
          | None -> (fun _ _ -> raise Nvector.OperationNotProvided)
          | Some f -> f

        let constrmask =
          match A.ops.constrmask_local with
          | None -> (fun _ _ _ -> raise Nvector.OperationNotProvided)
          | Some f -> f

        let minquotient =
          match A.ops.minquotient_local with
          | None -> (fun _ _ -> raise Nvector.OperationNotProvided)
          | Some f -> f

        let wsqrsum =
          match A.ops.wsqrsum_local with
          | None -> (fun _ _ -> raise Nvector.OperationNotProvided)
          | Some f -> f

        let wsqrsummask =
          match A.ops.wsqrsummask_local with
          | None -> (fun _ _ _ -> raise Nvector.OperationNotProvided)
          | Some f -> f
      end
    end (* }}} *)
  end

module Any = struct (* {{{ *)

  let option_map f = function
    | None -> None
    | Some x -> Some (f x)

  let convert_ops (type d)
        ~(inject:d -> Nvector.gdata)
        ~(project:Nvector.gdata -> d) ops
    =
    let { (* {{{ *)
        check;
        clone;
        space;
        getlength;
        print;
        linearsum;
        const;
        prod;
        div;
        scale;
        abs;
        inv;
        addconst;
        maxnorm;
        wrmsnorm;
        min;

        dotprod;
        compare;
        invtest;

        wl2norm;
        l1norm;
        wrmsnormmask;
        constrmask;
        minquotient;

        getcommunicator;

        linearcombination;
        scaleaddmulti;
        dotprodmulti;

        linearsumvectorarray;
        scalevectorarray;
        constvectorarray;
        wrmsnormvectorarray;
        wrmsnormmaskvectorarray;
        scaleaddmultivectorarray;
        linearcombinationvectorarray;

        dotprod_local;
        maxnorm_local;
        min_local;
        l1norm_local;
        invtest_local;
        constrmask_local;
        minquotient_local;
        wsqrsum_local;
        wsqrsummask_local;
      } = ops (* }}} *)
    in
    let projecta = Array.map project in
    let projectaa = Array.map (Array.map project) in
    let single f = (fun v1 -> f (project v1)) in
    let double f = (fun v1 v2 -> f (project v1) (project v2)) in
    let triple f = (fun v1 v2 v3 -> f (project v1) (project v2) (project v3)) in
    let lifto f f' = match f with None -> None | Some f -> Some (f' f) in
    {
        (* {{{ *)
        check        = double check;
        clone        = (fun v -> inject (clone (project v)));
        space        = option_map single space;
        getlength    = single getlength;
        print        = (match print with None -> None
                       | Some f -> Some (fun v lf -> f (project v) lf));
        linearsum    = (fun c1 v1 c2 v2 vr -> linearsum c1 (project v1)
                                                              c2 (project v2)
                                                              (project vr));
        const        = (fun c v -> const c (project v));
        prod         = triple prod;
        div          = triple div;
        scale        = (fun c v1 vr -> scale c (project v1) (project vr));
        abs          = double abs;
        inv          = double inv;
        addconst     = (fun v1 c vr -> addconst (project v1) c (project vr));
        maxnorm      = single maxnorm;
        wrmsnorm     = double wrmsnorm;
        min          = single min;

        dotprod      = double dotprod;
        compare      = (fun c vx vr -> compare c (project vx) (project vr));
        invtest      = double invtest;

        wl2norm      = lifto wl2norm double;
        l1norm       = lifto l1norm single;
        wrmsnormmask = lifto wrmsnormmask triple;
        constrmask   = lifto constrmask triple;
        minquotient  = lifto minquotient double;

        getcommunicator = lifto getcommunicator single;

        linearcombination = lifto linearcombination
          (fun f a vxs vr -> f a (projecta vxs) (project vr));
        scaleaddmulti = lifto scaleaddmulti
          (fun f a vx vxs vrs -> f a (project vx) (projecta vxs) (projecta vrs));
        dotprodmulti = lifto dotprodmulti
          (fun f vx vxs r -> f (project vx) (projecta vxs) r);

        linearsumvectorarray = lifto linearsumvectorarray
          (fun f c1 vxs1 c2 vxs2 vrs -> f c1 (projecta vxs1)
                                          c2 (projecta vxs2)
                                             (projecta vrs));
        scalevectorarray = lifto scalevectorarray
          (fun f a vxs vrs -> f a (projecta vxs) (projecta vrs));
        constvectorarray = lifto constvectorarray
          (fun f c vrs -> f c (projecta vrs));
        wrmsnormvectorarray = lifto wrmsnormvectorarray
          (fun f vxs vys rs -> f (projecta vxs) (projecta vys) rs);
        wrmsnormmaskvectorarray = lifto wrmsnormmaskvectorarray
          (fun f vxs vys vz ra -> f (projecta vxs) (projecta vys) (project vz) ra);
        scaleaddmultivectorarray = lifto scaleaddmultivectorarray
          (fun f a vxs vyss vrss -> f a (projecta vxs) (projectaa vyss) (projectaa vrss));
        linearcombinationvectorarray = lifto linearcombinationvectorarray
          (fun f a vxss vrs -> f a (projectaa vxss) (projecta vrs));

        dotprod_local     = lifto dotprod_local double;
        maxnorm_local     = lifto maxnorm_local single;
        min_local         = lifto min_local single;
        l1norm_local      = lifto l1norm_local single;
        invtest_local     = lifto invtest_local double;
        constrmask_local  = lifto constrmask_local triple;
        minquotient_local = lifto minquotient_local double;
        wsqrsum_local     = lifto wsqrsum_local double;
        wsqrsummask_local = lifto wsqrsummask_local triple;
        (* }}} *)
    }

  (* Same underlying call as c_make_wrap but with types declared differently. *)
  external c_make_any_wrap
    : Nvector.gdata nvector_ops
      -> Nvector.gdata
      -> (Nvector.any -> bool)
      -> (Nvector.any -> Nvector.any)
      -> Nvector.any
    = "sunml_nvec_wrap_custom"

  let do_enable f nv v = if not (f nv v) then raise Nvector.OperationNotProvided

  let rec make_wrap_injected ops
      ?(with_fused_ops=false)
      ?(with_linear_combination=false)
      ?(with_scale_add_multi=false)
      ?(with_dot_prod_multi=false)
      ?(with_linear_sum_vector_array=false)
      ?(with_scale_vector_array=false)
      ?(with_const_vector_array=false)
      ?(with_wrms_norm_vector_array=false)
      ?(with_wrms_norm_mask_vector_array=false)
      ?(with_scale_add_multi_vector_array=false)
      ?(with_linear_combination_vector_array=false)
      check v
    =
      let nv = c_make_any_wrap ops v check (clone ops check) in
      do_enable c_enablefusedops_custom nv
                with_fused_ops;
      do_enable c_enablefusedops_custom nv
                with_fused_ops;
      do_enable c_enablelinearcombination_custom nv
                with_linear_combination;
      do_enable c_enablescaleaddmulti_custom nv
                with_scale_add_multi;
      do_enable c_enabledotprodmulti_custom nv
                with_dot_prod_multi;
      do_enable c_enablelinearsumvectorarray_custom nv
                with_linear_sum_vector_array;
      do_enable c_enablescalevectorarray_custom nv
                with_scale_vector_array;
      do_enable c_enableconstvectorarray_custom nv
                with_const_vector_array;
      do_enable c_enablewrmsnormvectorarray_custom nv
                with_wrms_norm_vector_array;
      do_enable c_enablewrmsnormmaskvectorarray_custom nv
                with_wrms_norm_mask_vector_array;
      do_enable c_enablescaleaddmultivectorarray_custom nv
                with_scale_add_multi_vector_array;
      do_enable c_enablelinearcombinationvectorarray_custom nv
                with_linear_combination_vector_array;
      nv

  and make_wrap ops ~inject
      ?(with_fused_ops=false)
      ?(with_linear_combination=false)
      ?(with_scale_add_multi=false)
      ?(with_dot_prod_multi=false)
      ?(with_linear_sum_vector_array=false)
      ?(with_scale_vector_array=false)
      ?(with_const_vector_array=false)
      ?(with_wrms_norm_vector_array=false)
      ?(with_wrms_norm_mask_vector_array=false)
      ?(with_scale_add_multi_vector_array=false)
      ?(with_linear_combination_vector_array=false)
      rv
    =
      if not Sundials_impl.Version.has_nvector_get_id
        then raise Sundials.Config.NotImplementedBySundialsVersion;
      let v = inject rv in
      let check nv' =
        ops.check v (uv nv') && Nvector.get_id nv' = Nvector.Custom
      in
      make_wrap_injected ops
        ~with_fused_ops
        ~with_linear_combination
        ~with_scale_add_multi
        ~with_dot_prod_multi
        ~with_linear_sum_vector_array
        ~with_scale_vector_array
        ~with_const_vector_array
        ~with_wrms_norm_vector_array
        ~with_wrms_norm_mask_vector_array
        ~with_scale_add_multi_vector_array
        ~with_linear_combination_vector_array
        check v

  and clone ops check nv =
    let nv' = make_wrap_injected ops check (ops.clone (uv nv)) in
    ignore (c_enablelinearcombination_custom nv'
              (Nvector.Ops.has_linearcombination nv));
    ignore (c_enablescaleaddmulti_custom nv'
              (Nvector.Ops.has_scaleaddmulti nv));
    ignore (c_enabledotprodmulti_custom nv'
              (Nvector.Ops.has_dotprodmulti nv));
    ignore (c_enablelinearsumvectorarray_custom nv'
              (Nvector.Ops.has_linearsumvectorarray nv));
    ignore (c_enablescalevectorarray_custom nv'
              (Nvector.Ops.has_scalevectorarray nv));
    ignore (c_enableconstvectorarray_custom nv'
              (Nvector.Ops.has_constvectorarray nv));
    ignore (c_enablewrmsnormvectorarray_custom nv'
              (Nvector.Ops.has_wrmsnormvectorarray nv));
    ignore (c_enablewrmsnormmaskvectorarray_custom nv'
              (Nvector.Ops.has_wrmsnormmaskvectorarray nv));
    ignore (c_enablescaleaddmultivectorarray_custom nv'
              (Nvector.Ops.has_scaleaddmultivectorarray nv));
    ignore (c_enablelinearcombinationvectorarray_custom nv'
              (Nvector.Ops.has_linearcombinationvectorarray nv));
    nv'

end (* }}} *)

