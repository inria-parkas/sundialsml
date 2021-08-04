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

(* Must match with nvector_ml.h:nvector_ops_tag *)
type 'a nvector_ops = { (* {{{ *)
  n_vcheck           : 'a -> 'a -> bool;
  n_vclone           : 'a -> 'a;
  n_vspace           : ('a -> (int * int)) option;
  n_vgetlength       : 'a -> int;
  n_vlinearsum       : float -> 'a -> float -> 'a -> 'a -> unit;
  n_vconst           : float -> 'a -> unit;
  n_vprod            : 'a -> 'a -> 'a -> unit;
  n_vdiv             : 'a -> 'a -> 'a -> unit;
  n_vscale           : float -> 'a -> 'a -> unit;
  n_vabs             : 'a -> 'a -> unit;
  n_vinv             : 'a -> 'a -> unit;
  n_vaddconst        : 'a -> float -> 'a -> unit;
  n_vmaxnorm         : 'a -> float;
  n_vwrmsnorm        : 'a -> 'a -> float;
  n_vmin             : 'a -> float;

  n_vdotprod         : 'a -> 'a -> float;
  n_vcompare         : float -> 'a -> 'a -> unit;
  n_vinvtest         : 'a -> 'a -> bool;

  n_vwl2norm         : ('a -> 'a -> float) option;
  n_vl1norm          : ('a -> float) option;
  n_vwrmsnormmask    : ('a -> 'a -> 'a -> float) option;
  n_vconstrmask      : ('a -> 'a -> 'a -> bool) option;
  n_vminquotient     : ('a -> 'a -> float) option;

  n_vlinearcombination :
    (Sundials.RealArray.t -> 'a array -> 'a -> unit) option;
  n_vscaleaddmulti :
    (Sundials.RealArray.t -> 'a -> 'a array -> 'a array -> unit) option;
  n_vdotprodmulti : ('a -> 'a array -> Sundials.RealArray.t -> unit) option;

  n_vlinearsumvectorarray :
    (float -> 'a array -> float -> 'a array -> 'a array -> unit) option;
  n_vscalevectorarray :
    (Sundials.RealArray.t -> 'a array -> 'a array -> unit) option;
  n_vconstvectorarray : (float -> 'a array -> unit) option;
  n_vwrmsnormvectorarray :
    ('a array -> 'a array -> Sundials.RealArray.t -> unit) option;
  n_vwrmsnormmaskvectorarray :
    ('a array -> 'a array -> 'a -> Sundials.RealArray.t -> unit) option;
  n_vscaleaddmultivectorarray :
    (Sundials.RealArray.t -> 'a array -> 'a array array -> 'a array array ->
      unit) option;
  n_vlinearcombinationvectorarray :
    (Sundials.RealArray.t -> 'a array array -> 'a array -> unit) option;

  n_vdotprod_local      : ('a -> 'a -> float) option;
  n_vmaxnorm_local      : ('a -> float) option;
  n_vmin_local          : ('a -> float) option;
  n_vl1norm_local       : ('a -> float) option;
  n_vinvtest_local      : ('a -> 'a -> bool) option;
  n_vconstrmask_local   : ('a -> 'a -> 'a -> bool) option;
  n_vminquotient_local  : ('a -> 'a -> float) option;
  n_vwsqrsum_local      : ('a -> 'a -> float) option;
  n_vwsqrsummask_local  : ('a -> 'a -> 'a -> float) option;

} (* }}} *)

exception OperationNotSupported

(* Selectively enable and disable fused and array operations *)
external c_enablefusedops_custom                       : 'd t -> bool -> bool
  = "sunml_nvec_custom_enablefusedops"
external c_enablelinearcombination_custom              : 'd t -> bool -> bool
  = "sunml_nvec_custom_enablelinearcombination"
external c_enablescaleaddmulti_custom                  : 'd t -> bool -> bool
  = "sunml_nvec_custom_enablescaleaddmulti"
external c_enabledotprodmulti_custom                   : 'd t -> bool -> bool
  = "sunml_nvec_custom_enabledotprodmulti"
external c_enablelinearsumvectorarray_custom           : 'd t -> bool -> bool
  = "sunml_nvec_custom_enablelinearsumvectorarray"
external c_enablescalevectorarray_custom               : 'd t -> bool -> bool
  = "sunml_nvec_custom_enablescalevectorarray"
external c_enableconstvectorarray_custom               : 'd t -> bool -> bool
  = "sunml_nvec_custom_enableconstvectorarray"
external c_enablewrmsnormvectorarray_custom            : 'd t -> bool -> bool
  = "sunml_nvec_custom_enablewrmsnormvectorarray"
external c_enablewrmsnormmaskvectorarray_custom        : 'd t -> bool -> bool
  = "sunml_nvec_custom_enablewrmsnormmaskvectorarray"
external c_enablescaleaddmultivectorarray_custom       : 'd t -> bool -> bool
  = "sunml_nvec_custom_enablescaleaddmultivectorarray"
external c_enablelinearcombinationvectorarray_custom   : 'd t -> bool -> bool
  = "sunml_nvec_custom_enablelinearcombinationvectorarray"

external c_make_wrap : 'a nvector_ops -> 'a -> ('a -> bool) -> 'a t
    = "sunml_nvec_wrap_custom"

let do_enable f nv v =
  match v with
  | None -> ()
  | Some v -> if not (f nv v) then raise OperationNotSupported

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

let make_wrap ops ?(with_fused_ops=false) v =
  let nv = c_make_wrap ops v (ops.n_vcheck v) in
  if with_fused_ops && not (c_enablefusedops_custom nv true)
    then raise OperationNotSupported;
  nv

let add_tracing msg ops =
  let pr s = print_string msg; print_endline s in
  let { (* {{{ *)
      n_vcheck           = n_vcheck;
      n_vclone           = n_vclone;
      n_vspace           = n_vspace;
      n_vgetlength       = n_vgetlength;
      n_vlinearsum       = n_vlinearsum;
      n_vconst           = n_vconst;
      n_vprod            = n_vprod;
      n_vdiv             = n_vdiv;
      n_vscale           = n_vscale;
      n_vabs             = n_vabs;
      n_vinv             = n_vinv;
      n_vaddconst        = n_vaddconst;
      n_vmaxnorm         = n_vmaxnorm;
      n_vwrmsnorm        = n_vwrmsnorm;
      n_vmin             = n_vmin;

      n_vdotprod         = n_vdotprod;
      n_vcompare         = n_vcompare;
      n_vinvtest         = n_vinvtest;

      n_vwl2norm         = n_vwl2norm;
      n_vl1norm          = n_vl1norm;
      n_vwrmsnormmask    = n_vwrmsnormmask;
      n_vconstrmask      = n_vconstrmask;
      n_vminquotient     = n_vminquotient;

      n_vlinearcombination            = n_vlinearcombination;
      n_vscaleaddmulti                = n_vscaleaddmulti;
      n_vdotprodmulti                 = n_vdotprodmulti;

      n_vlinearsumvectorarray         = n_vlinearsumvectorarray;
      n_vscalevectorarray             = n_vscalevectorarray;
      n_vconstvectorarray             = n_vconstvectorarray;
      n_vwrmsnormvectorarray          = n_vwrmsnormvectorarray;
      n_vwrmsnormmaskvectorarray      = n_vwrmsnormmaskvectorarray;
      n_vscaleaddmultivectorarray     = n_vscaleaddmultivectorarray;
      n_vlinearcombinationvectorarray = n_vlinearcombinationvectorarray;

      n_vdotprod_local     = n_vdotprod_local;
      n_vmaxnorm_local     = n_vmaxnorm_local;
      n_vmin_local         = n_vmin_local;
      n_vl1norm_local      = n_vl1norm_local;
      n_vinvtest_local     = n_vinvtest_local;
      n_vconstrmask_local  = n_vconstrmask_local;
      n_vminquotient_local = n_vminquotient_local;
      n_vwsqrsum_local     = n_vwsqrsum_local;
      n_vwsqrsummask_local = n_vwsqrsummask_local;
    } = ops (* }}} *)
  in
  let fo f f' = match f with None -> None | Some f -> Some (f' f) in

  let tr_nvcheck x =
    pr "nvcheck-create";
    let check = n_vcheck x in
    function y ->
      pr "nvcheck-check";
      check y
  and tr_nvclone a = pr "nvclone"; n_vclone a
  (* ... {{{ *)
  and tr_nvspace = fo n_vspace (fun f -> fun a -> (pr "nvspace"; f a))
  and tr_nvgetlength a = pr "nvgetlength"; n_vgetlength a
  and tr_nvlinearsum a x b y z = pr "nvlinearsum"; n_vlinearsum a x b y z
  and tr_nvconst c z = pr "nvconst"; n_vconst c z
  and tr_nvprod x y z = pr "nvprod"; n_vprod x y z
  and tr_nvdiv x y z = pr "nvdiv"; n_vdiv x y z
  and tr_nvscale c x z = pr "nvscale"; n_vscale c x z
  and tr_nvabs x z = pr "nvabs"; n_vabs x z
  and tr_nvinv x z = pr "nvinv"; n_vinv x z
  and tr_nvaddconst x b z = pr "nvaddconst"; n_vaddconst x b z
  and tr_nvmaxnorm x = pr "nvmaxnorm"; n_vmaxnorm x
  and tr_nvwrmsnorm x w = pr "nvwrmsnorm"; n_vwrmsnorm x w
  and tr_nvmin x = pr "nvmin"; n_vmin x
  and tr_nvdotprod x y = pr "nvdotprod"; n_vdotprod x y
  and tr_nvcompare c x z = pr "nvcompare"; n_vcompare c x z
  and tr_nvinvtest x z = pr "nvinvtest"; n_vinvtest x z

  and tr_nvwl2norm = fo n_vwl2norm (fun f -> fun x w -> pr "nvwl2norm"; f x w)
  and tr_nvl1norm = fo n_vl1norm (fun f -> fun x -> pr "nvl1norm"; f x)
  and tr_nvwrmsnormmask =
    fo n_vwrmsnormmask (fun f -> fun x w id -> pr "nvwrmsnormmask"; f x w id)
  and tr_nvconstrmask =
    fo n_vconstrmask (fun f -> fun c x m -> pr "nvconstrmask"; f c x m)
  and tr_nvminquotient =
    fo n_vminquotient (fun f -> fun n d -> pr "nvminquotient"; f n d)

  and tr_nvlinearcombination =
    fo n_vlinearcombination
    (fun f -> fun c x z -> pr "nvlinearcombination"; f c x z)
  and tr_nvscaleaddmulti =
    fo n_vscaleaddmulti
    (fun f -> fun c x y z -> pr "nvscaleaddmulti"; f c x y z)
  and tr_nvdotprodmulti =
    fo n_vdotprodmulti
    (fun f -> fun x y d -> pr "nvdotprodmulti"; f x y d)

  and tr_nvlinearsumvectorarray =
    fo n_vlinearsumvectorarray
    (fun f -> fun a x b y z -> pr "nvlinearsumvectorarray"; f a x b y z)
  and tr_nvscalevectorarray =
    fo n_vscalevectorarray
    (fun f -> fun c x z -> pr "nvscalevectorarray"; f c x z)
  and tr_nvconstvectorarray =
    fo n_vconstvectorarray
    (fun f -> fun c x -> pr "nvconstvectorarray"; f c x)
  and tr_nvwrmsnormvectorarray =
    fo n_vwrmsnormvectorarray
    (fun f -> fun x w m -> pr "nvwrmsnormvectorarray"; f x w m)
  and tr_nvwrmsnormmaskvectorarray =
    fo n_vwrmsnormmaskvectorarray
    (fun f -> fun x w m -> pr "nvwrmsnormmaskvectorarray"; f x w m)
  and tr_nvscaleaddmultivectorarray =
    fo n_vscaleaddmultivectorarray
    (fun f -> fun a x yy zz -> pr "nvscaleaddmultivectorarray"; f a x yy zz)
  and tr_nvlinearcombinationvectorarray =
    fo n_vlinearcombinationvectorarray
    (fun f -> fun c xx z -> pr "nvlinearcombinationvectorarray"; f c xx z)

  and tr_nvdotprod_local     =
    fo n_vdotprod_local (fun f -> fun x t -> pr "nvdotprod_local"; f x t)
  and tr_nvmaxnorm_local     =
    fo n_vmaxnorm_local (fun f -> fun x -> pr "nvmaxnorm_local"; f x)
  and tr_nvmin_local         =
    fo n_vmin_local (fun f -> fun x -> pr "nvmin_local"; f x)
  and tr_nvl1norm_local      =
    fo n_vl1norm_local (fun f -> fun x -> pr "nvl1norm_local"; f x)
  and tr_nvinvtest_local     =
    fo n_vinvtest_local (fun f -> fun x z -> pr "nvinvtest_local"; f x z)
  and tr_nvconstrmask_local  =
    fo n_vconstrmask_local (fun f -> fun c x m -> pr "nvconstrmask_local"; f c x m)
  and tr_nvminquotient_local =
    fo n_vminquotient_local (fun f -> fun n d -> pr "nvminquotient_local"; f n d)
  and tr_nvwsqrsum_local     =
    fo n_vwsqrsum_local (fun f -> fun x w -> pr "nvwsqrsum_local"; f x w)
  and tr_nvwsqrsummask_local =
    fo n_vwsqrsummask_local (fun f -> fun x w id -> pr "nvwsqrsummask_local"; f x w id)
  (* }}} *)
  in
  {
      (* {{{ *)
      n_vcheck           = tr_nvcheck;
      n_vclone           = tr_nvclone;
      n_vspace           = tr_nvspace;
      n_vgetlength       = tr_nvgetlength;
      n_vlinearsum       = tr_nvlinearsum;
      n_vconst           = tr_nvconst;
      n_vprod            = tr_nvprod;
      n_vdiv             = tr_nvdiv;
      n_vscale           = tr_nvscale;
      n_vabs             = tr_nvabs;
      n_vinv             = tr_nvinv;
      n_vaddconst        = tr_nvaddconst;
      n_vmaxnorm         = tr_nvmaxnorm;
      n_vwrmsnorm        = tr_nvwrmsnorm;
      n_vmin             = tr_nvmin;

      n_vdotprod         = tr_nvdotprod;
      n_vcompare         = tr_nvcompare;
      n_vinvtest         = tr_nvinvtest;

      n_vwl2norm         = tr_nvwl2norm;
      n_vl1norm          = tr_nvl1norm;
      n_vwrmsnormmask    = tr_nvwrmsnormmask;
      n_vconstrmask      = tr_nvconstrmask;
      n_vminquotient     = tr_nvminquotient;

      n_vlinearcombination            = tr_nvlinearcombination;
      n_vscaleaddmulti                = tr_nvscaleaddmulti;
      n_vdotprodmulti                 = tr_nvdotprodmulti;

      n_vlinearsumvectorarray         = tr_nvlinearsumvectorarray;
      n_vscalevectorarray             = tr_nvscalevectorarray;
      n_vconstvectorarray             = tr_nvconstvectorarray;
      n_vwrmsnormvectorarray          = tr_nvwrmsnormvectorarray;
      n_vwrmsnormmaskvectorarray      = tr_nvwrmsnormmaskvectorarray;
      n_vscaleaddmultivectorarray     = tr_nvscaleaddmultivectorarray;
      n_vlinearcombinationvectorarray = tr_nvlinearcombinationvectorarray;

      n_vdotprod_local     = tr_nvdotprod_local;
      n_vmaxnorm_local     = tr_nvmaxnorm_local;
      n_vmin_local         = tr_nvmin_local;
      n_vl1norm_local      = tr_nvl1norm_local;
      n_vinvtest_local     = tr_nvinvtest_local;
      n_vconstrmask_local  = tr_nvconstrmask_local;
      n_vminquotient_local = tr_nvminquotient_local;
      n_vwsqrsum_local     = tr_nvwsqrsum_local;
      n_vwsqrsummask_local = tr_nvwsqrsummask_local;
      (* }}} *)
  }

let uv = Nvector.unwrap

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

      let n_vclone n = wrap (A.ops.n_vclone (uv n))
      let n_vlinearsum a x b y z = A.ops.n_vlinearsum a (uv x) b (uv y) (uv z)
      let n_vconst c z = A.ops.n_vconst c (uv z)
      let n_vprod x y z = A.ops.n_vprod (uv x) (uv y) (uv z)
      let n_vdiv x y z = A.ops.n_vdiv (uv x) (uv y) (uv z)
      let n_vscale c x z = A.ops.n_vscale c (uv x) (uv z)
      let n_vabs x z = A.ops.n_vabs (uv x) (uv z)
      let n_vinv x z = A.ops.n_vinv (uv x) (uv z)
      let n_vaddconst x b z = A.ops.n_vaddconst (uv x) b (uv z)
      let n_vdotprod x y = A.ops.n_vdotprod (uv x) (uv y)
      let n_vmaxnorm x = A.ops.n_vmaxnorm (uv x)
      let n_vwrmsnorm x w = A.ops.n_vwrmsnorm (uv x) (uv w)
      let n_vmin x = A.ops.n_vmin (uv x)
      let n_vcompare c x z = A.ops.n_vcompare c (uv x) (uv z)
      let n_vinvtest x z = A.ops.n_vinvtest (uv x) (uv z)

      let n_vwl2norm =
        match A.ops.n_vwl2norm with
        | None -> (fun x w -> raise OperationNotSupported)
        | Some f -> (fun x w -> f (uv x) (uv w))

      let n_vl1norm =
        match A.ops.n_vl1norm with
        | None -> (fun x -> raise OperationNotSupported)
        | Some f -> (fun x -> f (uv x))

      let n_vwrmsnormmask =
        match A.ops.n_vwrmsnormmask with
        | None -> (fun x w id -> raise OperationNotSupported)
        | Some f -> (fun x w id -> f (uv x) (uv w) (uv id))

      let n_vconstrmask =
        match A.ops.n_vconstrmask with
        | None -> (fun c x m -> raise OperationNotSupported)
        | Some f -> (fun c x m -> f (uv c) (uv x) (uv m))

      let n_vminquotient =
        match A.ops.n_vminquotient with
        | None -> (fun num denom -> raise OperationNotSupported)
        | Some f -> (fun num denom -> f (uv num) (uv denom))

      let n_vspace =
        match A.ops.n_vspace with
        | None -> (fun x -> raise OperationNotSupported)
        | Some f -> (fun x -> f (uv x))

      let n_vgetlength a = A.ops.n_vgetlength (uv a)

      let n_vlinearcombination =
        match A.ops.n_vlinearcombination with
        | None -> (fun c x z -> raise OperationNotSupported)
        | Some f -> (fun c x z -> f c (Array.map uv x) (uv z))

      let n_vscaleaddmulti =
        match A.ops.n_vscaleaddmulti with
        | None -> (fun c x y z -> raise OperationNotSupported)
        | Some f ->
            (fun c x y z -> f c (uv x) (Array.map uv y) (Array.map uv z))

      let n_vdotprodmulti =
        match A.ops.n_vdotprodmulti with
        | None -> (fun x y d -> raise OperationNotSupported)
        | Some f -> (fun x y d -> f (uv x) (Array.map uv y) d)

      let n_vlinearsumvectorarray =
        match A.ops.n_vlinearsumvectorarray with
        | None -> (fun a x b y z -> raise OperationNotSupported)
        | Some f -> (fun a x b y z -> f a (Array.map uv x) b (Array.map uv y)
                                          (Array.map uv z))

      let n_vscalevectorarray =
        match A.ops.n_vscalevectorarray with
        | None -> (fun c x z -> raise OperationNotSupported)
        | Some f -> (fun c x z -> f c (Array.map uv x) (Array.map uv z))

      let n_vconstvectorarray =
        match A.ops.n_vconstvectorarray with
        | None -> (fun c x -> raise OperationNotSupported)
        | Some f -> (fun c x -> f c (Array.map uv x))

      let n_vwrmsnormvectorarray =
        match A.ops.n_vwrmsnormvectorarray with
        | None -> (fun x w m -> raise OperationNotSupported)
        | Some f -> (fun x w m -> f (Array.map uv x) (Array.map uv w) m)

      let n_vwrmsnormmaskvectorarray =
        match A.ops.n_vwrmsnormmaskvectorarray with
        | None -> (fun x w id m -> raise OperationNotSupported)
        | Some f -> (fun x w id m -> f (Array.map uv x) (Array.map uv w)
                                       (uv id) m)

      let n_vscaleaddmultivectorarray =
        match A.ops.n_vscaleaddmultivectorarray with
        | None -> (fun a x yy zz -> raise OperationNotSupported)
        | Some f -> (fun a x yy zz -> f a (Array.map uv x)
                                        (Array.map (Array.map uv) yy)
                                        (Array.map (Array.map uv) zz))

      let n_vlinearcombinationvectorarray =
        match A.ops.n_vlinearcombinationvectorarray with
        | None -> (fun c xx z -> raise OperationNotSupported)
        | Some f -> (fun c xx z -> f c (Array.map (Array.map uv) xx)
                                       (Array.map uv z))

      module Local = struct
        let n_vdotprod =
          match A.ops.n_vdotprod_local with
          | None -> (fun _ _ -> raise OperationNotSupported)
          | Some f -> (fun x t -> f (uv x) (uv t))

        let n_vmaxnorm =
          match A.ops.n_vmaxnorm_local with
          | None -> (fun _ -> raise OperationNotSupported)
          | Some f -> (fun x -> f (uv x))

        let n_vmin =
          match A.ops.n_vmin_local with
          | None -> (fun _ -> raise OperationNotSupported)
          | Some f -> (fun x -> f (uv x))

        let n_vl1norm =
          match A.ops.n_vl1norm_local with
          | None -> (fun _ -> raise OperationNotSupported)
          | Some f -> (fun x -> f (uv x))

        let n_vinvtest =
          match A.ops.n_vinvtest_local with
          | None -> (fun _ _ -> raise OperationNotSupported)
          | Some f -> (fun x z -> f (uv x) (uv z))

        let n_vconstrmask =
          match A.ops.n_vconstrmask_local with
          | None -> (fun _ _ _ -> raise OperationNotSupported)
          | Some f -> (fun c x m -> f (uv c) (uv x) (uv m))

        let n_vminquotient =
          match A.ops.n_vminquotient_local with
          | None -> (fun _ _ -> raise OperationNotSupported)
          | Some f -> (fun n d -> f (uv n) (uv d))

        let n_vwsqrsum =
          match A.ops.n_vwsqrsum_local with
          | None -> (fun _ _ -> raise OperationNotSupported)
          | Some f -> (fun x w -> f (uv x) (uv w))

        let n_vwsqrsummask =
          match A.ops.n_vwsqrsummask_local with
          | None -> (fun _ _ _ -> raise OperationNotSupported)
          | Some f -> (fun x w id -> f (uv x) (uv w) (uv id))
      end
    end (* }}} *)

    module DataOps = struct (* {{{ *)
      type t = data

      let n_vclone        = A.ops.n_vclone
      let n_vlinearsum    = A.ops.n_vlinearsum
      let n_vconst        = A.ops.n_vconst
      let n_vprod         = A.ops.n_vprod
      let n_vdiv          = A.ops.n_vdiv
      let n_vscale        = A.ops.n_vscale
      let n_vabs          = A.ops.n_vabs
      let n_vinv          = A.ops.n_vinv
      let n_vaddconst     = A.ops.n_vaddconst
      let n_vdotprod      = A.ops.n_vdotprod
      let n_vmaxnorm      = A.ops.n_vmaxnorm
      let n_vwrmsnorm     = A.ops.n_vwrmsnorm
      let n_vmin          = A.ops.n_vmin
      let n_vcompare      = A.ops.n_vcompare
      let n_vinvtest      = A.ops.n_vinvtest

      let n_vwl2norm =
        match A.ops.n_vwl2norm with
        | None -> (fun x w -> raise OperationNotSupported)
        | Some f -> f

      let n_vl1norm =
        match A.ops.n_vl1norm with
        | None -> (fun x -> raise OperationNotSupported)
        | Some f -> f

      let n_vwrmsnormmask =
        match A.ops.n_vwrmsnormmask with
        | None -> (fun x w id -> raise OperationNotSupported)
        | Some f -> f

      let n_vconstrmask =
        match A.ops.n_vconstrmask with
        | None -> (fun c x m -> raise OperationNotSupported)
        | Some f -> f

      let n_vminquotient =
        match A.ops.n_vminquotient with
        | None -> (fun num denom -> raise OperationNotSupported)
        | Some f -> f

      let n_vspace =
        match A.ops.n_vspace with
        | None -> (fun x -> raise OperationNotSupported)
        | Some f -> f

      let n_vgetlength = A.ops.n_vgetlength

      let n_vlinearcombination =
        match A.ops.n_vlinearcombination with
        | None -> (fun c x z -> raise OperationNotSupported)
        | Some f -> f

      let n_vscaleaddmulti =
        match A.ops.n_vscaleaddmulti with
        | None -> (fun c x y z -> raise OperationNotSupported)
        | Some f -> f

      let n_vdotprodmulti =
        match A.ops.n_vdotprodmulti with
        | None -> (fun x y d -> raise OperationNotSupported)
        | Some f -> f

      let n_vlinearsumvectorarray =
        match A.ops.n_vlinearsumvectorarray with
        | None -> (fun a x b y z -> raise OperationNotSupported)
        | Some f -> f

      let n_vscalevectorarray =
        match A.ops.n_vscalevectorarray with
        | None -> (fun c x z -> raise OperationNotSupported)
        | Some f -> f

      let n_vconstvectorarray =
        match A.ops.n_vconstvectorarray with
        | None -> (fun c x -> raise OperationNotSupported)
        | Some f -> f

      let n_vwrmsnormvectorarray =
        match A.ops.n_vwrmsnormvectorarray with
        | None -> (fun x w m -> raise OperationNotSupported)
        | Some f -> f

      let n_vwrmsnormmaskvectorarray =
        match A.ops.n_vwrmsnormmaskvectorarray with
        | None -> (fun x w id m -> raise OperationNotSupported)
        | Some f -> f

      let n_vscaleaddmultivectorarray =
        match A.ops.n_vscaleaddmultivectorarray with
        | None -> (fun a x yy zz -> raise OperationNotSupported)
        | Some f -> f

      let n_vlinearcombinationvectorarray =
        match A.ops.n_vlinearcombinationvectorarray with
        | None -> (fun c xx z -> raise OperationNotSupported)
        | Some f -> f

      module Local = struct
        let n_vdotprod =
          match A.ops.n_vdotprod_local with
          | None -> (fun _ _ -> raise OperationNotSupported)
          | Some f -> f

        let n_vmaxnorm =
          match A.ops.n_vmaxnorm_local with
          | None -> (fun _ -> raise OperationNotSupported)
          | Some f -> f

        let n_vmin =
          match A.ops.n_vmin_local with
          | None -> (fun _ -> raise OperationNotSupported)
          | Some f -> f

        let n_vl1norm =
          match A.ops.n_vl1norm_local with
          | None -> (fun _ -> raise OperationNotSupported)
          | Some f -> f

        let n_vinvtest =
          match A.ops.n_vinvtest_local with
          | None -> (fun _ _ -> raise OperationNotSupported)
          | Some f -> f

        let n_vconstrmask =
          match A.ops.n_vconstrmask_local with
          | None -> (fun _ _ _ -> raise OperationNotSupported)
          | Some f -> f

        let n_vminquotient =
          match A.ops.n_vminquotient_local with
          | None -> (fun _ _ -> raise OperationNotSupported)
          | Some f -> f

        let n_vwsqrsum =
          match A.ops.n_vwsqrsum_local with
          | None -> (fun _ _ -> raise OperationNotSupported)
          | Some f -> f

        let n_vwsqrsummask =
          match A.ops.n_vwsqrsummask_local with
          | None -> (fun _ _ _ -> raise OperationNotSupported)
          | Some f -> f
      end
    end (* }}} *)
  end

