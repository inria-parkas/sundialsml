(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a BSD 2-Clause License, refer to the file LICENSE.           *)
(*                                                                     *)
(***********************************************************************)

type kind
type 'a t = ('a, kind) Nvector.t

type 'a nvector_ops = {
  n_vclone           : 'a -> 'a;
  n_vdestroy         : ('a -> unit) option;
  n_vspace           : ('a -> (int * int)) option;
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
}

external make_wrap : 'a nvector_ops -> 'a -> 'a t
    = "ml_nvec_wrap_custom"

let add_tracing msg ops =
  let pr s = print_string msg; print_endline s in
  let {
      n_vclone           = n_vclone;
      n_vdestroy         = n_vdestroy;
      n_vspace           = n_vspace;
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
    } = ops
  in
  let fo f f' = match f with None -> None | Some f -> Some (f' f) in

  let tr_nvclone a = pr "nvclone"; n_vclone a
  and tr_nvdestroy = fo n_vdestroy (fun f -> fun a -> (pr "nvdestroy"; f a))
  and tr_nvspace = fo n_vspace (fun f -> fun a -> (pr "nvspace"; f a))
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
  in
  {
      n_vclone           = tr_nvclone;
      n_vdestroy         = tr_nvdestroy;
      n_vspace           = tr_nvspace;
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
   }

exception OperationNotSupported

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

    module Ops = struct
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

      let n_vwl2norm = match A.ops.n_vwl2norm with
                       | None -> raise OperationNotSupported
                       | Some f -> (fun x w -> f (uv x) (uv w))

      let n_vl1norm = match A.ops.n_vl1norm with
                      | None -> raise OperationNotSupported
                      | Some f -> (fun x -> f (uv x))

      let n_vwrmsnormmask = match A.ops.n_vwrmsnormmask with
                            | None -> raise OperationNotSupported
                            | Some f -> (fun x w id -> f (uv x) (uv w) (uv id))

      let n_vconstrmask = match A.ops.n_vconstrmask with
                          | None -> raise OperationNotSupported
                          | Some f -> (fun c x m -> f (uv c) (uv x) (uv m))

      let n_vminquotient = match A.ops.n_vminquotient with
                           | None -> raise OperationNotSupported
                           | Some f -> (fun num denom -> f (uv num) (uv denom))
    end

    module DataOps = struct
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

      let n_vwl2norm      = match A.ops.n_vwl2norm with
                            | None -> raise OperationNotSupported
                            | Some f -> f
      let n_vl1norm       = match A.ops.n_vl1norm with
                            | None -> raise OperationNotSupported
                            | Some f -> f
      let n_vwrmsnormmask = match A.ops.n_vwrmsnormmask with
                            | None -> raise OperationNotSupported
                            | Some f -> f
      let n_vconstrmask   = match A.ops.n_vconstrmask with
                            | None -> raise OperationNotSupported
                            | Some f -> f
      let n_vminquotient  = match A.ops.n_vminquotient with
                            | None -> raise OperationNotSupported
                            | Some f -> f
    end
  end

