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
type 'a t = ('a, kind) Sundials.nvector

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

  n_vdotprod         : ('a -> 'a -> float) option;
  n_vcompare         : (float -> 'a -> 'a -> unit) option;
  n_vinvtest         : ('a -> 'a -> bool) option;

  n_vwl2norm         : ('a -> 'a -> float) option;
  n_vl1norm          : ('a -> float) option;
  n_vwrmsnormmask    : ('a -> 'a -> 'a -> float) option;
  n_vconstrmask      : ('a -> 'a -> 'a -> bool) option;
  n_vminquotient     : ('a -> 'a -> float) option;
}

external make : 'a nvector_ops -> 'a -> 'a t
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
  and tr_nvdotprod = fo n_vdotprod (fun f -> fun x y -> pr "nvdotprod"; f x y)
  and tr_nvcompare =
    fo n_vcompare (fun f -> fun c x z -> pr "nvcompare"; f c x z)
  and tr_nvinvtest = fo n_vinvtest (fun f -> fun x z -> pr "nvinvtest"; f x z)
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

