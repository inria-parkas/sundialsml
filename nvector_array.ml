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

module type ARRAY_NVECTOR =
  sig
    include Nvector.NVECTOR_OPS

    val array_nvec_ops  : t Nvector_custom.nvector_ops
    val make            : int -> float -> t Nvector_custom.t
    val wrap            : t -> t Nvector_custom.t
  end

module MakeOps =
  functor (A : sig
      type data
      val get       : data -> int -> float
      val set       : data -> int -> float -> unit
      val fill      : data -> float -> unit
      val make      : int -> float -> data
      val clone     : data -> data
      val length    : data -> int
      val fold_left : ('a -> float -> 'a) -> 'a -> data -> 'a
    end) ->
  struct
    type t = A.data

    let n_vclone = A.clone

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

    let n_vlinearsum a x b y z =
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

    let n_vconst c a = A.fill a c

    let n_vscale c x z =
      if c = 1.0 then
        lift_op (fun x -> x) x z
      else if c = -1.0 then
        lift_op (fun x -> -. x) x z
      else
        lift_op (fun x -> c *. x) x z

    let n_vaddconst x b = lift_op (fun x -> x +. b) x

    let n_vmaxnorm =
      let f max x =
        let ax = abs_float x in
        if ax > max then ax else max
      in A.fold_left f 0.0

    let n_vwrmsnorm x w =
      let f a x w = a +. ((x *. w) ** 2.0) in
      sqrt ((zip_fold_left f 0.0 x w) /. float (A.length x))

    let n_vwrmsnormmask x w id =
      let f a id x w =
        if id > 0.0 then a +. ((x *. w) ** 2.0) else a
      in
      sqrt ((triple_zip_fold_left f 0.0 id x w) /. float (A.length x))

    let n_vmin =
      let f min x = if x < min then x else min
      in A.fold_left f max_float

    let n_vdotprod =
      let f a x y = a +. x *. y in
      zip_fold_left f 0.0

    let n_vcompare c =
      lift_op (fun x -> if abs_float x >= c then 1.0 else 0.0)

    let n_vinvtest x z =
      let l = A.length x in
      let rec f r i =
        if i = l then r
        else if (A.get x i) = 0.0 then f false (i + 1)
        else (A.set z i (1.0 /. (A.get x i)); f r (i + 1))
      in f true 0

    let n_vwl2norm x w =
      let f a x w = a +. ((x *. w) ** 2.0) in
      sqrt (zip_fold_left f 0.0 x w)

    let n_vl1norm =
      let f a x = a +. abs_float x in
      A.fold_left f 0.0

    let n_vconstrmask c x m =
      let test = ref true in
      let check b = if b then 0.0 else (test := false; 1.0) in
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
      !test

    let n_vminquotient =
      let f m n d =
        if d = 0.0 then m
        else min m (n /. d)
      in
      zip_fold_left f Sundials.big_real

    let n_vprod = lift_bop ( *. )
    let n_vdiv  = lift_bop ( /. )
    let n_vabs  = lift_op abs_float
    let n_vinv  = lift_op (fun x -> 1.0 /. x)

    let array_nvec_ops = {
          Nvector_custom.n_vclone        = n_vclone;
          Nvector_custom.n_vdestroy      = None;
          Nvector_custom.n_vspace        = None;
          Nvector_custom.n_vlinearsum    = n_vlinearsum;
          Nvector_custom.n_vconst        = n_vconst;
          Nvector_custom.n_vprod         = n_vprod;
          Nvector_custom.n_vdiv          = n_vdiv;
          Nvector_custom.n_vscale        = n_vscale;
          Nvector_custom.n_vabs          = n_vabs;
          Nvector_custom.n_vinv          = n_vinv;
          Nvector_custom.n_vaddconst     = n_vaddconst;
          Nvector_custom.n_vmaxnorm      = n_vmaxnorm;
          Nvector_custom.n_vwrmsnorm     = n_vwrmsnorm;
          Nvector_custom.n_vmin          = n_vmin;
          Nvector_custom.n_vdotprod      = Some n_vdotprod;
          Nvector_custom.n_vcompare      = Some n_vcompare;
          Nvector_custom.n_vinvtest      = Some n_vinvtest;
          Nvector_custom.n_vwl2norm      = Some n_vwl2norm;
          Nvector_custom.n_vl1norm       = Some n_vl1norm;
          Nvector_custom.n_vwrmsnormmask = Some n_vwrmsnormmask;
          Nvector_custom.n_vconstrmask   = Some n_vconstrmask;
          Nvector_custom.n_vminquotient  = Some n_vminquotient;
    }

    let make n e =
      Nvector_custom.make array_nvec_ops (A.make n e)

    let wrap a =
      Nvector_custom.make array_nvec_ops a
      (* (Nvector.Mutable.add_tracing "::" array_nvec_ops) *)
  end

module Array = MakeOps (
  struct
    type data = float array
    include Array
    let clone = Array.copy
    let fill a c = fill a 0 (length a) c
  end)

include Array

