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
    type t

    val array_nvec_ops  : t Nvector_custom.nvector_ops
    val make            : int -> float -> t Nvector_custom.t
    val wrap            : t -> t Nvector_custom.t

    module Ops : Nvector.NVECTOR_OPS with type t = t Nvector_custom.t
    module DataOps : Nvector.NVECTOR_OPS with type t = t
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

    module DataOps = struct
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
    end

    let array_nvec_ops = {
          Nvector_custom.n_vclone        = DataOps.n_vclone;
          Nvector_custom.n_vdestroy      = None;
          Nvector_custom.n_vspace        = None;
          Nvector_custom.n_vlinearsum    = DataOps.n_vlinearsum;
          Nvector_custom.n_vconst        = DataOps.n_vconst;
          Nvector_custom.n_vprod         = DataOps.n_vprod;
          Nvector_custom.n_vdiv          = DataOps.n_vdiv;
          Nvector_custom.n_vscale        = DataOps.n_vscale;
          Nvector_custom.n_vabs          = DataOps.n_vabs;
          Nvector_custom.n_vinv          = DataOps.n_vinv;
          Nvector_custom.n_vaddconst     = DataOps.n_vaddconst;
          Nvector_custom.n_vmaxnorm      = DataOps.n_vmaxnorm;
          Nvector_custom.n_vwrmsnorm     = DataOps.n_vwrmsnorm;
          Nvector_custom.n_vmin          = DataOps.n_vmin;
          Nvector_custom.n_vdotprod      = Some DataOps.n_vdotprod;
          Nvector_custom.n_vcompare      = Some DataOps.n_vcompare;
          Nvector_custom.n_vinvtest      = Some DataOps.n_vinvtest;
          Nvector_custom.n_vwl2norm      = Some DataOps.n_vwl2norm;
          Nvector_custom.n_vl1norm       = Some DataOps.n_vl1norm;
          Nvector_custom.n_vwrmsnormmask = Some DataOps.n_vwrmsnormmask;
          Nvector_custom.n_vconstrmask   = Some DataOps.n_vconstrmask;
          Nvector_custom.n_vminquotient  = Some DataOps.n_vminquotient;
    }

    let make n e =
      Nvector_custom.make array_nvec_ops (A.make n e)

    let wrap a =
      Nvector_custom.make array_nvec_ops a
      (* (Nvector.Mutable.add_tracing "::" array_nvec_ops) *)

    module Ops = struct
      type t = A.data Nvector_custom.t

      let unvec = Sundials.unvec

      let n_vclone x =
        let xd = unvec x in
        wrap (DataOps.n_vclone xd)

      let n_vlinearsum a x b y z
            = DataOps.n_vlinearsum a (unvec x) b (unvec y) (unvec z)
      let n_vconst c a = DataOps.n_vconst c (unvec a)
      let n_vprod x y z = DataOps.n_vprod (unvec x) (unvec y) (unvec z)
      let n_vdiv x y z = DataOps.n_vdiv (unvec x) (unvec y) (unvec z)
      let n_vscale c x z = DataOps.n_vscale c (unvec x) (unvec z)
      let n_vabs x z = DataOps.n_vabs (unvec x) (unvec z)
      let n_vinv x z = DataOps.n_vinv (unvec x) (unvec z)
      let n_vaddconst x b z = DataOps.n_vaddconst (unvec x) b (unvec z)
      let n_vdotprod x y = DataOps.n_vdotprod (unvec x) (unvec y)
      let n_vmaxnorm x = DataOps.n_vmaxnorm (unvec x)
      let n_vwrmsnorm x w = DataOps.n_vwrmsnorm (unvec x) (unvec w)
      let n_vwrmsnormmask x w id
            = DataOps.n_vwrmsnormmask (unvec x) (unvec w) (unvec id)
      let n_vmin x = DataOps.n_vmin (unvec x)
      let n_vwl2norm x w = DataOps.n_vwl2norm (unvec x) (unvec w)
      let n_vl1norm x = DataOps.n_vl1norm (unvec x)
      let n_vcompare c x z = DataOps.n_vcompare c (unvec x) (unvec z)
      let n_vinvtest x z = DataOps.n_vinvtest (unvec x) (unvec z)
      let n_vconstrmask c x m
            = DataOps.n_vconstrmask (unvec c) (unvec x) (unvec m)
      let n_vminquotient n d = DataOps.n_vminquotient (unvec n) (unvec d)
    end
  end

module Array = MakeOps (
  struct
    type data = float array
    include Array
    let clone = Array.copy
    let fill a c = fill a 0 (length a) c
  end)

include Array

