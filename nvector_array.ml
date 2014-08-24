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
  end

module NvectorFn =
  functor (A : sig
      type data
      val get       : data -> int -> float
      val set       : data -> int -> float -> unit
      val fill      : data -> float -> unit
      val make      : int -> float -> data
      val copy      : data -> data
      val length    : data -> int
      val fold_left : ('a -> float -> 'a) -> 'a -> data -> 'a
    end) ->
  struct
    type t = A.data

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

    let arr_nvlinearsum a x b y z =
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

    let arr_nvconst c a = A.fill a c

    let arr_nvscale c x z =
      if c = 1.0 then
        lift_op (fun x -> x) x z
      else if c = -1.0 then
        lift_op (fun x -> -. x) x z
      else
        lift_op (fun x -> c *. x) x z

    let arr_nvaddconst x b = lift_op (fun x -> x +. b) x

    let arr_nvmaxnorm =
      let f max x =
        let ax = abs_float x in
        if ax > max then ax else max
      in A.fold_left f 0.0

    let arr_nvwrmsnorm x w =
      let f a x w = a +. ((x *. w) ** 2.0) in
      sqrt ((zip_fold_left f 0.0 x w) /. float (A.length x))

    let arr_nvwrmsnormmask x w id =
      let f a id x w =
        if id > 0.0 then a +. ((x *. w) ** 2.0) else a
      in
      sqrt ((triple_zip_fold_left f 0.0 id x w) /. float (A.length x))

    let arr_nvmin =
      let f min x = if x < min then x else min
      in A.fold_left f max_float

    let arr_nvdotprod =
      let f a x y = a +. x *. y in
      zip_fold_left f 0.0

    let arr_nvcompare c =
      lift_op (fun x -> if abs_float x >= c then 1.0 else 0.0)

    let arr_nvinvtest x z =
      let l = A.length x in
      let rec f i =
        if i = l then true
        else if (A.get x i) = 0.0 then false
        else (A.set z i (1.0 /. (A.get x i)); f (i + 1))
      in f 0

    let arr_nvwl2norm x w =
      let f a x w = a +. ((x *. w) ** 2.0) in
      sqrt (zip_fold_left f 0.0 x w)

    let arr_nvl1norm =
      let f a x = a +. abs_float x in
      A.fold_left f 0.0

    let arr_nvconstrmask c x m =
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

    let arr_nvminquotient =
      let f m n d =
        if d = 0.0 then m
        else min m (n /. d)
      in
      zip_fold_left f Sundials.big_real

    let array_nvec_ops = {
          Nvector_custom.nvclone        = A.copy;
          Nvector_custom.nvdestroy      = None;
          Nvector_custom.nvspace        = None;
          Nvector_custom.nvlinearsum    = arr_nvlinearsum;
          Nvector_custom.nvconst        = arr_nvconst;
          Nvector_custom.nvprod         = lift_bop ( *. );
          Nvector_custom.nvdiv          = lift_bop ( /. );
          Nvector_custom.nvscale        = arr_nvscale;
          Nvector_custom.nvabs          = lift_op abs_float;
          Nvector_custom.nvinv          = lift_op (fun x -> 1.0 /. x);
          Nvector_custom.nvaddconst     = arr_nvaddconst;
          Nvector_custom.nvmaxnorm      = arr_nvmaxnorm;
          Nvector_custom.nvwrmsnorm     = arr_nvwrmsnorm;
          Nvector_custom.nvmin          = arr_nvmin;
          Nvector_custom.nvdotprod      = Some arr_nvdotprod;
          Nvector_custom.nvcompare      = Some arr_nvcompare;
          Nvector_custom.nvinvtest      = Some arr_nvinvtest;
          Nvector_custom.nvwl2norm      = Some arr_nvwl2norm;
          Nvector_custom.nvl1norm       = Some arr_nvl1norm;
          Nvector_custom.nvwrmsnormmask = Some arr_nvwrmsnormmask;
          Nvector_custom.nvconstrmask   = Some arr_nvconstrmask;
          Nvector_custom.nvminquotient  = Some arr_nvminquotient;
    }

    let make n e =
      Nvector_custom.make array_nvec_ops (A.make n e)

    let wrap a =
      Nvector_custom.make array_nvec_ops a
      (* (Nvector.Mutable.add_tracing "::" array_nvec_ops) *)
  end

module Array = NvectorFn (
  struct
    type data = float array
    include Array
    let fill a c = fill a 0 (length a) c
  end)

module Bigarray = NvectorFn (
  struct
    include Bigarray.Array1

    type data = (float, Bigarray.float64_elt, Bigarray.c_layout) t

    let make n d =
      let arr = create Bigarray.float64 Bigarray.c_layout n in
      fill arr d;
      arr

    let length = dim

    let copy src =
      let dst = create Bigarray.float64 Bigarray.c_layout (length src) in
      blit src dst;
      dst

    let fold_left f a vs =
      let len = length vs - 1 in
      let a = ref a in
      for i = 0 to len do
        a := f !a vs.{i}
      done; !a
  end)

include Array

