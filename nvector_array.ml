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
    type data
    type kind = Nvector_custom.kind

    val array_nvec_ops  : data Nvector_custom.nvector_ops
    val make            : int -> float -> data Nvector_custom.t
    val wrap            : data -> data Nvector_custom.t
    val unwrap          : data Nvector_custom.t -> data

    module Ops : Nvector.NVECTOR_OPS with type t = data Nvector_custom.t
    module DataOps : Nvector.NVECTOR_OPS with type t = data
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
    end) ->
  struct
    type data = A.data
    type kind = Nvector_custom.kind

    let checkfn v1 =
      let l = A.length v1 in
      (fun v2 -> l = A.length v2)

    module DataOps = struct
      type t = A.data

      let n_vclone = A.clone

      let arr_vaxpy a x y =
        if a = 1.0 then
          for i = 0 to A.length x - 1 do
            A.set y i (A.get y i +. A.get x i)
          done
        else if a = -1.0 then
          for i = 0 to A.length x - 1 do
            A.set y i (A.get y i -. A.get x i)
          done
        else
          for i = 0 to A.length x - 1 do
            A.set y i (A.get y i +. a *. A.get x i)
          done

      let n_vlinearsum a x b y z =
        if b = 1.0 && z == y then
          arr_vaxpy a x y
        else if a = 1.0 && z == x then
          arr_vaxpy b y x
        else if a = 1.0 && b = 1.0 then
          for i = 0 to A.length x - 1 do
            A.set z i (A.get x i +. A.get y i)
          done
        else if (a = 1.0 && b = -1.0) || (a = -1.0 && b == 1.0) then
          let v1, v2 = if (a = 1.0 && b = -1.0) then y, x else x, y in
          for i = 0 to A.length v1 - 1 do
            A.set z i (A.get v1 i -. A.get v2 i)
          done
        else if a = 1.0 || b = 1.0 then
          let c, v1, v2 = if a = 1.0 then b, y, x else a, x, y in
          for i = 0 to A.length v1 - 1 do
            A.set z i (c *. A.get v1 i +. A.get v2 i)
          done
        else if a = -1.0 || b = -1.0 then
          let c, v1, v2 = if a = -1.0 then b, y, x else a, x, y in
          for i = 0 to A.length v1 - 1 do
            A.set z i (c *. A.get v1 i -. A.get v2 i)
          done
        else if a = b then
          for i = 0 to A.length x - 1 do
            A.set z i (a *. (A.get x i +. A.get y i))
          done
        else if a = -.b then
          for i = 0 to A.length x - 1 do
            A.set z i (a *. (A.get x i -. A.get y i))
          done
        else
          for i = 0 to A.length x - 1 do
            A.set z i (a *. A.get x i +. b *. A.get y i)
          done

      let n_vconst c a = A.fill a c

      let n_vscale c x z =
        if c = 1.0 then
          for i = 0 to A.length x - 1 do
            A.set z i (A.get x i)
          done
        else if c = -1.0 then
          for i = 0 to A.length x - 1 do
            A.set z i (-. A.get x i)
          done
        else
          for i = 0 to A.length x - 1 do
            A.set z i (c *. A.get x i)
          done

      let n_vaddconst x b z =
        for i = 0 to A.length x - 1 do
          A.set z i (A.get x i +. b)
        done

      let n_vmaxnorm x =
        let max = ref 0.0 in
        for i = 0 to A.length x - 1 do
          let ax = abs_float (A.get x i) in
          if ax > !max then max := ax
        done;
        !max

      let n_vwrmsnorm x w =
        let a = ref 0.0 in
        let lx = A.length x in
        for i = 0 to lx - 1 do
          a := !a +. ((A.get x i) *. (A.get w i) *. (A.get x i) *. (A.get w i))
        done;
        sqrt (!a /. float lx)

      let n_vwrmsnormmask x w id =
        let a = ref 0.0 in
        let lx = A.length x in
        for i = 0 to lx - 1 do
          if A.get id i > 0.0 then
            a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
        done;
        sqrt (!a /. float lx)

      let n_vmin x =
        let min = ref max_float in
        for i = 0 to A.length x - 1 do
          let xv = A.get x i in
          if xv < !min then min := xv
        done;
        !min

      let n_vdotprod x y =
        let a = ref 0.0 in
        for i = 0 to A.length x - 1 do
          a := !a +. (A.get x i *. A.get y i)
        done;
        !a

      let n_vcompare c x z =
        for i = 0 to A.length x - 1 do
          A.set z i (if abs_float (A.get x i) >= c then 1.0 else 0.0)
        done

      let n_vinvtest x z =
        let l = A.length x in
        let rec f r i =
          if i = l then r
          else if (A.get x i) = 0.0 then f false (i + 1)
          else (A.set z i (1.0 /. (A.get x i)); f r (i + 1))
        in f true 0

      let n_vwl2norm x w =
        let a = ref 0.0 in
        for i = 0 to A.length x - 1 do
          a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
        done;
        sqrt !a

      let n_vl1norm x =
        let a = ref 0.0 in
        for i = 0 to A.length x - 1 do
          a := !a +. abs_float (A.get x i)
        done;
        !a

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
        for i = 0 to A.length c - 1 do
          A.set m i (f (A.get c i) (A.get x i))
        done;
        !test

      let n_vminquotient n d =
        let m = ref Sundials.big_real in
        for i = 0 to A.length n - 1 do
          if (A.get d i) <> 0.0 then
            m := min !m (A.get n i /. A.get d i)
        done;
        !m

      let n_vprod x y z =
        for i = 0 to A.length x - 1 do
          A.set z i (A.get x i *. A.get y i)
        done

      let n_vdiv x y z =
        for i = 0 to A.length x - 1 do
          A.set z i (A.get x i /. A.get y i)
        done

      let n_vabs x z =
        for i = 0 to A.length x - 1 do
          A.set z i (abs_float (A.get x i))
        done

      let n_vinv x z =
        for i = 0 to A.length x - 1 do
          A.set z i (1.0 /. (A.get x i))
        done

    end

    let array_nvec_ops = {
          Nvector_custom.n_vcheck        = checkfn;
          Nvector_custom.n_vclone        = DataOps.n_vclone;
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
          Nvector_custom.n_vdotprod      = DataOps.n_vdotprod;
          Nvector_custom.n_vcompare      = DataOps.n_vcompare;
          Nvector_custom.n_vinvtest      = DataOps.n_vinvtest;

          Nvector_custom.n_vwl2norm      = Some DataOps.n_vwl2norm;
          Nvector_custom.n_vl1norm       = Some DataOps.n_vl1norm;
          Nvector_custom.n_vwrmsnormmask = Some DataOps.n_vwrmsnormmask;
          Nvector_custom.n_vconstrmask   = Some DataOps.n_vconstrmask;
          Nvector_custom.n_vminquotient  = Some DataOps.n_vminquotient;
    }

    let make n e = Nvector_custom.make_wrap array_nvec_ops (A.make n e)

    let wrap = Nvector_custom.make_wrap array_nvec_ops
      (* (Nvector.Mutable.add_tracing "::" array_nvec_ops) *)

    let unwrap = Nvector.unwrap

    module Ops = struct
      type t = A.data Nvector_custom.t

      let n_vclone x =
        let xd = unwrap x in
        wrap (DataOps.n_vclone xd)

      let n_vlinearsum a x b y z
            = DataOps.n_vlinearsum a (unwrap x) b (unwrap y) (unwrap z)
      let n_vconst c a = DataOps.n_vconst c (unwrap a)
      let n_vprod x y z = DataOps.n_vprod (unwrap x) (unwrap y) (unwrap z)
      let n_vdiv x y z = DataOps.n_vdiv (unwrap x) (unwrap y) (unwrap z)
      let n_vscale c x z = DataOps.n_vscale c (unwrap x) (unwrap z)
      let n_vabs x z = DataOps.n_vabs (unwrap x) (unwrap z)
      let n_vinv x z = DataOps.n_vinv (unwrap x) (unwrap z)
      let n_vaddconst x b z = DataOps.n_vaddconst (unwrap x) b (unwrap z)
      let n_vdotprod x y = DataOps.n_vdotprod (unwrap x) (unwrap y)
      let n_vmaxnorm x = DataOps.n_vmaxnorm (unwrap x)
      let n_vwrmsnorm x w = DataOps.n_vwrmsnorm (unwrap x) (unwrap w)
      let n_vwrmsnormmask x w id
            = DataOps.n_vwrmsnormmask (unwrap x) (unwrap w) (unwrap id)
      let n_vmin x = DataOps.n_vmin (unwrap x)
      let n_vwl2norm x w = DataOps.n_vwl2norm (unwrap x) (unwrap w)
      let n_vl1norm x = DataOps.n_vl1norm (unwrap x)
      let n_vcompare c x z = DataOps.n_vcompare c (unwrap x) (unwrap z)
      let n_vinvtest x z = DataOps.n_vinvtest (unwrap x) (unwrap z)
      let n_vconstrmask c x m
            = DataOps.n_vconstrmask (unwrap c) (unwrap x) (unwrap m)
      let n_vminquotient n d = DataOps.n_vminquotient (unwrap n) (unwrap d)
    end
  end

(* (* Too slow! *)
module SlowerArray = MakeOps (
  struct
    type data = float array
    include Array
    let clone = Array.copy
    let fill a c = fill a 0 (length a) c
  end)
*)

module Array =
  struct
    type data = float array
    type kind = Nvector_custom.kind

    let checkfn v1 =
      let l = Array.length v1 in
      (fun v2 -> l = Array.length v2)

    module DataOps = struct
      type t = float array
      module A = Array

      let fill a c = Array.fill a 0 (Array.length a) c

      let n_vclone = Array.copy

      let arr_vaxpy a x y =
        if a = 1.0 then
          for i = 0 to A.length x - 1 do
            A.set y i (A.get y i +. A.get x i)
          done
        else if a = -1.0 then
          for i = 0 to A.length x - 1 do
            A.set y i (A.get y i -. A.get x i)
          done
        else
          for i = 0 to A.length x - 1 do
            A.set y i (A.get y i +. a *. A.get x i)
          done

      let n_vlinearsum a x b y z =
        if b = 1.0 && z == y then
          arr_vaxpy a x y
        else if a = 1.0 && z == x then
          arr_vaxpy b y x
        else if a = 1.0 && b = 1.0 then
          for i = 0 to A.length x - 1 do
            A.set z i (A.get x i +. A.get y i)
          done
        else if (a = 1.0 && b = -1.0) || (a = -1.0 && b == 1.0) then
          let v1, v2 = if (a = 1.0 && b = -1.0) then y, x else x, y in
          for i = 0 to A.length v1 - 1 do
            A.set z i (A.get v1 i -. A.get v2 i)
          done
        else if a = 1.0 || b = 1.0 then
          let c, v1, v2 = if a = 1.0 then b, y, x else a, x, y in
          for i = 0 to A.length v1 - 1 do
            A.set z i (c *. A.get v1 i +. A.get v2 i)
          done
        else if a = -1.0 || b = -1.0 then
          let c, v1, v2 = if a = -1.0 then b, y, x else a, x, y in
          for i = 0 to A.length v1 - 1 do
            A.set z i (c *. A.get v1 i -. A.get v2 i)
          done
        else if a = b then
          for i = 0 to A.length x - 1 do
            A.set z i (a *. (A.get x i +. A.get y i))
          done
        else if a = -.b then
          for i = 0 to A.length x - 1 do
            A.set z i (a *. (A.get x i -. A.get y i))
          done
        else
          for i = 0 to A.length x - 1 do
            A.set z i (a *. A.get x i +. b *. A.get y i)
          done

      let n_vconst c a = fill a c

      let n_vscale c x z =
        if c = 1.0 then
          for i = 0 to A.length x - 1 do
            A.set z i (A.get x i)
          done
        else if c = -1.0 then
          for i = 0 to A.length x - 1 do
            A.set z i (-. A.get x i)
          done
        else
          for i = 0 to A.length x - 1 do
            A.set z i (c *. A.get x i)
          done

      let n_vaddconst x b z =
        for i = 0 to A.length x - 1 do
          A.set z i (A.get x i +. b)
        done

      let n_vmaxnorm x =
        let max = ref 0.0 in
        for i = 0 to A.length x - 1 do
          let ax = abs_float (A.get x i) in
          if ax > !max then max := ax
        done;
        !max

      let n_vwrmsnorm x w =
        let a = ref 0.0 in
        let lx = A.length x in
        for i = 0 to lx - 1 do
          a := !a +. ((A.get x i) *. (A.get w i) *. (A.get x i) *. (A.get w i))
        done;
        sqrt (!a /. float lx)

      let n_vwrmsnormmask x w id =
        let a = ref 0.0 in
        let lx = A.length x in
        for i = 0 to lx - 1 do
          if A.get id i > 0.0 then
            a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
        done;
        sqrt (!a /. float lx)

      let n_vmin x =
        let min = ref max_float in
        for i = 0 to A.length x - 1 do
          let xv = A.get x i in
          if xv < !min then min := xv
        done;
        !min

      let n_vdotprod x y =
        let a = ref 0.0 in
        for i = 0 to A.length x - 1 do
          a := !a +. (A.get x i *. A.get y i)
        done;
        !a

      let n_vcompare c x z =
        for i = 0 to A.length x - 1 do
          A.set z i (if abs_float (A.get x i) >= c then 1.0 else 0.0)
        done

      let n_vinvtest x z =
        let r = ref true in
        for i = 0 to A.length x do
          if A.get x i = 0.0 then r := false
          else A.set z i (1.0 /. (A.get x i))
        done;
        !r

      let n_vwl2norm x w =
        let a = ref 0.0 in
        for i = 0 to A.length x - 1 do
          a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
        done;
        sqrt !a

      let n_vl1norm x =
        let a = ref 0.0 in
        for i = 0 to A.length x - 1 do
          a := !a +. abs_float (A.get x i)
        done;
        !a

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
        for i = 0 to A.length c - 1 do
          A.set m i (f (A.get c i) (A.get x i))
        done;
        !test

      let n_vminquotient n d =
        let m = ref Sundials.big_real in
        for i = 0 to A.length n - 1 do
          if (A.get d i) <> 0.0 then
            m := min !m (A.get n i /. A.get d i)
        done;
        !m

      let n_vprod x y z =
        for i = 0 to A.length x - 1 do
          A.set z i (A.get x i *. A.get y i)
        done

      let n_vdiv x y z =
        for i = 0 to A.length x - 1 do
          A.set z i (A.get x i /. A.get y i)
        done

      let n_vabs x z =
        for i = 0 to A.length x - 1 do
          A.set z i (abs_float (A.get x i))
        done

      let n_vinv x z =
        for i = 0 to A.length x - 1 do
          A.set z i (1.0 /. (A.get x i))
        done
    end

    let array_nvec_ops = {
          Nvector_custom.n_vcheck        = checkfn;
          Nvector_custom.n_vclone        = DataOps.n_vclone;
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
          Nvector_custom.n_vdotprod      = DataOps.n_vdotprod;
          Nvector_custom.n_vcompare      = DataOps.n_vcompare;
          Nvector_custom.n_vinvtest      = DataOps.n_vinvtest;

          Nvector_custom.n_vwl2norm      = Some DataOps.n_vwl2norm;
          Nvector_custom.n_vl1norm       = Some DataOps.n_vl1norm;
          Nvector_custom.n_vwrmsnormmask = Some DataOps.n_vwrmsnormmask;
          Nvector_custom.n_vconstrmask   = Some DataOps.n_vconstrmask;
          Nvector_custom.n_vminquotient  = Some DataOps.n_vminquotient;
    }

    let make n e =
      Nvector_custom.make_wrap array_nvec_ops (Array.make n e)

    let wrap a =
      Nvector_custom.make_wrap array_nvec_ops a
      (* (Nvector.Mutable.add_tracing "::" array_nvec_ops) *)

    let unwrap = Nvector.unwrap

    module Ops = struct
      type t = float array Nvector_custom.t

      let n_vclone x =
        let xd = unwrap x in
        wrap (DataOps.n_vclone xd)

      let n_vlinearsum a x b y z
            = DataOps.n_vlinearsum a (unwrap x) b (unwrap y) (unwrap z)
      let n_vconst c a = DataOps.n_vconst c (unwrap a)
      let n_vprod x y z = DataOps.n_vprod (unwrap x) (unwrap y) (unwrap z)
      let n_vdiv x y z = DataOps.n_vdiv (unwrap x) (unwrap y) (unwrap z)
      let n_vscale c x z = DataOps.n_vscale c (unwrap x) (unwrap z)
      let n_vabs x z = DataOps.n_vabs (unwrap x) (unwrap z)
      let n_vinv x z = DataOps.n_vinv (unwrap x) (unwrap z)
      let n_vaddconst x b z = DataOps.n_vaddconst (unwrap x) b (unwrap z)
      let n_vdotprod x y = DataOps.n_vdotprod (unwrap x) (unwrap y)
      let n_vmaxnorm x = DataOps.n_vmaxnorm (unwrap x)
      let n_vwrmsnorm x w = DataOps.n_vwrmsnorm (unwrap x) (unwrap w)
      let n_vwrmsnormmask x w id
            = DataOps.n_vwrmsnormmask (unwrap x) (unwrap w) (unwrap id)
      let n_vmin x = DataOps.n_vmin (unwrap x)
      let n_vwl2norm x w = DataOps.n_vwl2norm (unwrap x) (unwrap w)
      let n_vl1norm x = DataOps.n_vl1norm (unwrap x)
      let n_vcompare c x z = DataOps.n_vcompare c (unwrap x) (unwrap z)
      let n_vinvtest x z = DataOps.n_vinvtest (unwrap x) (unwrap z)
      let n_vconstrmask c x m
            = DataOps.n_vconstrmask (unwrap c) (unwrap x) (unwrap m)
      let n_vminquotient n d = DataOps.n_vminquotient (unwrap n) (unwrap d)
    end
  end

include Array

