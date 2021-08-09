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

open Sundials

module type ARRAY_NVECTOR =
  sig (* {{{ *)
    type data
    type kind = Nvector_custom.kind
    type t = data Nvector_custom.t

    val array_nvec_ops  : data Nvector_custom.nvector_ops
    val make            : int -> float -> t
    val wrap            : ?with_fused_ops:bool -> data -> t
    val unwrap          : t -> data

    val enable :
         ?with_fused_ops                       : bool
      -> ?with_linear_combination              : bool
      -> ?with_scale_add_multi                 : bool
      -> ?with_dot_prod_multi                  : bool
      -> ?with_linear_sum_vector_array         : bool
      -> ?with_scale_vector_array              : bool
      -> ?with_const_vector_array              : bool
      -> ?with_wrms_norm_vector_array          : bool
      -> ?with_wrms_norm_mask_vector_array     : bool
      -> ?with_scale_add_multi_vector_array    : bool
      -> ?with_linear_combination_vector_array : bool
      -> t
      -> unit

    module Any : sig

      type Nvector.gdata += Arr of data

      val wrap :
           ?with_fused_ops                       : bool
        -> ?with_linear_combination              : bool
        -> ?with_scale_add_multi                 : bool
        -> ?with_dot_prod_multi                  : bool
        -> ?with_linear_sum_vector_array         : bool
        -> ?with_scale_vector_array              : bool
        -> ?with_const_vector_array              : bool
        -> ?with_wrms_norm_vector_array          : bool
        -> ?with_wrms_norm_mask_vector_array     : bool
        -> ?with_scale_add_multi_vector_array    : bool
        -> ?with_linear_combination_vector_array : bool
        -> data
        -> Nvector.any
    end

    module Ops : Nvector.NVECTOR_OPS with type t = t
    module DataOps : Nvector.NVECTOR_OPS with type t = data
  end (* }}} *)

module type ArrayOps = sig
  type data
  val get       : data -> int -> float
  val set       : data -> int -> float -> unit
  val fill      : data -> float -> unit
  val make      : int -> float -> data
  val clone     : data -> data
  val length    : data -> int
end

module MakeOps = functor
  (A : sig
         include ArrayOps

         type t
         val project : t -> data
         val inject  : data -> t
       end) ->
  struct (* {{{ *)
    type t = A.t

    let n_vclone wx = A.(inject (clone (project wx)))

    let uw = A.project

    let arr_vaxpy a wx wy =
      let x, y = uw wx, uw wy in
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

    let n_vlinearsum a wx b wy wz =
      if b = 1.0 && wz == wy then
        arr_vaxpy a wx wy
      else if a = 1.0 && wz == wx then
        arr_vaxpy b wy wx
      else
        let x, y, z = uw wx, uw wy, uw wz in
        if a = 1.0 && b = 1.0 then
          for i = 0 to A.length x - 1 do
            A.set z i (A.get x i +. A.get y i)
          done
        else if (a = 1.0 && b = -1.0) || (a = -1.0 && b = 1.0) then
          let v1, v2 = if (a = 1.0 && b = -1.0) then y, x else x, y in
          for i = 0 to A.length v1 - 1 do
            A.set z i (A.get v2 i -. A.get v1 i)
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

    let n_vconst c wa = A.fill (uw wa) c

    let n_vscale c wx wz =
      let x, z = uw wx, uw wz in
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

    let n_vaddconst wx b wz =
      let x, z = uw wx, uw wz in
      for i = 0 to A.length x - 1 do
        A.set z i (A.get x i +. b)
      done

    let n_vmaxnorm wx =
      let x = uw wx in
      let max = ref 0.0 in
      for i = 0 to A.length x - 1 do
        let ax = abs_float (A.get x i) in
        if ax > !max then max := ax
      done;
      !max

    let n_vwrmsnorm wx ww =
      let x, w = uw wx, uw ww in
      let a = ref 0.0 in
      let lx = A.length x in
      for i = 0 to lx - 1 do
        a := !a +. ((A.get x i) *. (A.get w i) *. (A.get x i) *. (A.get w i))
      done;
      sqrt (!a /. float lx)

    let n_vwrmsnormmask wx ww wid =
      let x, w, id = uw wx, uw ww, uw wid in
      let a = ref 0.0 in
      let lx = A.length x in
      for i = 0 to lx - 1 do
        if A.get id i > 0.0 then
          a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
      done;
      sqrt (!a /. float lx)

    let n_vmin wx =
      let x = uw wx in
      let min = ref max_float in
      for i = 0 to A.length x - 1 do
        let xv = A.get x i in
        if xv < !min then min := xv
      done;
      !min

    let n_vdotprod wx wy =
      let x, y = uw wx, uw wy in
      let a = ref 0.0 in
      for i = 0 to A.length x - 1 do
        a := !a +. (A.get x i *. A.get y i)
      done;
      !a

    let n_vcompare c wx wz =
      let x, z = uw wx, uw wz in
      for i = 0 to A.length x - 1 do
        A.set z i (if abs_float (A.get x i) >= c then 1.0 else 0.0)
      done

    let n_vinvtest wx wz =
      let x, z = uw wx, uw wz in
      let l = A.length x in
      let rec f r i =
        if i = l then r
        else if A.get x i = 0.0 then f false (i + 1)
        else (A.set z i (1.0 /. (A.get x i)); f r (i + 1))
      in f true 0

    let n_vwl2norm wx ww =
      let x, w = uw wx, uw ww in
      let a = ref 0.0 in
      for i = 0 to A.length x - 1 do
        a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
      done;
      sqrt !a

    let n_vl1norm wx =
      let x = uw wx in
      let a = ref 0.0 in
      for i = 0 to A.length x - 1 do
        a := !a +. abs_float (A.get x i)
      done;
      !a

    let n_vconstrmask wc wx wm =
      let c, x, m = uw wc, uw wx, uw wm in
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

    let n_vminquotient wn wd =
      let n, d = uw wn, uw wd in
      let m = ref Config.big_real in
      for i = 0 to A.length n - 1 do
        if (A.get d i) <> 0.0 then
          m := min !m (A.get n i /. A.get d i)
      done;
      !m

    let n_vprod wx wy wz =
      let x, y, z = uw wx, uw wy, uw wz in
      for i = 0 to A.length x - 1 do
        A.set z i (A.get x i *. A.get y i)
      done

    let n_vdiv wx wy wz =
      let x, y, z = uw wx, uw wy, uw wz in
      for i = 0 to A.length x - 1 do
        A.set z i (A.get x i /. A.get y i)
      done

    let n_vabs wx wz =
      let x, z = uw wx, uw wz in
      for i = 0 to A.length x - 1 do
        A.set z i (abs_float (A.get x i))
      done

    let n_vinv wx wz =
      let x, z = uw wx, uw wz in
      for i = 0 to A.length x - 1 do
        A.set z i (1.0 /. (A.get x i))
      done

    let n_vspace wx = (A.length (uw wx), 1)

    let n_vgetlength wx = A.length (uw wx)

    (* fused and array operations *)

    let n_vlinearcombination (ca : RealArray.t) (wxa : t array) (wz : t) =
      let nvec = Array.length wxa in
      if nvec = 1 then n_vscale ca.{0} wxa.(0) wz
      else if nvec = 2 then n_vlinearsum ca.{0} wxa.(0) ca.{1} wxa.(1) wz
      else
        let z = uw wz in
        let n = A.length z in
        if wxa.(0) == wz then begin
          let c0 = ca.{0} in
          if c0 <> 1.0 then
            for j = 0 to n - 1 do
              A.set z j (A.get z j *. c0)
            done;
          for i = 1 to nvec - 1 do
            let ci, x = ca.{i}, uw wxa.(i) in
            for j = 0 to n - 1 do
              A.set z j (A.get z j +. ci *. A.get x j)
            done
          done
        end
        else begin
          let c0, x = ca.{0}, uw wxa.(0) in
          for j = 0 to n - 1 do
            A.set z j (c0 *. A.get x j)
          done;
          for i = 1 to nvec - 1 do
            let ci, x = ca.{i}, uw wxa.(i) in
            for j = 0 to n - 1 do
              A.set z j (A.get z j +. ci *. A.get x j)
            done
          done
        end

    let n_vscaleaddmulti (aa : RealArray.t) (wx : t) (wya : t array)
                                                     (wza : t array) =
      let nvec = Array.length wya in
      if nvec = 1 then n_vlinearsum aa.{0} wx 1.0 wya.(0) wza.(0)
      else
        let x = uw wx in
        let n = A.length x in
        if wya == wza then
          for i = 0 to nvec - 1 do
            let a, y = aa.{i}, uw wya.(i) in
            for j = 0 to n - 1 do
              A.set y j (A.get y j +. a *. A.get x j)
            done
          done
        else
          for i = 0 to nvec - 1 do
            let ai, y, z = aa.{i}, uw wya.(i), uw wza.(i) in
            for j = 0 to n - 1 do
              A.set z j (ai *. A.get x j +. A.get y j)
            done
          done

    let n_vdotprodmulti (wx : t) (wya : t array) (dp : RealArray.t) =
      let nvec = Array.length wya in
      if nvec = 1 then dp.{0} <- n_vdotprod wx wya.(0)
      else
        let x = uw wx in
        let n = A.length x in
        for i = 0 to nvec - 1 do
          let y = uw wya.(i) in
          dp.{i} <- 0.0;
          for j = 0 to n - 1 do
            dp.{i} <- dp.{i} +. A.get x j *. A.get y j
          done
        done

    let arr_vaxpy_array a (wxa : t array) (wya : t array) =
      let nvec = Array.length wxa in
      let n = A.length (uw wxa.(0)) in
      if a = 1.0 then
        for i = 0 to nvec - 1 do
          let x, y = uw wxa.(i), uw wya.(i) in
          for j = 0 to n - 1 do
            A.set y j (A.get y j +. A.get x j)
          done
        done
      else if a = -1.0 then
        for i = 0 to nvec - 1 do
          let x, y = uw wxa.(i), uw wya.(i) in
          for j = 0 to n - 1 do
            A.set y j (A.get y j -. A.get x j)
          done
        done
      else
        for i = 0 to nvec - 1 do
          let x, y = uw wxa.(i), uw wya.(i) in
          for j = 0 to n - 1 do
            A.set y j (A.get y j +. a *. A.get x j)
          done
        done

    let v_sumvectorarray (wxa : t array) (wya : t array)
                                         (wza : t array) =
      let nvec = Array.length wxa in
      let n = A.length (uw wxa.(0)) in
      for i = 0 to nvec - 1 do
        let x, y, z = uw wxa.(i), uw wya.(i), uw wza.(i) in
        for j = 0 to n - 1 do
          A.set z j (A.get x j +. A.get y j)
        done
      done

    let v_diffvectorarray (wxa : t array) (wya : t array)
                                          (wza : t array) =
      let nvec = Array.length wxa in
      let n = A.length (uw wxa.(0)) in
      for i = 0 to nvec - 1 do
        let x, y, z = uw wxa.(i), uw wya.(i), uw wza.(i) in
        for j = 0 to n - 1 do
          A.set z j (A.get x j -. A.get y j)
        done
      done

    let v_lin1vectorarray a (wxa : t array) (wya : t array)
                                            (wza : t array) =
      let nvec = Array.length wxa in
      let n = A.length (uw wxa.(0)) in
      for i = 0 to nvec - 1 do
        let x, y, z = uw wxa.(i), uw wya.(i), uw wza.(i) in
        for j = 0 to n - 1 do
          A.set z j (a *. A.get x j +. A.get y j)
        done
      done

    let v_lin2vectorarray a (wxa : t array) (wya : t array)
                                            (wza : t array) =
      let nvec = Array.length wxa in
      let n = A.length (uw wxa.(0)) in
      for i = 0 to nvec - 1 do
        let x, y, z = uw wxa.(i), uw wya.(i), uw wza.(i) in
        for j = 0 to n - 1 do
          A.set z j (a *. A.get x j -. A.get y j)
        done
      done

    let v_scalesumvectorarray c (wxa : t array) (wya : t array)
                                                (wza : t array) =
      let nvec = Array.length wxa in
      let n = A.length (uw wxa.(0)) in
      for i = 0 to nvec - 1 do
        let x, y, z = uw wxa.(i), uw wya.(i), uw wza.(i) in
        for j = 0 to n - 1 do
          A.set z j (c *. (A.get x j +. A.get y j))
        done
      done

    let v_scalediffvectorarray c (wxa : t array) (wya : t array)
                                                 (wza : t array) =
      let nvec = Array.length wxa in
      let n = A.length (uw wxa.(0)) in
      for i = 0 to nvec - 1 do
        let x, y, z = uw wxa.(i), uw wya.(i), uw wza.(i) in
        for j = 0 to n - 1 do
          A.set z j (c *. (A.get x j -. A.get y j))
        done
      done

    let n_vlinearsumvectorarray a (wxa : t array) b (wya : t array)
                                                    (wza : t array) =
      let nvec = Array.length wya in
      if nvec = 1 then n_vlinearsum a wxa.(0) b wya.(0) wza.(0)
      else if b =  1.0 && (wza == wya) then arr_vaxpy_array a wxa wya
      else if a =  1.0 && (wza == wxa) then arr_vaxpy_array a wya wxa
      else if a =  1.0 && b =  1.0 then v_sumvectorarray wxa wya wza
      else if a =  1.0 && b = -1.0 then v_diffvectorarray wxa wya wza
      else if a = -1.0 && b =  1.0 then v_diffvectorarray wya wxa wza
      else if a =  1.0 then v_lin1vectorarray b wya wxa wza
      else if b =  1.0 then v_lin1vectorarray a wxa wya wza
      else if a = -1.0 then v_lin2vectorarray b wya wxa wza
      else if b = -1.0 then v_lin2vectorarray a wxa wya wza
      else if a = b then v_scalesumvectorarray a wxa wya wza
      else if a = -. b then v_scalediffvectorarray a wxa wya wza
      else
        let n = A.length (uw wxa.(0)) in
        for i = 0 to nvec - 1 do
          let x, y, z = uw wxa.(i), uw wya.(i), uw wza.(i) in
          for j = 0 to n - 1 do
            A.set z j (a *. A.get x j +. b *. A.get y j)
          done
        done

    let n_vscalevectorarray (c : RealArray.t) (wxa : t array)
                                              (wza : t array) =
      let nvec = Array.length wxa in
      if nvec = 1 then n_vscale c.{0} wxa.(0) wza.(0)
      else
        let n = A.length (uw wxa.(0)) in
        if wxa == wza then
          for i = 0 to nvec - 1 do
            let x, c = uw wxa.(i), c.{i} in
            for j = 0 to n - 1 do
              A.set x j (c *. A.get x j)
            done
          done
        else
          for i = 0 to nvec - 1 do
            let x, z, c = uw wxa.(i), uw wza.(i), c.{i} in
            for j = 0 to n - 1 do
              A.set z j (c *. A.get x j)
            done
          done

    let n_vconstvectorarray c (wza : t array) =
      let nvec = Array.length wza in
      if nvec = 1 then n_vconst c wza.(0)
      else
        let n = A.length (uw wza.(0)) in
        for i = 0 to nvec - 1 do
          let z = uw wza.(i) in
          for j = 0 to n - 1 do
            A.set z j c
          done
        done

    let n_vwrmsnormvectorarray (wxa : t array) (wwa : t array)
                                               (nrm : RealArray.t) =
      let nvec = Array.length wxa in
      if nvec = 1 then nrm.{0} <- n_vwrmsnorm wxa.(0) wwa.(0)
      else
        let n = A.length (uw wxa.(0)) in
        let nf = float n in
        let a = ref 0.0 in
        for i = 0 to nvec - 1 do
          let x, w = uw wxa.(i), uw wwa.(i) in
          a := 0.0;
          for j = 0 to n - 1 do
            let s = A.get x j *. A.get w j in
            a := !a +. s *. s
          done;
          nrm.{i} <- sqrt (!a /. nf)
        done

    let n_vwrmsnormmaskvectorarray (wxa : t array) (wwa : t array)
                                   (wid : t) (nrm : RealArray.t) =
      let nvec = Array.length wxa in
      if nvec = 1 then nrm.{0} <- n_vwrmsnormmask wxa.(0) wwa.(0) wid
      else
        let n = A.length (uw wxa.(0)) in
        let nf = float n in
        let id = uw wid in
        let a = ref 0.0 in
        for i = 0 to nvec - 1 do
          let x, w = uw wxa.(i), uw wwa.(i) in
          a := 0.0;
          for j = 0 to n - 1 do
            if A.get id j > 0.0 then begin
              let s = A.get x j *. A.get w j in
              a := !a +. s *. s
            end
          done;
          nrm.{i} <- sqrt (!a /. nf)
        done

    let n_vscaleaddmultivectorarray (ra : RealArray.t) (wxa : t array)
                                    (wyaa : t array array)
                                    (wzaa : t array array) =
      let nsum = Array.length wyaa in
      let nvec = Array.length wyaa.(0) in
      if nvec = 1 then begin
        if nsum = 1 then n_vlinearsum ra.{0} wxa.(0) 1.0 wyaa.(0).(0) wzaa.(0).(0)
        else
          let yya = Array.init nsum (fun j -> wyaa.(j).(0)) in
          let zza = Array.init nsum (fun j -> wzaa.(j).(0)) in
          n_vscaleaddmulti ra wxa.(0) yya zza
      end
      else if nsum = 1 then n_vlinearsumvectorarray ra.{0} wxa 1.0 wyaa.(0) wzaa.(0)
      else
        let n = A.length (uw wxa.(0)) in
        if (wyaa == wzaa) then
          for i = 0 to nvec - 1 do
            let x = uw wxa.(i) in
            for j = 0 to nsum - 1 do
              let a, y = ra.{j}, uw wyaa.(j).(i) in
              for k = 0 to n - 1 do
                A.set y k (A.get y k +. a *. A.get x k)
              done
            done
          done
        else
          for i = 0 to nvec - 1 do
            let x = uw wxa.(i) in
            for j = 0 to nsum - 1 do
              let a, y, z = ra.{j}, uw wyaa.(j).(i), uw wzaa.(j).(i) in
              for k = 0 to n - 1 do
                A.set z k (a *. A.get x k +. A.get y k)
              done
            done
          done

    let n_vlinearcombinationvectorarray (ca : RealArray.t)
                                        (wxaa : t array array)
                                        (wza : t array) =
      let nsum = Array.length wxaa in
      let nvec = Array.length wxaa.(0) in
      if nvec = 1 then begin
        if nsum = 1 then n_vscale ca.{0} wxaa.(0).(0) wza.(0)
        else if nsum = 2
             then n_vlinearsum ca.{0} wxaa.(0).(0) ca.{1} wxaa.(1).(0) wza.(0)
        else
          let wya = Array.init nsum (fun i -> wxaa.(i).(0)) in
          n_vlinearcombination ca wya wza.(0)
      end
      else
        if nsum = 1 then
          let ctmp = RealArray.make nvec ca.{0} in
          n_vscalevectorarray ctmp wxaa.(0) wza
        else if nsum = 2 then
          n_vlinearsumvectorarray ca.{0} wxaa.(0) ca.{1} wxaa.(1) wza
        else
          let n = A.length (uw wza.(0)) in
          if wxaa.(0) == wza then begin
            if ca.{0} = 1.0 then
              for j = 0 to nvec - 1 do
                let z = uw wza.(j) in
                for i = 1 to nsum - 1 do
                  let c, x = ca.{i}, uw wxaa.(i).(j) in
                  for k = 0 to n - 1 do
                    A.set z k (A.get z k +. c *. A.get x k)
                  done
                done
              done
            else
              let c0 = ca.{0} in
              for j = 0 to nvec - 1 do
                let z = uw wza.(j) in
                for k = 0 to n - 1 do
                  A.set z k (A.get z k *. c0)
                done;
                for i = 1 to nsum - 1 do
                  let c, x = ca.{i}, uw wxaa.(i).(j) in
                  for k = 0 to n - 1 do
                    A.set z k (A.get z k +. c *. A.get x k)
                  done
                done
              done
          end
          else
            let c0 = ca.{0} in
            for j = 0 to nvec - 1 do
              let x, z = uw wxaa.(0).(j), uw wza.(j) in
              for k = 0 to n - 1 do
                A.set z k (c0 *. A.get x k)
              done;
              for i = 1 to nsum - 1 do
                let c, x = ca.{i}, uw wxaa.(i).(j) in
                for k = 0 to n - 1 do
                  A.set z k (A.get z k +. c *. A.get x k)
                done
              done
            done

    module Local = struct
      let n_vdotprod     = n_vdotprod
      let n_vmaxnorm     = n_vmaxnorm
      let n_vmin         = n_vmin
      let n_vl1norm      = n_vl1norm
      let n_vinvtest     = n_vinvtest
      let n_vconstrmask  = n_vconstrmask
      let n_vminquotient = n_vminquotient

      let n_vwsqrsum wx ww =
        let x, w = uw wx, uw ww in
        let a = ref 0.0 in
        let lx = A.length x in
        for i = 0 to lx - 1 do
          a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
        done;
        !a

      let n_vwsqrsummask wx ww wid =
        let x, w, id = uw wx, uw ww, uw wid in
        let a = ref 0.0 in
        let lx = A.length x in
        for i = 0 to lx - 1 do
          if A.get id i > 0.0 then
            a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
        done;
        !a
    end

    let checkfn wv1 =
      let l = A.length (uw wv1) in
      (fun wv2 -> l = A.length (uw wv2))

    let array_nvec_ops = { (* {{{ *)
        Nvector_custom.n_vcheck        = checkfn;
        Nvector_custom.n_vclone        = n_vclone;
        Nvector_custom.n_vspace        = Some n_vspace;
        Nvector_custom.n_vgetlength    = n_vgetlength;
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
        Nvector_custom.n_vdotprod      = n_vdotprod;
        Nvector_custom.n_vcompare      = n_vcompare;
        Nvector_custom.n_vinvtest      = n_vinvtest;

        Nvector_custom.n_vwl2norm      = Some n_vwl2norm;
        Nvector_custom.n_vl1norm       = Some n_vl1norm;
        Nvector_custom.n_vwrmsnormmask = Some n_vwrmsnormmask;
        Nvector_custom.n_vconstrmask   = Some n_vconstrmask;
        Nvector_custom.n_vminquotient  = Some n_vminquotient;

        Nvector_custom.n_vgetcommunicator = None;

        Nvector_custom.n_vlinearcombination = Some n_vlinearcombination;
        Nvector_custom.n_vscaleaddmulti = Some n_vscaleaddmulti;
        Nvector_custom.n_vdotprodmulti = Some n_vdotprodmulti;

        Nvector_custom.n_vlinearsumvectorarray
          = Some n_vlinearsumvectorarray;
        Nvector_custom.n_vscalevectorarray
          = Some n_vscalevectorarray;
        Nvector_custom.n_vconstvectorarray
          = Some n_vconstvectorarray;
        Nvector_custom.n_vwrmsnormvectorarray
          = Some n_vwrmsnormvectorarray;
        Nvector_custom.n_vwrmsnormmaskvectorarray
          = Some n_vwrmsnormmaskvectorarray;
        Nvector_custom.n_vscaleaddmultivectorarray
          = Some n_vscaleaddmultivectorarray;
        Nvector_custom.n_vlinearcombinationvectorarray
          = Some n_vlinearcombinationvectorarray;

        Nvector_custom.n_vdotprod_local = Some Local.n_vdotprod;
        Nvector_custom.n_vmaxnorm_local = Some Local.n_vmaxnorm;
        Nvector_custom.n_vmin_local = Some Local.n_vmin;
        Nvector_custom.n_vl1norm_local = Some Local.n_vl1norm;
        Nvector_custom.n_vinvtest_local = Some Local.n_vinvtest;
        Nvector_custom.n_vconstrmask_local = Some Local.n_vconstrmask;
        Nvector_custom.n_vminquotient_local = Some Local.n_vminquotient;
        Nvector_custom.n_vwsqrsum_local = Some Local.n_vwsqrsum;
        Nvector_custom.n_vwsqrsummask_local = Some Local.n_vwsqrsummask;
      } (* }}} *)
  end (* }}} *)

module Make =
  functor (A : ArrayOps) ->
  struct (* {{{ *)
    type data = A.data
    type kind = Nvector_custom.kind
    type t = data Nvector_custom.t

    let checkfn v1 =
      let l = A.length v1 in
      (fun v2 -> l = A.length v2)

    module DataOps = MakeOps (struct
        include A
        type t = A.data
        let project x = x
        let inject x = x
      end)

    let array_nvec_ops = DataOps.array_nvec_ops

    let wrap = Nvector_custom.make_wrap DataOps.array_nvec_ops
      (* (Nvector.Mutable.add_tracing "::" array_nvec_ops) *)

    let make n e = wrap (A.make n e)

    let unwrap = Nvector.unwrap

    let enable = Nvector_custom.enable

    module Any = struct (* {{{ *)

      type Nvector.gdata += Arr of data

      module AnyOps = MakeOps (struct
          include A
          type t = Nvector.gdata
          let project = function Arr x -> x | _ -> raise Nvector.BadGenericType
          let inject x = Arr x
        end)

      let wrap =
        let ops = AnyOps.array_nvec_ops in
        Nvector_custom.Any.make_wrap ops ~inject:(fun x -> Arr x)

    end (* }}} *)

    module Ops = struct
      include MakeOps (struct
        include A
        type t = A.data Nvector_custom.t
        let project = unwrap
        let inject x = assert false (* n_vclone is overridden *)
      end)

      (* Cloning in this way ensures that the enabled status of fused and
         array operations is the same in the clone nvector. *)
      let n_vclone x = Nvector_custom.clone x
    end

  end (* }}} *)

(* (* The following Array module can be created like this SlowerArray,
      albeit with significant performance hits.  It may be useful for
      prototyping.  *)
module SlowerArray = Make (
  struct
    type data = float array
    include Array
    let clone = Array.copy
    let fill a c = fill a 0 (length a) c
  end)
*)

module Array =
  struct (* {{{ *)
    type data = float array
    type kind = Nvector_custom.kind
    type t = data Nvector_custom.t

    let checkfn v1 =
      let l = Array.length v1 in
      (fun v2 -> l = Array.length v2)

    module DataOps = struct (* {{{ *)
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
        else if (a = 1.0 && b = -1.0) || (a = -1.0 && b = 1.0) then
          let v1, v2 = if (a = 1.0 && b = -1.0) then y, x else x, y in
          for i = 0 to A.length v1 - 1 do
            A.set z i (A.get v2 i -. A.get v1 i)
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
        for i = 0 to A.length x - 1 do
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
        let m = ref Config.big_real in
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

      let n_vspace (x : float array) = (A.length x, 1)

      let n_vgetlength (x : float array) = A.length x

      (* fused and array operations *)

      let n_vlinearcombination (ca : RealArray.t)
                               (xa : (float array) array)
                               (z  : float array) =
        let nvec = Array.length xa in
        if nvec = 1 then n_vscale ca.{0} xa.(0) z
        else if nvec = 2 then n_vlinearsum ca.{0} xa.(0) ca.{1} xa.(1) z
        else
          let n = A.length z in
          if xa.(0) == z then begin
            let c0 = ca.{0} in
            if c0 <> 1.0 then
              for j = 0 to n - 1 do
                A.set z j (A.get z j *. c0)
              done;
            for i = 1 to nvec - 1 do
              let ci, x = ca.{i}, xa.(i) in
              for j = 0 to n - 1 do
                A.set z j (A.get z j +. ci *. A.get x j)
              done
            done
          end
          else begin
            let c0, x = ca.{0}, xa.(0) in
            for j = 0 to n - 1 do
              A.set z j (c0 *. A.get x j)
            done;
            for i = 1 to nvec - 1 do
              let ci, x = ca.{i}, xa.(i) in
              for j = 0 to n - 1 do
                A.set z j (A.get z j +. ci *. A.get x j)
              done
            done
          end

      let n_vscaleaddmulti (aa : RealArray.t)
                           (x : float array)
                           (ya : (float array) array)
                           (za : (float array) array) =
        let nvec = Array.length ya in
        if nvec = 1 then n_vlinearsum aa.{0} x 1.0 ya.(0) za.(0)
        else
          let n = A.length x in
          if ya == za then
            for i = 0 to nvec - 1 do
              let a, y = aa.{i}, ya.(i) in
              for j = 0 to n - 1 do
                A.set y j (A.get y j +. a *. A.get x j)
              done
            done
          else
            for i = 0 to nvec - 1 do
              let ai, y, z = aa.{i}, ya.(i), za.(i) in
              for j = 0 to n - 1 do
                A.set z j (ai *. A.get x j +. A.get y j)
              done
            done

      let n_vdotprodmulti (x : float array)
                          (ya : (float array) array)
                          (dp : RealArray.t) =
        let nvec = Array.length ya in
        if nvec = 1 then dp.{0} <- n_vdotprod x ya.(0)
        else
          let n = A.length x in
          for i = 0 to nvec - 1 do
            let y = ya.(i) in
            dp.{i} <- 0.0;
            for j = 0 to n - 1 do
              dp.{i} <- dp.{i} +. A.get x j *. A.get y j
            done
          done

      let arr_vaxpy_array a (xa : (float array) array)
                            (ya : (float array) array) =
        let nvec = Array.length xa in
        let n = A.length xa.(0) in
        if a = 1.0 then
          for i = 0 to nvec - 1 do
            let x, y = xa.(i), ya.(i) in
            for j = 0 to n - 1 do
              A.set y j (A.get y j +. A.get x j)
            done
          done
        else if a = -1.0 then
          for i = 0 to nvec - 1 do
            let x, y = xa.(i), ya.(i) in
            for j = 0 to n - 1 do
              A.set y j (A.get y j -. A.get x j)
            done
          done
        else
          for i = 0 to nvec - 1 do
            let x, y = xa.(i), ya.(i) in
            for j = 0 to n - 1 do
              A.set y j (A.get y j +. a *. A.get x j)
            done
          done

      let v_sumvectorarray (xa : (float array) array)
                           (ya : (float array) array)
                           (za : (float array) array) =
        let nvec = Array.length xa in
        let n = A.length xa.(0) in
        for i = 0 to nvec - 1 do
          let x, y, z = xa.(i), ya.(i), za.(i) in
          for j = 0 to n - 1 do
            A.set z j (A.get x j +. A.get y j)
          done
        done

      let v_diffvectorarray (xa : (float array) array)
                            (ya : (float array) array)
                            (za : (float array) array) =
        let nvec = Array.length xa in
        let n = A.length xa.(0) in
        for i = 0 to nvec - 1 do
          let x, y, z = xa.(i), ya.(i), za.(i) in
          for j = 0 to n - 1 do
            A.set z j (A.get x j -. A.get y j)
          done
        done

      let v_lin1vectorarray a (xa : (float array) array)
                              (ya : (float array) array)
                              (za : (float array) array) =
        let nvec = Array.length xa in
        let n = A.length xa.(0) in
        for i = 0 to nvec - 1 do
          let x, y, z = xa.(i), ya.(i), za.(i) in
          for j = 0 to n - 1 do
            A.set z j (a *. A.get x j +. A.get y j)
          done
        done

      let v_lin2vectorarray a (xa : (float array) array)
                              (ya : (float array) array)
                              (za : (float array) array) =
        let nvec = Array.length xa in
        let n = A.length xa.(0) in
        for i = 0 to nvec - 1 do
          let x, y, z = xa.(i), ya.(i), za.(i) in
          for j = 0 to n - 1 do
            A.set z j (a *. A.get x j -. A.get y j)
          done
        done

      let v_scalesumvectorarray c (xa : (float array) array)
                                  (ya : (float array) array)
                                  (za : (float array) array) =
        let nvec = Array.length xa in
        let n = A.length xa.(0) in
        for i = 0 to nvec - 1 do
          let x, y, z = xa.(i), ya.(i), za.(i) in
          for j = 0 to n - 1 do
            A.set z j (c *. (A.get x j +. A.get y j))
          done
        done

      let v_scalediffvectorarray c (xa : (float array) array)
                                   (ya : (float array) array)
                                   (za : (float array) array) =
        let nvec = Array.length xa in
        let n = A.length xa.(0) in
        for i = 0 to nvec - 1 do
          let x, y, z = xa.(i), ya.(i), za.(i) in
          for j = 0 to n - 1 do
            A.set z j (c *. (A.get x j -. A.get y j))
          done
        done

      let n_vlinearsumvectorarray a (xa : (float array) array)
                                  b (ya : (float array) array)
                                    (za : (float array) array) =
        let nvec = Array.length ya in
        if nvec = 1 then n_vlinearsum a xa.(0) b ya.(0) za.(0)
        else if b =  1.0 && (za == ya) then arr_vaxpy_array a xa ya
        else if a =  1.0 && (za == xa) then arr_vaxpy_array a ya xa
        else if a =  1.0 && b =  1.0 then v_sumvectorarray xa ya za
        else if a =  1.0 && b = -1.0 then v_diffvectorarray xa ya za
        else if a = -1.0 && b =  1.0 then v_diffvectorarray ya xa za
        else if a =  1.0 then v_lin1vectorarray b ya xa za
        else if b =  1.0 then v_lin1vectorarray a xa ya za
        else if a = -1.0 then v_lin2vectorarray b ya xa za
        else if b = -1.0 then v_lin2vectorarray a xa ya za
        else if a = b then v_scalesumvectorarray a xa ya za
        else if a = -. b then v_scalediffvectorarray a xa ya za
        else
          let n = A.length xa.(0) in
          for i = 0 to nvec - 1 do
            let x, y, z = xa.(i), ya.(i), za.(i) in
            for j = 0 to n - 1 do
              A.set z j (a *. A.get x j +. b *. A.get y j)
            done
          done

      let n_vscalevectorarray (c  : RealArray.t)
                              (xa : (float array) array)
                              (za : (float array) array) =
        let nvec = Array.length xa in
        if nvec = 1 then n_vscale c.{0} xa.(0) za.(0)
        else
          let n = A.length xa.(0) in
          if xa == za then
            for i = 0 to nvec - 1 do
              let x, c = xa.(i), c.{i} in
              for j = 0 to n - 1 do
                A.set x j (c *. A.get x j)
              done
            done
          else
            for i = 0 to nvec - 1 do
              let x, z, c = xa.(i), za.(i), c.{i} in
              for j = 0 to n - 1 do
                A.set z j (c *. A.get x j)
              done
            done

      let n_vconstvectorarray c (za : float array array) =
        let nvec = Array.length za in
        if nvec = 1 then n_vconst c za.(0)
        else
          let n = A.length za.(0) in
          for i = 0 to nvec - 1 do
            let z = za.(i) in
            for j = 0 to n - 1 do
              A.set z j c
            done
          done

      let n_vwrmsnormvectorarray (xa  : (float array) array)
                                 (wa  : (float array) array)
                                 (nrm : RealArray.t) =
        let nvec = Array.length xa in
        if nvec = 1 then nrm.{0} <- n_vwrmsnorm xa.(0) wa.(0)
        else
          let n = A.length xa.(0) in
          let nf = float n in
          let a = ref 0.0 in
          for i = 0 to nvec - 1 do
            let x, w = xa.(i), wa.(i) in
            a := 0.0;
            for j = 0 to n - 1 do
              let s = A.get x j *. A.get w j in
              a := !a +. s *. s
            done;
            nrm.{i} <- sqrt (!a /. nf)
          done

      let n_vwrmsnormmaskvectorarray (xa  : (float array) array)
                                     (wa  : (float array) array)
                                     (id  : float array)
                                     (nrm : RealArray.t) =
        let nvec = Array.length xa in
        if nvec = 1 then nrm.{0} <- n_vwrmsnormmask xa.(0) wa.(0) id
        else
          let n = A.length xa.(0) in
          let nf = float n in
          let a = ref 0.0 in
          for i = 0 to nvec - 1 do
            let x, w = xa.(i), wa.(i) in
            a := 0.0;
            for j = 0 to n - 1 do
              if A.get id j > 0.0 then begin
                let s = A.get x j *. A.get w j in
                a := !a +. s *. s
              end
            done;
            nrm.{i} <- sqrt (!a /. nf)
          done

      let n_vscaleaddmultivectorarray (ra  : RealArray.t)
                                      (xa  : (float array) array)
                                      (yaa : (float array) array array)
                                      (zaa : (float array) array array) =
        let nsum = Array.length yaa in
        let nvec = Array.length yaa.(0) in
        if nvec = 1 then begin
          if nsum = 1 then n_vlinearsum ra.{0} xa.(0) 1.0 yaa.(0).(0) zaa.(0).(0)
          else
            let yya = Array.init nsum (fun j -> yaa.(j).(0)) in
            let zza = Array.init nsum (fun j -> zaa.(j).(0)) in
            n_vscaleaddmulti ra xa.(0) yya zza
        end
        else if nsum = 1 then n_vlinearsumvectorarray ra.{0} xa 1.0 yaa.(0) zaa.(0)
        else
          let n = A.length xa.(0) in
          if (yaa == zaa) then
            for i = 0 to nvec - 1 do
              let x = xa.(i) in
              for j = 0 to nsum - 1 do
                let a, y = ra.{j}, yaa.(j).(i) in
                for k = 0 to n - 1 do
                  A.set y k (A.get y k +. a *. A.get x k)
                done
              done
            done
          else
            for i = 0 to nvec - 1 do
              let x = xa.(i) in
              for j = 0 to nsum - 1 do
                let a, y, z = ra.{j}, yaa.(j).(i), zaa.(j).(i) in
                for k = 0 to n - 1 do
                  A.set z k (a *. A.get x k +. A.get y k)
                done
              done
            done

      let n_vlinearcombinationvectorarray (ca  : RealArray.t)
                                          (xaa : (float array) array array)
                                          (za  : (float array) array) =
        let nsum = Array.length xaa in
        let nvec = Array.length xaa.(0) in
        if nvec = 1 then begin
          if nsum = 1 then n_vscale ca.{0} xaa.(0).(0) za.(0)
          else if nsum = 2
               then n_vlinearsum ca.{0} xaa.(0).(0) ca.{1} xaa.(1).(0) za.(0)
          else
            let ya = Array.init nsum (fun i -> xaa.(i).(0)) in
            n_vlinearcombination ca ya za.(0)
        end
        else
          if nsum = 1 then
            let ctmp = RealArray.make nvec ca.{0} in
            n_vscalevectorarray ctmp xaa.(0) za
          else if nsum = 2 then
            n_vlinearsumvectorarray ca.{0} xaa.(0) ca.{1} xaa.(1) za
          else
            let n = A.length za.(0) in
            if xaa.(0) == za then begin
              if ca.{0} = 1.0 then
                for j = 0 to nvec - 1 do
                  let z = za.(j) in
                  for i = 1 to nsum - 1 do
                    let c, x = ca.{i}, xaa.(i).(j) in
                    for k = 0 to n - 1 do
                      A.set z k (A.get z k +. c *. A.get x k)
                    done
                  done
                done
              else
                let c0 = ca.{0} in
                for j = 0 to nvec - 1 do
                  let z = za.(j) in
                  for k = 0 to n - 1 do
                    A.set z k (A.get z k *. c0)
                  done;
                  for i = 1 to nsum - 1 do
                    let c, x = ca.{i}, xaa.(i).(j) in
                    for k = 0 to n - 1 do
                      A.set z k (A.get z k +. c *. A.get x k)
                    done
                  done
                done
            end
            else
              let c0 = ca.{0} in
              for j = 0 to nvec - 1 do
                let x, z = xaa.(0).(j), za.(j) in
                for k = 0 to n - 1 do
                  A.set z k (c0 *. A.get x k)
                done;
                for i = 1 to nsum - 1 do
                  let c, x = ca.{i}, xaa.(i).(j) in
                  for k = 0 to n - 1 do
                    A.set z k (A.get z k +. c *. A.get x k)
                  done
                done
              done

      module Local = struct
        let n_vdotprod     = n_vdotprod
        let n_vmaxnorm     = n_vmaxnorm
        let n_vmin         = n_vmin
        let n_vl1norm      = n_vl1norm
        let n_vinvtest     = n_vinvtest
        let n_vconstrmask  = n_vconstrmask
        let n_vminquotient = n_vminquotient

        let n_vwsqrsum x w =
          let a = ref 0.0 in
          let lx = A.length x in
          for i = 0 to lx - 1 do
            a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
          done;
          !a

        let n_vwsqrsummask x w id =
          let a = ref 0.0 in
          let lx = A.length x in
          for i = 0 to lx - 1 do
            if A.get id i > 0.0 then
              a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
          done;
          !a
      end
    end (* }}} *)

    let array_nvec_ops = { (* {{{ *)
          Nvector_custom.n_vcheck        = checkfn;
          Nvector_custom.n_vclone        = DataOps.n_vclone;
          Nvector_custom.n_vspace        = Some DataOps.n_vspace;
          Nvector_custom.n_vgetlength    = DataOps.n_vgetlength;
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

          Nvector_custom.n_vgetcommunicator = None;

          Nvector_custom.n_vlinearcombination
            = Some DataOps.n_vlinearcombination;
          Nvector_custom.n_vscaleaddmulti
            = Some DataOps.n_vscaleaddmulti;
          Nvector_custom.n_vdotprodmulti
            = Some DataOps.n_vdotprodmulti;
          Nvector_custom.n_vlinearsumvectorarray
            = Some DataOps.n_vlinearsumvectorarray;
          Nvector_custom.n_vscalevectorarray
            = Some DataOps.n_vscalevectorarray;
          Nvector_custom.n_vconstvectorarray
            = Some DataOps.n_vconstvectorarray;
          Nvector_custom.n_vwrmsnormvectorarray
            = Some DataOps.n_vwrmsnormvectorarray;
          Nvector_custom.n_vwrmsnormmaskvectorarray
            = Some DataOps.n_vwrmsnormmaskvectorarray;
          Nvector_custom.n_vscaleaddmultivectorarray
            = Some DataOps.n_vscaleaddmultivectorarray;
          Nvector_custom.n_vlinearcombinationvectorarray
            = Some DataOps.n_vlinearcombinationvectorarray;

          Nvector_custom.n_vdotprod_local
            = Some DataOps.Local.n_vdotprod;

          Nvector_custom.n_vmaxnorm_local
            = Some DataOps.Local.n_vmaxnorm;

          Nvector_custom.n_vmin_local
            = Some DataOps.Local.n_vmin;

          Nvector_custom.n_vl1norm_local
            = Some DataOps.Local.n_vl1norm;

          Nvector_custom.n_vinvtest_local
            = Some DataOps.Local.n_vinvtest;

          Nvector_custom.n_vconstrmask_local
            = Some DataOps.Local.n_vconstrmask;

          Nvector_custom.n_vminquotient_local
            = Some DataOps.Local.n_vminquotient;

          Nvector_custom.n_vwsqrsum_local
            = Some DataOps.Local.n_vwsqrsum;

          Nvector_custom.n_vwsqrsummask_local
            = Some DataOps.Local.n_vwsqrsummask;
    } (* }}} *)

    let make n e =
      Nvector_custom.make_wrap array_nvec_ops (Array.make n e)

    let wrap = Nvector_custom.make_wrap array_nvec_ops
      (* (Nvector.Mutable.add_tracing "::" array_nvec_ops) *)

    let enable = Nvector_custom.enable

    let unwrap = Nvector.unwrap

    module Any = struct (* {{{ *)

      type Nvector.gdata += Arr of data

      let project = function Arr x -> x | _ -> raise Nvector.BadGenericType
      let inject x = Arr x

      module AnyOps = MakeOps (struct
          type data = float array

          let get (x : float array) = Array.get x
          let set (x : float array) = Array.set x
          let fill (x : float array) = Array.fill x 0 (Array.length x - 1)
          let make = Array.make
          let clone (x : float array) = Array.copy x
          let length (x : float array) = Array.length x

          type t = Nvector.gdata
          let project = project
          let inject = inject
        end)

      let wrap =
        let ops = {
          (* Use non-polymorphic versions for most operations... *)
          (Nvector_custom.Any.convert_ops ~inject ~project array_nvec_ops)
          with
            (* ...but avoid generating intermediate unwrapped arrays. *)
            Nvector_custom.n_vlinearcombination
              = Some AnyOps.n_vlinearcombination;
            Nvector_custom.n_vscaleaddmulti
              = Some AnyOps.n_vscaleaddmulti;
            Nvector_custom.n_vdotprodmulti
              = Some AnyOps.n_vdotprodmulti;
            Nvector_custom.n_vlinearsumvectorarray
              = Some AnyOps.n_vlinearsumvectorarray;
            Nvector_custom.n_vscalevectorarray
              = Some AnyOps.n_vscalevectorarray;
            Nvector_custom.n_vconstvectorarray
              = Some AnyOps.n_vconstvectorarray;
            Nvector_custom.n_vwrmsnormvectorarray
              = Some AnyOps.n_vwrmsnormvectorarray;
            Nvector_custom.n_vwrmsnormmaskvectorarray
              = Some AnyOps.n_vwrmsnormmaskvectorarray;
            Nvector_custom.n_vscaleaddmultivectorarray
              = Some AnyOps.n_vscaleaddmultivectorarray;
            Nvector_custom.n_vlinearcombinationvectorarray
              = Some AnyOps.n_vlinearcombinationvectorarray;
        }
        in
        Nvector_custom.Any.make_wrap ops ~inject

    end (* }}} *)

    module Ops = struct

      include MakeOps (struct
        type data = float array

        let get (x : float array) = Array.get x
        let set (x : float array) = Array.set x
        let fill (x : float array) = Array.fill x 0 (Array.length x - 1)
        let make = Array.make
        let clone (x : float array) = Array.copy x
        let length (x : float array) = Array.length x

        type t = float array Nvector_custom.t
        let project = unwrap
        let inject x = assert false (* n_vclone is overridden *)
      end)

      (* Cloning in this way ensures that the enabled status of fused and
         array operations is the same in the clone nvector. *)
      let n_vclone x = Nvector_custom.clone x
    end

  end (* }}} *)

include Array

