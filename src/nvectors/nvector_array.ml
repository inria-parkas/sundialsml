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

    let clone wx = A.inject (A.clone (A.project wx))

    let uw = A.project

    let floatmin (x : float) (y : float) = min x y

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

    let linearsum a wx b wy wz =
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

    let const c wa = A.fill (uw wa) c

    let scale c wx wz =
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

    let addconst wx b wz =
      let x, z = uw wx, uw wz in
      for i = 0 to A.length x - 1 do
        A.set z i (A.get x i +. b)
      done

    let maxnorm wx =
      let x = uw wx in
      let max = ref 0.0 in
      for i = 0 to A.length x - 1 do
        let ax = abs_float (A.get x i) in
        if ax > !max then max := ax
      done;
      !max

    let wrmsnorm wx ww =
      let x, w = uw wx, uw ww in
      let a = ref 0.0 in
      let lx = A.length x in
      for i = 0 to lx - 1 do
        a := !a +. ((A.get x i) *. (A.get w i) *. (A.get x i) *. (A.get w i))
      done;
      sqrt (!a /. float lx)

    let wrmsnormmask wx ww wid =
      let x, w, id = uw wx, uw ww, uw wid in
      let a = ref 0.0 in
      let lx = A.length x in
      for i = 0 to lx - 1 do
        if A.get id i > 0.0 then
          a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
      done;
      sqrt (!a /. float lx)

    let min wx =
      let x = uw wx in
      let min = ref max_float in
      for i = 0 to A.length x - 1 do
        let xv = A.get x i in
        if xv < !min then min := xv
      done;
      !min

    let dotprod wx wy =
      let x, y = uw wx, uw wy in
      let a = ref 0.0 in
      for i = 0 to A.length x - 1 do
        a := !a +. (A.get x i *. A.get y i)
      done;
      !a

    let compare c wx wz =
      let x, z = uw wx, uw wz in
      for i = 0 to A.length x - 1 do
        A.set z i (if abs_float (A.get x i) >= c then 1.0 else 0.0)
      done

    let invtest wx wz =
      let x, z = uw wx, uw wz in
      let l = A.length x in
      let rec f r i =
        if i = l then r
        else if A.get x i = 0.0 then f false (i + 1)
        else (A.set z i (1.0 /. (A.get x i)); f r (i + 1))
      in f true 0

    let wl2norm wx ww =
      let x, w = uw wx, uw ww in
      let a = ref 0.0 in
      for i = 0 to A.length x - 1 do
        a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
      done;
      sqrt !a

    let l1norm wx =
      let x = uw wx in
      let a = ref 0.0 in
      for i = 0 to A.length x - 1 do
        a := !a +. abs_float (A.get x i)
      done;
      !a

    let constrmask wc wx wm =
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

    let minquotient wn wd =
      let n, d = uw wn, uw wd in
      let m = ref Config.big_real in
      for i = 0 to A.length n - 1 do
        if (A.get d i) <> 0.0 then
          m := floatmin !m (A.get n i /. A.get d i)
      done;
      !m

    let prod wx wy wz =
      let x, y, z = uw wx, uw wy, uw wz in
      for i = 0 to A.length x - 1 do
        A.set z i (A.get x i *. A.get y i)
      done

    let div wx wy wz =
      let x, y, z = uw wx, uw wy, uw wz in
      for i = 0 to A.length x - 1 do
        A.set z i (A.get x i /. A.get y i)
      done

    let abs wx wz =
      let x, z = uw wx, uw wz in
      for i = 0 to A.length x - 1 do
        A.set z i (abs_float (A.get x i))
      done

    let inv wx wz =
      let x, z = uw wx, uw wz in
      for i = 0 to A.length x - 1 do
        A.set z i (1.0 /. (A.get x i))
      done

    let space wx = (A.length (uw wx), 1)

    let getlength wx = A.length (uw wx)

    let print ?(logfile=Logfile.stdout) wx =
      let x = uw wx in
      for i = 0 to A.length x - 1 do
        Logfile.output_string logfile (Printf.sprintf "%19.16g" (A.get x i));
        Logfile.output_string logfile "\n"
      done;
      Logfile.output_string logfile "\n"

    (* fused and array operations *)

    let linearcombination (ca : RealArray.t) (wxa : t array) (wz : t) =
      let nvec = Array.length wxa in
      if nvec = 1 then scale ca.{0} wxa.(0) wz
      else if nvec = 2 then linearsum ca.{0} wxa.(0) ca.{1} wxa.(1) wz
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

    let scaleaddmulti (aa : RealArray.t) (wx : t) (wya : t array)
                                                     (wza : t array) =
      let nvec = Array.length wya in
      if nvec = 1 then linearsum aa.{0} wx 1.0 wya.(0) wza.(0)
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

    let dotprodmulti (wx : t) (wya : t array) (dp : RealArray.t) =
      let nvec = Array.length wya in
      if nvec = 1 then dp.{0} <- dotprod wx wya.(0)
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

    let linearsumvectorarray a (wxa : t array) b (wya : t array)
                                                    (wza : t array) =
      let nvec = Array.length wya in
      if nvec = 1 then linearsum a wxa.(0) b wya.(0) wza.(0)
      else if b =  1.0 && (wza == wya) then arr_vaxpy_array a wxa wya
      else if a =  1.0 && (wza == wxa) then arr_vaxpy_array b wya wxa
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

    let scalevectorarray (c : RealArray.t) (wxa : t array)
                                              (wza : t array) =
      let nvec = Array.length wxa in
      if nvec = 1 then scale c.{0} wxa.(0) wza.(0)
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

    let constvectorarray c (wza : t array) =
      let nvec = Array.length wza in
      if nvec = 1 then const c wza.(0)
      else
        let n = A.length (uw wza.(0)) in
        for i = 0 to nvec - 1 do
          let z = uw wza.(i) in
          for j = 0 to n - 1 do
            A.set z j c
          done
        done

    let wrmsnormvectorarray (wxa : t array) (wwa : t array)
                                               (nrm : RealArray.t) =
      let nvec = Array.length wxa in
      if nvec = 1 then nrm.{0} <- wrmsnorm wxa.(0) wwa.(0)
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

    let wrmsnormmaskvectorarray (wxa : t array) (wwa : t array)
                                   (wid : t) (nrm : RealArray.t) =
      let nvec = Array.length wxa in
      if nvec = 1 then nrm.{0} <- wrmsnormmask wxa.(0) wwa.(0) wid
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

    let scaleaddmultivectorarray (ra : RealArray.t) (wxa : t array)
                                    (wyaa : t array array)
                                    (wzaa : t array array) =
      let nsum = Array.length wyaa in
      let nvec = Array.length wyaa.(0) in
      if nvec = 1 then begin
        if nsum = 1 then linearsum ra.{0} wxa.(0) 1.0 wyaa.(0).(0) wzaa.(0).(0)
        else
          let yya = Array.init nsum (fun j -> wyaa.(j).(0)) in
          let zza = Array.init nsum (fun j -> wzaa.(j).(0)) in
          scaleaddmulti ra wxa.(0) yya zza
      end
      else if nsum = 1 then linearsumvectorarray ra.{0} wxa 1.0 wyaa.(0) wzaa.(0)
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

    let linearcombinationvectorarray (ca : RealArray.t)
                                        (wxaa : t array array)
                                        (wza : t array) =
      let nsum = Array.length wxaa in
      let nvec = Array.length wxaa.(0) in
      if nvec = 1 then begin
        if nsum = 1 then scale ca.{0} wxaa.(0).(0) wza.(0)
        else if nsum = 2
             then linearsum ca.{0} wxaa.(0).(0) ca.{1} wxaa.(1).(0) wza.(0)
        else
          let wya = Array.init nsum (fun i -> wxaa.(i).(0)) in
          linearcombination ca wya wza.(0)
      end
      else
        if nsum = 1 then
          let ctmp = RealArray.make nvec ca.{0} in
          scalevectorarray ctmp wxaa.(0) wza
        else if nsum = 2 then
          linearsumvectorarray ca.{0} wxaa.(0) ca.{1} wxaa.(1) wza
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
      let dotprod     = dotprod
      let maxnorm     = maxnorm
      let min         = min
      let l1norm      = l1norm
      let invtest     = invtest
      let constrmask  = constrmask
      let minquotient = minquotient

      let wsqrsum wx ww =
        let x, w = uw wx, uw ww in
        let a = ref 0.0 in
        let lx = A.length x in
        for i = 0 to lx - 1 do
          a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
        done;
        !a

      let wsqrsummask wx ww wid =
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
        Nvector_custom.check        = checkfn;
        Nvector_custom.clone        = clone;
        Nvector_custom.space        = Some space;
        Nvector_custom.getlength    = getlength;
        Nvector_custom.print        = Some (fun x logfile -> print ?logfile x);
        Nvector_custom.linearsum    = linearsum;
        Nvector_custom.const        = const;
        Nvector_custom.prod         = prod;
        Nvector_custom.div          = div;
        Nvector_custom.scale        = scale;
        Nvector_custom.abs          = abs;
        Nvector_custom.inv          = inv;
        Nvector_custom.addconst     = addconst;
        Nvector_custom.maxnorm      = maxnorm;
        Nvector_custom.wrmsnorm     = wrmsnorm;
        Nvector_custom.min          = min;
        Nvector_custom.dotprod      = dotprod;
        Nvector_custom.compare      = compare;
        Nvector_custom.invtest      = invtest;

        Nvector_custom.wl2norm      = Some wl2norm;
        Nvector_custom.l1norm       = Some l1norm;
        Nvector_custom.wrmsnormmask = Some wrmsnormmask;
        Nvector_custom.constrmask   = Some constrmask;
        Nvector_custom.minquotient  = Some minquotient;

        Nvector_custom.getcommunicator = None;

        Nvector_custom.linearcombination = Some linearcombination;
        Nvector_custom.scaleaddmulti = Some scaleaddmulti;
        Nvector_custom.dotprodmulti = Some dotprodmulti;

        Nvector_custom.linearsumvectorarray
          = Some linearsumvectorarray;
        Nvector_custom.scalevectorarray
          = Some scalevectorarray;
        Nvector_custom.constvectorarray
          = Some constvectorarray;
        Nvector_custom.wrmsnormvectorarray
          = Some wrmsnormvectorarray;
        Nvector_custom.wrmsnormmaskvectorarray
          = Some wrmsnormmaskvectorarray;
        Nvector_custom.scaleaddmultivectorarray
          = Some scaleaddmultivectorarray;
        Nvector_custom.linearcombinationvectorarray
          = Some linearcombinationvectorarray;

        Nvector_custom.dotprod_local = Some Local.dotprod;
        Nvector_custom.maxnorm_local = Some Local.maxnorm;
        Nvector_custom.min_local = Some Local.min;
        Nvector_custom.l1norm_local = Some Local.l1norm;
        Nvector_custom.invtest_local = Some Local.invtest;
        Nvector_custom.constrmask_local = Some Local.constrmask;
        Nvector_custom.minquotient_local = Some Local.minquotient;
        Nvector_custom.wsqrsum_local = Some Local.wsqrsum;
        Nvector_custom.wsqrsummask_local = Some Local.wsqrsummask;
      } (* }}} *)
  end (* }}} *)

module Make =
  functor (A : ArrayOps) ->
  struct (* {{{ *)
    type data = A.data
    type kind = Nvector_custom.kind
    type t = data Nvector_custom.t

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
        let inject _ = assert false (* clone is overridden *)
      end)

      (* Cloning in this way ensures that the enabled status of fused and
         array operations is the same in the clone nvector. *)
      let clone x = Nvector.clone x
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

      let clone = Array.copy

      let floatmin (x : float) (y : float) = min x y

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

      let linearsum a x b y z =
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

      let const c a = fill a c

      let scale c x z =
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

      let addconst x b z =
        for i = 0 to A.length x - 1 do
          A.set z i (A.get x i +. b)
        done

      let maxnorm x =
        let max = ref 0.0 in
        for i = 0 to A.length x - 1 do
          let ax = abs_float (A.get x i) in
          if ax > !max then max := ax
        done;
        !max

      let wrmsnorm x w =
        let a = ref 0.0 in
        let lx = A.length x in
        for i = 0 to lx - 1 do
          a := !a +. ((A.get x i) *. (A.get w i) *. (A.get x i) *. (A.get w i))
        done;
        sqrt (!a /. float lx)

      let wrmsnormmask x w id =
        let a = ref 0.0 in
        let lx = A.length x in
        for i = 0 to lx - 1 do
          if A.get id i > 0.0 then
            a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
        done;
        sqrt (!a /. float lx)

      let min x =
        let min = ref max_float in
        for i = 0 to A.length x - 1 do
          let xv = A.get x i in
          if xv < !min then min := xv
        done;
        !min

      let dotprod x y =
        let a = ref 0.0 in
        for i = 0 to A.length x - 1 do
          a := !a +. (A.get x i *. A.get y i)
        done;
        !a

      let compare c x z =
        for i = 0 to A.length x - 1 do
          A.set z i (if abs_float (A.get x i) >= c then 1.0 else 0.0)
        done

      let invtest x z =
        let r = ref true in
        for i = 0 to A.length x - 1 do
          if A.get x i = 0.0 then r := false
          else A.set z i (1.0 /. (A.get x i))
        done;
        !r

      let wl2norm x w =
        let a = ref 0.0 in
        for i = 0 to A.length x - 1 do
          a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
        done;
        sqrt !a

      let l1norm x =
        let a = ref 0.0 in
        for i = 0 to A.length x - 1 do
          a := !a +. abs_float (A.get x i)
        done;
        !a

      let constrmask c x m =
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

      let minquotient n d =
        let m = ref Config.big_real in
        for i = 0 to A.length n - 1 do
          if (A.get d i) <> 0.0 then
            m := floatmin !m (A.get n i /. A.get d i)
        done;
        !m

      let prod x y z =
        for i = 0 to A.length x - 1 do
          A.set z i (A.get x i *. A.get y i)
        done

      let div x y z =
        for i = 0 to A.length x - 1 do
          A.set z i (A.get x i /. A.get y i)
        done

      let abs x z =
        for i = 0 to A.length x - 1 do
          A.set z i (abs_float (A.get x i))
        done

      let inv x z =
        for i = 0 to A.length x - 1 do
          A.set z i (1.0 /. (A.get x i))
        done

      let space (x : float array) = (A.length x, 1)

      let getlength (x : float array) = A.length x

      let print ?(logfile=Logfile.stdout) (x : float array) =
        for i = 0 to A.length x - 1 do
          Logfile.output_string logfile (Printf.sprintf "%19.16g" (A.get x i));
          Logfile.output_string logfile "\n"
        done;
        Logfile.output_string logfile "\n"

      (* fused and array operations *)

      let linearcombination (ca : RealArray.t)
                               (xa : (float array) array)
                               (z  : float array) =
        let nvec = Array.length xa in
        if nvec = 1 then scale ca.{0} xa.(0) z
        else if nvec = 2 then linearsum ca.{0} xa.(0) ca.{1} xa.(1) z
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

      let scaleaddmulti (aa : RealArray.t)
                           (x : float array)
                           (ya : (float array) array)
                           (za : (float array) array) =
        let nvec = Array.length ya in
        if nvec = 1 then linearsum aa.{0} x 1.0 ya.(0) za.(0)
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

      let dotprodmulti (x : float array)
                          (ya : (float array) array)
                          (dp : RealArray.t) =
        let nvec = Array.length ya in
        if nvec = 1 then dp.{0} <- dotprod x ya.(0)
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

      let linearsumvectorarray a (xa : (float array) array)
                                  b (ya : (float array) array)
                                    (za : (float array) array) =
        let nvec = Array.length ya in
        if nvec = 1 then linearsum a xa.(0) b ya.(0) za.(0)
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

      let scalevectorarray (c  : RealArray.t)
                              (xa : (float array) array)
                              (za : (float array) array) =
        let nvec = Array.length xa in
        if nvec = 1 then scale c.{0} xa.(0) za.(0)
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

      let constvectorarray c (za : float array array) =
        let nvec = Array.length za in
        if nvec = 1 then const c za.(0)
        else
          let n = A.length za.(0) in
          for i = 0 to nvec - 1 do
            let z = za.(i) in
            for j = 0 to n - 1 do
              A.set z j c
            done
          done

      let wrmsnormvectorarray (xa  : (float array) array)
                                 (wa  : (float array) array)
                                 (nrm : RealArray.t) =
        let nvec = Array.length xa in
        if nvec = 1 then nrm.{0} <- wrmsnorm xa.(0) wa.(0)
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

      let wrmsnormmaskvectorarray (xa  : (float array) array)
                                     (wa  : (float array) array)
                                     (id  : float array)
                                     (nrm : RealArray.t) =
        let nvec = Array.length xa in
        if nvec = 1 then nrm.{0} <- wrmsnormmask xa.(0) wa.(0) id
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

      let scaleaddmultivectorarray (ra  : RealArray.t)
                                      (xa  : (float array) array)
                                      (yaa : (float array) array array)
                                      (zaa : (float array) array array) =
        let nsum = Array.length yaa in
        let nvec = Array.length yaa.(0) in
        if nvec = 1 then begin
          if nsum = 1 then linearsum ra.{0} xa.(0) 1.0 yaa.(0).(0) zaa.(0).(0)
          else
            let yya = Array.init nsum (fun j -> yaa.(j).(0)) in
            let zza = Array.init nsum (fun j -> zaa.(j).(0)) in
            scaleaddmulti ra xa.(0) yya zza
        end
        else if nsum = 1 then linearsumvectorarray ra.{0} xa 1.0 yaa.(0) zaa.(0)
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

      let linearcombinationvectorarray (ca  : RealArray.t)
                                          (xaa : (float array) array array)
                                          (za  : (float array) array) =
        let nsum = Array.length xaa in
        let nvec = Array.length xaa.(0) in
        if nvec = 1 then begin
          if nsum = 1 then scale ca.{0} xaa.(0).(0) za.(0)
          else if nsum = 2
               then linearsum ca.{0} xaa.(0).(0) ca.{1} xaa.(1).(0) za.(0)
          else
            let ya = Array.init nsum (fun i -> xaa.(i).(0)) in
            linearcombination ca ya za.(0)
        end
        else
          if nsum = 1 then
            let ctmp = RealArray.make nvec ca.{0} in
            scalevectorarray ctmp xaa.(0) za
          else if nsum = 2 then
            linearsumvectorarray ca.{0} xaa.(0) ca.{1} xaa.(1) za
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
        let dotprod     = dotprod
        let maxnorm     = maxnorm
        let min         = min
        let l1norm      = l1norm
        let invtest     = invtest
        let constrmask  = constrmask
        let minquotient = minquotient

        let wsqrsum x w =
          let a = ref 0.0 in
          let lx = A.length x in
          for i = 0 to lx - 1 do
            a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
          done;
          !a

        let wsqrsummask x w id =
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
          Nvector_custom.check        = checkfn;
          Nvector_custom.clone        = DataOps.clone;
          Nvector_custom.space        = Some DataOps.space;
          Nvector_custom.getlength    = DataOps.getlength;
          Nvector_custom.print        = Some (fun x logfile -> DataOps.print ?logfile x);
          Nvector_custom.linearsum    = DataOps.linearsum;
          Nvector_custom.const        = DataOps.const;
          Nvector_custom.prod         = DataOps.prod;
          Nvector_custom.div          = DataOps.div;
          Nvector_custom.scale        = DataOps.scale;
          Nvector_custom.abs          = DataOps.abs;
          Nvector_custom.inv          = DataOps.inv;
          Nvector_custom.addconst     = DataOps.addconst;
          Nvector_custom.maxnorm      = DataOps.maxnorm;
          Nvector_custom.wrmsnorm     = DataOps.wrmsnorm;
          Nvector_custom.min          = DataOps.min;
          Nvector_custom.dotprod      = DataOps.dotprod;
          Nvector_custom.compare      = DataOps.compare;
          Nvector_custom.invtest      = DataOps.invtest;

          Nvector_custom.wl2norm      = Some DataOps.wl2norm;
          Nvector_custom.l1norm       = Some DataOps.l1norm;
          Nvector_custom.wrmsnormmask = Some DataOps.wrmsnormmask;
          Nvector_custom.constrmask   = Some DataOps.constrmask;
          Nvector_custom.minquotient  = Some DataOps.minquotient;

          Nvector_custom.getcommunicator = None;

          Nvector_custom.linearcombination
            = Some DataOps.linearcombination;
          Nvector_custom.scaleaddmulti
            = Some DataOps.scaleaddmulti;
          Nvector_custom.dotprodmulti
            = Some DataOps.dotprodmulti;
          Nvector_custom.linearsumvectorarray
            = Some DataOps.linearsumvectorarray;
          Nvector_custom.scalevectorarray
            = Some DataOps.scalevectorarray;
          Nvector_custom.constvectorarray
            = Some DataOps.constvectorarray;
          Nvector_custom.wrmsnormvectorarray
            = Some DataOps.wrmsnormvectorarray;
          Nvector_custom.wrmsnormmaskvectorarray
            = Some DataOps.wrmsnormmaskvectorarray;
          Nvector_custom.scaleaddmultivectorarray
            = Some DataOps.scaleaddmultivectorarray;
          Nvector_custom.linearcombinationvectorarray
            = Some DataOps.linearcombinationvectorarray;

          Nvector_custom.dotprod_local
            = Some DataOps.Local.dotprod;

          Nvector_custom.maxnorm_local
            = Some DataOps.Local.maxnorm;

          Nvector_custom.min_local
            = Some DataOps.Local.min;

          Nvector_custom.l1norm_local
            = Some DataOps.Local.l1norm;

          Nvector_custom.invtest_local
            = Some DataOps.Local.invtest;

          Nvector_custom.constrmask_local
            = Some DataOps.Local.constrmask;

          Nvector_custom.minquotient_local
            = Some DataOps.Local.minquotient;

          Nvector_custom.wsqrsum_local
            = Some DataOps.Local.wsqrsum;

          Nvector_custom.wsqrsummask_local
            = Some DataOps.Local.wsqrsummask;
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
            Nvector_custom.linearcombination
              = Some AnyOps.linearcombination;
            Nvector_custom.scaleaddmulti
              = Some AnyOps.scaleaddmulti;
            Nvector_custom.dotprodmulti
              = Some AnyOps.dotprodmulti;
            Nvector_custom.linearsumvectorarray
              = Some AnyOps.linearsumvectorarray;
            Nvector_custom.scalevectorarray
              = Some AnyOps.scalevectorarray;
            Nvector_custom.constvectorarray
              = Some AnyOps.constvectorarray;
            Nvector_custom.wrmsnormvectorarray
              = Some AnyOps.wrmsnormvectorarray;
            Nvector_custom.wrmsnormmaskvectorarray
              = Some AnyOps.wrmsnormmaskvectorarray;
            Nvector_custom.scaleaddmultivectorarray
              = Some AnyOps.scaleaddmultivectorarray;
            Nvector_custom.linearcombinationvectorarray
              = Some AnyOps.linearcombinationvectorarray;
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
        let inject _ = assert false (* clone is overridden *)
      end)

      (* Cloning in this way ensures that the enabled status of fused and
         array operations is the same in the clone nvector. *)
      let clone x = Nvector.clone x
    end

  end (* }}} *)

include Array

