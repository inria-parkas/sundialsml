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

module Make =
  functor (A : ArrayOps) ->
  struct (* {{{ *)
    type data = A.data
    type kind = Nvector_custom.kind
    type t = data Nvector_custom.t

    let checkfn v1 =
      let l = A.length v1 in
      (fun v2 -> l = A.length v2)

    module DataOps = struct (* {{{ *)
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
          else if A.get x i = 0.0 then f false (i + 1)
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

      let n_vspace x = (A.length x, 1)

      let n_vgetlength x = A.length x

      (* fused and array operations *)

      let n_vlinearcombination (ca : RealArray.t) (xa : A.data array) (z : A.data) =
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

      let n_vscaleaddmulti (aa : RealArray.t) (x : A.data) (ya : A.data array)
                                                           (za : A.data array) =
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

      let n_vdotprodmulti (x : A.data) (ya : A.data array) (dp : RealArray.t) =
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

      let arr_vaxpy_array a (xa : A.data array) (ya : A.data array) =
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

      let v_sumvectorarray (xa : A.data array) (ya : A.data array)
                                               (za : A.data array) =
        let nvec = Array.length xa in
        let n = A.length xa.(0) in
        for i = 0 to nvec - 1 do
          let x, y, z = xa.(i), ya.(i), za.(i) in
          for j = 0 to n - 1 do
            A.set z j (A.get x j +. A.get y j)
          done
        done

      let v_diffvectorarray (xa : A.data array) (ya : A.data array)
                                                (za : A.data array) =
        let nvec = Array.length xa in
        let n = A.length xa.(0) in
        for i = 0 to nvec - 1 do
          let x, y, z = xa.(i), ya.(i), za.(i) in
          for j = 0 to n - 1 do
            A.set z j (A.get x j -. A.get y j)
          done
        done

      let v_lin1vectorarray a (xa : A.data array) (ya : A.data array)
                                                  (za : A.data array) =
        let nvec = Array.length xa in
        let n = A.length xa.(0) in
        for i = 0 to nvec - 1 do
          let x, y, z = xa.(i), ya.(i), za.(i) in
          for j = 0 to n - 1 do
            A.set z j (a *. A.get x j +. A.get y j)
          done
        done

      let v_lin2vectorarray a (xa : A.data array) (ya : A.data array)
                                                  (za : A.data array) =
        let nvec = Array.length xa in
        let n = A.length xa.(0) in
        for i = 0 to nvec - 1 do
          let x, y, z = xa.(i), ya.(i), za.(i) in
          for j = 0 to n - 1 do
            A.set z j (a *. A.get x j -. A.get y j)
          done
        done

      let v_scalesumvectorarray c (xa : A.data array) (ya : A.data array)
                                                      (za : A.data array) =
        let nvec = Array.length xa in
        let n = A.length xa.(0) in
        for i = 0 to nvec - 1 do
          let x, y, z = xa.(i), ya.(i), za.(i) in
          for j = 0 to n - 1 do
            A.set z j (c *. (A.get x j +. A.get y j))
          done
        done

      let v_scalediffvectorarray c (xa : A.data array) (ya : A.data array)
                                                       (za : A.data array) =
        let nvec = Array.length xa in
        let n = A.length xa.(0) in
        for i = 0 to nvec - 1 do
          let x, y, z = xa.(i), ya.(i), za.(i) in
          for j = 0 to n - 1 do
            A.set z j (c *. (A.get x j -. A.get y j))
          done
        done

      let n_vlinearsumvectorarray a (xa : A.data array) b (ya : A.data array)
                                                          (za : A.data array) =
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

      let n_vscalevectorarray (c : RealArray.t) (xa : A.data array)
                                                (za : A.data array) =
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

      let n_vconstvectorarray c (za : A.data array) =
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

      let n_vwrmsnormvectorarray (xa : A.data array) (wa : A.data array)
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

      let n_vwrmsnormmaskvectorarray (xa : A.data array) (wa : A.data array)
                                     (id : A.data) (nrm : RealArray.t) =
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

      let n_vscaleaddmultivectorarray (ra : RealArray.t) (xa : A.data array)
                                      (yaa : A.data array array)
                                      (zaa : A.data array array) =
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

      let n_vlinearcombinationvectorarray (ca : RealArray.t)
                                          (xaa : A.data array array)
                                          (za : A.data array) =
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

          Nvector_custom.n_vlinearcombination = Some DataOps.n_vlinearcombination;
          Nvector_custom.n_vscaleaddmulti = Some DataOps.n_vscaleaddmulti;
          Nvector_custom.n_vdotprodmulti = Some DataOps.n_vdotprodmulti;

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

          Nvector_custom.n_vdotprod_local = Some DataOps.Local.n_vdotprod;
          Nvector_custom.n_vmaxnorm_local = Some DataOps.Local.n_vmaxnorm;
          Nvector_custom.n_vmin_local = Some DataOps.Local.n_vmin;
          Nvector_custom.n_vl1norm_local = Some DataOps.Local.n_vl1norm;
          Nvector_custom.n_vinvtest_local = Some DataOps.Local.n_vinvtest;
          Nvector_custom.n_vconstrmask_local = Some DataOps.Local.n_vconstrmask;
          Nvector_custom.n_vminquotient_local = Some DataOps.Local.n_vminquotient;
          Nvector_custom.n_vwsqrsum_local = Some DataOps.Local.n_vwsqrsum;
          Nvector_custom.n_vwsqrsummask_local = Some DataOps.Local.n_vwsqrsummask;
    } (* }}} *)

    let make n e = Nvector_custom.make_wrap array_nvec_ops (A.make n e)

    let wrap = Nvector_custom.make_wrap array_nvec_ops
      (* (Nvector.Mutable.add_tracing "::" array_nvec_ops) *)

    let unwrap = Nvector.unwrap

    let enable = Nvector_custom.enable

    module Any = struct (* {{{ *)

      type Nvector.gdata += Arr of data

      let wrap =
        let inject x = Arr x in
        let project = function Arr x -> x | _ -> raise Nvector.BadGenericType in
        let ops = Nvector_custom.Any.convert_ops ~inject ~project array_nvec_ops in
        (* TODO: override array ops? *)
        Nvector_custom.Any.make_wrap ops ~inject

    end (* }}} *)

    module Ops = struct (* {{{ *)
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
      let n_vspace x = DataOps.n_vspace (unwrap x)
      let n_vgetlength x = DataOps.n_vgetlength (unwrap x)

      let n_vlinearcombination c x z
        = DataOps.n_vlinearcombination c (Array.map unwrap x) (unwrap z)

      let n_vscaleaddmulti c x y z
        = DataOps.n_vscaleaddmulti c (unwrap x)
                                     (Array.map unwrap y)
                                     (Array.map unwrap z)

      let n_vdotprodmulti x y d
        = DataOps.n_vdotprodmulti (unwrap x) (Array.map unwrap y) d

      let n_vlinearsumvectorarray a x b y z
        = DataOps.n_vlinearsumvectorarray a (Array.map unwrap x)
                                          b (Array.map unwrap y)
                                          (Array.map unwrap z)

      let n_vscalevectorarray c x z
        = DataOps.n_vscalevectorarray c (Array.map unwrap x)
                                        (Array.map unwrap z)

      let n_vconstvectorarray c x
        = DataOps.n_vconstvectorarray c (Array.map unwrap x)

      let n_vwrmsnormvectorarray x w m
        = DataOps.n_vwrmsnormvectorarray (Array.map unwrap x)
                                         (Array.map unwrap w) m

      let n_vwrmsnormmaskvectorarray x w id m
        = DataOps.n_vwrmsnormmaskvectorarray (Array.map unwrap x)
                                             (Array.map unwrap w)
                                             (unwrap id) m

      let n_vscaleaddmultivectorarray a x yy zz
        = DataOps.n_vscaleaddmultivectorarray a (Array.map unwrap x)
                                              (Array.map (Array.map unwrap) yy)
                                              (Array.map (Array.map unwrap) zz)

      let n_vlinearcombinationvectorarray c xx z
        = DataOps.n_vlinearcombinationvectorarray c
                                            (Array.map (Array.map unwrap) xx)
                                            (Array.map unwrap z)

      module Local = struct
        let n_vdotprod     = n_vdotprod
        let n_vmaxnorm     = n_vmaxnorm
        let n_vmin         = n_vmin
        let n_vl1norm      = n_vl1norm
        let n_vinvtest     = n_vinvtest
        let n_vconstrmask  = n_vconstrmask
        let n_vminquotient = n_vminquotient

        let n_vwsqrsum x w = DataOps.Local.n_vwsqrsum (unwrap x) (unwrap w)

        let n_vwsqrsummask x w id =
          DataOps.Local.n_vwsqrsummask (unwrap x) (unwrap w) (unwrap id)
      end
    end (* }}} *)
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

      let wrap =
        let inject x = Arr x in
        let project = function Arr x -> x | _ -> raise Nvector.BadGenericType in
        let ops = Nvector_custom.Any.convert_ops ~inject ~project array_nvec_ops in
        (* TODO: override array ops? *)
        Nvector_custom.Any.make_wrap ops ~inject

    end (* }}} *)

    module Ops = struct (* {{{ *)
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
      let n_vspace x = DataOps.n_vspace (unwrap x)
      let n_vgetlength x = DataOps.n_vgetlength (unwrap x)

      let n_vlinearcombination c x z
        = DataOps.n_vlinearcombination c (Array.map unwrap x) (unwrap z)

      let n_vscaleaddmulti c x y z
        = DataOps.n_vscaleaddmulti c (unwrap x)
                                     (Array.map unwrap y)
                                     (Array.map unwrap z)

      let n_vdotprodmulti x y d
        = DataOps.n_vdotprodmulti (unwrap x) (Array.map unwrap y) d

      let n_vlinearsumvectorarray a x b y z
        = DataOps.n_vlinearsumvectorarray a (Array.map unwrap x)
                                          b (Array.map unwrap y)
                                          (Array.map unwrap z)

      let n_vscalevectorarray c x z
        = DataOps.n_vscalevectorarray c (Array.map unwrap x)
                                        (Array.map unwrap z)

      let n_vconstvectorarray c x
        = DataOps.n_vconstvectorarray c (Array.map unwrap x)

      let n_vwrmsnormvectorarray x w m
        = DataOps.n_vwrmsnormvectorarray (Array.map unwrap x)
                                         (Array.map unwrap w) m

      let n_vwrmsnormmaskvectorarray x w id m
        = DataOps.n_vwrmsnormmaskvectorarray (Array.map unwrap x)
                                             (Array.map unwrap w)
                                             (unwrap id) m

      let n_vscaleaddmultivectorarray a x yy zz
        = DataOps.n_vscaleaddmultivectorarray a (Array.map unwrap x)
                                              (Array.map (Array.map unwrap) yy)
                                              (Array.map (Array.map unwrap) zz)

      let n_vlinearcombinationvectorarray c xx z
        = DataOps.n_vlinearcombinationvectorarray c
                                            (Array.map (Array.map unwrap) xx)
                                            (Array.map unwrap z)

      module Local = struct
        let n_vdotprod     = n_vdotprod
        let n_vmaxnorm     = n_vmaxnorm
        let n_vmin         = n_vmin
        let n_vl1norm      = n_vl1norm
        let n_vinvtest     = n_vinvtest
        let n_vconstrmask  = n_vconstrmask
        let n_vminquotient = n_vminquotient

        let n_vwsqrsum x w = DataOps.Local.n_vwsqrsum (unwrap x) (unwrap w)

        let n_vwsqrsummask x w id =
          DataOps.Local.n_vwsqrsummask (unwrap x) (unwrap w) (unwrap id)
      end
    end (* }}} *)
  end (* }}} *)

include Array

