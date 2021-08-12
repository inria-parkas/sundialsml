
open Sundials

type data = RealArray.t
type kind = [`Serial]
type t = (data, kind) Nvector.t
type 'k any = (data, [>kind] as 'k) Nvector.t

(* Selectively enable and disable fused and array operations *)
external c_enablefusedops_serial                     : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_ser_enablefusedops"
external c_enablelinearcombination_serial            : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_ser_enablelinearcombination"
external c_enablescaleaddmulti_serial                : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_ser_enablescaleaddmulti"
external c_enabledotprodmulti_serial                 : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_ser_enabledotprodmulti"
external c_enablelinearsumvectorarray_serial         : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_ser_enablelinearsumvectorarray"
external c_enablescalevectorarray_serial             : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_ser_enablescalevectorarray"
external c_enableconstvectorarray_serial             : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_ser_enableconstvectorarray"
external c_enablewrmsnormvectorarray_serial          : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_ser_enablewrmsnormvectorarray"
external c_enablewrmsnormmaskvectorarray_serial      : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_ser_enablewrmsnormmaskvectorarray"
external c_enablescaleaddmultivectorarray_serial     : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_ser_enablescaleaddmultivectorarray"
external c_enablelinearcombinationvectorarray_serial : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_ser_enablelinearcombinationvectorarray"

let unwrap = Nvector.unwrap

external c_wrap : RealArray.t -> (t -> bool) -> (t -> t) -> t
  = "sunml_nvec_wrap_serial"

let rec wrap ?(with_fused_ops=false) v =
  let len = RealArray.length v in
  let nv = c_wrap v (fun nv' -> len = RealArray.length (unwrap nv')) clone in
  if with_fused_ops then c_enablefusedops_serial nv true;
  nv

and clone nv =
  let nv' = wrap (RealArray.copy (unwrap nv)) in
  c_enablelinearcombination_serial nv'
    (Nvector.has_n_vlinearcombination nv);
  c_enablescaleaddmulti_serial nv'
    (Nvector.has_n_vscaleaddmulti nv);
  c_enabledotprodmulti_serial nv'
    (Nvector.has_n_vdotprodmulti nv);
  c_enablelinearsumvectorarray_serial nv'
    (Nvector.has_n_vlinearsumvectorarray nv);
  c_enablescalevectorarray_serial nv'
    (Nvector.has_n_vscalevectorarray nv);
  c_enableconstvectorarray_serial nv'
    (Nvector.has_n_vconstvectorarray nv);
  c_enablewrmsnormvectorarray_serial nv'
    (Nvector.has_n_vwrmsnormvectorarray nv);
  c_enablewrmsnormmaskvectorarray_serial nv'
    (Nvector.has_n_vwrmsnormmaskvectorarray nv);
  c_enablescaleaddmultivectorarray_serial nv'
    (Nvector.has_n_vscaleaddmultivectorarray nv);
  c_enablelinearcombinationvectorarray_serial nv'
    (Nvector.has_n_vlinearcombinationvectorarray nv);
  nv'

let make ?with_fused_ops n iv = wrap ?with_fused_ops (RealArray.make n iv)

let pp fmt v = RealArray.pp fmt (unwrap v)

let do_enable f nv v =
  match v with
  | None -> ()
  | Some v -> f nv v

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
  = do_enable c_enablefusedops_serial nv
              with_fused_ops;
    do_enable c_enablelinearcombination_serial nv
              with_linear_combination;
    do_enable c_enablescaleaddmulti_serial nv
              with_scale_add_multi;
    do_enable c_enabledotprodmulti_serial nv
              with_dot_prod_multi;
    do_enable c_enablelinearsumvectorarray_serial nv
              with_linear_sum_vector_array;
    do_enable c_enablescalevectorarray_serial nv
              with_scale_vector_array;
    do_enable c_enableconstvectorarray_serial nv
              with_const_vector_array;
    do_enable c_enablewrmsnormvectorarray_serial nv
              with_wrms_norm_vector_array;
    do_enable c_enablewrmsnormmaskvectorarray_serial nv
              with_wrms_norm_mask_vector_array;
    do_enable c_enablescaleaddmultivectorarray_serial nv
              with_scale_add_multi_vector_array;
    do_enable c_enablelinearcombinationvectorarray_serial nv
              with_linear_combination_vector_array

module Any = struct (* {{{ *)

  external c_any_wrap
    : extension_constructor
      -> RealArray.t
      -> (Nvector.any -> bool)
      -> (Nvector.any -> Nvector.any)
      -> Nvector.any
    = "sunml_nvec_anywrap_serial"

  let rec wrap
      ?(with_fused_ops=false)
      ?(with_linear_combination=false)
      ?(with_scale_add_multi=false)
      ?(with_dot_prod_multi=false)
      ?(with_linear_sum_vector_array=false)
      ?(with_scale_vector_array=false)
      ?(with_const_vector_array=false)
      ?(with_wrms_norm_vector_array=false)
      ?(with_wrms_norm_mask_vector_array=false)
      ?(with_scale_add_multi_vector_array=false)
      ?(with_linear_combination_vector_array=false)
      v
    =
      if not Sundials_impl.Versions.has_nvector_get_id
        then raise Config.NotImplementedBySundialsVersion;
      let len = RealArray.length v in
      let check nv' =
        match unwrap nv' with
        | Nvector.RA ra ->
            len = RealArray.length ra && Nvector.get_id nv' = Nvector.Serial
        | _ -> false
      in
      let nv = c_any_wrap [%extension_constructor Nvector.RA] v check clone in
      if with_fused_ops
        then c_enablefusedops_serial nv true;
      if with_fused_ops
        then c_enablefusedops_serial nv true;
      if with_linear_combination
        then c_enablelinearcombination_serial nv true;
      if with_scale_add_multi
        then c_enablescaleaddmulti_serial nv true;
      if with_dot_prod_multi
        then c_enabledotprodmulti_serial nv true;
      if with_linear_sum_vector_array
        then c_enablelinearsumvectorarray_serial nv true;
      if with_scale_vector_array
        then c_enablescalevectorarray_serial nv true;
      if with_const_vector_array
        then c_enableconstvectorarray_serial nv true;
      if with_wrms_norm_vector_array
        then c_enablewrmsnormvectorarray_serial nv true;
      if with_wrms_norm_mask_vector_array
        then c_enablewrmsnormmaskvectorarray_serial nv true;
      if with_scale_add_multi_vector_array
        then c_enablescaleaddmultivectorarray_serial nv true;
      if with_linear_combination_vector_array
        then c_enablelinearcombinationvectorarray_serial nv true;
      nv

  and clone nv =
    let v = match unwrap nv with
            | Nvector.RA v -> v
            | _ -> assert false
    in
    let nv' = wrap (RealArray.copy v) in
    c_enablelinearcombination_serial nv'
      (Nvector.has_n_vlinearcombination nv);
    c_enablescaleaddmulti_serial nv'
      (Nvector.has_n_vscaleaddmulti nv);
    c_enabledotprodmulti_serial nv'
      (Nvector.has_n_vdotprodmulti nv);
    c_enablelinearsumvectorarray_serial nv'
      (Nvector.has_n_vlinearsumvectorarray nv);
    c_enablescalevectorarray_serial nv'
      (Nvector.has_n_vscalevectorarray nv);
    c_enableconstvectorarray_serial nv'
      (Nvector.has_n_vconstvectorarray nv);
    c_enablewrmsnormvectorarray_serial nv'
      (Nvector.has_n_vwrmsnormvectorarray nv);
    c_enablewrmsnormmaskvectorarray_serial nv'
      (Nvector.has_n_vwrmsnormmaskvectorarray nv);
    c_enablescaleaddmultivectorarray_serial nv'
      (Nvector.has_n_vscaleaddmultivectorarray nv);
    c_enablelinearcombinationvectorarray_serial nv'
      (Nvector.has_n_vlinearcombinationvectorarray nv);
    nv'

  let make
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
      n iv
    = wrap ?with_fused_ops
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
           (RealArray.make n iv)

end (* }}} *)

module Ops : Nvector.NVECTOR_OPS with type t = t =
struct (* {{{ *)
  type t = (data, kind) Nvector.t
  let check = Nvector.check

  let n_vclone = Nvector.clone

  external c_n_vlinearsum    : float -> t -> float -> t -> t -> unit
    = "sunml_nvec_ser_n_vlinearsum" [@@noalloc]

  let n_vlinearsum a (x : t) b (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_n_vlinearsum a x b y z

  external n_vconst          : float -> t -> unit
    = "sunml_nvec_ser_n_vconst" [@@noalloc]

  external c_n_vprod         : t -> t -> t -> unit
    = "sunml_nvec_ser_n_vprod" [@@noalloc]

  let n_vprod (x : t) (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_n_vprod x y z

  external c_n_vdiv          : t -> t -> t -> unit
    = "sunml_nvec_ser_n_vdiv" [@@noalloc]

  let n_vdiv (x : t) (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_n_vdiv x y z

  external c_n_vscale        : float -> t -> t -> unit
    = "sunml_nvec_ser_n_vscale" [@@noalloc]

  let n_vscale c (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_n_vscale c x z

  external c_n_vabs          : t -> t -> unit
    = "sunml_nvec_ser_n_vabs" [@@noalloc]

  let n_vabs (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_n_vabs x z

  external c_n_vinv          : t -> t -> unit
    = "sunml_nvec_ser_n_vinv" [@@noalloc]

  let n_vinv (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_n_vinv x z

  external c_n_vaddconst     : t -> float -> t -> unit
    = "sunml_nvec_ser_n_vaddconst" [@@noalloc]

  let n_vaddconst (x : t) b (z : t) =
    if Sundials_configuration.safe then check x z;
    c_n_vaddconst x b z

  external c_n_vdotprod      : t -> t -> float
    = "sunml_nvec_ser_n_vdotprod"

  let n_vdotprod (x : t) (y : t) =
    if Sundials_configuration.safe then check x y;
    c_n_vdotprod x y

  external n_vmaxnorm        : t -> float
    = "sunml_nvec_ser_n_vmaxnorm"

  external c_n_vwrmsnorm     : t -> t -> float
    = "sunml_nvec_ser_n_vwrmsnorm"

  let n_vwrmsnorm (x : t) (w : t) =
    if Sundials_configuration.safe then check x w;
    c_n_vwrmsnorm x w

  external c_n_vwrmsnormmask : t -> t -> t -> float
    = "sunml_nvec_ser_n_vwrmsnormmask"

  let n_vwrmsnormmask (x : t) (w : t) (id : t) =
    if Sundials_configuration.safe then (check x w; check x id);
    c_n_vwrmsnormmask x w id

  external n_vmin            : t -> float
    = "sunml_nvec_ser_n_vmin"

  external c_n_vwl2norm      : t -> t -> float
    = "sunml_nvec_ser_n_vwl2norm"

  let n_vwl2norm (x : t) (w : t) =
    if Sundials_configuration.safe then check x w;
    c_n_vwl2norm x w

  external n_vl1norm         : t -> float
    = "sunml_nvec_ser_n_vl1norm"

  external c_n_vcompare      : float -> t -> t -> unit
    = "sunml_nvec_ser_n_vcompare" [@@noalloc]

  let n_vcompare c (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_n_vcompare c x z

  external c_n_vinvtest      : t -> t -> bool
    = "sunml_nvec_ser_n_vinvtest" [@@noalloc]

  let n_vinvtest (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_n_vinvtest x z

  external c_n_vconstrmask   : t -> t -> t -> bool
    = "sunml_nvec_ser_n_vconstrmask" [@@noalloc]

  let n_vconstrmask (c : t) (x : t) (m : t) =
    if Sundials_configuration.safe then (check c x; check c m);
    c_n_vconstrmask c x m

  external c_n_vminquotient  : t -> t -> float
    = "sunml_nvec_ser_n_vminquotient"

  let n_vminquotient (n : t) (d : t) =
    if Sundials_configuration.safe then check n d;
    c_n_vminquotient n d

  external n_vspace          : t -> int * int
    = "sunml_nvec_ser_n_vspace" [@@noalloc]

  external n_vgetlength      : t -> int
    = "sunml_nvec_ser_n_vgetlength"

  external c_n_vlinearcombination : RealArray.t -> t array -> t -> unit
    = "sunml_nvec_ser_n_vlinearcombination"

  let n_vlinearcombination ca (xa : t array) (z : t) =
    if Sundials_impl.Versions.sundials_lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe then Array.iter (check z) xa;
    c_n_vlinearcombination ca xa z

  external c_n_vscaleaddmulti : RealArray.t -> t -> t array -> t array -> unit
    = "sunml_nvec_ser_n_vscaleaddmulti"

  let n_vscaleaddmulti aa (x : t) (ya : t array) (za : t array) =
    if Sundials_impl.Versions.sundials_lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe then
      (Array.iter (check x) ya; Array.iter (check x) za);
    c_n_vscaleaddmulti aa x ya za

  external c_n_vdotprodmulti  : t -> t array -> RealArray.t -> unit
    = "sunml_nvec_ser_n_vdotprodmulti"

  let n_vdotprodmulti (x : t) (ya : t array) dp =
    if Sundials_impl.Versions.sundials_lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe then Array.iter (check x) ya;
    c_n_vdotprodmulti x ya dp

  external c_n_vlinearsumvectorarray
    : float -> t array -> float -> t array -> t array -> unit
    = "sunml_nvec_ser_n_vlinearsumvectorarray"

  let n_vlinearsumvectorarray a (xa : t array) b (ya : t array) (za : t array) =
    if Sundials_impl.Versions.sundials_lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (check x) ya;
          Array.iter (check x) za);
    c_n_vlinearsumvectorarray a xa b ya za

  external c_n_vscalevectorarray
    : RealArray.t -> t array -> t array -> unit
    = "sunml_nvec_ser_n_vscalevectorarray"

  let n_vscalevectorarray c (xa : t array) (za : t array) =
    if Sundials_impl.Versions.sundials_lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (check x) za);
    c_n_vscalevectorarray c xa za

  external c_n_vconstvectorarray
    : float -> t array -> unit
    = "sunml_nvec_ser_n_vconstvectorarray"

  let n_vconstvectorarray c (za : t array) =
    if Sundials_impl.Versions.sundials_lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe
    then (let z = Array.get za 0 in
          Array.iter (check z) za);
    c_n_vconstvectorarray c za

  external c_n_vwrmsnormvectorarray
    : t array -> t array -> RealArray.t -> unit
    = "sunml_nvec_ser_n_vwrmsnormvectorarray"

  let n_vwrmsnormvectorarray (xa : t array) (wa : t array) nrm =
    if Sundials_impl.Versions.sundials_lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (check x) wa);
    c_n_vwrmsnormvectorarray xa wa nrm

  external c_n_vwrmsnormmaskvectorarray
    : t array -> t array -> t -> RealArray.t -> unit
    = "sunml_nvec_ser_n_vwrmsnormmaskvectorarray"

  let n_vwrmsnormmaskvectorarray (xa : t array) (wa : t array) (id : t) nrm =
    if Sundials_impl.Versions.sundials_lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe
    then (Array.iter (check id) xa;
          Array.iter (check id) wa);
    c_n_vwrmsnormmaskvectorarray xa wa id nrm

  external c_n_vscaleaddmultivectorarray
    : RealArray.t -> t array -> t array array -> t array array -> unit
    = "sunml_nvec_ser_n_vscaleaddmultivectorarray"

  let n_vscaleaddmultivectorarray ra (xa : t array) (yaa : t array array)
                                     (zaa : t array array) =
    if Sundials_impl.Versions.sundials_lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (Array.iter (check x)) yaa;
          Array.iter (Array.iter (check x)) zaa);
    c_n_vscaleaddmultivectorarray ra xa yaa zaa

  external c_n_vlinearcombinationvectorarray
    : RealArray.t -> t array array -> t array -> unit
    = "sunml_nvec_ser_n_vlinearcombinationvectorarray"

  let n_vlinearcombinationvectorarray ca (xaa : t array array) (za : t array) =
    if Sundials_impl.Versions.sundials_lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe
    then (let z = Array.get za 0 in
          Array.iter (check z) za;
          Array.iter (Array.iter (check z)) xaa);
    c_n_vlinearcombinationvectorarray ca xaa za

  module Local = struct
    let n_vdotprod     = n_vdotprod
    let n_vmaxnorm     = n_vmaxnorm
    let n_vmin         = n_vmin
    let n_vl1norm      = n_vl1norm
    let n_vinvtest     = n_vinvtest
    let n_vconstrmask  = n_vconstrmask
    let n_vminquotient = n_vminquotient

    external c_n_vwsqrsum
      : t -> t -> float
      = "sunml_nvec_ser_n_vwsqrsumlocal"

    let n_vwsqrsum (x : t) (w : t) =
      if Sundials_configuration.safe then check x w;
      c_n_vwsqrsum x w

    external c_n_vwsqrsummask
      : t -> t -> t -> float
      = "sunml_nvec_ser_n_vwsqrsummasklocal"

    let n_vwsqrsummask (x : t) (w : t) (id : t) =
      if Sundials_configuration.safe then (check x w; check x id);
      c_n_vwsqrsummask x w id
  end
end (* }}} *)

(* (* Too slow! *)
module ArrayOps = Nvector_array.MakeOps (struct
    type data = RealArray.t

    let get       = Bigarray.Array1.get
    let set       = Bigarray.Array1.set
    let fill      = Bigarray.Array1.fill

    let make      = RealArray.make
    let length    = RealArray.length
    let clone     = RealArray.clone
  end)
module DataOps = ArrayOps.DataOps
*)

module DataOps =
  struct (* {{{ *)

    module A = Bigarray.Array1
    type t = RealArray.t

    let n_vclone     = RealArray.copy

    let arr_vaxpy a (x : t) (y : t) =
      if a = 1.0 then
        for i = 0 to A.dim x - 1 do
          A.set y i (A.get y i +. A.get x i)
        done
      else if a = -1.0 then
        for i = 0 to A.dim x - 1 do
          A.set y i (A.get y i -. A.get x i)
        done
      else
        for i = 0 to A.dim x - 1 do
          A.set y i (A.get y i +. a *. A.get x i)
        done

    let n_vlinearsum a (x : t) b (y : t) (z : t) =
      if b = 1.0 && z == y then
        arr_vaxpy a x y
      else if a = 1.0 && z == x then
        arr_vaxpy b y x
      else if a = 1.0 && b = 1.0 then
        for i = 0 to A.dim x - 1 do
          A.set z i (A.get x i +. A.get y i)
        done
      else if (a = 1.0 && b = -1.0) || (a = -1.0 && b == 1.0) then
        let v1, v2 = if (a = 1.0 && b = -1.0) then y, x else x, y in
        for i = 0 to A.dim v1 - 1 do
          A.set z i (A.get v2 i -. A.get v1 i)
        done
      else if a = 1.0 || b = 1.0 then
        let c, v1, v2 = if a = 1.0 then b, y, x else a, x, y in
        for i = 0 to A.dim v1 - 1 do
          A.set z i (c *. A.get v1 i +. A.get v2 i)
        done
      else if a = -1.0 || b = -1.0 then
        let c, v1, v2 = if a = -1.0 then b, y, x else a, x, y in
        for i = 0 to A.dim v1 - 1 do
          A.set z i (c *. A.get v1 i -. A.get v2 i)
        done
      else if a = b then
        for i = 0 to A.dim x - 1 do
          A.set z i (a *. (A.get x i +. A.get y i))
        done
      else if a = -.b then
        for i = 0 to A.dim x - 1 do
          A.set z i (a *. (A.get x i -. A.get y i))
        done
      else
        for i = 0 to A.dim x - 1 do
          A.set z i (a *. A.get x i +. b *. A.get y i)
        done

    let n_vconst c a = A.fill a c

    let n_vscale c (x : t) (z : t) =
      if c = 1.0 then
        for i = 0 to A.dim x - 1 do
          A.set z i (A.get x i)
        done
      else if c = -1.0 then
        for i = 0 to A.dim x - 1 do
          A.set z i (-. A.get x i)
        done
      else
        for i = 0 to A.dim x - 1 do
          A.set z i (c *. A.get x i)
        done

    let n_vaddconst (x : t) b (z : t) =
      for i = 0 to A.dim x - 1 do
        A.set z i (A.get x i +. b)
      done

    let n_vmaxnorm (x : t) =
      let max = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        let ax = abs_float (A.get x i) in
        if ax > !max then max := ax
      done;
      !max

    let n_vwrmsnorm (x : t) (w : t) =
      let a = ref 0.0 in
      let lx = A.dim x in
      for i = 0 to lx - 1 do
        a := !a +. ((A.get x i) *. (A.get w i) *. (A.get x i) *. (A.get w i))
      done;
      sqrt (!a /. float lx)

    let n_vwrmsnormmask (x : t) (w : t) (id : t) =
      let a = ref 0.0 in
      let lx = A.dim x in
      for i = 0 to lx - 1 do
        if A.get id i > 0.0 then
          a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
      done;
      sqrt (!a /. float lx)

    let n_vmin (x : t) =
      let min = ref max_float in
      for i = 0 to A.dim x - 1 do
        let xv = A.get x i in
        if xv < !min then min := xv
      done;
      !min

    let n_vdotprod (x : t) (y : t) =
      let a = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        a := !a +. (A.get x i *. A.get y i)
      done;
      !a

    let n_vcompare c (x : t) (z : t) =
      for i = 0 to A.dim x - 1 do
        A.set z i (if abs_float (A.get x i) >= c then 1.0 else 0.0)
      done

    let n_vinvtest (x : t) (z : t) =
      let r = ref true in
      for i = 0 to A.dim x - 1 do
        if A.get x i = 0.0 then r := false
        else A.set z i (1.0 /. (A.get x i))
      done;
      !r

    let n_vwl2norm (x : t) (w : t) =
      let a = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
      done;
      sqrt !a

    let n_vl1norm (x : t) =
      let a = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        a := !a +. abs_float (A.get x i)
      done;
      !a

    let n_vconstrmask c (x : t) (m : t) =
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
      for i = 0 to A.dim c - 1 do
        A.set m i (f (A.get c i) (A.get x i))
      done;
      !test

    let n_vminquotient (n : t) (d : t) =
      let m = ref Config.big_real in
      for i = 0 to A.dim n - 1 do
        if (A.get d i) <> 0.0 then
          m := min !m (A.get n i /. A.get d i)
      done;
      !m

    let n_vprod (x : t) (y : t) (z : t) =
      for i = 0 to A.dim x - 1 do
        A.set z i (A.get x i *. A.get y i)
      done

    let n_vdiv (x : t) (y : t) (z : t) =
      for i = 0 to A.dim x - 1 do
        A.set z i (A.get x i /. A.get y i)
      done

    let n_vabs (x : t) (z : t) =
      for i = 0 to A.dim x - 1 do
        A.set z i (abs_float (A.get x i))
      done

    let n_vinv (x : t) (z : t) =
      for i = 0 to A.dim x - 1 do
        A.set z i (1.0 /. (A.get x i))
      done

    let n_vspace (x : t) = (A.dim x, 1)

    let n_vgetlength (x : t) = A.dim x

    (* fused and array operations *)

    let n_vlinearcombination (ca : RealArray.t) (xa : t array) (z : t) =
      let nvec = Array.length xa in
      if nvec = 1 then n_vscale ca.{0} xa.(0) z
      else if nvec = 2 then n_vlinearsum ca.{0} xa.(0) ca.{1} xa.(1) z
      else
        let n = RealArray.length z in
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

    let n_vscaleaddmulti (aa : RealArray.t) (x : t) (ya : t array) (za : t array) =
      let nvec = Array.length ya in
      if nvec = 1 then n_vlinearsum aa.{0} x 1.0 ya.(0) za.(0)
      else
        let n = RealArray.length x in
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

    let n_vdotprodmulti (x : t) (ya : t array) (dp : RealArray.t) =
      let nvec = Array.length ya in
      if nvec = 1 then dp.{0} <- n_vdotprod x ya.(0)
      else
        let n = RealArray.length x in
        for i = 0 to nvec - 1 do
          let y = ya.(i) in
          dp.{i} <- 0.0;
          for j = 0 to n - 1 do
            dp.{i} <- dp.{i} +. A.get x j *. A.get y j
          done
        done

    let arr_vaxpy_array a (xa : t array) (ya : t array) =
      let nvec = Array.length xa in
      let n = A.dim xa.(0) in
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

    let v_sumvectorarray (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let n = A.dim xa.(0) in
      for i = 0 to nvec - 1 do
        let x, y, z = xa.(i), ya.(i), za.(i) in
        for j = 0 to n - 1 do
          A.set z j (A.get x j +. A.get y j)
        done
      done

    let v_diffvectorarray (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let n = A.dim xa.(0) in
      for i = 0 to nvec - 1 do
        let x, y, z = xa.(i), ya.(i), za.(i) in
        for j = 0 to n - 1 do
          A.set z j (A.get x j -. A.get y j)
        done
      done

    let v_lin1vectorarray a (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let n = A.dim xa.(0) in
      for i = 0 to nvec - 1 do
        let x, y, z = xa.(i), ya.(i), za.(i) in
        for j = 0 to n - 1 do
          A.set z j (a *. A.get x j +. A.get y j)
        done
      done

    let v_lin2vectorarray a (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let n = A.dim xa.(0) in
      for i = 0 to nvec - 1 do
        let x, y, z = xa.(i), ya.(i), za.(i) in
        for j = 0 to n - 1 do
          A.set z j (a *. A.get x j -. A.get y j)
        done
      done

    let v_scalesumvectorarray c (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let n = A.dim xa.(0) in
      for i = 0 to nvec - 1 do
        let x, y, z = xa.(i), ya.(i), za.(i) in
        for j = 0 to n - 1 do
          A.set z j (c *. (A.get x j +. A.get y j))
        done
      done

    let v_scalediffvectorarray c (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let n = A.dim xa.(0) in
      for i = 0 to nvec - 1 do
        let x, y, z = xa.(i), ya.(i), za.(i) in
        for j = 0 to n - 1 do
          A.set z j (c *. (A.get x j -. A.get y j))
        done
      done

    let n_vlinearsumvectorarray a (xa : t array) b (ya : t array) (za : t array) =
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
        let n = A.dim xa.(0) in
        for i = 0 to nvec - 1 do
          let x, y, z = xa.(i), ya.(i), za.(i) in
          for j = 0 to n - 1 do
            A.set z j (a *. A.get x j +. b *. A.get y j)
          done
        done

    let n_vscalevectorarray (c : RealArray.t) (xa : t array) (za : t array) =
      let nvec = Array.length xa in
      if nvec = 1 then n_vscale c.{0} xa.(0) za.(0)
      else
        let n = A.dim xa.(0) in
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

    let n_vconstvectorarray c (za : t array) =
      let nvec = Array.length za in
      if nvec = 1 then n_vconst c za.(0)
      else
        let n = A.dim za.(0) in
        for i = 0 to nvec - 1 do
          let z = za.(i) in
          for j = 0 to n - 1 do
            A.set z j c
          done
        done

    let n_vwrmsnormvectorarray (xa : t array) (wa : t array) (nrm : RealArray.t) =
      let nvec = Array.length xa in
      if nvec = 1 then nrm.{0} <- n_vwrmsnorm xa.(0) wa.(0)
      else
        let n = A.dim xa.(0) in
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

    let n_vwrmsnormmaskvectorarray (xa : t array) (wa : t array) (id : t)
                                   (nrm : RealArray.t) =
      let nvec = Array.length xa in
      if nvec = 1 then nrm.{0} <- n_vwrmsnormmask xa.(0) wa.(0) id
      else
        let n = A.dim xa.(0) in
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

    let n_vscaleaddmultivectorarray (ra : RealArray.t) (xa : t array)
                                    (yaa : t array array) (zaa : t array array) =
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
        let n = A.dim xa.(0) in
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

    let n_vlinearcombinationvectorarray (ca : RealArray.t) (xaa : t array array)
                                        (za : t array) =
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
          let n = A.dim za.(0) in
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

      let n_vwsqrsum (x : t) (w : t) =
        let a = ref 0.0 in
        let lx = A.dim x in
        for i = 0 to lx - 1 do
          a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
        done;
        !a

      let n_vwsqrsummask (x : t) (w : t) (id : t) =
        let a = ref 0.0 in
        let lx = A.dim x in
        for i = 0 to lx - 1 do
          if A.get id i > 0.0 then
            a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
        done;
        !a
    end
  end (* }}} *)

