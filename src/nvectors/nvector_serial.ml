
open Sundials

type data = RealArray.t
type kind = [`Serial]
type t = (data, kind) Nvector.t

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

external c_wrap : RealArray.t -> (t -> bool) -> (t -> t) -> Context.t -> t
  = "sunml_nvec_wrap_serial"

let rec wrap ?context ?(with_fused_ops=false) v =
  let len = RealArray.length v in
  let ctx = Sundials_impl.Context.get context in
  let nv =
    c_wrap v (fun nv' -> len = RealArray.length (unwrap nv')) clone ctx
  in
  if with_fused_ops then c_enablefusedops_serial nv true;
  nv

and clone nv =
  let nv' = wrap ~context:(Nvector.context nv) (RealArray.copy (unwrap nv)) in
  if Sundials_impl.Version.lt400 then ()
  else begin
    c_enablelinearcombination_serial nv'
      (Nvector.Ops.has_linearcombination nv);
    c_enablescaleaddmulti_serial nv'
      (Nvector.Ops.has_scaleaddmulti nv);
    c_enabledotprodmulti_serial nv'
      (Nvector.Ops.has_dotprodmulti nv);
    c_enablelinearsumvectorarray_serial nv'
      (Nvector.Ops.has_linearsumvectorarray nv);
    c_enablescalevectorarray_serial nv'
      (Nvector.Ops.has_scalevectorarray nv);
    c_enableconstvectorarray_serial nv'
      (Nvector.Ops.has_constvectorarray nv);
    c_enablewrmsnormvectorarray_serial nv'
      (Nvector.Ops.has_wrmsnormvectorarray nv);
    c_enablewrmsnormmaskvectorarray_serial nv'
      (Nvector.Ops.has_wrmsnormmaskvectorarray nv);
    c_enablescaleaddmultivectorarray_serial nv'
      (Nvector.Ops.has_scaleaddmultivectorarray nv);
    c_enablelinearcombinationvectorarray_serial nv'
      (Nvector.Ops.has_linearcombinationvectorarray nv)
  end;
  nv'

let make ?context ?with_fused_ops n iv =
  wrap ?context ?with_fused_ops (RealArray.make n iv)

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
      -> Context.t
      -> Nvector.any
    = "sunml_nvec_anywrap_serial"

  let rec wrap
      ?context
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
      if not Sundials_impl.Version.has_nvector_get_id
        then raise Config.NotImplementedBySundialsVersion;
      let len = RealArray.length v in
      let check nv' =
        match unwrap nv' with
        | Nvector.RA ra ->
            len = RealArray.length ra && Nvector.get_id nv' = Nvector.Serial
        | _ -> false
      in
      let ctx = Sundials_impl.Context.get context in
      let nv =
        c_any_wrap [%extension_constructor Nvector.RA] v check clone ctx
      in
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
    let nv' = wrap ~context:(Nvector.context nv) (RealArray.copy v) in
    c_enablelinearcombination_serial nv'
      (Nvector.Ops.has_linearcombination nv);
    c_enablescaleaddmulti_serial nv'
      (Nvector.Ops.has_scaleaddmulti nv);
    c_enabledotprodmulti_serial nv'
      (Nvector.Ops.has_dotprodmulti nv);
    c_enablelinearsumvectorarray_serial nv'
      (Nvector.Ops.has_linearsumvectorarray nv);
    c_enablescalevectorarray_serial nv'
      (Nvector.Ops.has_scalevectorarray nv);
    c_enableconstvectorarray_serial nv'
      (Nvector.Ops.has_constvectorarray nv);
    c_enablewrmsnormvectorarray_serial nv'
      (Nvector.Ops.has_wrmsnormvectorarray nv);
    c_enablewrmsnormmaskvectorarray_serial nv'
      (Nvector.Ops.has_wrmsnormmaskvectorarray nv);
    c_enablescaleaddmultivectorarray_serial nv'
      (Nvector.Ops.has_scaleaddmultivectorarray nv);
    c_enablelinearcombinationvectorarray_serial nv'
      (Nvector.Ops.has_linearcombinationvectorarray nv);
    nv'

  let make
      ?context
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
    = wrap ?context
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
           (RealArray.make n iv)

  let unwrap nv =
    match Nvector.unwrap nv with
    | Nvector.RA a -> a
    | _ -> raise Nvector.BadGenericType

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
    = if Sundials_impl.Version.lt400
        then raise Config.NotImplementedBySundialsVersion;
      if Nvector.get_id nv <> Nvector.Serial then raise Nvector.BadGenericType;
      do_enable c_enablefusedops_serial nv
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
end (* }}} *)

module Ops : Nvector.NVECTOR_OPS with type t = t =
struct (* {{{ *)
  type t = (data, kind) Nvector.t
  let check = Nvector.check

  let clone = Nvector.clone

  external c_linearsum    : float -> t -> float -> t -> t -> unit
    = "sunml_nvec_ser_linearsum" [@@noalloc]

  let linearsum a (x : t) b (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_linearsum a x b y z

  external const          : float -> t -> unit
    = "sunml_nvec_ser_const" [@@noalloc]

  external c_prod         : t -> t -> t -> unit
    = "sunml_nvec_ser_prod" [@@noalloc]

  let prod (x : t) (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_prod x y z

  external c_div          : t -> t -> t -> unit
    = "sunml_nvec_ser_div" [@@noalloc]

  let div (x : t) (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_div x y z

  external c_scale        : float -> t -> t -> unit
    = "sunml_nvec_ser_scale" [@@noalloc]

  let scale c (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_scale c x z

  external c_abs          : t -> t -> unit
    = "sunml_nvec_ser_abs" [@@noalloc]

  let abs (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_abs x z

  external c_inv          : t -> t -> unit
    = "sunml_nvec_ser_inv" [@@noalloc]

  let inv (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_inv x z

  external c_addconst     : t -> float -> t -> unit
    = "sunml_nvec_ser_addconst" [@@noalloc]

  let addconst (x : t) b (z : t) =
    if Sundials_configuration.safe then check x z;
    c_addconst x b z

  external c_dotprod      : t -> t -> float
    = "sunml_nvec_ser_dotprod"

  let dotprod (x : t) (y : t) =
    if Sundials_configuration.safe then check x y;
    c_dotprod x y

  external maxnorm        : t -> float
    = "sunml_nvec_ser_maxnorm"

  external c_wrmsnorm     : t -> t -> float
    = "sunml_nvec_ser_wrmsnorm"

  let wrmsnorm (x : t) (w : t) =
    if Sundials_configuration.safe then check x w;
    c_wrmsnorm x w

  external c_wrmsnormmask : t -> t -> t -> float
    = "sunml_nvec_ser_wrmsnormmask"

  let wrmsnormmask (x : t) (w : t) (id : t) =
    if Sundials_configuration.safe then (check x w; check x id);
    c_wrmsnormmask x w id

  external min            : t -> float
    = "sunml_nvec_ser_min"

  external c_wl2norm      : t -> t -> float
    = "sunml_nvec_ser_wl2norm"

  let wl2norm (x : t) (w : t) =
    if Sundials_configuration.safe then check x w;
    c_wl2norm x w

  external l1norm         : t -> float
    = "sunml_nvec_ser_l1norm"

  external c_compare      : float -> t -> t -> unit
    = "sunml_nvec_ser_compare" [@@noalloc]

  let compare c (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_compare c x z

  external c_invtest      : t -> t -> bool
    = "sunml_nvec_ser_invtest" [@@noalloc]

  let invtest (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_invtest x z

  external c_constrmask   : t -> t -> t -> bool
    = "sunml_nvec_ser_constrmask" [@@noalloc]

  let constrmask (c : t) (x : t) (m : t) =
    if Sundials_configuration.safe then (check c x; check c m);
    c_constrmask c x m

  external c_minquotient  : t -> t -> float
    = "sunml_nvec_ser_minquotient"

  let minquotient (n : t) (d : t) =
    if Sundials_configuration.safe then check n d;
    c_minquotient n d

  external space          : t -> int * int
    = "sunml_nvec_ser_space"

  external getlength      : t -> int
    = "sunml_nvec_ser_getlength"

  external c_print_file : t -> Logfile.t option -> unit
    = "sunml_nvec_ser_print_file"

  let print ?logfile nv = c_print_file nv logfile

  external c_linearcombination : RealArray.t -> t array -> t -> unit
    = "sunml_nvec_ser_linearcombination"

  let linearcombination ca (xa : t array) (z : t) =
    if Sundials_impl.Version.lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe then Array.iter (check z) xa;
    c_linearcombination ca xa z

  let same_len' n ya =
    if n <> Array.length ya then invalid_arg "arrays of unequal length"
  let same_len xa ya = same_len' (Array.length xa) ya

  external c_scaleaddmulti : RealArray.t -> t -> t array -> t array -> unit
    = "sunml_nvec_ser_scaleaddmulti"

  let scaleaddmulti aa (x : t) (ya : t array) (za : t array) =
    if Sundials_impl.Version.lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe then
      (Array.iter (check x) ya; Array.iter (check x) za;
       let nv = RealArray.length aa in
       same_len' nv ya; same_len' nv za);
    c_scaleaddmulti aa x ya za

  external c_dotprodmulti  : t -> t array -> RealArray.t -> unit
    = "sunml_nvec_ser_dotprodmulti"

  let dotprodmulti (x : t) (ya : t array) dp =
    if Sundials_impl.Version.lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe then
      (let nv = RealArray.length dp in
       same_len' nv ya;
       Array.iter (check x) ya);
    c_dotprodmulti x ya dp

  external c_linearsumvectorarray
    : float -> t array -> float -> t array -> t array -> unit
    = "sunml_nvec_ser_linearsumvectorarray"

  let linearsumvectorarray a (xa : t array) b (ya : t array) (za : t array) =
    if Sundials_impl.Version.lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (check x) ya;
          Array.iter (check x) za;
          same_len xa ya; same_len xa za);
    c_linearsumvectorarray a xa b ya za

  external c_scalevectorarray
    : RealArray.t -> t array -> t array -> unit
    = "sunml_nvec_ser_scalevectorarray"

  let scalevectorarray c (xa : t array) (za : t array) =
    if Sundials_impl.Version.lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (check x) za;
          same_len xa za);
    c_scalevectorarray c xa za

  external c_constvectorarray
    : float -> t array -> unit
    = "sunml_nvec_ser_constvectorarray"

  let constvectorarray c (za : t array) =
    if Sundials_impl.Version.lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe
    then (let z = Array.get za 0 in
          Array.iter (check z) za);
    c_constvectorarray c za

  external c_wrmsnormvectorarray
    : t array -> t array -> RealArray.t -> unit
    = "sunml_nvec_ser_wrmsnormvectorarray"

  let wrmsnormvectorarray (xa : t array) (wa : t array) nrm =
    if Sundials_impl.Version.lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (check x) wa;
         same_len xa wa);
    c_wrmsnormvectorarray xa wa nrm

  external c_wrmsnormmaskvectorarray
    : t array -> t array -> t -> RealArray.t -> unit
    = "sunml_nvec_ser_wrmsnormmaskvectorarray"

  let wrmsnormmaskvectorarray (xa : t array) (wa : t array) (id : t) nrm =
    if Sundials_impl.Version.lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe
    then (Array.iter (check id) xa;
          Array.iter (check id) wa;
          same_len xa wa);
    c_wrmsnormmaskvectorarray xa wa id nrm

  external c_scaleaddmultivectorarray
    : RealArray.t -> t array -> t array array -> t array array -> unit
    = "sunml_nvec_ser_scaleaddmultivectorarray"

  let scaleaddmultivectorarray ra (xa : t array) (yaa : t array array)
                                     (zaa : t array array) =
    if Sundials_impl.Version.lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          let ns = RealArray.length ra in
          let nv = Array.length xa in
          same_len' ns yaa;
          same_len' ns zaa;
          Array.iter (check x) xa;
          Array.iter (fun ya -> same_len' nv ya; Array.iter (check x) ya) yaa;
          Array.iter (fun za -> same_len' nv za; Array.iter (check x) za) zaa;
          same_len yaa zaa);
    c_scaleaddmultivectorarray ra xa yaa zaa

  external c_linearcombinationvectorarray
    : RealArray.t -> t array array -> t array -> unit
    = "sunml_nvec_ser_linearcombinationvectorarray"

  let linearcombinationvectorarray ca (xaa : t array array) (za : t array) =
    if Sundials_impl.Version.lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe
    then (let z = Array.get za 0 in
          let ns = RealArray.length ca in
          let nv = Array.length za in
          same_len' ns xaa;
          Array.iter (check z) za;
          Array.iter (fun xa -> same_len' nv xa; Array.iter (check z) xa) xaa);
    c_linearcombinationvectorarray ca xaa za

  module Local = struct
    let dotprod     = dotprod
    let maxnorm     = maxnorm
    let min         = min
    let l1norm      = l1norm
    let invtest     = invtest
    let constrmask  = constrmask
    let minquotient = minquotient

    external c_wsqrsum
      : t -> t -> float
      = "sunml_nvec_ser_wsqrsumlocal"

    let wsqrsum (x : t) (w : t) =
      if Sundials_impl.Version.lt500
        then raise Config.NotImplementedBySundialsVersion;
      if Sundials_configuration.safe then check x w;
      c_wsqrsum x w

    external c_wsqrsummask
      : t -> t -> t -> float
      = "sunml_nvec_ser_wsqrsummasklocal"

    let wsqrsummask (x : t) (w : t) (id : t) =
      if Sundials_impl.Version.lt500
        then raise Config.NotImplementedBySundialsVersion;
      if Sundials_configuration.safe then (check x w; check x id);
      c_wsqrsummask x w id

    let dotprodmulti = dotprodmulti
    let dotprodmulti_allreduce _ _ = raise Nvector.OperationNotProvided
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

    let clone     = RealArray.copy

    let floatmin (x : float) (y : float) = min x y

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

    let linearsum a (x : t) b (y : t) (z : t) =
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

    let const c a = A.fill a c

    let scale c (x : t) (z : t) =
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

    let addconst (x : t) b (z : t) =
      for i = 0 to A.dim x - 1 do
        A.set z i (A.get x i +. b)
      done

    let maxnorm (x : t) =
      let max = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        let ax = abs_float (A.get x i) in
        if ax > !max then max := ax
      done;
      !max

    let wrmsnorm (x : t) (w : t) =
      let a = ref 0.0 in
      let lx = A.dim x in
      for i = 0 to lx - 1 do
        a := !a +. ((A.get x i) *. (A.get w i) *. (A.get x i) *. (A.get w i))
      done;
      sqrt (!a /. float lx)

    let wrmsnormmask (x : t) (w : t) (id : t) =
      let a = ref 0.0 in
      let lx = A.dim x in
      for i = 0 to lx - 1 do
        if A.get id i > 0.0 then
          a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
      done;
      sqrt (!a /. float lx)

    let min (x : t) =
      let min = ref max_float in
      for i = 0 to A.dim x - 1 do
        let xv = A.get x i in
        if xv < !min then min := xv
      done;
      !min

    let dotprod (x : t) (y : t) =
      let a = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        a := !a +. (A.get x i *. A.get y i)
      done;
      !a

    let compare c (x : t) (z : t) =
      for i = 0 to A.dim x - 1 do
        A.set z i (if abs_float (A.get x i) >= c then 1.0 else 0.0)
      done

    let invtest (x : t) (z : t) =
      let r = ref true in
      for i = 0 to A.dim x - 1 do
        if A.get x i = 0.0 then r := false
        else A.set z i (1.0 /. (A.get x i))
      done;
      !r

    let wl2norm (x : t) (w : t) =
      let a = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
      done;
      sqrt !a

    let l1norm (x : t) =
      let a = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        a := !a +. abs_float (A.get x i)
      done;
      !a

    let constrmask c (x : t) (m : t) =
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

    let minquotient (n : t) (d : t) =
      let m = ref Config.big_real in
      for i = 0 to A.dim n - 1 do
        if (A.get d i) <> 0.0 then
          m := floatmin !m (A.get n i /. A.get d i)
      done;
      !m

    let prod (x : t) (y : t) (z : t) =
      for i = 0 to A.dim x - 1 do
        A.set z i (A.get x i *. A.get y i)
      done

    let div (x : t) (y : t) (z : t) =
      for i = 0 to A.dim x - 1 do
        A.set z i (A.get x i /. A.get y i)
      done

    let abs (x : t) (z : t) =
      for i = 0 to A.dim x - 1 do
        A.set z i (abs_float (A.get x i))
      done

    let inv (x : t) (z : t) =
      for i = 0 to A.dim x - 1 do
        A.set z i (1.0 /. (A.get x i))
      done

    let space (x : t) = (A.dim x, 1)

    let getlength (x : t) = A.dim x

    let print ?(logfile=Logfile.stdout) (x : t) =
      for i = 0 to A.dim x - 1 do
        Logfile.output_string logfile (Printf.sprintf "%19.16g" (A.get x i));
        Logfile.output_string logfile "\n"
      done;
      Logfile.output_string logfile "\n"

    (* fused and array operations *)

    let linearcombination (ca : RealArray.t) (xa : t array) (z : t) =
      let nvec = Array.length xa in
      if nvec = 1 then scale ca.{0} xa.(0) z
      else if nvec = 2 then linearsum ca.{0} xa.(0) ca.{1} xa.(1) z
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

    let scaleaddmulti (aa : RealArray.t) (x : t) (ya : t array) (za : t array) =
      let nvec = Array.length ya in
      if nvec = 1 then linearsum aa.{0} x 1.0 ya.(0) za.(0)
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

    let dotprodmulti (x : t) (ya : t array) (dp : RealArray.t) =
      let nvec = Array.length ya in
      if nvec = 1 then dp.{0} <- dotprod x ya.(0)
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

    let linearsumvectorarray a (xa : t array) b (ya : t array) (za : t array) =
      let nvec = Array.length ya in
      if nvec = 1 then linearsum a xa.(0) b ya.(0) za.(0)
      else if b =  1.0 && (za == ya) then arr_vaxpy_array a xa ya
      else if a =  1.0 && (za == xa) then arr_vaxpy_array b ya xa
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

    let scalevectorarray (c : RealArray.t) (xa : t array) (za : t array) =
      let nvec = Array.length xa in
      if nvec = 1 then scale c.{0} xa.(0) za.(0)
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

    let constvectorarray c (za : t array) =
      let nvec = Array.length za in
      if nvec = 1 then const c za.(0)
      else
        let n = A.dim za.(0) in
        for i = 0 to nvec - 1 do
          let z = za.(i) in
          for j = 0 to n - 1 do
            A.set z j c
          done
        done

    let wrmsnormvectorarray (xa : t array) (wa : t array) (nrm : RealArray.t) =
      let nvec = Array.length xa in
      if nvec = 1 then nrm.{0} <- wrmsnorm xa.(0) wa.(0)
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

    let wrmsnormmaskvectorarray (xa : t array) (wa : t array) (id : t)
                                   (nrm : RealArray.t) =
      let nvec = Array.length xa in
      if nvec = 1 then nrm.{0} <- wrmsnormmask xa.(0) wa.(0) id
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

    let scaleaddmultivectorarray (ra : RealArray.t) (xa : t array)
                                    (yaa : t array array) (zaa : t array array) =
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

    let linearcombinationvectorarray (ca : RealArray.t) (xaa : t array array)
                                        (za : t array) =
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
      let dotprod     = dotprod
      let maxnorm     = maxnorm
      let min         = min
      let l1norm      = l1norm
      let invtest     = invtest
      let constrmask  = constrmask
      let minquotient = minquotient

      let wsqrsum (x : t) (w : t) =
        let a = ref 0.0 in
        let lx = A.dim x in
        for i = 0 to lx - 1 do
          a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
        done;
        !a

      let wsqrsummask (x : t) (w : t) (id : t) =
        let a = ref 0.0 in
        let lx = A.dim x in
        for i = 0 to lx - 1 do
          if A.get id i > 0.0 then
            a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
        done;
        !a

      let dotprodmulti = dotprodmulti
      let dotprodmulti_allreduce _ _ = raise Nvector.OperationNotProvided
    end
  end (* }}} *)

