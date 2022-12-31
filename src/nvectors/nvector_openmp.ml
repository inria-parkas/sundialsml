open Sundials
type data = RealArray.t
type kind = [`OpenMP|Nvector_serial.kind]
type t = (data, kind) Nvector.t

(* Selectively enable and disable fused and array operations *)
external c_enablefusedops_openmp                     : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_openmp_enablefusedops"
external c_enablelinearcombination_openmp            : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_openmp_enablelinearcombination"
external c_enablescaleaddmulti_openmp                : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_openmp_enablescaleaddmulti"
external c_enabledotprodmulti_openmp                 : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_openmp_enabledotprodmulti"
external c_enablelinearsumvectorarray_openmp         : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_openmp_enablelinearsumvectorarray"
external c_enablescalevectorarray_openmp             : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_openmp_enablescalevectorarray"
external c_enableconstvectorarray_openmp             : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_openmp_enableconstvectorarray"
external c_enablewrmsnormvectorarray_openmp          : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_openmp_enablewrmsnormvectorarray"
external c_enablewrmsnormmaskvectorarray_openmp      : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_openmp_enablewrmsnormmaskvectorarray"
external c_enablescaleaddmultivectorarray_openmp     : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_openmp_enablescaleaddmultivectorarray"
external c_enablelinearcombinationvectorarray_openmp : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_openmp_enablelinearcombinationvectorarray"

let unwrap = Nvector.unwrap

external c_wrap :
     int
  -> RealArray.t
  -> (t -> bool)
  -> (t -> t)
  -> Context.t
  -> t
  = "sunml_nvec_wrap_openmp"

let rec wrap ?context ?(with_fused_ops=false) nthreads v =
  let len = RealArray.length v in
  let check nv' = (len = RealArray.length (unwrap nv')) in
  let ctx = Sundials_impl.Context.get context in
  let nv = c_wrap nthreads v check (clone nthreads) ctx in
  if with_fused_ops then c_enablefusedops_openmp nv true;
  nv

and clone nthreads nv =
  let nv' = wrap ~context:(Nvector.context nv) nthreads
                 (RealArray.copy (unwrap nv))
  in
  if Sundials_impl.Version.lt400 then ()
  else begin
    c_enablelinearcombination_openmp nv'
      (Nvector.Ops.has_linearcombination nv);
    c_enablescaleaddmulti_openmp nv'
      (Nvector.Ops.has_scaleaddmulti nv);
    c_enabledotprodmulti_openmp nv'
      (Nvector.Ops.has_dotprodmulti nv);
    c_enablelinearsumvectorarray_openmp nv'
      (Nvector.Ops.has_linearsumvectorarray nv);
    c_enablescalevectorarray_openmp nv'
      (Nvector.Ops.has_scalevectorarray nv);
    c_enableconstvectorarray_openmp nv'
      (Nvector.Ops.has_constvectorarray nv);
    c_enablewrmsnormvectorarray_openmp nv'
      (Nvector.Ops.has_wrmsnormvectorarray nv);
    c_enablewrmsnormmaskvectorarray_openmp nv'
      (Nvector.Ops.has_wrmsnormmaskvectorarray nv);
    c_enablescaleaddmultivectorarray_openmp nv'
      (Nvector.Ops.has_scaleaddmultivectorarray nv);
    c_enablelinearcombinationvectorarray_openmp nv'
      (Nvector.Ops.has_linearcombinationvectorarray nv)
  end;
  nv'

let pp fmt v = RealArray.pp fmt (unwrap v)

let make ?context ?with_fused_ops nthreads n iv =
  wrap ?context ?with_fused_ops nthreads (RealArray.make n iv)

external num_threads : t -> int
  = "sunml_nvec_openmp_num_threads"

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
  = do_enable c_enablefusedops_openmp nv
              with_fused_ops;
    do_enable c_enablelinearcombination_openmp nv
              with_linear_combination;
    do_enable c_enablescaleaddmulti_openmp nv
              with_scale_add_multi;
    do_enable c_enabledotprodmulti_openmp nv
              with_dot_prod_multi;
    do_enable c_enablelinearsumvectorarray_openmp nv
              with_linear_sum_vector_array;
    do_enable c_enablescalevectorarray_openmp nv
              with_scale_vector_array;
    do_enable c_enableconstvectorarray_openmp nv
              with_const_vector_array;
    do_enable c_enablewrmsnormvectorarray_openmp nv
              with_wrms_norm_vector_array;
    do_enable c_enablewrmsnormmaskvectorarray_openmp nv
              with_wrms_norm_mask_vector_array;
    do_enable c_enablescaleaddmultivectorarray_openmp nv
              with_scale_add_multi_vector_array;
    do_enable c_enablelinearcombinationvectorarray_openmp nv
              with_linear_combination_vector_array

module Any = struct (* {{{ *)

  external c_any_wrap
    : extension_constructor
      -> int
      -> RealArray.t
      -> (Nvector.any -> bool)
      -> (Nvector.any -> Nvector.any)
      -> Context.t
      -> Nvector.any
    = "sunml_nvec_anywrap_openmp_byte"
      "sunml_nvec_anywrap_openmp"

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
      nthreads v
    =
      if not Sundials_impl.Version.has_nvector_get_id
        then raise Config.NotImplementedBySundialsVersion;
      let len = RealArray.length v in
      let check nv =
        match unwrap nv with
        | Nvector.RA ra ->
            len = RealArray.length ra && Nvector.get_id nv = Nvector.OpenMP
        | _ -> false
      in
      let ctx = Sundials_impl.Context.get context in
      let nv = c_any_wrap [%extension_constructor Nvector.RA]
                          nthreads v
                          check (clone nthreads) ctx
      in
      if with_fused_ops
        then c_enablefusedops_openmp nv true;
      if with_fused_ops
        then c_enablefusedops_openmp nv true;
      if with_linear_combination
        then c_enablelinearcombination_openmp nv true;
      if with_scale_add_multi
        then c_enablescaleaddmulti_openmp nv true;
      if with_dot_prod_multi
        then c_enabledotprodmulti_openmp nv true;
      if with_linear_sum_vector_array
        then c_enablelinearsumvectorarray_openmp nv true;
      if with_scale_vector_array
        then c_enablescalevectorarray_openmp nv true;
      if with_const_vector_array
        then c_enableconstvectorarray_openmp nv true;
      if with_wrms_norm_vector_array
        then c_enablewrmsnormvectorarray_openmp nv true;
      if with_wrms_norm_mask_vector_array
        then c_enablewrmsnormmaskvectorarray_openmp nv true;
      if with_scale_add_multi_vector_array
        then c_enablescaleaddmultivectorarray_openmp nv true;
      if with_linear_combination_vector_array
        then c_enablelinearcombinationvectorarray_openmp nv true;
      nv

  and clone nthreads nv =
    let v = match unwrap nv with
            | Nvector.RA v -> v
            | _ -> assert false
    in
    let nv' = wrap ~context:(Nvector.context nv) nthreads (RealArray.copy v) in
    c_enablelinearcombination_openmp nv'
      (Nvector.Ops.has_linearcombination nv);
    c_enablescaleaddmulti_openmp nv'
      (Nvector.Ops.has_scaleaddmulti nv);
    c_enabledotprodmulti_openmp nv'
      (Nvector.Ops.has_dotprodmulti nv);
    c_enablelinearsumvectorarray_openmp nv'
      (Nvector.Ops.has_linearsumvectorarray nv);
    c_enablescalevectorarray_openmp nv'
      (Nvector.Ops.has_scalevectorarray nv);
    c_enableconstvectorarray_openmp nv'
      (Nvector.Ops.has_constvectorarray nv);
    c_enablewrmsnormvectorarray_openmp nv'
      (Nvector.Ops.has_wrmsnormvectorarray nv);
    c_enablewrmsnormmaskvectorarray_openmp nv'
      (Nvector.Ops.has_wrmsnormmaskvectorarray nv);
    c_enablescaleaddmultivectorarray_openmp nv'
      (Nvector.Ops.has_scaleaddmultivectorarray nv);
    c_enablelinearcombinationvectorarray_openmp nv'
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
      nthreads n iv
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
           nthreads (RealArray.make n iv)

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
      do_enable c_enablefusedops_openmp nv
                with_fused_ops;
      do_enable c_enablelinearcombination_openmp nv
                with_linear_combination;
      do_enable c_enablescaleaddmulti_openmp nv
                with_scale_add_multi;
      do_enable c_enabledotprodmulti_openmp nv
                with_dot_prod_multi;
      do_enable c_enablelinearsumvectorarray_openmp nv
                with_linear_sum_vector_array;
      do_enable c_enablescalevectorarray_openmp nv
                with_scale_vector_array;
      do_enable c_enableconstvectorarray_openmp nv
                with_const_vector_array;
      do_enable c_enablewrmsnormvectorarray_openmp nv
                with_wrms_norm_vector_array;
      do_enable c_enablewrmsnormmaskvectorarray_openmp nv
                with_wrms_norm_mask_vector_array;
      do_enable c_enablescaleaddmultivectorarray_openmp nv
                with_scale_add_multi_vector_array;
      do_enable c_enablelinearcombinationvectorarray_openmp nv
                with_linear_combination_vector_array
end (* }}} *)

module Ops = struct (* {{{ *)
  type t = (RealArray.t, kind) Nvector.t
  let check = Nvector.check

  let clone nv =
    let data = Nvector.unwrap nv in
    wrap (num_threads nv) (RealArray.copy data)

  external c_linearsum    : float -> t -> float -> t -> t -> unit
    = "sunml_nvec_openmp_linearsum"

  let linearsum a (x : t) b (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_linearsum a x b y z

  external const        : float -> t -> unit
    = "sunml_nvec_openmp_const" [@@noalloc]

  external c_prod         : t -> t -> t -> unit
    = "sunml_nvec_openmp_prod" [@@noalloc]

  let prod (x : t) (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_prod x y z

  external c_div          : t -> t -> t -> unit
    = "sunml_nvec_openmp_div" [@@noalloc]

  let div (x : t) (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_div x y z

  external c_scale        : float -> t -> t -> unit
    = "sunml_nvec_openmp_scale" [@@noalloc]

  let scale c (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_scale c x z

  external c_abs          : t -> t -> unit
    = "sunml_nvec_openmp_abs" [@@noalloc]

  let abs (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_abs x z

  external c_inv          : t -> t -> unit
    = "sunml_nvec_openmp_inv" [@@noalloc]

  let inv (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_inv x z

  external c_addconst     : t -> float -> t -> unit
    = "sunml_nvec_openmp_addconst" [@@noalloc]

  let addconst (x : t) b (z : t) =
    if Sundials_configuration.safe then check x z;
    c_addconst x b z

  external c_dotprod      : t -> t -> float
    = "sunml_nvec_openmp_dotprod"

  let dotprod (x : t) (y : t) =
    if Sundials_configuration.safe then check x y;
    c_dotprod x y

  external maxnorm      : t -> float
    = "sunml_nvec_openmp_maxnorm"

  external c_wrmsnorm     : t -> t -> float
    = "sunml_nvec_openmp_wrmsnorm"

  let wrmsnorm (x : t) (w : t) =
    if Sundials_configuration.safe then check x w;
    c_wrmsnorm x w

  external c_wrmsnormmask : t -> t -> t -> float
    = "sunml_nvec_openmp_wrmsnormmask"

  let wrmsnormmask (x : t) (w : t) (id : t) =
    if Sundials_configuration.safe then (check x w; check x id);
    c_wrmsnormmask x w id

  external min          : t -> float
    = "sunml_nvec_openmp_min"

  external c_wl2norm      : t -> t -> float
    = "sunml_nvec_openmp_wl2norm"

  let wl2norm (x : t) (w : t) =
    if Sundials_configuration.safe then check x w;
    c_wl2norm x w

  external l1norm       : t -> float
    = "sunml_nvec_openmp_l1norm"

  external c_compare      : float -> t -> t -> unit
    = "sunml_nvec_openmp_compare" [@@noalloc]

  let compare c (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_compare c x z

  external c_invtest      : t -> t -> bool
    = "sunml_nvec_openmp_invtest" [@@noalloc]

  let invtest (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_invtest x z

  external c_constrmask   : t -> t -> t -> bool
    = "sunml_nvec_openmp_constrmask" [@@noalloc]

  let constrmask (c : t) (x : t) (m : t) =
    if Sundials_configuration.safe then (check c x; check c m);
    c_constrmask c x m

  external c_minquotient  : t -> t -> float
    = "sunml_nvec_openmp_minquotient"

  let minquotient (n : t) (d : t) =
    if Sundials_configuration.safe then check n d;
    c_minquotient n d

  external space  : t -> int * int
    = "sunml_nvec_openmp_space"

  external getlength  : t -> int
    = "sunml_nvec_openmp_getlength"

  let getlocallength = getlength

  external c_print_file : t -> Logfile.t option -> unit
    = "sunml_nvec_openmp_print_file"

  let print ?logfile nv = c_print_file nv logfile

  external c_linearcombination : RealArray.t -> t array -> t -> unit
    = "sunml_nvec_openmp_linearcombination"

  let linearcombination ca (xa : t array) (z : t) =
    if Sundials_impl.Version.lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe then Array.iter (check z) xa;
    c_linearcombination ca xa z

  let same_len' n ya =
    if n <> Array.length ya then invalid_arg "arrays of unequal length"
  let same_len xa ya = same_len' (Array.length xa) ya

  external c_scaleaddmulti : RealArray.t -> t -> t array -> t array -> unit
    = "sunml_nvec_openmp_scaleaddmulti"

  let scaleaddmulti aa (x : t) (ya : t array) (za : t array) =
    if Sundials_configuration.safe then
      (Array.iter (check x) ya; Array.iter (check x) za;
       let nv = RealArray.length aa in
       same_len' nv ya; same_len' nv za);
    c_scaleaddmulti aa x ya za

  external c_dotprodmulti : t -> t array -> RealArray.t -> unit
    = "sunml_nvec_openmp_dotprodmulti"

  let dotprodmulti (x : t) (ya : t array) (dp : RealArray.t) =
    if Sundials_configuration.safe then
      (let nv = RealArray.length dp in
       same_len' nv ya;
       Array.iter (check x) ya);
    c_dotprodmulti x ya dp

  external c_linearsumvectorarray
    : float -> t array -> float -> t array -> t array -> unit
    = "sunml_nvec_openmp_linearsumvectorarray"

  let linearsumvectorarray a (xa : t array) b (ya : t array) (za : t array) =
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (check x) ya;
          Array.iter (check x) za;
          same_len xa ya; same_len xa za);
    c_linearsumvectorarray a xa b ya za

  external c_scalevectorarray
    : RealArray.t -> t array -> t array -> unit
    = "sunml_nvec_openmp_scalevectorarray"

  let scalevectorarray c (xa : t array) (za : t array) =
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (check x) za;
          same_len xa za);
    c_scalevectorarray c xa za

  external c_constvectorarray
    : float -> t array -> unit
    = "sunml_nvec_openmp_constvectorarray"

  let constvectorarray c (za : t array) =
    if Sundials_configuration.safe
    then (let z = Array.get za 0 in
          Array.iter (check z) za);
    c_constvectorarray c za

  external c_wrmsnormvectorarray
    : t array -> t array -> RealArray.t -> unit
    = "sunml_nvec_openmp_wrmsnormvectorarray"

  let wrmsnormvectorarray (xa : t array) (wa : t array) nrm =
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (check x) wa;
         same_len xa wa);
    c_wrmsnormvectorarray xa wa nrm

  external c_wrmsnormmaskvectorarray
    : t array -> t array -> t -> RealArray.t -> unit
    = "sunml_nvec_openmp_wrmsnormmaskvectorarray"

  let wrmsnormmaskvectorarray (xa : t array) (wa : t array) (id : t) nrm =
    if Sundials_configuration.safe
    then (Array.iter (check id) xa;
          Array.iter (check id) wa;
          same_len xa wa);
    c_wrmsnormmaskvectorarray xa wa id nrm

  external c_scaleaddmultivectorarray
    : RealArray.t -> t array -> t array array -> t array array -> unit
    = "sunml_nvec_openmp_scaleaddmultivectorarray"

  let scaleaddmultivectorarray ra (xa : t array) (yaa : t array array)
                                  (zaa : t array array) =
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
    = "sunml_nvec_openmp_linearcombinationvectorarray"

  let linearcombinationvectorarray ca (xaa : t array array) (za : t array) =
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
      = "sunml_nvec_openmp_wsqrsumlocal"

    let wsqrsum (x : t) (w : t) =
      if Sundials_configuration.safe then check x w;
      c_wsqrsum x w

    external c_wsqrsummask
      : t -> t -> t -> float
      = "sunml_nvec_openmp_wsqrsummasklocal"

    let wsqrsummask (x : t) (w : t) (id : t) =
      if Sundials_configuration.safe then (check x w; check x id);
      c_wsqrsummask x w id

    let dotprodmulti = dotprodmulti
    let dotprodmulti_allreduce _ _ = raise Nvector.OperationNotProvided
  end
end (* }}} *)

