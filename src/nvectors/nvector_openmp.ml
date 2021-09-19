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

external c_wrap : int -> RealArray.t -> (t -> bool) -> (t -> t) -> t
  = "sunml_nvec_wrap_openmp"

let rec wrap ?(with_fused_ops=false) nthreads v =
  let len = RealArray.length v in
  let check nv' = (len = RealArray.length (unwrap nv')) in
  let nv = c_wrap nthreads v check (clone nthreads) in
  if with_fused_ops then c_enablefusedops_openmp nv true;
  nv

and clone nthreads nv =
  let nv' = wrap nthreads (RealArray.copy (unwrap nv)) in
  c_enablelinearcombination_openmp nv'
    (Nvector.has_linearcombination nv);
  c_enablescaleaddmulti_openmp nv'
    (Nvector.has_scaleaddmulti nv);
  c_enabledotprodmulti_openmp nv'
    (Nvector.has_dotprodmulti nv);
  c_enablelinearsumvectorarray_openmp nv'
    (Nvector.has_linearsumvectorarray nv);
  c_enablescalevectorarray_openmp nv'
    (Nvector.has_scalevectorarray nv);
  c_enableconstvectorarray_openmp nv'
    (Nvector.has_constvectorarray nv);
  c_enablewrmsnormvectorarray_openmp nv'
    (Nvector.has_wrmsnormvectorarray nv);
  c_enablewrmsnormmaskvectorarray_openmp nv'
    (Nvector.has_wrmsnormmaskvectorarray nv);
  c_enablescaleaddmultivectorarray_openmp nv'
    (Nvector.has_scaleaddmultivectorarray nv);
  c_enablelinearcombinationvectorarray_openmp nv'
    (Nvector.has_linearcombinationvectorarray nv);
  nv'

let pp fmt v = RealArray.pp fmt (unwrap v)

let make ?with_fused_ops nthreads n iv =
  wrap ?with_fused_ops nthreads (RealArray.make n iv)

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
      -> Nvector.any
    = "sunml_nvec_anywrap_openmp"

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
      nthreads v
    =
      if not Sundials_impl.Versions.has_nvector_get_id
        then raise Config.NotImplementedBySundialsVersion;
      let len = RealArray.length v in
      let check nv =
        match unwrap nv with
        | Nvector.RA ra ->
            len = RealArray.length ra && Nvector.get_id nv = Nvector.OpenMP
        | _ -> false
      in
      let nv = c_any_wrap [%extension_constructor Nvector.RA]
                          nthreads v
                          check (clone nthreads)
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
    let nv' = wrap nthreads (RealArray.copy v) in
    c_enablelinearcombination_openmp nv'
      (Nvector.has_linearcombination nv);
    c_enablescaleaddmulti_openmp nv'
      (Nvector.has_scaleaddmulti nv);
    c_enabledotprodmulti_openmp nv'
      (Nvector.has_dotprodmulti nv);
    c_enablelinearsumvectorarray_openmp nv'
      (Nvector.has_linearsumvectorarray nv);
    c_enablescalevectorarray_openmp nv'
      (Nvector.has_scalevectorarray nv);
    c_enableconstvectorarray_openmp nv'
      (Nvector.has_constvectorarray nv);
    c_enablewrmsnormvectorarray_openmp nv'
      (Nvector.has_wrmsnormvectorarray nv);
    c_enablewrmsnormmaskvectorarray_openmp nv'
      (Nvector.has_wrmsnormmaskvectorarray nv);
    c_enablescaleaddmultivectorarray_openmp nv'
      (Nvector.has_scaleaddmultivectorarray nv);
    c_enablelinearcombinationvectorarray_openmp nv'
      (Nvector.has_linearcombinationvectorarray nv);
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
      nthreads n iv
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
           nthreads (RealArray.make n iv)

  let unwrap nv =
    match Nvector.unwrap nv with
    | Nvector.RA a -> a
    | _ -> raise Nvector.BadGenericType

end (* }}} *)

module Ops = struct (* {{{ *)
  type t = (RealArray.t, kind) Nvector.t

  let clone nv =
    let data = Nvector.unwrap nv in
    wrap (num_threads nv) (RealArray.copy data)

  external linearsum    : float -> t -> float -> t -> t -> unit
    = "sunml_nvec_openmp_linearsum"

  external const        : float -> t -> unit
    = "sunml_nvec_openmp_const"

  external prod         : t -> t -> t -> unit
    = "sunml_nvec_openmp_prod"

  external div          : t -> t -> t -> unit
    = "sunml_nvec_openmp_div"

  external scale        : float -> t -> t -> unit
    = "sunml_nvec_openmp_scale"

  external abs          : t -> t -> unit
    = "sunml_nvec_openmp_abs"

  external inv          : t -> t -> unit
    = "sunml_nvec_openmp_inv"

  external addconst     : t -> float -> t -> unit
    = "sunml_nvec_openmp_addconst"

  external dotprod      : t -> t -> float
    = "sunml_nvec_openmp_dotprod"

  external maxnorm      : t -> float
    = "sunml_nvec_openmp_maxnorm"

  external wrmsnorm     : t -> t -> float
    = "sunml_nvec_openmp_wrmsnorm"

  external wrmsnormmask : t -> t -> t -> float
    = "sunml_nvec_openmp_wrmsnormmask"

  external min          : t -> float
    = "sunml_nvec_openmp_min"

  external wl2norm      : t -> t -> float
    = "sunml_nvec_openmp_wl2norm"

  external l1norm       : t -> float
    = "sunml_nvec_openmp_l1norm"

  external compare      : float -> t -> t -> unit
    = "sunml_nvec_openmp_compare"

  external invtest      : t -> t -> bool
    = "sunml_nvec_openmp_invtest"

  external constrmask   : t -> t -> t -> bool
    = "sunml_nvec_openmp_constrmask"

  external minquotient  : t -> t -> float
    = "sunml_nvec_openmp_minquotient"

  external space  : t -> int * int
    = "sunml_nvec_openmp_space"

  external getlength  : t -> int
    = "sunml_nvec_openmp_getlength"

  external linearcombination : RealArray.t -> t array -> t -> unit
    = "sunml_nvec_openmp_linearcombination"

  external scaleaddmulti : RealArray.t -> t -> t array -> t array -> unit
    = "sunml_nvec_openmp_scaleaddmulti"

  external dotprodmulti : t -> t array -> RealArray.t -> unit
    = "sunml_nvec_openmp_dotprodmulti"

  external linearsumvectorarray
    : float -> t array -> float -> t array -> t array -> unit
    = "sunml_nvec_openmp_linearsumvectorarray"

  external scalevectorarray
    : RealArray.t -> t array -> t array -> unit
    = "sunml_nvec_openmp_scalevectorarray"

  external constvectorarray
    : float -> t array -> unit
    = "sunml_nvec_openmp_constvectorarray"

  external wrmsnormvectorarray
    : t array -> t array -> RealArray.t -> unit
    = "sunml_nvec_openmp_wrmsnormvectorarray"

  external wrmsnormmaskvectorarray
    : t array -> t array -> t -> RealArray.t -> unit
    = "sunml_nvec_openmp_wrmsnormmaskvectorarray"

  external scaleaddmultivectorarray
    : RealArray.t -> t array -> t array array -> t array array -> unit
    = "sunml_nvec_openmp_scaleaddmultivectorarray"

  external linearcombinationvectorarray
    : RealArray.t -> t array array -> t array -> unit
    = "sunml_nvec_openmp_linearcombinationvectorarray"

  module Local = struct
    let dotprod     = dotprod
    let maxnorm     = maxnorm
    let min         = min
    let l1norm      = l1norm
    let invtest     = invtest
    let constrmask  = constrmask
    let minquotient = minquotient

    external wsqrsum
      : t -> t -> float
      = "sunml_nvec_openmp_wsqrsumlocal"

    external wsqrsummask
      : t -> t -> t -> float
      = "sunml_nvec_openmp_wsqrsummasklocal"
  end
end (* }}} *)

