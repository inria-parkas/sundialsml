open Sundials
type data = RealArray.t
type kind = [`Pthreads|Nvector_serial.kind]
type t = (data, kind) Nvector.t

(* Selectively enable and disable fused and array operations *)
external c_enablefusedops_pthreads                     : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_pthreads_enablefusedops"
external c_enablelinearcombination_pthreads            : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_pthreads_enablelinearcombination"
external c_enablescaleaddmulti_pthreads                : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_pthreads_enablescaleaddmulti"
external c_enabledotprodmulti_pthreads                 : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_pthreads_enabledotprodmulti"
external c_enablelinearsumvectorarray_pthreads         : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_pthreads_enablelinearsumvectorarray"
external c_enablescalevectorarray_pthreads             : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_pthreads_enablescalevectorarray"
external c_enableconstvectorarray_pthreads             : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_pthreads_enableconstvectorarray"
external c_enablewrmsnormvectorarray_pthreads          : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_pthreads_enablewrmsnormvectorarray"
external c_enablewrmsnormmaskvectorarray_pthreads      : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_pthreads_enablewrmsnormmaskvectorarray"
external c_enablescaleaddmultivectorarray_pthreads     : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_pthreads_enablescaleaddmultivectorarray"
external c_enablelinearcombinationvectorarray_pthreads : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_pthreads_enablelinearcombinationvectorarray"

let unwrap = Nvector.unwrap

external c_wrap : int -> RealArray.t -> (t -> bool) -> (t -> t) -> t
  = "sunml_nvec_wrap_pthreads"

let rec wrap ?(with_fused_ops=false) nthreads v =
  let len = RealArray.length v in
  let check nv' = (len = RealArray.length (unwrap nv')) in
  let nv = c_wrap nthreads v check (clone nthreads) in
  if with_fused_ops then c_enablefusedops_pthreads nv true;
  nv

and clone nthreads nv =
  let nv' = wrap nthreads (RealArray.copy (unwrap nv)) in
  c_enablelinearcombination_pthreads nv'
    (Nvector.has_n_vlinearcombination nv);
  c_enablescaleaddmulti_pthreads nv'
    (Nvector.has_n_vscaleaddmulti nv);
  c_enabledotprodmulti_pthreads nv'
    (Nvector.has_n_vdotprodmulti nv);
  c_enablelinearsumvectorarray_pthreads nv'
    (Nvector.has_n_vlinearsumvectorarray nv);
  c_enablescalevectorarray_pthreads nv'
    (Nvector.has_n_vscalevectorarray nv);
  c_enableconstvectorarray_pthreads nv'
    (Nvector.has_n_vconstvectorarray nv);
  c_enablewrmsnormvectorarray_pthreads nv'
    (Nvector.has_n_vwrmsnormvectorarray nv);
  c_enablewrmsnormmaskvectorarray_pthreads nv'
    (Nvector.has_n_vwrmsnormmaskvectorarray nv);
  c_enablescaleaddmultivectorarray_pthreads nv'
    (Nvector.has_n_vscaleaddmultivectorarray nv);
  c_enablelinearcombinationvectorarray_pthreads nv'
    (Nvector.has_n_vlinearcombinationvectorarray nv);
  nv'

let pp fmt v = RealArray.pp fmt (unwrap v)

let make ?with_fused_ops nthreads n iv =
  wrap ?with_fused_ops nthreads (RealArray.make n iv)

external num_threads : t -> int
  = "sunml_nvec_pthreads_num_threads"

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
  = do_enable c_enablefusedops_pthreads nv
              with_fused_ops;
    do_enable c_enablelinearcombination_pthreads nv
              with_linear_combination;
    do_enable c_enablescaleaddmulti_pthreads nv
              with_scale_add_multi;
    do_enable c_enabledotprodmulti_pthreads nv
              with_dot_prod_multi;
    do_enable c_enablelinearsumvectorarray_pthreads nv
              with_linear_sum_vector_array;
    do_enable c_enablescalevectorarray_pthreads nv
              with_scale_vector_array;
    do_enable c_enableconstvectorarray_pthreads nv
              with_const_vector_array;
    do_enable c_enablewrmsnormvectorarray_pthreads nv
              with_wrms_norm_vector_array;
    do_enable c_enablewrmsnormmaskvectorarray_pthreads nv
              with_wrms_norm_mask_vector_array;
    do_enable c_enablescaleaddmultivectorarray_pthreads nv
              with_scale_add_multi_vector_array;
    do_enable c_enablelinearcombinationvectorarray_pthreads nv
              with_linear_combination_vector_array

module Any = struct (* {{{ *)

  external c_any_wrap
    : extension_constructor
      -> int
      -> RealArray.t
      -> (Nvector.any -> bool)
      -> (Nvector.any -> Nvector.any)
      -> Nvector.any
    = "sunml_nvec_anywrap_pthreads"

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
            len = RealArray.length ra && Nvector.get_id nv = Nvector.Pthreads
        | _ -> false
      in
      let nv = c_any_wrap [%extension_constructor Nvector.RA]
                          nthreads v
                          check (clone nthreads)
      in
      if with_fused_ops
        then c_enablefusedops_pthreads nv true;
      if with_fused_ops
        then c_enablefusedops_pthreads nv true;
      if with_linear_combination
        then c_enablelinearcombination_pthreads nv true;
      if with_scale_add_multi
        then c_enablescaleaddmulti_pthreads nv true;
      if with_dot_prod_multi
        then c_enabledotprodmulti_pthreads nv true;
      if with_linear_sum_vector_array
        then c_enablelinearsumvectorarray_pthreads nv true;
      if with_scale_vector_array
        then c_enablescalevectorarray_pthreads nv true;
      if with_const_vector_array
        then c_enableconstvectorarray_pthreads nv true;
      if with_wrms_norm_vector_array
        then c_enablewrmsnormvectorarray_pthreads nv true;
      if with_wrms_norm_mask_vector_array
        then c_enablewrmsnormmaskvectorarray_pthreads nv true;
      if with_scale_add_multi_vector_array
        then c_enablescaleaddmultivectorarray_pthreads nv true;
      if with_linear_combination_vector_array
        then c_enablelinearcombinationvectorarray_pthreads nv true;
      nv

  and clone nthreads nv =
    let v = match unwrap nv with
            | Nvector.RA v -> v
            | _ -> assert false
    in
    let nv' = wrap nthreads (RealArray.copy v) in
    c_enablelinearcombination_pthreads nv'
      (Nvector.has_n_vlinearcombination nv);
    c_enablescaleaddmulti_pthreads nv'
      (Nvector.has_n_vscaleaddmulti nv);
    c_enabledotprodmulti_pthreads nv'
      (Nvector.has_n_vdotprodmulti nv);
    c_enablelinearsumvectorarray_pthreads nv'
      (Nvector.has_n_vlinearsumvectorarray nv);
    c_enablescalevectorarray_pthreads nv'
      (Nvector.has_n_vscalevectorarray nv);
    c_enableconstvectorarray_pthreads nv'
      (Nvector.has_n_vconstvectorarray nv);
    c_enablewrmsnormvectorarray_pthreads nv'
      (Nvector.has_n_vwrmsnormvectorarray nv);
    c_enablewrmsnormmaskvectorarray_pthreads nv'
      (Nvector.has_n_vwrmsnormmaskvectorarray nv);
    c_enablescaleaddmultivectorarray_pthreads nv'
      (Nvector.has_n_vscaleaddmultivectorarray nv);
    c_enablelinearcombinationvectorarray_pthreads nv'
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

end (* }}} *)

module Ops = struct (* {{{ *)
  type t = (RealArray.t, kind) Nvector.t

  let n_vclone nv =
    let data = Nvector.unwrap nv in
    wrap (num_threads nv) (RealArray.copy data)

  external n_vlinearsum    : float -> t -> float -> t -> t -> unit
    = "sunml_nvec_pthreads_n_vlinearsum"

  external n_vconst        : float -> t -> unit
    = "sunml_nvec_pthreads_n_vconst"

  external n_vprod         : t -> t -> t -> unit
    = "sunml_nvec_pthreads_n_vprod"

  external n_vdiv          : t -> t -> t -> unit
    = "sunml_nvec_pthreads_n_vdiv"

  external n_vscale        : float -> t -> t -> unit
    = "sunml_nvec_pthreads_n_vscale"

  external n_vabs          : t -> t -> unit
    = "sunml_nvec_pthreads_n_vabs"

  external n_vinv          : t -> t -> unit
    = "sunml_nvec_pthreads_n_vinv"

  external n_vaddconst     : t -> float -> t -> unit
    = "sunml_nvec_pthreads_n_vaddconst"

  external n_vdotprod      : t -> t -> float
    = "sunml_nvec_pthreads_n_vdotprod"

  external n_vmaxnorm      : t -> float
    = "sunml_nvec_pthreads_n_vmaxnorm"

  external n_vwrmsnorm     : t -> t -> float
    = "sunml_nvec_pthreads_n_vwrmsnorm"

  external n_vwrmsnormmask : t -> t -> t -> float
    = "sunml_nvec_pthreads_n_vwrmsnormmask"

  external n_vmin          : t -> float
    = "sunml_nvec_pthreads_n_vmin"

  external n_vwl2norm      : t -> t -> float
    = "sunml_nvec_pthreads_n_vwl2norm"

  external n_vl1norm       : t -> float
    = "sunml_nvec_pthreads_n_vl1norm"

  external n_vcompare      : float -> t -> t -> unit
    = "sunml_nvec_pthreads_n_vcompare"

  external n_vinvtest      : t -> t -> bool
    = "sunml_nvec_pthreads_n_vinvtest"

  external n_vconstrmask   : t -> t -> t -> bool
    = "sunml_nvec_pthreads_n_vconstrmask"

  external n_vminquotient  : t -> t -> float
    = "sunml_nvec_pthreads_n_vminquotient"

  external n_vspace  : t -> int * int
    = "sunml_nvec_pthreads_n_vspace"

  external n_vgetlength : t -> int
    = "sunml_nvec_pthreads_n_vgetlength"

  external n_vlinearcombination : RealArray.t -> t array -> t -> unit
    = "sunml_nvec_pthreads_n_vlinearcombination"

  external n_vscaleaddmulti : RealArray.t -> t -> t array -> t array -> unit
    = "sunml_nvec_pthreads_n_vscaleaddmulti"

  external n_vdotprodmulti : t -> t array -> RealArray.t -> unit
    = "sunml_nvec_pthreads_n_vdotprodmulti"

  external n_vlinearsumvectorarray
    : float -> t array -> float -> t array -> t array -> unit
    = "sunml_nvec_pthreads_n_vlinearsumvectorarray"

  external n_vscalevectorarray
    : RealArray.t -> t array -> t array -> unit
    = "sunml_nvec_pthreads_n_vscalevectorarray"

  external n_vconstvectorarray
    : float -> t array -> unit
    = "sunml_nvec_pthreads_n_vconstvectorarray"

  external n_vwrmsnormvectorarray
    : t array -> t array -> RealArray.t -> unit
    = "sunml_nvec_pthreads_n_vwrmsnormvectorarray"

  external n_vwrmsnormmaskvectorarray
    : t array -> t array -> t -> RealArray.t -> unit
    = "sunml_nvec_pthreads_n_vwrmsnormmaskvectorarray"

  external n_vscaleaddmultivectorarray
    : RealArray.t -> t array -> t array array -> t array array -> unit
    = "sunml_nvec_pthreads_n_vscaleaddmultivectorarray"

  external n_vlinearcombinationvectorarray
    : RealArray.t -> t array array -> t array -> unit
    = "sunml_nvec_pthreads_n_vlinearcombinationvectorarray"

  module Local = struct
    let n_vdotprod     = n_vdotprod
    let n_vmaxnorm     = n_vmaxnorm
    let n_vmin         = n_vmin
    let n_vl1norm      = n_vl1norm
    let n_vinvtest     = n_vinvtest
    let n_vconstrmask  = n_vconstrmask
    let n_vminquotient = n_vminquotient

    external n_vwsqrsum
      : t -> t -> float
      = "sunml_nvec_pthreads_n_vwsqrsumlocal"

    external n_vwsqrsummask
      : t -> t -> t -> float
      = "sunml_nvec_pthreads_n_vwsqrsummasklocal"
  end
end (* }}} *)

