open Sundials

type kind
type data = RealArray.t * int * Mpi.communicator
type t = (data, kind) Nvector.t

type Nvector.gdata += Par of data

exception IncorrectGlobalSize

external c_wrap
  : (RealArray.t * int * Mpi.communicator) -> (t -> bool) -> (t -> t) -> t
  = "sunml_nvec_wrap_parallel"

(* Selectively enable and disable fused and array operations *)
external c_enablefusedops_parallel                     : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_par_enablefusedops"
external c_enablelinearcombination_parallel            : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_par_enablelinearcombination"
external c_enablescaleaddmulti_parallel                : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_par_enablescaleaddmulti"
external c_enabledotprodmulti_parallel                 : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_par_enabledotprodmulti"
external c_enablelinearsumvectorarray_parallel         : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_par_enablelinearsumvectorarray"
external c_enablescalevectorarray_parallel             : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_par_enablescalevectorarray"
external c_enableconstvectorarray_parallel             : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_par_enableconstvectorarray"
external c_enablewrmsnormvectorarray_parallel          : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_par_enablewrmsnormvectorarray"
external c_enablewrmsnormmaskvectorarray_parallel      : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_par_enablewrmsnormmaskvectorarray"
external c_enablescaleaddmultivectorarray_parallel     : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_par_enablescaleaddmultivectorarray"
external c_enablelinearcombinationvectorarray_parallel : ('d, 'k) Nvector.t -> bool -> unit
  = "sunml_nvec_par_enablelinearcombinationvectorarray"

let unwrap = Nvector.unwrap

let copydata (a, ng, comm) = (RealArray.copy a, ng, comm)

let rec wrap ?(with_fused_ops=false) ((nl, ng, comm) as v) =
  let nl_len = RealArray.length nl in
  let check nv =
    let (nl', ng', comm') = unwrap nv in
    (nl_len = RealArray.length nl') && (ng <= ng') && (comm == comm')
  in
  let nv = c_wrap v check clone in
  if with_fused_ops then c_enablefusedops_parallel nv true;
  nv

and clone nv =
  let nv' = wrap (copydata (unwrap nv)) in
  c_enablelinearcombination_parallel nv'
    (Nvector.has_linearcombination nv);
  c_enablescaleaddmulti_parallel nv'
    (Nvector.has_scaleaddmulti nv);
  c_enabledotprodmulti_parallel nv'
    (Nvector.has_dotprodmulti nv);
  c_enablelinearsumvectorarray_parallel nv'
    (Nvector.has_linearsumvectorarray nv);
  c_enablescalevectorarray_parallel nv'
    (Nvector.has_scalevectorarray nv);
  c_enableconstvectorarray_parallel nv'
    (Nvector.has_constvectorarray nv);
  c_enablewrmsnormvectorarray_parallel nv'
    (Nvector.has_wrmsnormvectorarray nv);
  c_enablewrmsnormmaskvectorarray_parallel nv'
    (Nvector.has_wrmsnormmaskvectorarray nv);
  c_enablescaleaddmultivectorarray_parallel nv'
    (Nvector.has_scaleaddmultivectorarray nv);
  c_enablelinearcombinationvectorarray_parallel nv'
    (Nvector.has_linearcombinationvectorarray nv);
  nv'

let make ?with_fused_ops nl ng comm iv =
  wrap ?with_fused_ops (RealArray.make nl iv, ng, comm)

let clone nv =
  let loc, glen, comm = Nvector.unwrap nv in
  wrap (RealArray.copy loc, glen, comm)

let pp fmt nv =
  let data, _, _ = Nvector.unwrap nv in
  RealArray.pp fmt data

let local_array nv =
  let data, _, _ = Nvector.unwrap nv in
  data

let local_length nv = RealArray.length (local_array nv)

let global_length nv =
  let _, gl, _ = Nvector.unwrap nv in
  gl

let communicator nv =
  let _, _, comm = Nvector.unwrap nv in
  comm

external get_communicator : ('d, 'k) Nvector.t -> Mpi.communicator option
  = "sunml_nvec_par_getcommunicator"

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
  = do_enable c_enablefusedops_parallel nv
              with_fused_ops;
    do_enable c_enablelinearcombination_parallel nv
              with_linear_combination;
    do_enable c_enablescaleaddmulti_parallel nv
              with_scale_add_multi;
    do_enable c_enabledotprodmulti_parallel nv
              with_dot_prod_multi;
    do_enable c_enablelinearsumvectorarray_parallel nv
              with_linear_sum_vector_array;
    do_enable c_enablescalevectorarray_parallel nv
              with_scale_vector_array;
    do_enable c_enableconstvectorarray_parallel nv
              with_const_vector_array;
    do_enable c_enablewrmsnormvectorarray_parallel nv
              with_wrms_norm_vector_array;
    do_enable c_enablewrmsnormmaskvectorarray_parallel nv
              with_wrms_norm_mask_vector_array;
    do_enable c_enablescaleaddmultivectorarray_parallel nv
              with_scale_add_multi_vector_array;
    do_enable c_enablelinearcombinationvectorarray_parallel nv
              with_linear_combination_vector_array

external hide_communicator
  : Mpi.communicator -> Nvector_custom.communicator
  = "%identity"

module Any = struct (* {{{ *)

  external c_any_wrap
    : extension_constructor
      -> data
      -> (Nvector.any -> bool)
      -> (Nvector.any -> Nvector.any)
      -> Nvector.any
    = "sunml_nvec_anywrap_parallel"

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
      ((nl, ng, comm) as v)
    =
      if not Sundials_impl.Version.has_nvector_get_id
        then raise Config.NotImplementedBySundialsVersion;
      let len = RealArray.length nl in
      let check nv' =
        match unwrap nv' with
        | Par (nl', ng', comm') ->
            (len = RealArray.length nl') && (ng <= ng') && (comm == comm')
            && (Nvector.get_id nv' = Nvector.Parallel)
        | _ -> false
      in
      let nv = c_any_wrap [%extension_constructor Par] v check clone in
      if with_fused_ops
        then c_enablefusedops_parallel nv true;
      if with_fused_ops
        then c_enablefusedops_parallel nv true;
      if with_linear_combination
        then c_enablelinearcombination_parallel nv true;
      if with_scale_add_multi
        then c_enablescaleaddmulti_parallel nv true;
      if with_dot_prod_multi
        then c_enabledotprodmulti_parallel nv true;
      if with_linear_sum_vector_array
        then c_enablelinearsumvectorarray_parallel nv true;
      if with_scale_vector_array
        then c_enablescalevectorarray_parallel nv true;
      if with_const_vector_array
        then c_enableconstvectorarray_parallel nv true;
      if with_wrms_norm_vector_array
        then c_enablewrmsnormvectorarray_parallel nv true;
      if with_wrms_norm_mask_vector_array
        then c_enablewrmsnormmaskvectorarray_parallel nv true;
      if with_scale_add_multi_vector_array
        then c_enablescaleaddmultivectorarray_parallel nv true;
      if with_linear_combination_vector_array
        then c_enablelinearcombinationvectorarray_parallel nv true;
      nv

  and clone nv =
    let v = match unwrap nv with
            | Par v -> v
            | _ -> assert false
    in
    let nv' = wrap (copydata v) in
    c_enablelinearcombination_parallel nv'
      (Nvector.has_linearcombination nv);
    c_enablescaleaddmulti_parallel nv'
      (Nvector.has_scaleaddmulti nv);
    c_enabledotprodmulti_parallel nv'
      (Nvector.has_dotprodmulti nv);
    c_enablelinearsumvectorarray_parallel nv'
      (Nvector.has_linearsumvectorarray nv);
    c_enablescalevectorarray_parallel nv'
      (Nvector.has_scalevectorarray nv);
    c_enableconstvectorarray_parallel nv'
      (Nvector.has_constvectorarray nv);
    c_enablewrmsnormvectorarray_parallel nv'
      (Nvector.has_wrmsnormvectorarray nv);
    c_enablewrmsnormmaskvectorarray_parallel nv'
      (Nvector.has_wrmsnormmaskvectorarray nv);
    c_enablescaleaddmultivectorarray_parallel nv'
      (Nvector.has_scaleaddmultivectorarray nv);
    c_enablelinearcombinationvectorarray_parallel nv'
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
      nl ng comm iv
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
           (RealArray.make nl iv, ng, comm)

  let unwrap nv =
    match Nvector.unwrap nv with
    | Par a -> a
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
      if Nvector.get_id nv <> Nvector.Parallel then raise Nvector.BadGenericType;
      do_enable c_enablefusedops_parallel nv
                with_fused_ops;
      do_enable c_enablelinearcombination_parallel nv
                with_linear_combination;
      do_enable c_enablescaleaddmulti_parallel nv
                with_scale_add_multi;
      do_enable c_enabledotprodmulti_parallel nv
                with_dot_prod_multi;
      do_enable c_enablelinearsumvectorarray_parallel nv
                with_linear_sum_vector_array;
      do_enable c_enablescalevectorarray_parallel nv
                with_scale_vector_array;
      do_enable c_enableconstvectorarray_parallel nv
                with_const_vector_array;
      do_enable c_enablewrmsnormvectorarray_parallel nv
                with_wrms_norm_vector_array;
      do_enable c_enablewrmsnormmaskvectorarray_parallel nv
                with_wrms_norm_mask_vector_array;
      do_enable c_enablescaleaddmultivectorarray_parallel nv
                with_scale_add_multi_vector_array;
      do_enable c_enablelinearcombinationvectorarray_parallel nv
                with_linear_combination_vector_array

end (* }}} *)

module Ops = struct (* {{{ *)
  type t = (data, kind) Nvector.t
  let check = Nvector.check

  let clone = clone

  external c_linearsum    : float -> t -> float -> t -> t -> unit
    = "sunml_nvec_par_linearsum" [@@noalloc]

  let linearsum a (x : t) b (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_linearsum a x b y z

  external const        : float -> t -> unit
    = "sunml_nvec_par_const" [@@noalloc]

  external c_prod         : t -> t -> t -> unit
    = "sunml_nvec_par_prod" [@@noalloc]

  let prod (x : t) (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_prod x y z

  external c_div          : t -> t -> t -> unit
    = "sunml_nvec_par_div" [@@noalloc]

  let div (x : t) (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_div x y z

  external c_scale        : float -> t -> t -> unit
    = "sunml_nvec_par_scale" [@@noalloc]

  let scale c (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_scale c x z

  external c_abs          : t -> t -> unit
    = "sunml_nvec_par_abs" [@@noalloc]

  let abs (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_abs x z

  external c_inv          : t -> t -> unit
    = "sunml_nvec_par_inv" [@@noalloc]

  let inv (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_inv x z

  external c_addconst     : t -> float -> t -> unit
    = "sunml_nvec_par_addconst" [@@noalloc]

  let addconst (x : t) b (z : t) =
    if Sundials_configuration.safe then check x z;
    c_addconst x b z

  external c_dotprod      : t -> t -> float
    = "sunml_nvec_par_dotprod"

  let dotprod (x : t) (y : t) =
    if Sundials_configuration.safe then check x y;
    c_dotprod x y

  external maxnorm      : t -> float
    = "sunml_nvec_par_maxnorm"

  external c_wrmsnorm     : t -> t -> float
    = "sunml_nvec_par_wrmsnorm"

  let wrmsnorm (x : t) (w : t) =
    if Sundials_configuration.safe then check x w;
    c_wrmsnorm x w

  external c_wrmsnormmask : t -> t -> t -> float
    = "sunml_nvec_par_wrmsnormmask"

  let wrmsnormmask (x : t) (w : t) (id : t) =
    if Sundials_configuration.safe then (check x w; check x id);
    c_wrmsnormmask x w id

  external min          : t -> float
    = "sunml_nvec_par_min"

  external c_wl2norm      : t -> t -> float
    = "sunml_nvec_par_wl2norm"

  let wl2norm (x : t) (w : t) =
    if Sundials_configuration.safe then check x w;
    c_wl2norm x w

  external l1norm       : t -> float
    = "sunml_nvec_par_l1norm"

  external c_compare      : float -> t -> t -> unit
    = "sunml_nvec_par_compare"

  let compare c (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_compare c x z

  external c_invtest      : t -> t -> bool
    = "sunml_nvec_par_invtest"

  let invtest (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_invtest x z

  external c_constrmask   : t -> t -> t -> bool
    = "sunml_nvec_par_constrmask"

  let constrmask (c : t) (x : t) (m : t) =
    if Sundials_configuration.safe then (check c x; check c m);
    c_constrmask c x m

  external c_minquotient  : t -> t -> float
    = "sunml_nvec_par_minquotient"

  let minquotient (n : t) (d : t) =
    if Sundials_configuration.safe then check n d;
    c_minquotient n d

  external space  : t -> int * int
    = "sunml_nvec_par_space"

  external getlength  : t -> int
    = "sunml_nvec_par_getlength"

  external c_print_file : t -> Logfile.t option -> unit
    = "sunml_nvec_par_print_file"

  let print ?logfile nv = c_print_file nv logfile

  external c_linearcombination : RealArray.t -> t array -> t -> unit
    = "sunml_nvec_par_linearcombination"

  let linearcombination ca (xa : t array) (z : t) =
    if Sundials_impl.Version.lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe then Array.iter (check z) xa;
    c_linearcombination ca xa z

  let same_len' n ya =
    if n <> Array.length ya then invalid_arg "arrays of unequal length"
  let same_len xa ya = same_len' (Array.length xa) ya

  external c_scaleaddmulti : RealArray.t -> t -> t array -> t array -> unit
    = "sunml_nvec_par_scaleaddmulti"

  let scaleaddmulti aa (x : t) (ya : t array) (za : t array) =
    if Sundials_impl.Version.lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe then
      (Array.iter (check x) ya; Array.iter (check x) za;
       let nv = RealArray.length aa in
       same_len' nv ya; same_len' nv za);
    c_scaleaddmulti aa x ya za

  external c_dotprodmulti : t -> t array -> RealArray.t -> unit
    = "sunml_nvec_par_dotprodmulti"

  let dotprodmulti (x : t) (ya : t array) (dp : RealArray.t) =
    if Sundials_impl.Version.lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe then
      (let nv = RealArray.length dp in
       same_len' nv ya;
       Array.iter (check x) ya);
    c_dotprodmulti x ya dp

  external c_linearsumvectorarray
    : float -> t array -> float -> t array -> t array -> unit
    = "sunml_nvec_par_linearsumvectorarray"

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
    = "sunml_nvec_par_scalevectorarray"

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
    = "sunml_nvec_par_constvectorarray"

  let constvectorarray c (za : t array) =
    if Sundials_impl.Version.lt400
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_configuration.safe
    then (let z = Array.get za 0 in
          Array.iter (check z) za);
    c_constvectorarray c za

  external c_wrmsnormvectorarray
    : t array -> t array -> RealArray.t -> unit
    = "sunml_nvec_par_wrmsnormvectorarray"

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
    = "sunml_nvec_par_wrmsnormmaskvectorarray"

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
    = "sunml_nvec_par_scaleaddmultivectorarray"

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
    = "sunml_nvec_par_linearcombinationvectorarray"

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

    external c_dotprod : t -> t -> float
      = "sunml_nvec_par_dotprodlocal"

    let dotprod (x : t) (y : t) =
      if Sundials_impl.Version.lt500
        then raise Config.NotImplementedBySundialsVersion;
      if Sundials_configuration.safe then check x y;
      c_dotprod x y

    external c_maxnorm : t -> float
      = "sunml_nvec_par_maxnormlocal"

    let maxnorm (x : t) =
      if Sundials_impl.Version.lt500
        then raise Config.NotImplementedBySundialsVersion;
      c_maxnorm x

    external c_min : t -> float
      = "sunml_nvec_par_minlocal"

    let min (x : t) =
      if Sundials_impl.Version.lt500
        then raise Config.NotImplementedBySundialsVersion;
      c_min x

    external c_l1norm : t -> float
      = "sunml_nvec_par_l1normlocal"

    let l1norm (x : t) =
      if Sundials_impl.Version.lt500
        then raise Config.NotImplementedBySundialsVersion;
      c_l1norm x

    external c_invtest : t -> t -> bool
      = "sunml_nvec_par_invtestlocal" [@@noalloc]

    let invtest (x : t) (z : t) =
      if Sundials_impl.Version.lt500
        then raise Config.NotImplementedBySundialsVersion;
      if Sundials_configuration.safe then check x z;
      c_invtest x z

    external c_constrmask : t -> t -> t -> bool
      = "sunml_nvec_par_constrmasklocal" [@@noalloc]

    let constrmask (c : t) (x : t) (m : t) =
      if Sundials_impl.Version.lt500
        then raise Config.NotImplementedBySundialsVersion;
      if Sundials_configuration.safe then (check c x; check c m);
      c_constrmask c x m

    external c_minquotient : t -> t -> float
      = "sunml_nvec_par_minquotientlocal"

    let minquotient (n : t) (d : t) =
      if Sundials_impl.Version.lt500
        then raise Config.NotImplementedBySundialsVersion;
      if Sundials_configuration.safe then check n d;
      c_minquotient n d

    external c_wsqrsum : t -> t -> float
      = "sunml_nvec_par_wsqrsumlocal"

    let wsqrsum (x : t) (w : t) =
      if Sundials_impl.Version.lt500
        then raise Config.NotImplementedBySundialsVersion;
      if Sundials_configuration.safe then check x w;
      c_wsqrsum x w

    external c_wsqrsummask
      : t -> t -> t -> float
      = "sunml_nvec_par_wsqrsummasklocal"

    let wsqrsummask (x : t) (w : t) (id : t) =
      if Sundials_impl.Version.lt500
        then raise Config.NotImplementedBySundialsVersion;
      if Sundials_configuration.safe then (check x w; check x id);
      c_wsqrsummask x w id
  end
end (* }}} *)

module MakeOps =
  functor (A : sig
      type local_data
      val get       : local_data -> int -> float
      val set       : local_data -> int -> float -> unit
      val fill      : local_data -> float -> unit
      val make      : int -> float -> local_data
      val clone     : local_data -> local_data
      val length    : local_data -> int
    end) ->
  struct (* {{{ *)
    type t = A.local_data * int * Mpi.communicator

    (* let clone (d, gl, comm) = (A.clone d, gl, comm) *)
    let clone (d, gl, comm) = (A.clone d, gl, comm)

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

    let linearsum a (x, _, _) b (y, _, _) (z, _, _) =
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

    let const c (a, _, _) = A.fill a c

    let scale c (x, _, _) (z, _, _) =
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

    let addconst (x, _, _) b (z, _, _) =
      for i = 0 to A.length x - 1 do
        A.set z i (A.get x i +. b)
      done

    let maxnorm (x, _, comm) =
      let lmax = ref 0.0 in
      for i = 0 to A.length x - 1 do
        let ax = abs_float (A.get x i) in
        if ax > !lmax then lmax := ax
      done;
      Mpi.allreduce_float !lmax Mpi.Max comm

    let wrmsnorm (x, n_global, comm) (w, _, _) =
      let lsum = ref 0.0 in
      for i = 0 to A.length x - 1 do
        lsum := !lsum
                  +. ((A.get x i) *. (A.get w i) *. (A.get x i) *. (A.get w i))
      done;
      let gsum = Mpi.allreduce_float !lsum Mpi.Sum comm in
      sqrt (gsum /. float n_global)

    let wrmsnormmask (x, n_global, comm) (w, _, _) (id, _, _) =
      let lsum = ref 0.0 in
      for i = 0 to A.length x - 1 do
        if A.get id i > 0.0 then
          lsum := !lsum +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
      done;
      let gsum = Mpi.allreduce_float !lsum Mpi.Sum comm in
      sqrt (gsum /. float n_global)

    let min (x, _, comm) =
      let lmin = ref max_float in
      for i = 0 to A.length x - 1 do
        let xv = A.get x i in
        if xv < !lmin then lmin := xv
      done;
      Mpi.allreduce_float !lmin Mpi.Min comm

    let dotprod (x, _, comm) (y, _, _) =
      let lsum = ref 0.0 in
      for i = 0 to A.length x - 1 do
        lsum := !lsum +. (A.get x i *. A.get y i)
      done;
      Mpi.allreduce_float !lsum Mpi.Sum comm

    let compare c (x, _, _) (z, _, _) =
      for i = 0 to A.length x - 1 do
        A.set z i (if abs_float (A.get x i) >= c then 1.0 else 0.0)
      done

    let invtest (x, _, comm) (z, _, _) =
      let lval = ref 1 in
      for i = 0 to A.length x - 1 do
        if A.get x i = 0.0 then lval := 0
        else A.set z i (1.0 /. (A.get x i))
      done;
      (Mpi.allreduce_int !lval Mpi.Min comm <> 0)

    let wl2norm (x, _, comm) (w, _, _) =
      let lsum = ref 0.0 in
      for i = 0 to A.length x - 1 do
        lsum := !lsum +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
      done;
      let gsum = Mpi.allreduce_float !lsum Mpi.Sum comm in
      sqrt gsum

    let l1norm (x, _, comm) =
      let lsum = ref 0.0 in
      for i = 0 to A.length x - 1 do
        lsum := !lsum +. abs_float (A.get x i)
      done;
      Mpi.allreduce_float !lsum Mpi.Sum comm

    let constrmask (c, _, _) (x, _, comm) (m, _, _) =
      let test = ref 1.0 in
      let check b = if b then 0.0 else (test := 0.0; 1.0) in
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
      (Mpi.allreduce_float !test Mpi.Min comm = 1.0)

    let minquotient (num, _, comm) (denom, _, _) =
      let lmin = ref Config.big_real in
      for i = 0 to A.length num - 1 do
        if (A.get denom i) <> 0.0 then
          lmin := floatmin !lmin (A.get num i /. A.get denom i)
      done;
      Mpi.allreduce_float !lmin Mpi.Min comm

    let prod (x, _, _) (y, _, _) (z, _, _) =
      for i = 0 to A.length x - 1 do
        A.set z i (A.get x i *. A.get y i)
      done

    let div (x, _, _) (y, _, _) (z, _, _) =
      for i = 0 to A.length x - 1 do
        A.set z i (A.get x i /. A.get y i)
      done

    let abs (x, _, _) (z, _, _) =
      for i = 0 to A.length x - 1 do
        A.set z i (abs_float (A.get x i))
      done

    let inv (x, _, _) (z, _, _) =
      for i = 0 to A.length x - 1 do
        A.set z i (1.0 /. (A.get x i))
      done

    let space (_, ng, comm) = (ng, 2 * Mpi.comm_size comm)

    let getlength (_, ng, _) = ng

    let print ?(logfile=Logfile.stdout) (x, _, _) =
      for i = 0 to A.length x - 1 do
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
        let zd, _, _ = z in
        let n = A.length zd in
        if xa.(0) == z then begin
          let c0 = ca.{0} in
          if c0 <> 1.0 then
            for j = 0 to n - 1 do
              A.set zd j (A.get zd j *. c0)
            done;
          for i = 1 to nvec - 1 do
            let ci, (xd, _, _) = ca.{i}, xa.(i) in
            for j = 0 to n - 1 do
              A.set zd j (A.get zd j +. ci *. A.get xd j)
            done
          done
        end
        else begin
          let c0, (xd, _, _) = ca.{0}, xa.(0) in
          for j = 0 to n - 1 do
            A.set zd j (c0 *. A.get xd j)
          done;
          for i = 1 to nvec - 1 do
            let ci, (xd, _, _) = ca.{i}, xa.(i) in
            for j = 0 to n - 1 do
              A.set zd j (A.get zd j +. ci *. A.get xd j)
            done
          done
        end

    let scaleaddmulti (aa : RealArray.t) (x : t) (ya : t array) (za : t array) =
      let nvec = Array.length ya in
      if nvec = 1 then linearsum aa.{0} x 1.0 ya.(0) za.(0)
      else
        let xd, _, _ = x in
        let n = A.length xd in
        if ya == za then
          for i = 0 to nvec - 1 do
            let a, (yd, _, _) = aa.{i}, ya.(i) in
            for j = 0 to n - 1 do
              A.set yd j (A.get yd j +. a *. A.get xd j)
            done
          done
        else
          for i = 0 to nvec - 1 do
            let ai, (yd, _, _), (zd, _, _) = aa.{i}, ya.(i), za.(i) in
            for j = 0 to n - 1 do
              A.set zd j (ai *. A.get xd j +. A.get yd j)
            done
          done

    let dotprodmulti (x : t) (ya : t array) (dp : RealArray.t) =
      let nvec = Array.length ya in
      if nvec = 1 then dp.{0} <- dotprod x ya.(0)
      else
        let xd, _, comm = x in
        let n = A.length xd in
        for i = 0 to nvec - 1 do
          let yd, _, _ = ya.(i) in
          dp.{i} <- 0.0;
          for j = 0 to n - 1 do
            dp.{i} <- dp.{i} +. A.get xd j *. A.get yd j
          done
        done;
        (* Note: ocamlmpi does not provide MPI_IN_PLACE *)
        Mpi.(allreduce_bigarray1 dp dp Sum comm)

    let arr_vaxpy_array a (xa : t array) (ya : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.length xad in
      if a = 1.0 then
        for i = 0 to nvec - 1 do
          let xd, _, _ = xa.(i) in
          let yd, _, _ = ya.(i) in
          for j = 0 to n - 1 do
            A.set yd j (A.get yd j +. A.get xd j)
          done
        done
      else if a = -1.0 then
        for i = 0 to nvec - 1 do
          let xd, _, _ = xa.(i) in
          let yd, _, _ = ya.(i) in
          for j = 0 to n - 1 do
            A.set yd j (A.get yd j -. A.get xd j)
          done
        done
      else
        for i = 0 to nvec - 1 do
          let xd, _, _ = xa.(i) in
          let yd, _, _ = ya.(i) in
          for j = 0 to n - 1 do
            A.set yd j (A.get yd j +. a *. A.get xd j)
          done
        done

    let v_sumvectorarray (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.length xad in
      for i = 0 to nvec - 1 do
        let xd, _, _ = xa.(i) in
        let yd, _, _ = ya.(i) in
        let zd, _, _ = za.(i) in
        for j = 0 to n - 1 do
          A.set zd j (A.get xd j +. A.get yd j)
        done
      done

    let v_diffvectorarray (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.length xad in
      for i = 0 to nvec - 1 do
        let xd, _, _ = xa.(i) in
        let yd, _, _ = ya.(i) in
        let zd, _, _ = za.(i) in
        for j = 0 to n - 1 do
          A.set zd j (A.get xd j -. A.get yd j)
        done
      done

    let v_lin1vectorarray a (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.length xad in
      for i = 0 to nvec - 1 do
        let xd, _, _ = xa.(i) in
        let yd, _, _ = ya.(i) in
        let zd, _, _ = za.(i) in
        for j = 0 to n - 1 do
          A.set zd j (a *. A.get xd j +. A.get yd j)
        done
      done

    let v_lin2vectorarray a (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.length xad in
      for i = 0 to nvec - 1 do
        let xd, _, _ = xa.(i) in
        let yd, _, _ = ya.(i) in
        let zd, _, _ = za.(i) in
        for j = 0 to n - 1 do
          A.set zd j (a *. A.get xd j -. A.get yd j)
        done
      done

    let v_scalesumvectorarray c (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.length xad in
      for i = 0 to nvec - 1 do
        let xd, _, _ = xa.(i) in
        let yd, _, _ = ya.(i) in
        let zd, _, _ = za.(i) in
        for j = 0 to n - 1 do
          A.set zd j (c *. (A.get xd j +. A.get yd j))
        done
      done

    let v_scalediffvectorarray c (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.length xad in
      for i = 0 to nvec - 1 do
        let xd, _, _ = xa.(i) in
        let yd, _, _ = ya.(i) in
        let zd, _, _ = za.(i) in
        for j = 0 to n - 1 do
          A.set zd j (c *. (A.get xd j -. A.get yd j))
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
        let xad, _, _ = xa.(0) in
        let n = A.length xad in
        for i = 0 to nvec - 1 do
          let xd, _, _ = xa.(i) in
          let yd, _, _ = ya.(i) in
          let zd, _, _ = za.(i) in
          for j = 0 to n - 1 do
            A.set zd j (a *. A.get xd j +. b *. A.get yd j)
          done
        done

    let scalevectorarray (c : RealArray.t) (xa : t array) (za : t array) =
      let nvec = Array.length xa in
      if nvec = 1 then scale c.{0} xa.(0) za.(0)
      else
        let xad, _, _ = xa.(0) in
        let n = A.length xad in
        if xa == za then
          for i = 0 to nvec - 1 do
            let (xd, _, _), c = xa.(i), c.{i} in
            for j = 0 to n - 1 do
              A.set xd j (c *. A.get xd j)
            done
          done
        else
          for i = 0 to nvec - 1 do
            let c = c.{i} in
            let xd, _, _ = xa.(i) in
            let zd, _, _ = za.(i) in
            for j = 0 to n - 1 do
              A.set zd j (c *. A.get xd j)
            done
          done

    let constvectorarray c (za : t array) =
      let nvec = Array.length za in
      if nvec = 1 then const c za.(0)
      else
        let zad, _, _ = za.(0) in
        let n = A.length zad in
        for i = 0 to nvec - 1 do
          let zd, _, _ = za.(i) in
          for j = 0 to n - 1 do
            A.set zd j c
          done
        done

    let wrmsnormvectorarray (xa : t array) (wa : t array) (nrm : RealArray.t) =
      let nvec = Array.length xa in
      if nvec = 1 then nrm.{0} <- wrmsnorm xa.(0) wa.(0)
      else
        let xad, ng, comm = xa.(0) in
        let nl = A.length xad in
        let nf = float ng in
        for i = 0 to nvec - 1 do
          let xd, _, _ = xa.(i) in
          let wd, _, _ = wa.(i) in
          nrm.{i} <- 0.0;
          for j = 0 to nl - 1 do
            let s = A.get xd j *. A.get wd j in
            nrm.{i} <- nrm.{i} +. s *. s
          done
        done;
        (* Note: ocamlmpi does not provide MPI_IN_PLACE *)
        Mpi.(allreduce_bigarray1 nrm nrm Sum comm);
        for i = 0 to nvec - 1 do
          nrm.{i} <- sqrt (nrm.{i}/.nf)
        done

    let wrmsnormmaskvectorarray (xa : t array) (wa : t array) (id : t)
                                   (nrm : RealArray.t) =
      let nvec = Array.length xa in
      if nvec = 1 then nrm.{0} <- wrmsnormmask xa.(0) wa.(0) id
      else
        let xad, ng, comm = xa.(0) in
        let idd, _, _ = id in
        let nl = A.length xad in
        let nf = float ng in
        for i = 0 to nvec - 1 do
          let xd, _, _ = xa.(i) in
          let wd, _, _ = wa.(i) in
          nrm.{i} <- 0.0;
          for j = 0 to nl - 1 do
            if A.get idd j > 0.0 then begin
              let s = A.get xd j *. A.get wd j in
              nrm.{i} <- nrm.{i} +. s *. s
            end
          done;
        done;
        (* Note: ocamlmpi does not provide MPI_IN_PLACE *)
        Mpi.(allreduce_bigarray1 nrm nrm Sum comm);
        for i = 0 to nvec - 1 do
          nrm.{i} <- sqrt (nrm.{i}/.nf)
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
        let xad, _, _ = xa.(0) in
        let n = A.length xad in
        if (yaa == zaa) then
          for i = 0 to nvec - 1 do
            let xd, _, _ = xa.(i) in
            for j = 0 to nsum - 1 do
              let a, (yd, _, _) = ra.{j}, yaa.(j).(i) in
              for k = 0 to n - 1 do
                A.set yd k (A.get yd k +. a *. A.get xd k)
              done
            done
          done
        else
          for i = 0 to nvec - 1 do
            let xd, _, _ = xa.(i) in
            for j = 0 to nsum - 1 do
              let a = ra.{j} in
              let yd, _, _ = yaa.(j).(i) in
              let zd, _, _ = zaa.(j).(i) in
              for k = 0 to n - 1 do
                A.set zd k (a *. A.get xd k +. A.get yd k)
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
          let zad, _, _ = za.(0) in
          let n = A.length zad in
          if xaa.(0) == za then begin
            if ca.{0} = 1.0 then
              for j = 0 to nvec - 1 do
                let zd, _, _ = za.(j) in
                for i = 1 to nsum - 1 do
                  let c, (xd, _, _) = ca.{i}, xaa.(i).(j) in
                  for k = 0 to n - 1 do
                    A.set zd k (A.get zd k +. c *. A.get xd k)
                  done
                done
              done
            else
              let c0 = ca.{0} in
              for j = 0 to nvec - 1 do
                let zd, _, _ = za.(j) in
                for k = 0 to n - 1 do
                  A.set zd k (A.get zd k *. c0)
                done;
                for i = 1 to nsum - 1 do
                  let c, (xd, _, _) = ca.{i}, xaa.(i).(j) in
                  for k = 0 to n - 1 do
                    A.set zd k (A.get zd k +. c *. A.get xd k)
                  done
                done
              done
          end
          else
            let c0 = ca.{0} in
            for j = 0 to nvec - 1 do
              let xd, _, _ = xaa.(0).(j) in
              let zd, _, _ = za.(j) in
              for k = 0 to n - 1 do
                A.set zd k (c0 *. A.get xd k)
              done;
              for i = 1 to nsum - 1 do
                let c, (xd, _, _) = ca.{i}, xaa.(i).(j) in
                for k = 0 to n - 1 do
                  A.set zd k (A.get zd k +. c *. A.get xd k)
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

      let wsqrsum (x, _, _) (w, _, _) =
        let a = ref 0.0 in
        let lx = A.length x in
        for i = 0 to lx - 1 do
          a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
        done;
        !a

      let wsqrsummask (x, _, _) (w, _, _) (id, _, _) =
        let a = ref 0.0 in
        let lx = A.length x in
        for i = 0 to lx - 1 do
          if A.get id i > 0.0 then
            a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
        done;
        !a
    end
  end (* }}} *)

(* (* Too slow *)
module SlowerDataOps = MakeOps (struct
    type local_data = RealArray.t

    let get       = Bigarray.Array1.get
    let set       = Bigarray.Array1.set
    let fill      = Bigarray.Array1.fill

    let make      = RealArray.make
    let length    = RealArray.length
    let clone     = RealArray.clone
  end)
*)

module DataOps =
  struct (* {{{ *)
    module A = Bigarray.Array1

    let make      = RealArray.make
    let clone     = RealArray.copy

    let floatmin (x : float) (y : float) = min x y

    type t = RealArray.t * int * Mpi.communicator
    type d = RealArray.t

    let clone (d, gl, comm) = (clone d, gl, comm)

    let arr_vaxpy a (x : d) (y : d) =
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

    let linearsum a ((x : d), _, _) b ((y : d), _, _) ((z : d), _, _) =
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

    let const c ((a : d), _, _) = A.fill a c

    let scale c ((x : d), _, _) ((z : d), _, _) =
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

    let addconst ((x : d), _, _) b ((z : d), _, _) =
      for i = 0 to A.dim x - 1 do
        A.set z i (A.get x i +. b)
      done

    let maxnorm ((x : d), _, comm) =
      let lmax = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        let ax = abs_float (A.get x i) in
        if ax > !lmax then lmax := ax
      done;
      Mpi.allreduce_float !lmax Mpi.Max comm

    let wrmsnorm ((x : d), n_global, comm) ((w : d), _, _) =
      let lsum = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        lsum := !lsum
                  +. ((A.get x i) *. (A.get w i) *. (A.get x i) *. (A.get w i))
      done;
      let gsum = Mpi.allreduce_float !lsum Mpi.Sum comm in
      sqrt (gsum /. float n_global)

    let wrmsnormmask ((x : d), n_global, comm) ((w : d), _, _) ((id : d), _, _) =
      let lsum = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        if A.get id i > 0.0 then
          lsum := !lsum +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
      done;
      let gsum = Mpi.allreduce_float !lsum Mpi.Sum comm in
      sqrt (gsum /. float n_global)

    let min ((x : d), _, comm) =
      let lmin = ref max_float in
      for i = 0 to A.dim x - 1 do
        let xv = A.get x i in
        if xv < !lmin then lmin := xv
      done;
      Mpi.allreduce_float !lmin Mpi.Min comm

    let dotprod ((x : d), _, comm) ((y : d), _, _) =
      let lsum = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        lsum := !lsum +. (A.get x i *. A.get y i)
      done;
      Mpi.allreduce_float !lsum Mpi.Sum comm

    let compare c ((x : d), _, _) ((z : d), _, _) =
      for i = 0 to A.dim x - 1 do
        A.set z i (if abs_float (A.get x i) >= c then 1.0 else 0.0)
      done

    let invtest ((x : d), _, comm) ((z : d), _, _) =
      let lval = ref 1 in
      for i = 0 to A.dim x - 1 do
        if A.get x i = 0.0 then lval := 0
        else A.set z i (1.0 /. (A.get x i))
      done;
      (Mpi.allreduce_int !lval Mpi.Min comm <> 0)

    let wl2norm ((x : d), _, comm) ((w : d), _, _) =
      let lsum = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        lsum := !lsum +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
      done;
      let gsum = Mpi.allreduce_float !lsum Mpi.Sum comm in
      sqrt gsum

    let l1norm ((x : d), _, comm) =
      let lsum = ref 0.0 in
      for i = 0 to A.dim x - 1 do
        lsum := !lsum +. abs_float (A.get x i)
      done;
      Mpi.allreduce_float !lsum Mpi.Sum comm

    let constrmask ((c : d), _, _) ((x : d), _, comm) ((m : d), _, _) =
      let test = ref 1.0 in
      let check b = if b then 0.0 else (test := 0.0; 1.0) in
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
      (Mpi.allreduce_float !test Mpi.Min comm = 1.0)

    let minquotient ((num : d), _, comm) ((denom : d), _, _) =
      let lmin = ref Config.big_real in
      for i = 0 to A.dim num - 1 do
        if (A.get denom i) <> 0.0 then
          lmin := floatmin !lmin (A.get num i /. A.get denom i)
      done;
      Mpi.allreduce_float !lmin Mpi.Min comm

    let prod ((x : d), _, _) ((y : d), _, _) ((z : d), _, _) =
      for i = 0 to A.dim x - 1 do
        A.set z i (A.get x i *. A.get y i)
      done

    let div ((x : d), _, _) ((y : d), _, _) ((z : d), _, _) =
      for i = 0 to A.dim x - 1 do
        A.set z i (A.get x i /. A.get y i)
      done

    let abs ((x : d), _, _) ((z : d), _, _) =
      for i = 0 to A.dim x - 1 do
        A.set z i (abs_float (A.get x i))
      done

    let inv ((x : d), _, _) ((z : d), _, _) =
      for i = 0 to A.dim x - 1 do
        A.set z i (1.0 /. (A.get x i))
      done

    let space (_, ng, comm) = (ng, 2 * Mpi.comm_size comm)

    let getlength (_, ng, _) = ng

    let print ?(logfile=Logfile.stdout) (x, _, _) =
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
        let zd, _, _ = z in
        let n = RealArray.length zd in
        if xa.(0) == z then begin
          let c0 = ca.{0} in
          if c0 <> 1.0 then
            for j = 0 to n - 1 do
              A.set zd j (A.get zd j *. c0)
            done;
          for i = 1 to nvec - 1 do
            let ci, (xd, _, _) = ca.{i}, xa.(i) in
            for j = 0 to n - 1 do
              A.set zd j (A.get zd j +. ci *. A.get xd j)
            done
          done
        end
        else begin
          let c0, (xd, _, _) = ca.{0}, xa.(0) in
          for j = 0 to n - 1 do
            A.set zd j (c0 *. A.get xd j)
          done;
          for i = 1 to nvec - 1 do
            let ci, (xd, _, _) = ca.{i}, xa.(i) in
            for j = 0 to n - 1 do
              A.set zd j (A.get zd j +. ci *. A.get xd j)
            done
          done
        end

    let scaleaddmulti (aa : RealArray.t) (x : t) (ya : t array) (za : t array) =
      let nvec = Array.length ya in
      if nvec = 1 then linearsum aa.{0} x 1.0 ya.(0) za.(0)
      else
        let xd, _, _ = x in
        let n = RealArray.length xd in
        if ya == za then
          for i = 0 to nvec - 1 do
            let a, (yd, _, _) = aa.{i}, ya.(i) in
            for j = 0 to n - 1 do
              A.set yd j (A.get yd j +. a *. A.get xd j)
            done
          done
        else
          for i = 0 to nvec - 1 do
            let ai, (yd, _, _), (zd, _, _) = aa.{i}, ya.(i), za.(i) in
            for j = 0 to n - 1 do
              A.set zd j (ai *. A.get xd j +. A.get yd j)
            done
          done

    let dotprodmulti (x : t) (ya : t array) (dp : RealArray.t) =
      let nvec = Array.length ya in
      if nvec = 1 then dp.{0} <- dotprod x ya.(0)
      else
        let xd, _, comm = x in
        let n = RealArray.length xd in
        for i = 0 to nvec - 1 do
          let yd, _, _ = ya.(i) in
          dp.{i} <- 0.0;
          for j = 0 to n - 1 do
            dp.{i} <- dp.{i} +. A.get xd j *. A.get yd j
          done
        done;
        (* Note: ocamlmpi does not provide MPI_IN_PLACE *)
        Mpi.(allreduce_bigarray1 dp dp Sum comm)

    let arr_vaxpy_array a (xa : t array) (ya : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.dim xad in
      if a = 1.0 then
        for i = 0 to nvec - 1 do
          let xd, _, _ = xa.(i) in
          let yd, _, _ = ya.(i) in
          for j = 0 to n - 1 do
            A.set yd j (A.get yd j +. A.get xd j)
          done
        done
      else if a = -1.0 then
        for i = 0 to nvec - 1 do
          let xd, _, _ = xa.(i) in
          let yd, _, _ = ya.(i) in
          for j = 0 to n - 1 do
            A.set yd j (A.get yd j -. A.get xd j)
          done
        done
      else
        for i = 0 to nvec - 1 do
          let xd, _, _ = xa.(i) in
          let yd, _, _ = ya.(i) in
          for j = 0 to n - 1 do
            A.set yd j (A.get yd j +. a *. A.get xd j)
          done
        done

    let v_sumvectorarray (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.dim xad in
      for i = 0 to nvec - 1 do
        let xd, _, _ = xa.(i) in
        let yd, _, _ = ya.(i) in
        let zd, _, _ = za.(i) in
        for j = 0 to n - 1 do
          A.set zd j (A.get xd j +. A.get yd j)
        done
      done

    let v_diffvectorarray (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.dim xad in
      for i = 0 to nvec - 1 do
        let xd, _, _ = xa.(i) in
        let yd, _, _ = ya.(i) in
        let zd, _, _ = za.(i) in
        for j = 0 to n - 1 do
          A.set zd j (A.get xd j -. A.get yd j)
        done
      done

    let v_lin1vectorarray a (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.dim xad in
      for i = 0 to nvec - 1 do
        let xd, _, _ = xa.(i) in
        let yd, _, _ = ya.(i) in
        let zd, _, _ = za.(i) in
        for j = 0 to n - 1 do
          A.set zd j (a *. A.get xd j +. A.get yd j)
        done
      done

    let v_lin2vectorarray a (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.dim xad in
      for i = 0 to nvec - 1 do
        let xd, _, _ = xa.(i) in
        let yd, _, _ = ya.(i) in
        let zd, _, _ = za.(i) in
        for j = 0 to n - 1 do
          A.set zd j (a *. A.get xd j -. A.get yd j)
        done
      done

    let v_scalesumvectorarray c (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.dim xad in
      for i = 0 to nvec - 1 do
        let xd, _, _ = xa.(i) in
        let yd, _, _ = ya.(i) in
        let zd, _, _ = za.(i) in
        for j = 0 to n - 1 do
          A.set zd j (c *. (A.get xd j +. A.get yd j))
        done
      done

    let v_scalediffvectorarray c (xa : t array) (ya : t array) (za : t array) =
      let nvec = Array.length xa in
      let xad, _, _ = xa.(0) in
      let n = A.dim xad in
      for i = 0 to nvec - 1 do
        let xd, _, _ = xa.(i) in
        let yd, _, _ = ya.(i) in
        let zd, _, _ = za.(i) in
        for j = 0 to n - 1 do
          A.set zd j (c *. (A.get xd j -. A.get yd j))
        done
      done

    let linearsumvectorarray a (xa : t array) b (ya : t array) (za : t array) =
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
        let xad, _, _ = xa.(0) in
        let n = A.dim xad in
        for i = 0 to nvec - 1 do
          let xd, _, _ = xa.(i) in
          let yd, _, _ = ya.(i) in
          let zd, _, _ = za.(i) in
          for j = 0 to n - 1 do
            A.set zd j (a *. A.get xd j +. b *. A.get yd j)
          done
        done

    let scalevectorarray (c : RealArray.t) (xa : t array) (za : t array) =
      let nvec = Array.length xa in
      if nvec = 1 then scale c.{0} xa.(0) za.(0)
      else
        let xad, _, _ = xa.(0) in
        let n = A.dim xad in
        if xa == za then
          for i = 0 to nvec - 1 do
            let (xd, _, _), c = xa.(i), c.{i} in
            for j = 0 to n - 1 do
              A.set xd j (c *. A.get xd j)
            done
          done
        else
          for i = 0 to nvec - 1 do
            let c = c.{i} in
            let xd, _, _ = xa.(i) in
            let zd, _, _ = za.(i) in
            for j = 0 to n - 1 do
              A.set zd j (c *. A.get xd j)
            done
          done

    let constvectorarray c (za : t array) =
      let nvec = Array.length za in
      if nvec = 1 then const c za.(0)
      else
        let zad, _, _ = za.(0) in
        let n = A.dim zad in
        for i = 0 to nvec - 1 do
          let zd, _, _ = za.(i) in
          for j = 0 to n - 1 do
            A.set zd j c
          done
        done

    let wrmsnormvectorarray (xa : t array) (wa : t array) (nrm : RealArray.t) =
      let nvec = Array.length xa in
      if nvec = 1 then nrm.{0} <- wrmsnorm xa.(0) wa.(0)
      else
        let xad, ng, comm = xa.(0) in
        let nl = A.dim xad in
        let nf = float ng in
        for i = 0 to nvec - 1 do
          let xd, _, _ = xa.(i) in
          let wd, _, _ = wa.(i) in
          nrm.{i} <- 0.0;
          for j = 0 to nl - 1 do
            let s = A.get xd j *. A.get wd j in
            nrm.{i} <- nrm.{i} +. s *. s
          done
        done;
        (* Note: ocamlmpi does not provide MPI_IN_PLACE *)
        Mpi.(allreduce_bigarray1 nrm nrm Sum comm);
        for i = 0 to nvec - 1 do
          nrm.{i} <- sqrt (nrm.{i}/.nf)
        done

    let wrmsnormmaskvectorarray (xa : t array) (wa : t array) (id : t)
                                   (nrm : RealArray.t) =
      let nvec = Array.length xa in
      if nvec = 1 then nrm.{0} <- wrmsnormmask xa.(0) wa.(0) id
      else
        let xad, ng, comm = xa.(0) in
        let idd, _, _ = id in
        let nl = A.dim xad in
        let nf = float ng in
        for i = 0 to nvec - 1 do
          let xd, _, _ = xa.(i) in
          let wd, _, _ = wa.(i) in
          nrm.{i} <- 0.0;
          for j = 0 to nl - 1 do
            if A.get idd j > 0.0 then begin
              let s = A.get xd j *. A.get wd j in
              nrm.{i} <- nrm.{i} +. s *. s
            end
          done;
        done;
        (* Note: ocamlmpi does not provide MPI_IN_PLACE *)
        Mpi.(allreduce_bigarray1 nrm nrm Sum comm);
        for i = 0 to nvec - 1 do
          nrm.{i} <- sqrt (nrm.{i}/.nf)
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
        let xad, _, _ = xa.(0) in
        let n = A.dim xad in
        if (yaa == zaa) then
          for i = 0 to nvec - 1 do
            let xd, _, _ = xa.(i) in
            for j = 0 to nsum - 1 do
              let a, (yd, _, _) = ra.{j}, yaa.(j).(i) in
              for k = 0 to n - 1 do
                A.set yd k (A.get yd k +. a *. A.get xd k)
              done
            done
          done
        else
          for i = 0 to nvec - 1 do
            let xd, _, _ = xa.(i) in
            for j = 0 to nsum - 1 do
              let a = ra.{j} in
              let yd, _, _ = yaa.(j).(i) in
              let zd, _, _ = zaa.(j).(i) in
              for k = 0 to n - 1 do
                A.set zd k (a *. A.get xd k +. A.get yd k)
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
          let zad, _, _ = za.(0) in
          let n = A.dim zad in
          if xaa.(0) == za then begin
            if ca.{0} = 1.0 then
              for j = 0 to nvec - 1 do
                let zd, _, _ = za.(j) in
                for i = 1 to nsum - 1 do
                  let c, (xd, _, _) = ca.{i}, xaa.(i).(j) in
                  for k = 0 to n - 1 do
                    A.set zd k (A.get zd k +. c *. A.get xd k)
                  done
                done
              done
            else
              let c0 = ca.{0} in
              for j = 0 to nvec - 1 do
                let zd, _, _ = za.(j) in
                for k = 0 to n - 1 do
                  A.set zd k (A.get zd k *. c0)
                done;
                for i = 1 to nsum - 1 do
                  let c, (xd, _, _) = ca.{i}, xaa.(i).(j) in
                  for k = 0 to n - 1 do
                    A.set zd k (A.get zd k +. c *. A.get xd k)
                  done
                done
              done
          end
          else
            let c0 = ca.{0} in
            for j = 0 to nvec - 1 do
              let xd, _, _ = xaa.(0).(j) in
              let zd, _, _ = za.(j) in
              for k = 0 to n - 1 do
                A.set zd k (c0 *. A.get xd k)
              done;
              for i = 1 to nsum - 1 do
                let c, (xd, _, _) = ca.{i}, xaa.(i).(j) in
                for k = 0 to n - 1 do
                  A.set zd k (A.get zd k +. c *. A.get xd k)
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

      let wsqrsum ((x : d), _, _) ((w : d), _, _) =
        let a = ref 0.0 in
        let lx = A.dim x in
        for i = 0 to lx - 1 do
          a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
        done;
        !a

      let wsqrsummask ((x : d), _, _) ((w : d), _, _) ((id : d), _, _) =
        let a = ref 0.0 in
        let lx = A.dim x in
        for i = 0 to lx - 1 do
          if A.get id i > 0.0 then
            a := !a +. (A.get x i *. A.get w i *. A.get x i *. A.get w i)
        done;
        !a
    end
  end (* }}} *)


(* Let C code know about some of the values in this module.  *)
external c_init_module : exn array -> unit =
  "sunml_nvector_parallel_init_module"

let _ =
  c_init_module
    (* Exceptions must be listed in the same order as
       nvector_parallel_exn_index.  *)
    [|IncorrectGlobalSize|]
