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

type cnvec
type ('data, 'kind) nvector =
  NV of { payload: 'data;
          cptr: cnvec;
          check: (('data, 'kind) nvector -> bool);
          clone: ('data, 'kind) t -> ('data, 'kind) t;
        }
and ('data, 'kind) t = ('data, 'kind) nvector

type 'k serial = (Sundials.RealArray.t, [>`Serial] as 'k) t

let unwrap (NV { payload; _ }) = payload

exception IncompatibleNvector

let check (NV { check; _ }) nv2 =
  if not (check nv2) then raise IncompatibleNvector

let clone (NV { clone; _ } as nv) = clone nv

let _ = Callback.register "Nvector.clone" clone

type nvector_id =
    Serial
  | Parallel
  | OpenMP
  | Pthreads
  | ParHyp
  | PETSc
  | CUDA
  | RAJA
  | OpenMPdev
  | Trilinos
  | ManyVector
  | MpiManyVector
  | MpiPlusX
  | Custom

external get_id : ('data, 'kind) t -> nvector_id
  = "sunml_nvec_get_id"

module type NVECTOR_OPS =
  sig (* {{{ *)
    type t

    val clone        : t -> t
    val linearsum    : float -> t -> float -> t -> t -> unit
    val const        : float -> t -> unit
    val prod         : t -> t -> t -> unit
    val div          : t -> t -> t -> unit
    val scale        : float -> t -> t -> unit
    val abs          : t -> t -> unit
    val inv          : t -> t -> unit
    val addconst     : t -> float -> t -> unit
    val dotprod      : t -> t -> float
    val maxnorm      : t -> float
    val wrmsnorm     : t -> t -> float
    val min          : t -> float
    val compare      : float -> t -> t -> unit
    val invtest      : t -> t -> bool

    val wl2norm      : t -> t -> float
    val l1norm       : t -> float
    val wrmsnormmask : t -> t -> t -> float
    val constrmask   : t -> t -> t -> bool
    val minquotient  : t -> t -> float

    val space        : t -> int * int
    val getlength    : t -> int
    val print        : ?logfile:Sundials.Logfile.t -> t -> unit

    val linearcombination
      : Sundials.RealArray.t -> t array -> t -> unit
    val scaleaddmulti
      : Sundials.RealArray.t -> t -> t array -> t array -> unit
    val dotprodmulti
      : t -> t array -> Sundials.RealArray.t -> unit

    val linearsumvectorarray
      : float -> t array -> float -> t array -> t array -> unit
    val scalevectorarray
      : Sundials.RealArray.t -> t array -> t array -> unit
    val constvectorarray
      : float -> t array -> unit
    val wrmsnormvectorarray
      : t array -> t array -> Sundials.RealArray.t -> unit
    val wrmsnormmaskvectorarray
      : t array -> t array -> t -> Sundials.RealArray.t -> unit
    val scaleaddmultivectorarray
      : Sundials.RealArray.t -> t array -> t array array -> t array array -> unit
    val linearcombinationvectorarray
      : Sundials.RealArray.t -> t array array -> t array -> unit

    module Local : sig
      val dotprod     : t -> t -> float
      val maxnorm     : t -> float
      val min         : t -> float
      val l1norm      : t -> float
      val invtest     : t -> t -> bool
      val constrmask  : t -> t -> t -> bool
      val minquotient : t -> t -> float
      val wsqrsum     : t -> t -> float
      val wsqrsummask : t -> t -> t -> float
    end
  end (* }}} *)

module type NVECTOR =
  sig (* {{{ *)
    type kind
    type data
    type t = (data, kind) nvector
    val wrap : ?with_fused_ops:bool -> data -> t
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
    module Ops : NVECTOR_OPS with type t = t
    module DataOps : NVECTOR_OPS with type t = data
  end (* }}} *)

(* Hack to ensure that Sundials.c_init_module is executed so that the global
   exceptions are properly registered. *)
let e = Sundials.RecoverableFailure

type gdata = ..
type gdata += RA of Sundials.RealArray.t
type gkind
type any = (gdata, gkind) t
exception BadGenericType

exception OperationNotProvided

module MakeDataOps (Ops : sig
    include NVECTOR_OPS

    val unwrap : gdata -> t
    val wrap : t -> gdata
  end) : NVECTOR_OPS with type t = gdata
=
struct (* {{{ *)
  type t = gdata

  let uw = Ops.unwrap
  let uwa = Array.map uw
  let uwaa = Array.map uwa

  let clone a = Ops.wrap (Ops.clone (uw a))
  let linearsum a x b y z = Ops.linearsum a (uw x) b (uw y) (uw z)
  let const c x = Ops.const c (uw x)
  let prod x y z = Ops.prod (uw x) (uw y) (uw z)
  let div x y z = Ops.div (uw x) (uw y) (uw z)
  let scale a x z = Ops.scale a (uw x) (uw z)
  let abs x z = Ops.abs (uw x ) (uw z)
  let inv x z = Ops.inv (uw x) (uw z)
  let addconst x c z = Ops.addconst (uw x) c (uw z)
  let dotprod x y = Ops.dotprod (uw x) (uw y)
  let maxnorm x = Ops.maxnorm (uw x)
  let wrmsnorm x y = Ops.wrmsnorm (uw x) (uw y)
  let min x = Ops.min (uw x)
  let compare c x y = Ops.compare c (uw x) (uw y)
  let invtest x z = Ops.invtest (uw x) (uw z)

  let wl2norm x y = Ops.wl2norm (uw x) (uw y)
  let l1norm x = Ops.l1norm (uw x)
  let wrmsnormmask x y z = Ops.wrmsnormmask (uw x) (uw y) (uw z)
  let constrmask x y z = Ops.constrmask (uw x) (uw y) (uw z)
  let minquotient x y = Ops.minquotient (uw x) (uw y)

  let space x = Ops.space (uw x)
  let getlength x = Ops.getlength (uw x)
  let print ?logfile x = Ops.print ?logfile (uw x)

  let linearcombination a xa y = Ops.linearcombination a (uwa xa) (uw y)
  let scaleaddmulti a x ya za = Ops.scaleaddmulti a (uw x) (uwa ya) (uwa za)
  let dotprodmulti x ya a = Ops.dotprodmulti (uw x) (uwa ya) a

  let linearsumvectorarray a xa b ya za
    = Ops.linearsumvectorarray a (uwa xa) b (uwa ya) (uwa za)
  let scalevectorarray ra xa ya = Ops.scalevectorarray ra (uwa xa) (uwa ya)
  let constvectorarray c xa = Ops.constvectorarray c (uwa xa)
  let wrmsnormvectorarray xa ya ra
    = Ops.wrmsnormvectorarray (uwa xa) (uwa ya) ra
  let wrmsnormmaskvectorarray xa ya z ra
    = Ops.wrmsnormmaskvectorarray (uwa xa) (uwa ya) (uw z) ra
  let scaleaddmultivectorarray ra xa yaa zaa
    = Ops.scaleaddmultivectorarray ra (uwa xa) (uwaa yaa) (uwaa zaa)
  let linearcombinationvectorarray ra xaa ya
    = Ops.linearcombinationvectorarray ra (uwaa xaa) (uwa ya)

  module Local = struct
    let dotprod x y = Ops.Local.dotprod (uw x) (uw y)
    let maxnorm x = Ops.Local.maxnorm (uw x)
    let min x = Ops.Local.min (uw x)
    let l1norm x = Ops.Local.l1norm (uw x)
    let invtest x y = Ops.Local.invtest (uw x) (uw y)
    let constrmask x y z = Ops.Local.constrmask (uw x) (uw y) (uw z)
    let minquotient x y = Ops.Local.minquotient (uw x) (uw y)
    let wsqrsum x y = Ops.Local.wsqrsum (uw x) (uw y)
    let wsqrsummask x y z = Ops.Local.wsqrsummask (uw x) (uw y) (uw z)
  end
end (* }}} *)

module Ops =
  struct (* {{{ *)
    let clone = clone

    external c_linearsum
      : float -> ('d, 'k) t -> float -> ('d, 'k) t -> ('d, 'k) t -> unit
      = "sunml_nvec_linearsum" [@@noalloc]

    let linearsum a x b y z =
      if Sundials_configuration.safe then (check x y; check x z);
      c_linearsum a x b y z

    external const
      : float -> ('d, 'k) t -> unit
      = "sunml_nvec_const" [@@noalloc]

    external c_prod
      : ('d, 'k) t -> ('d, 'k) t -> ('d, 'k) t -> unit
      = "sunml_nvec_prod" [@@noalloc]

    let prod x y z =
      if Sundials_configuration.safe then (check x y; check x z);
      c_prod x y z

    external c_div
      : ('d, 'k) t -> ('d, 'k) t -> ('d, 'k) t -> unit
      = "sunml_nvec_div" [@@noalloc]

    let div x y z =
      if Sundials_configuration.safe then (check x y; check x z);
      c_div x y z

    external c_scale
      : float -> ('d, 'k) t -> ('d, 'k) t -> unit
      = "sunml_nvec_scale" [@@noalloc]

    let scale c x z =
      if Sundials_configuration.safe then check x z;
      c_scale c x z

    external c_abs
      : ('d, 'k) t -> ('d, 'k) t -> unit
      = "sunml_nvec_abs" [@@noalloc]

    let abs x z =
      if Sundials_configuration.safe then check x z;
      c_abs x z

    external c_inv
      : ('d, 'k) t -> ('d, 'k) t -> unit
      = "sunml_nvec_inv" [@@noalloc]

    let inv x z =
      if Sundials_configuration.safe then check x z;
      c_inv x z

    external c_addconst
      : ('d, 'k) t -> float -> ('d, 'k) t -> unit
      = "sunml_nvec_addconst" [@@noalloc]

    let addconst x b z =
      if Sundials_configuration.safe then check x z;
      c_addconst x b z

    external c_dotprod
      : ('d, 'k) t -> ('d, 'k) t -> float
      = "sunml_nvec_dotprod"

    let dotprod x y =
      if Sundials_configuration.safe then check x y;
      c_dotprod x y

    external maxnorm
      : ('d, 'k) t -> float
      = "sunml_nvec_maxnorm"

    external c_wrmsnorm
      : ('d, 'k) t -> ('d, 'k) t -> float
      = "sunml_nvec_wrmsnorm"

    let wrmsnorm x w =
      if Sundials_configuration.safe then check x w;
      c_wrmsnorm x w

    external c_wrmsnormmask
      : ('d, 'k) t -> ('d, 'k) t -> ('d, 'k) t -> float
      = "sunml_nvec_wrmsnormmask"

    let wrmsnormmask x w id =
      if Sundials_configuration.safe then (check x w; check x id);
      c_wrmsnormmask x w id

    external min
      : ('d, 'k) t -> float
      = "sunml_nvec_min"

    external c_wl2norm
      : ('d, 'k) t -> ('d, 'k) t -> float
      = "sunml_nvec_wl2norm"

    let wl2norm x w =
      if Sundials_configuration.safe then check x w;
      c_wl2norm x w

    external l1norm
      : ('d, 'k) t -> float
      = "sunml_nvec_l1norm"

    external c_compare
      : float -> ('d, 'k) t -> ('d, 'k) t -> unit
      = "sunml_nvec_compare" [@@noalloc]

    let compare c x z =
      if Sundials_configuration.safe then check x z;
      c_compare c x z

    external c_invtest
      : ('d, 'k) t -> ('d, 'k) t -> bool
      = "sunml_nvec_invtest" [@@noalloc]

    let invtest x z =
      if Sundials_configuration.safe then check x z;
      c_invtest x z

    external c_constrmask
      : ('d, 'k) t -> ('d, 'k) t -> ('d, 'k) t -> bool
      = "sunml_nvec_constrmask"

    let constrmask c x m =
      if Sundials_configuration.safe then (check c x; check c m);
      c_constrmask c x m

    external c_minquotient
      : ('d, 'k) t -> ('d, 'k) t -> float
      = "sunml_nvec_minquotient"

    let minquotient n d =
      if Sundials_configuration.safe then check n d;
      c_minquotient n d

    external space
      : ('d, 'k) t -> int * int
      = "sunml_nvec_space"

    external getlength
      : ('d, 'k) t -> int
      = "sunml_nvec_getlength"

    external c_print_file
      : ('d, 'k) t -> Sundials.Logfile.t option -> unit
      = "sunml_nvec_print_file"

    let print ?logfile nv = c_print_file nv logfile

    external c_linearcombination
      : Sundials.RealArray.t -> ('d, 'k) t array -> ('d, 'k) t -> unit
      = "sunml_nvec_linearcombination"

    let linearcombination ca xa z =
      if Sundials_impl.Version.lt400
        then raise Sundials.Config.NotImplementedBySundialsVersion;
      if Sundials_configuration.safe then Array.iter (check z) xa;
      c_linearcombination ca xa z

    let same_len' n ya =
      if n <> Array.length ya then invalid_arg "arrays of unequal length"
    let same_len xa ya = same_len' (Array.length xa) ya

    external c_scaleaddmulti
      : Sundials.RealArray.t -> ('d, 'k) t -> ('d, 'k) t array -> ('d, 'k) t array -> unit
      = "sunml_nvec_scaleaddmulti"

    let scaleaddmulti aa x ya za =
      if Sundials_impl.Version.lt400
        then raise Sundials.Config.NotImplementedBySundialsVersion;
      if Sundials_configuration.safe then
        (Array.iter (check x) ya; Array.iter (check x) za;
         let nv = Sundials.RealArray.length aa in
         same_len' nv ya; same_len' nv za);
      c_scaleaddmulti aa x ya za

    external c_dotprodmulti
      : ('d, 'k) t -> ('d, 'k) t array -> Sundials.RealArray.t -> unit
      = "sunml_nvec_dotprodmulti"

    let dotprodmulti x ya dp =
      if Sundials_impl.Version.lt400
        then raise Sundials.Config.NotImplementedBySundialsVersion;
      if Sundials_configuration.safe then
        (let nv = Sundials.RealArray.length dp in
         same_len' nv ya;
         Array.iter (check x) ya);
      c_dotprodmulti x ya dp

    external c_linearsumvectorarray
      : float -> ('d, 'k) t array -> float -> ('d, 'k) t array -> ('d, 'k) t array -> unit
      = "sunml_nvec_linearsumvectorarray"

    let linearsumvectorarray a xa b ya za =
      if Sundials_impl.Version.lt400
        then raise Sundials.Config.NotImplementedBySundialsVersion;
      if Sundials_configuration.safe
      then (let x = Array.get xa 0 in
            Array.iter (check x) xa;
            Array.iter (check x) ya;
            Array.iter (check x) za;
            same_len xa ya; same_len xa za);
      c_linearsumvectorarray a xa b ya za

    external c_scalevectorarray
      : Sundials.RealArray.t -> ('d, 'k) t array -> ('d, 'k) t array -> unit
      = "sunml_nvec_scalevectorarray"

    let scalevectorarray c xa za =
      if Sundials_impl.Version.lt400
        then raise Sundials.Config.NotImplementedBySundialsVersion;
      if Sundials_configuration.safe
      then (let x = Array.get xa 0 in
            Array.iter (check x) xa;
            Array.iter (check x) za;
            same_len xa za);
      c_scalevectorarray c xa za

    external c_constvectorarray
      : float -> ('d, 'k) t array -> unit
      = "sunml_nvec_constvectorarray"

    let constvectorarray c za =
      if Sundials_impl.Version.lt400
        then raise Sundials.Config.NotImplementedBySundialsVersion;
      if Sundials_configuration.safe
      then (let z = Array.get za 0 in
            Array.iter (check z) za);
      c_constvectorarray c za

    external c_wrmsnormvectorarray
      : ('d, 'k) t array -> ('d, 'k) t array -> Sundials.RealArray.t -> unit
      = "sunml_nvec_wrmsnormvectorarray"

    let wrmsnormvectorarray xa wa nrm =
      if Sundials_impl.Version.lt400
        then raise Sundials.Config.NotImplementedBySundialsVersion;
      if Sundials_configuration.safe
      then (let x = Array.get xa 0 in
            Array.iter (check x) xa;
            Array.iter (check x) wa;
           same_len xa wa);
      c_wrmsnormvectorarray xa wa nrm

    external c_wrmsnormmaskvectorarray
      : ('d, 'k) t array -> ('d, 'k) t array -> ('d, 'k) t -> Sundials.RealArray.t -> unit
      = "sunml_nvec_wrmsnormmaskvectorarray"

    let wrmsnormmaskvectorarray xa wa id nrm =
      if Sundials_impl.Version.lt400
        then raise Sundials.Config.NotImplementedBySundialsVersion;
      if Sundials_configuration.safe
      then (Array.iter (check id) xa;
            Array.iter (check id) wa;
            same_len xa wa);
      c_wrmsnormmaskvectorarray xa wa id nrm

    external c_scaleaddmultivectorarray
      : Sundials.RealArray.t -> ('d, 'k) t array -> ('d, 'k) t array array -> ('d, 'k) t array array -> unit
      = "sunml_nvec_scaleaddmultivectorarray"

    let scaleaddmultivectorarray ra xa yaa zaa =
      if Sundials_impl.Version.lt400
        then raise Sundials.Config.NotImplementedBySundialsVersion;
      if Sundials_configuration.safe
      then (let x = Array.get xa 0 in
            let ns = Sundials.RealArray.length ra in
            let nv = Array.length xa in
            same_len' ns yaa;
            same_len' ns zaa;
            Array.iter (check x) xa;
            Array.iter (fun ya -> same_len' nv ya; Array.iter (check x) ya) yaa;
            Array.iter (fun za -> same_len' nv za; Array.iter (check x) za) zaa;
            same_len yaa zaa);
      c_scaleaddmultivectorarray ra xa yaa zaa

    external c_linearcombinationvectorarray
      : Sundials.RealArray.t -> ('d, 'k) t array array -> ('d, 'k) t array -> unit
      = "sunml_nvec_linearcombinationvectorarray"

    let linearcombinationvectorarray ca xaa za =
      if Sundials_impl.Version.lt400
        then raise Sundials.Config.NotImplementedBySundialsVersion;
      if Sundials_configuration.safe
      then (let z = Array.get za 0 in
            let ns = Sundials.RealArray.length ca in
            let nv = Array.length za in
            same_len' ns xaa;
            Array.iter (check z) za;
            Array.iter (fun xa -> same_len' nv xa; Array.iter (check z) xa) xaa);
      c_linearcombinationvectorarray ca xaa za

    external has_linearcombination            : ('d, 'k) t -> bool
        = "sunml_nvec_has_linearcombination" [@@noalloc]
    external has_scaleaddmulti                : ('d, 'k) t -> bool
        = "sunml_nvec_has_scaleaddmulti" [@@noalloc]
    external has_dotprodmulti                 : ('d, 'k) t -> bool
        = "sunml_nvec_has_dotprodmulti" [@@noalloc]
    external has_linearsumvectorarray         : ('d, 'k) t -> bool
        = "sunml_nvec_has_linearsumvectorarray" [@@noalloc]
    external has_scalevectorarray             : ('d, 'k) t -> bool
        = "sunml_nvec_has_scalevectorarray" [@@noalloc]
    external has_constvectorarray             : ('d, 'k) t -> bool
        = "sunml_nvec_has_constvectorarray" [@@noalloc]
    external has_wrmsnormvectorarray          : ('d, 'k) t -> bool
        = "sunml_nvec_has_wrmsnormvectorarray" [@@noalloc]
    external has_wrmsnormmaskvectorarray      : ('d, 'k) t -> bool
        = "sunml_nvec_has_wrmsnormmaskvectorarray" [@@noalloc]
    external has_scaleaddmultivectorarray     : ('d, 'k) t -> bool
        = "sunml_nvec_has_scaleaddmultivectorarray" [@@noalloc]
    external has_linearcombinationvectorarray : ('d, 'k) t -> bool
        = "sunml_nvec_has_linearcombinationvectorarray" [@@noalloc]

    module Local = struct

      external c_dotprod : ('d, 'k) t -> ('d, 'k) t -> float
        = "sunml_nvec_dotprodlocal"

      let dotprod x y =
        if Sundials_impl.Version.lt500
          then raise Sundials.Config.NotImplementedBySundialsVersion;
        if Sundials_configuration.safe then check x y;
        c_dotprod x y

      external c_maxnorm : ('d, 'k) t -> float
        = "sunml_nvec_maxnormlocal"

      let maxnorm x =
        if Sundials_impl.Version.lt500
          then raise Sundials.Config.NotImplementedBySundialsVersion;
        c_maxnorm x

      external c_min : ('d, 'k) t -> float
        = "sunml_nvec_minlocal"

      let min x =
        if Sundials_impl.Version.lt500
          then raise Sundials.Config.NotImplementedBySundialsVersion;
        c_min x

      external c_l1norm : ('d, 'k) t -> float
        = "sunml_nvec_l1normlocal"

      let l1norm x =
        if Sundials_impl.Version.lt500
          then raise Sundials.Config.NotImplementedBySundialsVersion;
        c_l1norm x

      external c_invtest : ('d, 'k) t -> ('d, 'k) t -> bool
        = "sunml_nvec_invtestlocal"

      let invtest x z =
        if Sundials_impl.Version.lt500
          then raise Sundials.Config.NotImplementedBySundialsVersion;
        if Sundials_configuration.safe then check x z;
        c_invtest x z

      external c_constrmask : ('d, 'k) t -> ('d, 'k) t -> ('d, 'k) t -> bool
        = "sunml_nvec_constrmasklocal"

      let constrmask c x m =
        if Sundials_impl.Version.lt500
          then raise Sundials.Config.NotImplementedBySundialsVersion;
        if Sundials_configuration.safe then (check c x; check c m);
        c_constrmask c x m

      external c_minquotient : ('d, 'k) t -> ('d, 'k) t -> float
        = "sunml_nvec_minquotientlocal"

      let minquotient n d =
        if Sundials_impl.Version.lt500
          then raise Sundials.Config.NotImplementedBySundialsVersion;
        if Sundials_configuration.safe then check n d;
        c_minquotient n d

      external c_wsqrsum
        : ('d, 'k) t -> ('d, 'k) t -> float
        = "sunml_nvec_wsqrsumlocal"

      let wsqrsum x w =
        if Sundials_impl.Version.lt500
          then raise Sundials.Config.NotImplementedBySundialsVersion;
        if Sundials_configuration.safe then check x w;
        c_wsqrsum x w

      external c_wsqrsummask
        : ('d, 'k) t -> ('d, 'k) t -> ('d, 'k) t -> float
        = "sunml_nvec_wsqrsummasklocal"

      let wsqrsummask x w id =
        if Sundials_impl.Version.lt500
          then raise Sundials.Config.NotImplementedBySundialsVersion;
        if Sundials_configuration.safe then (check x w; check x id);
        c_wsqrsummask x w id

      external has_dotprod      : ('d, 'k) t -> bool
        = "sunml_nvec_has_dotprodlocal" [@@noalloc]
      external has_maxnorm      : ('d, 'k) t -> bool
        = "sunml_nvec_has_maxnormlocal" [@@noalloc]
      external has_min          : ('d, 'k) t -> bool
        = "sunml_nvec_has_minlocal" [@@noalloc]
      external has_l1norm       : ('d, 'k) t -> bool
        = "sunml_nvec_has_l1normlocal" [@@noalloc]
      external has_invtest      : ('d, 'k) t -> bool
        = "sunml_nvec_has_invtestlocal" [@@noalloc]
      external has_constrmask   : ('d, 'k) t -> bool
        = "sunml_nvec_has_constrmasklocal" [@@noalloc]
      external has_minquotient  : ('d, 'k) t -> bool
        = "sunml_nvec_has_minquotientlocal" [@@noalloc]
      external has_wsqrsum      : ('d, 'k) t -> bool
        = "sunml_nvec_has_wsqrsumlocal" [@@noalloc]
      external has_wsqrsummask  : ('d, 'k) t -> bool
        = "sunml_nvec_has_wsqrsummasklocal" [@@noalloc]
    end
  end (* }}} *)

(* Let C code know about some of the values in this module.  *)
external c_init_module : exn array -> unit =
  "sunml_nvec_init_module"

let _ =
  c_init_module
    (* Exceptions must be listed in the same order as nv_exn_index.  *)
    [|OperationNotProvided;
    |]

