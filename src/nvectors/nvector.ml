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

let unwrap (NV { payload; _ }) = payload

exception IncompatibleNvector

let check (NV { check; _ }) nv2 =
  if not (check nv2) then raise IncompatibleNvector

let clone (NV { clone; _ } as nv) = clone nv

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

    val n_vclone        : t -> t
    val n_vlinearsum    : float -> t -> float -> t -> t -> unit
    val n_vconst        : float -> t -> unit
    val n_vprod         : t -> t -> t -> unit
    val n_vdiv          : t -> t -> t -> unit
    val n_vscale        : float -> t -> t -> unit
    val n_vabs          : t -> t -> unit
    val n_vinv          : t -> t -> unit
    val n_vaddconst     : t -> float -> t -> unit
    val n_vdotprod      : t -> t -> float
    val n_vmaxnorm      : t -> float
    val n_vwrmsnorm     : t -> t -> float
    val n_vmin          : t -> float
    val n_vcompare      : float -> t -> t -> unit
    val n_vinvtest      : t -> t -> bool

    val n_vwl2norm      : t -> t -> float
    val n_vl1norm       : t -> float
    val n_vwrmsnormmask : t -> t -> t -> float
    val n_vconstrmask   : t -> t -> t -> bool
    val n_vminquotient  : t -> t -> float

    val n_vspace        : t -> int * int
    val n_vgetlength    : t -> int

    val n_vlinearcombination
      : Sundials.RealArray.t -> t array -> t -> unit
    val n_vscaleaddmulti
      : Sundials.RealArray.t -> t -> t array -> t array -> unit
    val n_vdotprodmulti
      : t -> t array -> Sundials.RealArray.t -> unit

    val n_vlinearsumvectorarray
      : float -> t array -> float -> t array -> t array -> unit
    val n_vscalevectorarray
      : Sundials.RealArray.t -> t array -> t array -> unit
    val n_vconstvectorarray
      : float -> t array -> unit
    val n_vwrmsnormvectorarray
      : t array -> t array -> Sundials.RealArray.t -> unit
    val n_vwrmsnormmaskvectorarray
      : t array -> t array -> t -> Sundials.RealArray.t -> unit
    val n_vscaleaddmultivectorarray
      : Sundials.RealArray.t -> t array -> t array array -> t array array -> unit
    val n_vlinearcombinationvectorarray
      : Sundials.RealArray.t -> t array array -> t array -> unit

    module Local : sig
      val n_vdotprod     : t -> t -> float
      val n_vmaxnorm     : t -> float
      val n_vmin         : t -> float
      val n_vl1norm      : t -> float
      val n_vinvtest     : t -> t -> bool
      val n_vconstrmask  : t -> t -> t -> bool
      val n_vminquotient : t -> t -> float
      val n_vwsqrsum     : t -> t -> float
      val n_vwsqrsummask : t -> t -> t -> float
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

(* {{{ *)

external has_n_vlinearcombination            : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vlinearcombination" [@@noalloc]
external has_n_vscaleaddmulti                : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vscaleaddmulti" [@@noalloc]
external has_n_vdotprodmulti                 : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vdotprodmulti" [@@noalloc]
external has_n_vlinearsumvectorarray         : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vlinearsumvectorarray" [@@noalloc]
external has_n_vscalevectorarray             : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vscalevectorarray" [@@noalloc]
external has_n_vconstvectorarray             : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vconstvectorarray" [@@noalloc]
external has_n_vwrmsnormvectorarray          : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vwrmsnormvectorarray" [@@noalloc]
external has_n_vwrmsnormmaskvectorarray      : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vwrmsnormmaskvectorarray" [@@noalloc]
external has_n_vscaleaddmultivectorarray     : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vscaleaddmultivectorarray" [@@noalloc]
external has_n_vlinearcombinationvectorarray : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vlinearcombinationvectorarray" [@@noalloc]

module Local = struct
  external has_n_vdotprod      : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vdotprodlocal" [@@noalloc]
  external has_n_vmaxnorm      : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vmaxnormlocal" [@@noalloc]
  external has_n_vmin          : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vminlocal" [@@noalloc]
  external has_n_vl1norm       : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vl1normlocal" [@@noalloc]
  external has_n_vinvtest      : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vinvtestlocal" [@@noalloc]
  external has_n_vconstrmask   : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vconstrmasklocal" [@@noalloc]
  external has_n_vminquotient  : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vminquotientlocal" [@@noalloc]
  external has_n_vwsqrsum      : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vwsqrsumlocal" [@@noalloc]
  external has_n_vwsqrsummask  : ('d, 'k) t -> bool
    = "sunml_nvec_has_n_vwsqrsummasklocal" [@@noalloc]
end

(* }}} *)

type gdata = ..
type gdata += RA of Sundials.RealArray.t
type gkind
type any = (gdata, gkind) t
exception BadGenericType

module Any = struct (* {{{ *)
  type t = any

  external has_n_vlinearcombination            : t -> bool
      = "sunml_nvec_has_n_vlinearcombination" [@@noalloc]
  external has_n_vscaleaddmulti                : t -> bool
      = "sunml_nvec_has_n_vscaleaddmulti" [@@noalloc]
  external has_n_vdotprodmulti                 : t -> bool
      = "sunml_nvec_has_n_vdotprodmulti" [@@noalloc]
  external has_n_vlinearsumvectorarray         : t -> bool
      = "sunml_nvec_has_n_vlinearsumvectorarray" [@@noalloc]
  external has_n_vscalevectorarray             : t -> bool
      = "sunml_nvec_has_n_vscalevectorarray" [@@noalloc]
  external has_n_vconstvectorarray             : t -> bool
      = "sunml_nvec_has_n_vconstvectorarray" [@@noalloc]
  external has_n_vwrmsnormvectorarray          : t -> bool
      = "sunml_nvec_has_n_vwrmsnormvectorarray" [@@noalloc]
  external has_n_vwrmsnormmaskvectorarray      : t -> bool
      = "sunml_nvec_has_n_vwrmsnormmaskvectorarray" [@@noalloc]
  external has_n_vscaleaddmultivectorarray     : t -> bool
      = "sunml_nvec_has_n_vscaleaddmultivectorarray" [@@noalloc]
  external has_n_vlinearcombinationvectorarray : t -> bool
      = "sunml_nvec_has_n_vlinearcombinationvectorarray" [@@noalloc]

  module Local = struct
    external has_n_vdotprod      : t -> bool
      = "sunml_nvec_has_n_vdotprodlocal" [@@noalloc]
    external has_n_vmaxnorm      : t -> bool
      = "sunml_nvec_has_n_vmaxnormlocal" [@@noalloc]
    external has_n_vmin          : t -> bool
      = "sunml_nvec_has_n_vminlocal" [@@noalloc]
    external has_n_vl1norm       : t -> bool
      = "sunml_nvec_has_n_vl1normlocal" [@@noalloc]
    external has_n_vinvtest      : t -> bool
      = "sunml_nvec_has_n_vinvtestlocal" [@@noalloc]
    external has_n_vconstrmask   : t -> bool
      = "sunml_nvec_has_n_vconstrmasklocal" [@@noalloc]
    external has_n_vminquotient  : t -> bool
      = "sunml_nvec_has_n_vminquotientlocal" [@@noalloc]
    external has_n_vwsqrsum      : t -> bool
      = "sunml_nvec_has_n_vwsqrsumlocal" [@@noalloc]
    external has_n_vwsqrsummask  : t -> bool
      = "sunml_nvec_has_n_vwsqrsummasklocal" [@@noalloc]
  end
end (* }}} *)

module Ops = struct (* {{{ *)
  type t = any

  let n_vclone = clone

  external c_n_vlinearsum    : float -> t -> float -> t -> t -> unit
    = "sunml_nvec_any_n_vlinearsum"

  external n_vconst          : float -> t -> unit
    = "sunml_nvec_any_n_vconst"

  external c_n_vprod         : t -> t -> t -> unit
    = "sunml_nvec_any_n_vprod"

  external c_n_vdiv          : t -> t -> t -> unit
    = "sunml_nvec_any_n_vdiv"

  external c_n_vscale        : float -> t -> t -> unit
    = "sunml_nvec_any_n_vscale"

  external c_n_vabs          : t -> t -> unit
    = "sunml_nvec_any_n_vabs"

  external c_n_vinv          : t -> t -> unit
    = "sunml_nvec_any_n_vinv"

  external c_n_vaddconst     : t -> float -> t -> unit
    = "sunml_nvec_any_n_vaddconst"

  external c_n_vdotprod      : t -> t -> float
    = "sunml_nvec_any_n_vdotprod"

  external n_vmaxnorm        : t -> float
    = "sunml_nvec_any_n_vmaxnorm"

  external c_n_vwrmsnorm     : t -> t -> float
    = "sunml_nvec_any_n_vwrmsnorm"

  external c_n_vwrmsnormmask : t -> t -> t -> float
    = "sunml_nvec_any_n_vwrmsnormmask"

  external n_vmin            : t -> float
    = "sunml_nvec_any_n_vmin"

  external c_n_vwl2norm      : t -> t -> float
    = "sunml_nvec_any_n_vwl2norm"

  external n_vl1norm         : t -> float
    = "sunml_nvec_any_n_vl1norm"

  external c_n_vcompare      : float -> t -> t -> unit
    = "sunml_nvec_any_n_vcompare"

  external c_n_vinvtest      : t -> t -> bool
    = "sunml_nvec_any_n_vinvtest"

  external c_n_vconstrmask   : t -> t -> t -> bool
    = "sunml_nvec_any_n_vconstrmask"

  external c_n_vminquotient  : t -> t -> float
    = "sunml_nvec_any_n_vminquotient"

  external n_vspace  : t -> int * int
    = "sunml_nvec_any_n_vspace"

  external n_vgetlength : t -> int
    = "sunml_nvec_any_n_vgetlength"

  external c_n_vlinearcombination
    : Sundials.RealArray.t -> t array -> t -> unit
    = "sunml_nvec_any_n_vlinearcombination"

  external c_n_vscaleaddmulti
    : Sundials.RealArray.t -> t -> t array -> t array -> unit
    = "sunml_nvec_any_n_vscaleaddmulti"

  external c_n_vdotprodmulti
    : t -> t array -> Sundials.RealArray.t -> unit
    = "sunml_nvec_any_n_vdotprodmulti"

  external c_n_vlinearsumvectorarray
    : float -> t array -> float -> t array -> t array -> unit
    = "sunml_nvec_any_n_vlinearsumvectorarray"

  external c_n_vscalevectorarray
    : Sundials.RealArray.t -> t array -> t array -> unit
    = "sunml_nvec_any_n_vscalevectorarray"

  external c_n_vconstvectorarray
    : float -> t array -> unit
    = "sunml_nvec_any_n_vconstvectorarray"

  external c_n_vwrmsnormvectorarray
    : t array -> t array -> Sundials.RealArray.t -> unit
    = "sunml_nvec_any_n_vwrmsnormvectorarray"

  external c_n_vwrmsnormmaskvectorarray
    : t array -> t array -> t -> Sundials.RealArray.t -> unit
    = "sunml_nvec_any_n_vwrmsnormmaskvectorarray"

  external c_n_vscaleaddmultivectorarray
    : Sundials.RealArray.t -> t array -> t array array -> t array array -> unit
    = "sunml_nvec_any_n_vscaleaddmultivectorarray"

  external c_n_vlinearcombinationvectorarray
    : Sundials.RealArray.t -> t array array -> t array -> unit
    = "sunml_nvec_any_n_vlinearcombinationvectorarray"

  let n_vlinearsum a (x : t) b (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_n_vlinearsum a x b y z

  let n_vscale c (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_n_vscale c x z

  let n_vaddconst (x : t) b (z : t) =
    if Sundials_configuration.safe then check x z;
    c_n_vaddconst x b z

  let n_vwrmsnorm (x : t) (w : t) =
    if Sundials_configuration.safe then check x w;
    c_n_vwrmsnorm x w

  let n_vwrmsnormmask (x : t) (w : t) (id : t) =
    if Sundials_configuration.safe then (check x w; check x id);
    c_n_vwrmsnormmask x w id

  let n_vdotprod (x : t) (y : t) =
    if Sundials_configuration.safe then check x y;
    c_n_vdotprod x y

  let n_vcompare c (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_n_vcompare c x z

  let n_vinvtest (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_n_vinvtest x z

  let n_vwl2norm (x : t) (w : t) =
    if Sundials_configuration.safe then check x w;
    c_n_vwl2norm x w

  let n_vconstrmask (c : t) (x : t) (m : t) =
    if Sundials_configuration.safe then (check c x; check c m);
    c_n_vconstrmask c x m

  let n_vminquotient (n : t) (d : t) =
    if Sundials_configuration.safe then check n d;
    c_n_vminquotient n d

  let n_vprod (x : t) (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_n_vprod x y z

  let n_vdiv (x : t) (y : t) (z : t) =
    if Sundials_configuration.safe then (check x y; check x z);
    c_n_vdiv x y z

  let n_vabs (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_n_vabs x z

  let n_vinv (x : t) (z : t) =
    if Sundials_configuration.safe then check x z;
    c_n_vinv x z

  let n_vlinearcombination ca (xa : t array) (z : t) =
    if Sundials_configuration.safe then Array.iter (check z) xa;
    c_n_vlinearcombination ca xa z

  let n_vscaleaddmulti aa (x : t) (ya : t array) (za : t array) =
    if Sundials_configuration.safe then
      (Array.iter (check x) ya; Array.iter (check x) za);
    c_n_vscaleaddmulti aa x ya za

  let n_vdotprodmulti (x : t) (ya : t array) dp =
    if Sundials_configuration.safe then Array.iter (check x) ya;
    c_n_vdotprodmulti x ya dp

  let n_vlinearsumvectorarray a (xa : t array) b (ya : t array) (za : t array) =
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (check x) ya;
          Array.iter (check x) za);
    c_n_vlinearsumvectorarray a xa b ya za

  let n_vscalevectorarray c (xa : t array) (za : t array) =
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (check x) za);
    c_n_vscalevectorarray c xa za

  let n_vconstvectorarray c (za : t array) =
    if Sundials_configuration.safe
    then (let z = Array.get za 0 in
          Array.iter (check z) za);
    c_n_vconstvectorarray c za

  let n_vwrmsnormvectorarray (xa : t array) (wa : t array) nrm =
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (check x) wa);
    c_n_vwrmsnormvectorarray xa wa nrm

  let n_vwrmsnormmaskvectorarray (xa : t array) (wa : t array) (id : t) nrm =
    if Sundials_configuration.safe
    then (Array.iter (check id) xa;
          Array.iter (check id) wa);
    c_n_vwrmsnormmaskvectorarray xa wa id nrm

  let n_vscaleaddmultivectorarray ra (xa : t array) (yaa : t array array)
                                     (zaa : t array array) =
    if Sundials_configuration.safe
    then (let x = Array.get xa 0 in
          Array.iter (check x) xa;
          Array.iter (Array.iter (check x)) yaa;
          Array.iter (Array.iter (check x)) zaa);
    c_n_vscaleaddmultivectorarray ra xa yaa zaa

  let n_vlinearcombinationvectorarray ca (xaa : t array array) (za : t array) =
    if Sundials_configuration.safe
    then (let z = Array.get za 0 in
          Array.iter (check z) za;
          Array.iter (Array.iter (check z)) xaa);
    c_n_vlinearcombinationvectorarray ca xaa za

  module Local = struct

    external c_n_vdotprod      : t -> t -> float
      = "sunml_nvec_any_n_vdotprodlocal"

    let n_vdotprod (x : t) (y : t) =
      if Sundials_configuration.safe then check x y;
      c_n_vdotprod x y

    external n_vmaxnorm        : t -> float
      = "sunml_nvec_any_n_vmaxnormlocal"

    external n_vmin            : t -> float
      = "sunml_nvec_any_n_vminlocal"

    external n_vl1norm         : t -> float
      = "sunml_nvec_any_n_vl1normlocal"

    external c_n_vinvtest      : t -> t -> bool
      = "sunml_nvec_any_n_vinvtestlocal"

    let n_vinvtest (x : t) (z : t) =
      if Sundials_configuration.safe then check x z;
      c_n_vinvtest x z

    external c_n_vconstrmask   : t -> t -> t -> bool
      = "sunml_nvec_any_n_vconstrmasklocal"

    let n_vconstrmask (c : t) (x : t) (m : t) =
      if Sundials_configuration.safe then (check c x; check c m);
      c_n_vconstrmask c x m

    external c_n_vminquotient  : t -> t -> float
      = "sunml_nvec_any_n_vminquotientlocal"

    let n_vminquotient (n : t) (d : t) =
      if Sundials_configuration.safe then check n d;
      c_n_vminquotient n d

    external c_n_vwsqrsum      : t -> t -> float
      = "sunml_nvec_any_n_vwsqrsumlocal"

    let n_vwsqrsum (x : t) (w : t) =
      if Sundials_configuration.safe then check x w;
      c_n_vwsqrsum x w

    external c_n_vwsqrsummask  : t -> t -> t -> float
      = "sunml_nvec_any_n_vwsqrsummasklocal"

    let n_vwsqrsummask (x : t) (w : t) (id : t) =
      if Sundials_configuration.safe then (check x w; check x id);
      c_n_vwsqrsummask x w id
  end
end (* }}} *)

