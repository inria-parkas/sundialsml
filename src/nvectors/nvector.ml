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
type ('data, 'kind) nvector = 'data * cnvec * ('data -> bool)
and ('data, 'kind) t = ('data, 'kind) nvector

let unwrap (payload, _, _) = payload

exception IncompatibleNvector

let check (_, _, checkfn) = (function (payload, _, _) ->
                              if not (checkfn payload) then
                                raise IncompatibleNvector)

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
  | Custom

external get_id : ('data, 'kind) t -> nvector_id
  = "sunml_nvec_get_id"

module type NVECTOR_OPS =
  sig
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
  end

module type NVECTOR =
  sig
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
  end

(* Hack to ensure that Sundials.c_init_module is executed so that the global
   exceptions are properly registered. *)
let e = Sundials.RecoverableFailure

