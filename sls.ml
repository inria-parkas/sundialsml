(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2015 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

type int_array =
  (int32, Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array1.t
type real_array = Sundials.RealArray.t

(* sparse linear solvers functions *)

(* note: uses DENSE_ELEM rather than the more efficient DENSE_COL. *)
module SparseMatrix =
  struct
    (* Must correspond with sls_ml.h:sls_sparsematrix_index *)
    type t = {
      colptrs : int_array;
      rowvals : int_array;
      data    : real_array;
      slsmat  : Obj.t;
      mutable valid : bool;
    }

    exception Invalidated

    external c_create : int -> int -> int -> t
        = "c_sparsematrix_new_sparse_mat"

    let create m n nnz =
      if Sundials_config.safe && (m <= 0 || n <= 0 || nnz <= 0)
      then failwith "M, N and NNZ must all be positive";
      c_create m n nnz

    (* Allowing direct access is not safe because the underlying data may
       have been invalidated. Invalidated by setting the dims of the
       bigarray to 0 is tempting but error prone and the mechanism can
       be circumvented using Bigarray.Array2.sub (see:
         https://groups.google.com/d/msg/fa.caml/ROr_PifT_44/aqQ8Z0TWzH8J).
       Worse, the underlying arrays may be reallocated by calls to
       add_identity, blit, or add. *)

    let invalidate v = v.valid <- false

    external c_size : Obj.t -> (int * int * int)
        = "c_sparsematrix_size"

    let size { slsmat; valid } =
      if Sundials_config.safe && not valid then raise Invalidated;
      c_size slsmat

    external c_print        : Obj.t -> unit
        = "c_sparsematrix_print_mat"

    let print { slsmat; valid } =
      if Sundials_config.safe && not valid then raise Invalidated;
      c_print slsmat

    external c_set_to_zero  : Obj.t -> unit
        = "c_sparsematrix_set_to_zero"

    let set_to_zero { slsmat; valid } =
      if Sundials_config.safe && not valid then raise Invalidated;
      c_set_to_zero slsmat

    external c_convert_dls : Obj.t -> t
        = "c_sparsematrix_convert_dls"

    let from_dense { Dls_impl.DenseTypes.valid = valid;
                     Dls_impl.DenseTypes.dlsmat = dlsmat } =
      if Sundials_config.safe && not valid then raise Invalidated;
      c_convert_dls dlsmat

    let from_band { Dls_impl.BandTypes.valid = valid;
                    Dls_impl.BandTypes.dlsmat = dlsmat } =
      if Sundials_config.safe && not valid then raise Invalidated;
      c_convert_dls dlsmat

    let set_col { colptrs; valid } j idx =
      if Sundials_config.safe && not valid then raise Invalidated;
      colptrs.{j} <- Int32.of_int idx

    let get_col { colptrs; valid } j =
      if Sundials_config.safe && not valid then raise Invalidated;
      Int32.to_int colptrs.{j}

    let set { rowvals; data; valid } idx i v =
      if Sundials_config.safe && not valid then raise Invalidated;
      rowvals.{idx} <- Int32.of_int i;
      data.{idx} <- v

    let get { rowvals; data; valid } idx =
      if Sundials_config.safe && not valid then raise Invalidated;
      Int32.to_int rowvals.{idx}, data.{idx}

    external c_add_identity : t -> unit
        = "c_sparsematrix_add_identity"

    let add_identity ({ valid } as d) =
      if Sundials_config.safe && not valid then raise Invalidated;
      c_add_identity d

    external c_copy : Obj.t -> t -> unit
        = "c_sparsematrix_copy"

    let blit ({ slsmat=a; valid=valid1 })
             ({ valid=valid2 } as bm) =
      if Sundials_config.safe && not (valid1 && valid2) then raise Invalidated;
      c_copy a bm

    external c_scale  : float -> Obj.t -> unit
        = "c_sparsematrix_scale"

    let scale a { slsmat; valid } =
      if Sundials_config.safe && not valid then raise Invalidated;
      c_scale a slsmat

    external c_add  : t -> Obj.t -> unit
        = "c_sparsematrix_add"

    let add ({ valid=valid1 } as am) ({ slsmat=b; valid=valid2 }) =
      if Sundials_config.safe && not (valid1 && valid2) then raise Invalidated;
      c_add am b

    external c_matvec  : Obj.t -> real_array -> real_array -> unit
        = "c_sparsematrix_matvec"

    let matvec ({ slsmat; valid }) x y =
      if Sundials_config.safe && not valid then raise Invalidated;
      c_matvec slsmat x y
    
    let make m n nnz =
      let r = create m n nnz in
      set_to_zero r;
      r

  end

