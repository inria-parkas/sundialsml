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

(* sparse linear solvers functions *)

(* note: uses DENSE_ELEM rather than the more efficient DENSE_COL. *)
module SparseMatrix =
  struct
    (* Must correspond with sls_ml.h:sls_sparsematrix_index *)
    include Sls_impl

    type csc
    type csr

    type t_csc = csc t
    type t_csr = csr t

    exception Invalidated

    external c_create : int -> int -> int -> sformat -> 'f t
        = "c_sparsematrix_new_sparse_mat"

    let create_csc m n nnz =
      if Sundials_config.safe && (m <= 0 || n <= 0 || nnz <= 0)
      then failwith "M, N and NNZ must all be positive";
      c_create m n nnz CSC_MAT

    let create_csr m n nnz =
      if Sundials_config.safe && (m <= 0 || n <= 0 || nnz <= 0)
      then failwith "M, N and NNZ must all be positive";
      if Sundials_config.safe then
        (match Sundials.sundials_version with
         | 2,v,_ when v < 7 -> raise Sundials.NotImplementedBySundialsVersion
         | _ -> ());
      c_create m n nnz CSR_MAT

    (* Allowing direct access is not safe because the underlying data may
       have been invalidated. Invalidated by setting the dims of the
       bigarray to 0 is tempting but error prone and the mechanism can
       be circumvented using Bigarray.Array2.sub (see:
         https://groups.google.com/d/msg/fa.caml/ROr_PifT_44/aqQ8Z0TWzH8J).
       Worse, the underlying arrays may be reallocated by calls to
       add_identity, blit, or add. *)

    external c_size : slsmat -> (int * int * int)
        = "c_sparsematrix_size"

    let size { slsmat; valid } =
      if Sundials_config.safe && not valid then raise Invalidated;
      c_size slsmat

    external c_print        : Sundials.Logfile.t -> slsmat -> unit
        = "c_sparsematrix_print_mat"

    let print logfile { slsmat; valid } =
      if Sundials_config.safe && not valid then raise Invalidated;
      c_print logfile slsmat

    external c_set_to_zero  : slsmat -> unit
        = "c_sparsematrix_set_to_zero"

    let set_to_zero { slsmat; valid } =
      if Sundials_config.safe && not valid then raise Invalidated;
      c_set_to_zero slsmat

    external c_realloc      : 'f t -> int -> unit
        = "c_sparsematrix_realloc"

    let realloc ({ valid } as d) =
      if Sundials_config.safe && not valid then raise Invalidated;
      c_realloc d

    external c_convert_dls : sformat -> Dls_impl.dlsmat -> 'f t
        = "c_sparsematrix_convert_dls"

    let csc_from_dense { Dls_impl.DenseTypes.valid = valid;
                         Dls_impl.DenseTypes.dlsmat = dlsmat } =
      if Sundials_config.safe && not valid then raise Invalidated;
      c_convert_dls CSC_MAT dlsmat

    let csr_from_dense { Dls_impl.DenseTypes.valid = valid;
                         Dls_impl.DenseTypes.dlsmat = dlsmat } =
      if Sundials_config.safe && not valid then raise Invalidated;
      c_convert_dls CSR_MAT dlsmat

    let csc_from_band { Dls_impl.BandTypes.valid = valid;
                    Dls_impl.BandTypes.dlsmat = dlsmat } =
      if Sundials_config.safe && not valid then raise Invalidated;
      c_convert_dls CSC_MAT dlsmat

    let csr_from_band { Dls_impl.BandTypes.valid = valid;
                    Dls_impl.BandTypes.dlsmat = dlsmat } =
      if Sundials_config.safe && not valid then raise Invalidated;
      c_convert_dls CSR_MAT dlsmat

    let set_col { idxptrs; valid } j idx =
      if Sundials_config.safe && not valid then raise Invalidated;
      idxptrs.{j} <- Int32.of_int idx

    let get_col { idxptrs; valid } j =
      if Sundials_config.safe && not valid then raise Invalidated;
      Int32.to_int idxptrs.{j}

    let set_row = set_col
    let get_row = get_col

    let set { idxvals; data; valid } idx i v =
      if Sundials_config.safe && not valid then raise Invalidated;
      idxvals.{idx} <- Int32.of_int i;
      data.{idx} <- v

    let set_data { data; valid } idx v =
      if Sundials_config.safe && not valid then raise Invalidated;
      data.{idx} <- v

    let set_rowval { idxvals; valid } idx i =
      if Sundials_config.safe && not valid then raise Invalidated;
      idxvals.{idx} <- Int32.of_int i

    let set_colval = set_rowval

    let get { idxvals; data; valid } idx =
      if Sundials_config.safe && not valid then raise Invalidated;
      Int32.to_int idxvals.{idx}, data.{idx}

    let get_rowval { idxvals; valid } idx =
      if Sundials_config.safe && not valid then raise Invalidated;
      Int32.to_int idxvals.{idx}

    let get_colval = get_rowval

    let get_data { data; valid } idx =
      if Sundials_config.safe && not valid then raise Invalidated;
      data.{idx}

    let pp fmt mat =
      if Sundials_config.safe && not mat.valid then raise Invalidated;
      let m, n, nnz = size mat in
      let end_row = ref false in
      let idxrole = match mat.sformat with
                    | CSC_MAT -> "col"
                    | CSR_MAT -> "row" in

      Format.pp_print_string fmt "[";
      Format.pp_open_vbox fmt 0;
      for j = 0 to n - 1 do
        let p, np = get_col mat j, get_col mat (j + 1) in
        if p < np then begin
          if !end_row then begin
            Format.pp_print_string fmt ";";
            Format.pp_print_cut fmt ();
            Format.pp_close_box fmt ()
          end;
          Format.pp_open_hovbox fmt 4;
          Format.fprintf fmt "%s %2d: " idxrole j;
          Format.pp_print_space fmt ();

          for i = p to np - 1 do
            if i > p then Format.pp_print_space fmt ();
            let r, v = get mat i in
            Format.fprintf fmt "%2d=% -15e" r v
          done;
          end_row := true
        end;
      done;
      if !end_row then Format.pp_close_box fmt ();
      Format.pp_close_box fmt ();
      Format.pp_print_string fmt "]"

    external c_add_identity : 'f t -> unit
        = "c_sparsematrix_add_identity"

    let add_identity ({ valid } as d) =
      if Sundials_config.safe && not valid then raise Invalidated;
      c_add_identity d

    external c_copy : slsmat -> 'f t -> unit
        = "c_sparsematrix_copy"

    let blit ({ slsmat=a; valid=valid1 })
             ({ valid=valid2 } as bm) =
      if Sundials_config.safe && not (valid1 && valid2) then raise Invalidated;
      c_copy a bm

    external c_scale  : float -> slsmat -> unit
        = "c_sparsematrix_scale"

    let scale a { slsmat; valid } =
      if Sundials_config.safe && not valid then raise Invalidated;
      c_scale a slsmat

    external c_add  : 'f t -> slsmat -> unit
        = "c_sparsematrix_add"

    let add ({ valid=valid1 } as am) ({ slsmat=b; valid=valid2 }) =
      if Sundials_config.safe && not (valid1 && valid2) then raise Invalidated;
      c_add am b

    external c_matvec  : slsmat -> real_array -> real_array -> unit
        = "c_sparsematrix_matvec"

    let matvec ({ slsmat; valid }) x y =
      if Sundials_config.safe && not valid then raise Invalidated;
      c_matvec slsmat x y
    
    let make_csc m n nnz =
      let r = create_csc m n nnz in
      set_to_zero r;
      r

    let make_csr m n nnz =
      let r = create_csr m n nnz in
      set_to_zero r;
      r

  end

