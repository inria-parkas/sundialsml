(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a BSD 2-Clause License, refer to the file LICENSE.           *)
(*                                                                     *)
(***********************************************************************)

type lint_array = Sundials.LintArray.t
type real_array = Sundials.RealArray.t

(* direct linear solvers functions *)

exception ZeroDiagonalElement of int

(* note: uses DENSE_ELEM rather than the more efficient DENSE_COL. *)
module DenseMatrix =
  struct
    type data = (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t

    (* (Bigarray into data, DlsMat with finalizer) *)
    type t = data * Obj.t

    exception Relinquished
    let assert_valid ba =
      if Bigarray.Array2.dim1 ba = 0 then raise Relinquished

    external c_create : int -> int -> t
        = "c_densematrix_new_dense_mat"

    let create i j =
      if i <= 0 || j <= 0 then failwith "Both M and N must be positive";
      c_create i j

    let unwrap = fst

    external c_relinquish : Obj.t -> unit
        = "c_dls_relinquish"

    let relinquish (_, v) = c_relinquish v

    external c_size : Obj.t -> (int * int)
        = "c_densematrix_size"

    let size (ba, v) =
      assert_valid ba;
      c_size v

    external c_print        : Obj.t -> unit
        = "c_densematrix_print_mat"

    let print (ba, v) =
      assert_valid ba;
      c_print v

    external c_set_to_zero  : Obj.t -> unit
        = "c_densematrix_set_to_zero"

    let set_to_zero (ba, v) =
      assert_valid ba;
      c_set_to_zero v

    external c_add_identity : Obj.t -> unit
        = "c_densematrix_add_identity"

    let add_identity (ba, v) =
      assert_valid ba;
      c_add_identity v

    external c_copy     : Obj.t -> Obj.t -> unit
        = "c_densematrix_copy"

    let copy (ba1, v1) (ba2, v2) =
      assert_valid ba1;
      assert_valid ba2;
      c_copy v1 v2

    external c_scale  : float -> Obj.t -> unit
        = "c_densematrix_scale"

    let scale a (ba, v) =
      assert_valid ba;
      c_scale a v

    external c_getrf  : Obj.t -> lint_array -> unit
        = "c_densematrix_getrf"

    let getrf (ba, v) la =
      assert_valid ba;
      c_getrf v la

    external c_getrs  : Obj.t -> lint_array -> real_array -> unit
        = "c_densematrix_getrs"

    let getrs (ba, v) la ra =
      assert_valid ba;
      c_getrs v la ra

    external c_potrf  : Obj.t -> unit
        = "c_densematrix_potrf"

    let potrf (ba, v) =
      assert_valid ba;
      c_potrf v

    external c_potrs  : Obj.t -> real_array -> unit
        = "c_densematrix_potrs"

    let potrs (ba, v) ra =
      assert_valid ba;
      c_potrs v ra

    external c_geqrf  : Obj.t -> real_array -> real_array -> unit
        = "c_densematrix_geqrf"

    let geqrf (ba, v) ra1 ra2 =
      assert_valid ba;
      c_geqrf v ra1 ra2

    external c_ormqr
        : Obj.t -> (real_array * real_array * real_array * real_array) -> unit
        = "c_densematrix_ormqr"

    let ormqr ~a ~beta ~v ~w ~work =
      assert_valid (fst a);
      c_ormqr (snd a) (beta, v, w, work)

    (*
    external c_get : Obj.t -> int -> int -> float
        = "c_densematrix_get"

    let get (ba, v) i j =
      assert_valid ba;
      c_get v i j
    *)
    let get ((ba : data), v) i j = ba.{j, i}

    (*
    external c_set : Obj.t -> int -> int -> float -> unit
        = "c_densematrix_set"

    let set (ba, v) i j e =
      assert_valid ba;
      c_set v i j e
    *)
    let set ((ba : data), v) i j v = ba.{j, i} <- v

    let make m n v =
      let r = create m n in
      for i = 0 to m - 1 do
        for j = 0 to n - 1 do
          set r i j v
        done
      done;
      r

  end

module ArrayDenseMatrix =
  struct
    type t = Sundials.RealArray2.t

    let make = Sundials.RealArray2.make
    let create = Sundials.RealArray2.create
    let get = Sundials.RealArray2.get
    let set = Sundials.RealArray2.set

    let set_to_zero x = Bigarray.Array2.fill (Sundials.RealArray2.unwrap x) 0.0

    let copy = Sundials.RealArray2.copyinto

    external scale : float -> t -> unit
        = "c_arraydensematrix_scale"

    external add_identity : t -> unit
        = "c_arraydensematrix_add_identity"

    external getrf : t -> lint_array -> unit
        = "c_arraydensematrix_getrf"

    external getrs : t -> lint_array -> real_array -> unit
        = "c_arraydensematrix_getrs"

    external potrf : t -> unit
        = "c_arraydensematrix_potrf"

    external potrs : t -> real_array -> unit
        = "c_arraydensematrix_potrs"

    external geqrf : t -> real_array -> real_array -> unit
        = "c_arraydensematrix_geqrf"

    external ormqr'
        : t -> (real_array * real_array * real_array * real_array) -> unit
        = "c_arraydensematrix_ormqr"

    let ormqr ~a ~beta ~v ~w ~work = ormqr' a (beta, v, w, work)
  end

module BandMatrix =
  struct
    type data = (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t

    (* (Bigarray into data, DlsMat with finalizer, smu) *)
    type t = data * Obj.t * int

    exception Relinquished
    let assert_valid ba =
      if Bigarray.Array2.dim1 ba = 0 then raise Relinquished

    external c_create : int -> int -> int -> int -> t
        = "c_bandmatrix_new_band_mat"

    let create i j =
      if i <= 0 || j <= 0 then failwith "Both M and N must be positive";
      c_create i j

    let unwrap (ba, _, _) = ba

    external c_relinquish : Obj.t -> unit
        = "c_dls_relinquish"

    let relinquish (_, v, _) = c_relinquish v

    external c_size : Obj.t -> (int * int * int * int)
        = "c_bandmatrix_size"

    let size (ba, v, _) =
      assert_valid ba;
      c_size v

    external c_print          : Obj.t -> unit
        = "c_densematrix_print_mat"
          (* NB: same as densematrix *)

    let print (ba, v, _) =
      assert_valid ba;
      c_print v

    external c_set_to_zero    : Obj.t -> unit
        = "c_densematrix_set_to_zero"
          (* NB: same as densematrix *)

    let set_to_zero (ba, v, _) =
      assert_valid ba;
      c_set_to_zero v

    external c_add_identity : Obj.t -> unit
        = "c_densematrix_add_identity"
          (* NB: same as densematrix *)

    let add_identity (ba, v, _) =
      assert_valid ba;
      c_add_identity v

    external c_copy : Obj.t -> Obj.t -> int -> int -> unit
        = "c_bandmatrix_copy"

    let copy (ba1, v1, _) (ba2, v2, _) copymu copyml =
      assert_valid ba1;
      assert_valid ba2;
      c_copy v1 v2 copymu copyml

    external c_scale : float -> Obj.t -> unit
        = "c_bandmatrix_scale"

    let scale a (ba, v, _) =
      assert_valid ba;
      c_scale a v

    external c_gbtrf : Obj.t -> lint_array -> unit
        = "c_bandmatrix_gbtrf"

    let gbtrf (ba, v, _) la =
      assert_valid ba;
      c_gbtrf v la

    external c_gbtrs : Obj.t -> lint_array -> real_array -> unit
        = "c_bandmatrix_gbtrs"

    let gbtrs (ba, v, _) la ra =
      assert_valid ba;
      c_gbtrs v la ra

    (*
    external c_get : Obj.t -> int -> int -> float
        = "c_bandmatrix_get"

    let get (ba, v, _) i j =
      assert_valid ba;
      c_get v i j
    *)
    let get ((ba : data), _, s_mu) i j = ba.{j, i - j + s_mu}

    (*
    external c_set : Obj.t -> int -> int -> float -> unit
        = "c_bandmatrix_set"

    let set (ba, v, _) i j e =
      assert_valid ba;
      c_set v i j e
    *)
    let set ((ba : data), _, s_mu) i j v = ba.{j, i - j + s_mu} <- v

    let make n mu ml smu v =
      let r = create n mu ml smu in
      for i = 0 to n - 1 do
        for j = (max 0 (i - 1)) to (min n (i + 1)) - 1 do
          set r i j v
        done
      done;
      r
  end

module ArrayBandMatrix =
  struct
    type t = Sundials.RealArray2.t

    let make n smu ml v =
      Sundials.RealArray2.make (smu + ml + 1) n v

    let create n smu ml =
      Sundials.RealArray2.create (smu + ml + 1) n

    let get a smu i j =
      Sundials.RealArray2.get a (i - j + smu) j

    let set a smu i j v =
      Sundials.RealArray2.set a (i - j + smu) j v

    external copy' : t -> t -> int * int * int * int -> unit
        = "c_arraybandmatrix_copy"

    let copy a b a_smu b_smu copymu copyml
        = copy' a b (a_smu, b_smu, copymu, copyml)

    external scale' : float -> t -> int * int * int -> unit
        = "c_arraybandmatrix_scale"

    let scale c a mu ml smu = scale' c a (mu, ml, smu)

    external add_identity : t -> int -> unit
        = "c_arraybandmatrix_add_identity"

    external gbtrf' : t -> int * int * int -> lint_array -> unit
        = "c_arraybandmatrix_gbtrf"

    let gbtrf a mu ml smu p = gbtrf' a (mu, ml, smu) p

    external gbtrs'
        : t -> int * int -> lint_array -> real_array -> unit
        = "c_arraybandmatrix_gbtrs"

    let gbtrs a smu ml p b = gbtrs' a (smu, ml) p b
  end


(* Let C code know about some of the values in this module.  *)
external c_init_module : exn array -> unit =
  "c_dls_init_module"

let _ =
  c_init_module
    (* Exceptions must be listed in the same order as
       dls_exn_index.  *)
    [|ZeroDiagonalElement 0|]
