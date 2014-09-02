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

let _ =
  List.iter (fun (nm, ex) -> Callback.register_exception nm ex)
  [
    ("dls_ZeroDiagonalElement",     ZeroDiagonalElement 0);
  ]

(* note: uses DENSE_ELEM rather than the more efficient DENSE_COL. *)
module DenseMatrix =
  struct
    type t

    external create : int -> int -> t
        = "c_densematrix_new_dense_mat"

    external size : t -> (int * int)
        = "c_densematrix_size"

    external print          : t -> unit
        = "c_densematrix_print_mat"

    external set_to_zero    : t -> unit
        = "c_densematrix_set_to_zero"

    external add_identity   : t -> unit
        = "c_densematrix_add_identity"

    external copy     : t -> t -> unit
        = "c_densematrix_copy"

    external scale    : float -> t -> unit
        = "c_densematrix_scale"

    external getrf    : t -> lint_array -> unit
        = "c_densematrix_getrf"

    external getrs    : t -> lint_array -> real_array -> unit
        = "c_densematrix_getrs"

    external potrf    : t -> unit
        = "c_densematrix_potrf"

    external potrs    : t -> real_array -> unit
        = "c_densematrix_potrs"

    external geqrf    : t -> real_array -> real_array -> unit
        = "c_densematrix_geqrf"

    external ormqr'
        : t -> (real_array * real_array * real_array * real_array) -> unit
        = "c_densematrix_ormqr"

    let ormqr ~a ~beta ~v ~w ~work = ormqr' a (beta, v, w, work)

    external get : t -> int -> int -> float
        = "c_densematrix_get"

    external set : t -> int -> int -> float -> unit
        = "c_densematrix_set"

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
    type t

    external create : int -> int -> int -> int -> t
        = "c_bandmatrix_new_band_mat"

    external size : t -> (int * int * int * int)
        = "c_bandmatrix_size"

    external print          : t -> unit
        = "c_densematrix_print_mat"
          (* NB: same as densematrix *)

    external set_to_zero    : t -> unit
        = "c_densematrix_set_to_zero"
          (* NB: same as densematrix *)

    external add_identity : t -> unit
        = "c_densematrix_add_identity"
          (* NB: same as densematrix *)

    external copy : t -> t -> int -> int -> unit
        = "c_bandmatrix_copy"

    external scale : float -> t -> unit
        = "c_bandmatrix_scale"

    external gbtrf : t -> lint_array -> unit
        = "c_bandmatrix_gbtrf"

    external gbtrs : t -> lint_array -> real_array -> unit
        = "c_bandmatrix_gbtrs"

    external get : t -> int -> int -> float
        = "c_bandmatrix_get"

    external set : t -> int -> int -> float -> unit
        = "c_bandmatrix_set"

    let make n mu ml smu v =
      let r = create n mu ml smu in
      for i = 0 to n - 1 do
        for j = (max 0 (i - 1)) to (min n (i + 1)) - 1 do
          set r i j v
        done
      done;
      r

    module Col =
      struct
        type c

        external get_col : t -> int -> c
            = "c_bandmatrix_col_get_col"

        external get : c -> int -> int -> float
            = "c_bandmatrix_col_get"

        external set : c -> int -> int -> float -> unit
            = "c_bandmatrix_col_set"
      end

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

