(***********************************************************************)
(*                                                                     *)
(*     OCaml interface to Sundials (serial) CVODE and IDA solvers      *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2013 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a BSD 2-Clause License, refer to the file LICENSE.           *)
(*                                                                     *)
(***********************************************************************)

type lint_array = Sundials.lint_array
type real_array = Sundials.real_array

(* direct linear solvers functions *)

(* note: uses DENSE_ELEM rather than the more efficient DENSE_COL. *)
module Densematrix =
  struct
    type t

    external new_dense_mat  : int * int -> t
        = "c_densematrix_new_dense_mat"

    external print_mat      : t -> unit
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

    external get : t -> (int * int) -> float
        = "c_densematrix_get"

    external set : t -> (int * int) -> float -> unit
        = "c_densematrix_set"
  end

module Directdensematrix =
  struct
    type t

    external new_dense_mat  : int * int -> t
        = "c_densematrix_direct_new_dense_mat"

    external get : t -> (int * int) -> float
        = "c_densematrix_direct_get"

    external set : t -> (int * int) -> float -> unit
        = "c_densematrix_direct_set"

    external copy  : t -> t -> int * int -> unit
        = "c_densematrix_direct_copy"

    external scale : float -> t -> int * int -> unit
        = "c_densematrix_direct_scale"

    external add_identity : t -> int -> unit
        = "c_densematrix_direct_add_identity"

    external getrf : t -> int * int -> lint_array -> unit
        = "c_densematrix_direct_getrf"

    external getrs : t -> int -> lint_array -> real_array -> unit
        = "c_densematrix_direct_getrs"

    external potrf : t -> int -> unit
        = "c_densematrix_direct_potrf"

    external potrs : t -> int -> real_array -> unit
        = "c_densematrix_direct_potrs"

    external geqrf : t -> int * int -> real_array -> real_array -> unit
        = "c_densematrix_direct_geqrf"

    external ormqr'
        : t -> int * int
          -> (real_array * real_array * real_array * real_array)
          -> unit
        = "c_densematrix_direct_ormqr"

    let ormqr ~a ~mn ~beta ~v ~w ~work = ormqr' a mn (beta, v, w, work)
  end

module Bandmatrix =
  struct
    type t

    external new_band_mat : int * int * int * int -> t
        = "c_bandmatrix_new_band_mat"

    external print_mat : t -> unit
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

    external get : t -> (int * int) -> float
        = "c_bandmatrix_get"

    external set : t -> (int * int) -> float -> unit
        = "c_bandmatrix_set"

    module Col =
      struct
        type c

        external get_col : t -> int -> c
            = "c_bandmatrix_col_get_col"

        external get : c -> (int * int) -> float
            = "c_bandmatrix_col_get"

        external set : c -> (int * int) -> float -> unit
            = "c_bandmatrix_col_set"
      end

  end

module Directbandmatrix =
  struct
    type t

    external new_band_mat : int * int * int -> t
        = "c_bandmatrix_direct_new_band_mat"

    external get : t -> (int * int) -> float
        = "c_densematrix_direct_get"
        (* NB: same as densematrix_direct *)

    external set : t -> (int * int) -> float -> unit
        = "c_densematrix_direct_set"
        (* NB: same as densematrix_direct *)

    external copy' : t -> t -> int * int * int * int * int -> unit
        = "c_bandmatrix_direct_copy"

    let copy a b n a_smu b_smu copymu copyml
        = copy' a b (n, a_smu, b_smu, copymu, copyml)

    external scale' : float -> t -> int * int * int * int -> unit
        = "c_bandmatrix_direct_scale"

    let scale c a n mu ml smu = scale' c a (n, mu, ml, smu)

    external add_identity : t -> int -> int -> unit
        = "c_bandmatrix_direct_add_identity"

    external gbtrf' : t -> int * int * int * int -> lint_array -> unit
        = "c_bandmatrix_direct_gbtrf"

    let gbtrf a n mu ml smu p = gbtrf' a (n, mu, ml, smu) p

    external gbtrs'
        : t -> int * int * int -> lint_array -> real_array -> unit
        = "c_bandmatrix_direct_gbtrs"

    let gbtrs a n smu ml p b = gbtrs' a (n, smu, ml) p b
  end

