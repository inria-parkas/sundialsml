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

type lint_array = Sundials.LintArray.t
type real_array = Sundials.RealArray.t

module ArrayDenseMatrix = struct (* {{{ *)
  type t = Sundials.RealArray2.t

  let make = Sundials.RealArray2.make
  let create = Sundials.RealArray2.create
  let get = Sundials.RealArray2.get
  let set = Sundials.RealArray2.set

  let update a i j f = set a i j (f (get a i j))

  let set_to_zero x = Bigarray.Array2.fill (Sundials.RealArray2.unwrap x) 0.0

  let blit = Sundials.RealArray2.blit

  external scale : float -> t -> unit
      = "c_arraydensematrix_scale"

  external add_identity : t -> unit
      = "c_arraydensematrix_add_identity"

  external matvec' : t -> real_array -> real_array -> unit
      = "c_arraydensematrix_matvec"

  let matvec a ~x ~y = matvec' a x y

  external getrf : t -> lint_array -> unit
      = "c_arraydensematrix_getrf"

  external getrs : t -> lint_array -> real_array -> unit
      = "c_arraydensematrix_getrs"

  external getrs' : t -> lint_array -> real_array -> int -> unit
      = "c_arraydensematrix_getrs_off"

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
end (* }}} *)

module ArrayBandMatrix = struct (* {{{ *)
  type t = Sundials.RealArray2.t

  type smu = int
  type mu = int
  type ml = int

  let make n smu ml v =
    Sundials.RealArray2.make (smu + ml + 1) n v

  let create n smu ml =
    Sundials.RealArray2.create (smu + ml + 1) n

  let get a smu i j =
    Sundials.RealArray2.get a (i - j + smu) j

  let set a smu i j v =
    Sundials.RealArray2.set a (i - j + smu) j v

  let update a smu i j f =
    let k = i - j + smu in
    Sundials.RealArray2.set a k j (f (Sundials.RealArray2.get a k j))

  external copy' : t -> t -> int * int * int * int -> unit
      = "c_arraybandmatrix_copy"

  let blit a b a_smu b_smu copymu copyml
      = copy' a b (a_smu, b_smu, copymu, copyml)

  external scale' : float -> t -> int * int * int -> unit
      = "c_arraybandmatrix_scale"

  let scale c a smu mu ml = scale' c a (mu, ml, smu)

  external add_identity : t -> int -> unit
      = "c_arraybandmatrix_add_identity"

  external matvec' : t -> int * int * int -> real_array -> real_array -> unit
      = "c_arraybandmatrix_matvec"

  let matvec a smu mu ml ~x ~y = matvec' a (mu, ml, smu) x y

  external gbtrf' : t -> int * int * int -> lint_array -> unit
      = "c_arraybandmatrix_gbtrf"

  let gbtrf a smu mu ml p = gbtrf' a (mu, ml, smu) p

  external gbtrs'
      : t -> int * int -> lint_array -> real_array -> unit
      = "c_arraybandmatrix_gbtrs"

  let gbtrs a smu ml p b = gbtrs' a (smu, ml) p b

end (* }}} *)

