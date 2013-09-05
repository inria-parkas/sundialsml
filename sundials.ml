(***********************************************************************)
(*                                                                     *)
(*     OCaml interface to Sundials (serial) CVODE and IDA solvers      *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2013 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under the terms of the GNU Library General Public License, with    *)
(*  the special exception on linking described in file LICENSE.        *)
(*                                                                     *)
(***********************************************************************)

let extra_precision = ref false

let print_time (s1, s2) t =
  if !extra_precision
  then Printf.printf "%s%.15e%s" s1 t s2
  else Printf.printf "%s%e%s" s1 t s2

external format_float : string -> float -> string
    = "caml_format_float"

let floata = format_float "%a"

external get_big_real : unit -> float
    = "sundials_ml_big_real"

let big_real = get_big_real ()

external get_unit_roundoff : unit -> float
    = "sundials_ml_unit_roundoff"

let unit_roundoff = get_unit_roundoff ()

type real_array =
  (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t

let make_real_array =
  Bigarray.Array1.create Bigarray.float64 Bigarray.c_layout

type real_array2 =
  (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t

let make_real_array2 =
  Bigarray.Array2.create Bigarray.float64 Bigarray.c_layout

module Carray =
  struct
    type t = real_array

    let kind = Bigarray.float64
    let layout = Bigarray.c_layout
    let empty = Bigarray.Array1.create kind layout 0

    let create = Bigarray.Array1.create kind layout
    let of_array = Bigarray.Array1.of_array kind layout

    let fill = Bigarray.Array1.fill

    let init size x =
      let a = create size in
      fill a x;
      a

    let length = Bigarray.Array1.dim

    let blit = Bigarray.Array1.blit

    let of_carray src =
      let dst = create (length src) in
      blit src dst;
      dst

    let app f v =
      for i = 0 to (length v - 1) do
        f v.{i}
      done

    let map f v =
      for i = 0 to (length v - 1) do
        v.{i} <- f v.{i}
      done

    let appi f v =
      for i = 0 to (length v - 1) do
        f i v.{i}
      done

    let mapi f v =
      for i = 0 to (length v - 1) do
        v.{i} <- f i v.{i}
      done

    let print_with_time t v =
      print_time ("", "") t;
      if !extra_precision
      then app (Printf.printf "\t% .15e") v
      else app (Printf.printf "\t% e") v;
      print_newline ()
  end

(* root arrays *)

type lint_array = (int, Bigarray.int_elt, Bigarray.c_layout) Bigarray.Array1.t
let make_lint_array = Bigarray.Array1.create Bigarray.int Carray.layout

module Roots =
  struct
    type t = (int32, Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array1.t
    type val_array = Carray.t

    type root_event =
      | NoRoot
      | Rising
      | Falling

    let string_of_root_event e =
      match e with
      | NoRoot -> "NoRoot"
      | Rising -> "Rising"
      | Falling -> "Falling"

    let from_int32 x =
      if x = 1l then Rising else if x = -1l then Falling else NoRoot

    let from_int x =
      if x = 1 then Rising else if x = -1 then Falling else NoRoot

    let to_int32 x =
      match x with
      | NoRoot -> 0l
      | Rising -> 1l
      | Falling -> -1l

    let to_int x =
      match x with
      | NoRoot -> 0
      | Rising -> 1
      | Falling -> -1

    let reset v = Bigarray.Array1.fill v 0l

    let create n =
      let a = Bigarray.Array1.create Bigarray.int32 Carray.layout n in
      reset a;
      a

    let init n x =
      let a = Bigarray.Array1.create Bigarray.int32 Carray.layout n in
      Bigarray.Array1.fill a (to_int32 x);
      a

    let empty = create 0

    let length = Bigarray.Array1.dim

    let get roots i = roots.{i} <> 0l
    let get' roots i = from_int32 (roots.{i})
    let set a i v =
      Bigarray.Array1.set a i (to_int32 v)

    let rising  roots i = roots.{i} = 1l
    let falling roots i = roots.{i} = -1l

    let set_rising a i v = set a i (if v then Rising else NoRoot)
    let set_falling a i v = set a i (if v then Falling else NoRoot)

    let fold_left f a vs =
      let len = Bigarray.Array1.dim vs - 1 in
      let a = ref a in
      for i = 0 to len do
        a := f !a (Int32.to_int vs.{i})
      done;
      !a

    let appi f v =
      for i = 0 to (length v - 1) do
        f i (from_int32 v.{i})
      done
    let exists = fold_left (fun a x -> a || x <> 0) false

    let app f v =
      for i = 0 to (length v - 1) do
        f (from_int32 v.{i})
      done

    let print vs =
      app (fun x -> Printf.printf "\t%s" (string_of_root_event x)) vs;
      print_newline ()

  end

module RootDirs =
  struct
    type t = (int32, Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array1.t

    type root_direction =
      | Increasing
      | Decreasing
      | IncreasingOrDecreasing

    let string_of_root_direction d =
      match d with
      | Increasing -> "Increasing"
      | Decreasing -> "Decreasing"
      | IncreasingOrDecreasing -> "IncreasingOrDecreasing"

    let int32_of_root_direction x =
      match x with
      | Increasing -> 1l
      | Decreasing -> -1l
      | IncreasingOrDecreasing -> 0l
    let root_direction_of_int32 x =
      match x with
      | 1l -> Increasing
      | -1l -> Decreasing
      | _ -> IncreasingOrDecreasing

    let make n x =
      let a = Bigarray.Array1.create Bigarray.int32 Carray.layout n in
      Bigarray.Array1.fill a (int32_of_root_direction x);
      a

    let length a = Bigarray.Array1.dim a

    let create n = make n IncreasingOrDecreasing

    let create' n src =
      let a = Bigarray.Array1.create Bigarray.int32 Carray.layout n in
      if n > Array.length src
      then Bigarray.Array1.fill a
            (int32_of_root_direction IncreasingOrDecreasing);
      Array.iteri (fun i v -> a.{i} <- int32_of_root_direction v) src;
      a

    let set a i v = a.{i} <- int32_of_root_direction v
    let get a i = root_direction_of_int32 a.{i}

    let of_array src =
      let a = Bigarray.Array1.create Bigarray.int32 Carray.layout
                (Array.length src)
      in
      Array.iteri (fun i v -> a.{i} <- int32_of_root_direction v) src;
      a

  end

