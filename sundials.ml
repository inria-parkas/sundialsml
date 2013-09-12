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

    let of_list src = of_array (Array.of_list src)

    let to_list v =
      let rec go ls i =
        if i < 0 then ls
        else go (v.{i}::ls) (i-1)
      in go [] (length v - 1)

    let to_array v =
      let n = length v in
      let a = Array.make n v.{0} in
      for i = 1 to n-1 do
        a.(i) <- v.{i}
      done;
      a

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

(* Opaque arrays *)

module type ArrayBaseOps =
  sig
    type t
    type elt
    val create : int -> t
    val get : t -> int -> elt
    val set : t -> int -> elt -> unit
    val length : t -> int
  end

module ArrayLike (A : ArrayBaseOps) =
  struct
    include A

    let make n x =
      let a = create n in
      for i = 0 to n-1 do
        set a i x
      done;
      a

    let copy a =
      let n = length a in
      let b = create n in
      for i = 0 to n-1 do
        set b i (get a i)
      done;
      b

    let fold_left f x a =
      let n = length a in
      let rec go x i =
        if i >= n then x
        else go (f x (get a i)) (i+1)
      in go x 0

    let fold_right f a x =
      let rec go x i =
        if i < 0 then x
        else go (f (get a i) x) (i-1)
      in go x (length a - 1)

    let iteri f a =
      for i = 0 to length a - 1 do
        f i (get a i)
      done

    let iter f a = iteri (fun _ x -> f x) a

    let mapi f a =
      let n = length a in
      let b = create n in
      for i = 0 to n-1 do
        set b i (f i (get a i))
      done;
      b

    let map f a = mapi (fun _ x -> f x) a

    let mapi_overwrite f a =
      for i = 0 to length a - 1 do
        set a i (f i (get a i))
      done;
      a

    let map_overwrite f a = mapi_overwrite (fun _ x -> f x) a

    let of_array a =
      let n = Array.length a in
      let b = create n in
      for i = 0 to n-1 do
        set b i a.(i)
      done;
      b

    let of_list a =
      let n = List.length a in
      let b = create n in
      ignore (List.fold_left (fun i x -> set b i x; (i+1)) 0 a);
      b

    let to_array a =
      let n = length a in
      if n = 0 then [||]
      else
        let b = Array.make n (get a 0) in
        iteri (fun i ai -> b.(i) <- ai) a;
        b

    let to_list a = fold_right (fun x xs -> x::xs) a []

    let fill a offs len x =
      let n = length a in
      if len < 0 || offs < 0 || n < offs+len
      then invalid_arg "Sundials.ArrayLike.fill"
      else
        for i = offs to offs + len - 1 do
          set a i x
        done

    let fill_all a x = fill a 0 (length a) x

    let blit a oa b ob len =
      let na = length a
      and nb = length b in
      if len < 0 || oa < 0 || na < oa+len || ob < 0 || nb < ob+len
      then invalid_arg "Sundials.ArrayLike.blit"
      else
        for i = 0 to len-1 do
          set b (ob + i) (get a (oa + i))
        done

    let blit_all a b =
      for i = 0 to min (length a) (length b) - 1 do
        set b i (get a i)
      done

    let init n f =
      let a = create n in
      for i = 0 to n-1 do
        set a i (f i)
      done;
      a
  end

(* root arrays *)

type lint_array = (int, Bigarray.int_elt, Bigarray.c_layout) Bigarray.Array1.t
let make_lint_array = Bigarray.Array1.create Bigarray.int Carray.layout

module Roots =
  struct
    open Bigarray
    type t = (int32, int32_elt, c_layout) Array1.t
    type val_array = Carray.t

    type root_event =
      | NoRoot
      | Rising
      | Falling

    let root_event_of_int32 = function
      | 1l -> Rising
      | -1l -> Falling
      | 0l -> NoRoot
      | n ->
        failwith
          (Printf.sprintf
             "Sundials.Roots.root_event_of_int32: invalid root event %ld" n)

    let root_event_of_int = function
      | 1 -> Rising
      | -1 -> Falling
      | 0 -> NoRoot
      | n ->
        failwith
          ("Sundials.Roots.root_event_of_int: invalid root event "
           ^ string_of_int n)

    let int32_of_root_event x =
      match x with
      | NoRoot -> 0l
      | Rising -> 1l
      | Falling -> -1l

    let int_of_root_event x =
      match x with
      | NoRoot -> 0
      | Rising -> 1
      | Falling -> -1

    let reset v = Array1.fill v 0l

    let detected roots i = roots.{i} <> 0l
    let get roots i = root_event_of_int32 (roots.{i})
    let set a i v = a.{i} <- int32_of_root_event v

    let create n =
      let a = Array1.create int32 c_layout n in
      reset a;
      a

    let length = Array1.dim

    module A = ArrayLike (struct
      type t = (int32, int32_elt, c_layout) Array1.t
      and elt = root_event
      let get = get
      let set = set
      let create = create
      let length = length
    end)

    let empty = create 0
    let make n x = A.make n x
    let copy = A.copy
    let fold_left = A.fold_left
    let fold_right = A.fold_right

    let of_array = A.of_array
    let of_list = A.of_list
    let to_array = A.to_array
    let to_list = A.to_list

    let fill = A.fill
    let fill_all = A.fill_all
    let blit = A.blit
    let blit_all = A.blit_all

    let rising  roots i = roots.{i} = 1l
    let falling roots i = roots.{i} = -1l

    let set_rising a i v = set a i (if v then Rising else NoRoot)
    let set_falling a i v = set a i (if v then Falling else NoRoot)

    let appi = A.iteri
    let app = A.iter

    let exists a =
      let n = length a in
      let rec go i = a.{i} <> 0l || (i < n-1 && go (i+1)) in
      go 0

    let string_of_root_event e =
      match e with
      | NoRoot -> "NoRoot"
      | Rising -> "Rising"
      | Falling -> "Falling"

    let print vs =
      app (fun x -> Printf.printf "\t%s" (string_of_root_event x)) vs;
      print_newline ()

  end

module RootDirs =
  struct
    open Bigarray
    type t = (int32, int32_elt, c_layout) Array1.t

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
      | 0l -> IncreasingOrDecreasing
      | n ->
        failwith
          (Printf.sprintf
             "Sundials.Roots.root_direction_of_int32: \
              invalid root direction code %ld" n)

    let make n x =
      let a = Array1.create int32 c_layout n in
      Array1.fill a (int32_of_root_direction x);
      a

    let create n = make n IncreasingOrDecreasing

    let length a = Array1.dim a

    let create' n src =
      let a = Array1.create int32 c_layout n in
      if n > Array.length src
      then Array1.fill a (int32_of_root_direction IncreasingOrDecreasing);
      Array.iteri (fun i v -> a.{i} <- int32_of_root_direction v) src;
      a

    let set a i v = a.{i} <- int32_of_root_direction v
    let get a i = root_direction_of_int32 a.{i}

    module A = ArrayLike (
      struct
        type t = (int32, Bigarray.int32_elt, Bigarray.c_layout)
                  Bigarray.Array1.t
        and elt = root_direction
        let create = create
        let set = set
        let get = get
        let length = length
      end)

    let of_array = A.of_array
    let of_list  = A.of_list
    let to_array = A.to_array
    let to_list  = A.to_list
    let fill     = A.fill
    let fill_all = A.fill_all
    let blit     = A.blit
    let blit_all = A.blit_all

    let init = A.init

  end

