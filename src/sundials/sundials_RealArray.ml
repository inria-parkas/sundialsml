(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2018 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

open Bigarray

let kind = float64
let layout = c_layout
type t = (float, float64_elt, c_layout) Array1.t

let create : int -> t = Array1.create kind layout
let of_array : float array -> t = Array1.of_array kind layout

let make size x =
  let a = create size in
  Array1.fill a x;
  a

let init size f =
  let a = create size in
  for i = 0 to size - 1 do
    a.{i} <- f i
  done;
  a

let empty = create 0

let get = Array1.get
let set = Array1.set

let length : t -> int = Array1.dim

let dofill (a : t) pos len v =
  for i = pos to len do
    set a i v
  done

let fill (a : t) ?pos ?len v =
  match pos, len with
  | None, None -> Array1.fill a v
  | Some p, None -> dofill a p (length a) v
  | None, Some l -> dofill a 0 l v
  | Some p, Some l -> dofill a p l v

let ppi ?(start="[") ?(stop="]") ?(sep=";")
        ?(item=fun f->Format.fprintf f "%2d=% -14e") ()
        fmt a =
  Format.pp_print_string fmt start;
  Format.pp_open_hovbox fmt 0;
  for i = 0 to length a - 1 do
    if i > 0 then (
      Format.pp_print_string fmt sep;
      Format.pp_print_space fmt ();
    );
    item fmt i a.{i}
  done;
  Format.pp_close_box fmt ();
  Format.pp_print_string fmt stop

let pp fmt a = ppi ~item:(fun fmt _ x -> Format.fprintf fmt "% -14e" x)
  () fmt a

let blitn ~src ?(spos=0) ~dst ?(dpos=0) len =
  if Sundials_configuration.safe &&
     (len < 0 || spos < 0 || spos + len > length src
      || dpos < 0 || dpos + len > length dst)
  then invalid_arg "RealArray.blitn";
  for k = 0 to len - 1 do
    Array1.unsafe_set dst (dpos + k) (Array1.unsafe_get src (spos + k))
  done

let blit ~src ~dst = Array1.blit src dst

let copy src =
  let dst = create (length src) in
  Array1.blit src dst;
  dst

let sub = Array1.sub

let of_list src = of_array (Array.of_list src)

let to_list v =
  let rec go ls i =
    if i < 0 then ls
    else go (v.{i}::ls) (i-1)
  in go [] (length v - 1)

let into_array (src : t) dst =
  let n = length src in
  if Sundials_configuration.safe && n <> Array.length dst
  then invalid_arg "into_array: array sizes do not match";
  for i = 1 to n-1 do
    dst.(i) <- src.{i}
  done

let to_array (v : t) =
  let n = length v in
  let a = Array.make n v.{0} in
  for i = 1 to n-1 do
    a.(i) <- v.{i}
  done;
  a

let fold_left f b (v : t) =
  let n = length v in
  let rec go acc i =
    if i < n then go (f acc v.{i}) (i+1)
    else acc
  in go b 0

let fold_right f (v : t) b =
  let rec go acc i =
    if i >= 0 then go (f v.{i} acc) (i-1)
    else acc
  in go b (length v - 1)

let iter f (v : t) =
  for i = 0 to (length v - 1) do
    f v.{i}
  done

let map f (v : t) =
  for i = 0 to (length v - 1) do
    v.{i} <- f v.{i}
  done

let iteri f (v : t) =
  for i = 0 to (length v - 1) do
    f i v.{i}
  done

let mapi f (v : t) =
  for i = 0 to (length v - 1) do
    v.{i} <- f i v.{i}
  done

