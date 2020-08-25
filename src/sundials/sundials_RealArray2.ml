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

let make_data = Array2.create float64 c_layout

type data =
  (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array2.t

type t = data * Obj.t

external wrap : data -> t
  = "sunml_sundials_realarray2_wrap"

let unwrap = fst

let create nr nc =
  let d = Array2.create float64 c_layout nc nr
  in wrap d

let make nr nc v =
  let d = Array2.create float64 c_layout nc nr
  in
  Array2.fill d v;
  wrap d

let empty = create 0 0

let size a =
  let d = unwrap a in
  (Array2.dim2 d, Array2.dim1 d)

let fill (a, _) v = Array2.fill a v

let ppi ?(start="[") ?(rowstart="[") ?(stop="]") ?(rowstop="]")
        ?(sep=";") ?(rowsep=";")
        ?(item=fun f -> Format.fprintf f "(%2d,%2d)=% -14e") ()
        fmt a =
  let d = unwrap a in
  let ni, nj = Array2.dim2 d - 1, Array2.dim1 d - 1 in
  Format.pp_print_string fmt start;
  Format.pp_open_vbox fmt 0;
  for i = 0 to ni do
    if i > 0 then (
      Format.pp_print_string fmt rowsep;
      Format.pp_print_cut fmt ()
    );

    Format.pp_print_string fmt rowstart;
    Format.pp_open_hovbox fmt 0;
    for j = 0 to nj do
      if j > 0 then (
        Format.pp_print_string fmt sep;
        Format.pp_print_space fmt ();
      );
      item fmt i j d.{j, i}
    done;
    Format.pp_close_box fmt ();
    Format.pp_print_string fmt rowstop

  done;
  Format.pp_close_box fmt ();
  Format.pp_print_string fmt stop

let pp fmt a = ppi ~item:(fun fmt _ _ x -> Format.fprintf fmt "% -14e" x) ()
    fmt a

let get x i j = Array2.get (unwrap x) j i
let set x i j = Array2.set (unwrap x) j i

let col x j = Array2.slice_left (unwrap x) j

let copy a =
  let d = unwrap a in
  let c = Array2.dim1 d in
  let r = Array2.dim2 d in
  let d' = Array2.create float64 c_layout c r in
  Array2.blit d d';
  wrap d'

let blit ~src ~dst =
  Array2.blit (unwrap src) (unwrap dst)

