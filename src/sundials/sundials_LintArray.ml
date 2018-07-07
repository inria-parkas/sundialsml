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

type t = (int, int_elt, c_layout) Array1.t

let create = Array1.create int c_layout

let make n x =
  let v = create n in
  Array1.fill v x;
  v

let empty = create 0

let pp fmt a =
  Format.pp_print_string fmt "[";
  Format.pp_open_hovbox fmt 0;
  for i = 0 to Array1.dim a - 1 do
    if i > 0 then (
      Format.pp_print_string fmt ";";
      Format.pp_print_space fmt ();
    );
    Format.fprintf fmt "% 6d" a.{i}
  done;
  Format.pp_close_box fmt ();
  Format.pp_print_string fmt "]"

let ppi ?(start="[") ?(stop="]") ?(sep=";")
        ?(item=fun fmt ->Format.fprintf fmt "%2d=% 6d") ()
        fmt a =
  Format.pp_print_string fmt start;
  Format.pp_open_hovbox fmt 0;
  for i = 0 to Array1.dim a - 1 do
    if i > 0 then (
      Format.pp_print_string fmt sep;
      Format.pp_print_space fmt ();
    );
    item fmt i a.{i}
  done;
  Format.pp_close_box fmt ();
  Format.pp_print_string fmt stop

