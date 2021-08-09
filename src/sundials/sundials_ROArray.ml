(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2021 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

include Array
type 'a t = 'a array

let from_array = Array.copy
let to_array = Array.copy

let iteri2 f a b =
  let la = length a in
  if la <> length b then
    invalid_arg "ROArray.iter2i: arrays must have the same length"
  else
    for i = 0 to la - 1 do
      f i (unsafe_get a i) (unsafe_get b i)
    done

let iter3 f a b c =
  let la = length a in
  if la <> length b || la <> length c then
    invalid_arg "ROArray.iter3: arrays must have the same length"
  else
    for i = 0 to la - 1 do
      f (unsafe_get a i) (unsafe_get b i) (unsafe_get c i)
    done

let iteri3 f a b c =
  let la = length a in
  if la <> length b || la <> length c then
    invalid_arg "ROArray.iter3i: arrays must have the same length"
  else
    for i = 0 to la - 1 do
      f i (unsafe_get a i) (unsafe_get b i) (unsafe_get c i)
    done

let map3 f a b c =
  let la = length a in
  if la <> length b || la <> length c then
    invalid_arg "ROArray.map3: arrays must have the same length"
  else begin
    if la = 0 then [||] else begin
      let r = make la (f (unsafe_get a 0) (unsafe_get b 0) (unsafe_get c 0)) in
      for i = 1 to la - 1 do
        unsafe_set r i (f (unsafe_get a i) (unsafe_get b i) (unsafe_get c i))
      done;
      r
    end
  end

let fold_left2 f x a b =
  let la = length a in
  if la <> length b then
    invalid_arg "ROArray.fold_left2: arrays must have the same length"
  else
    let r = ref x in
    for i = 0 to la - 1 do
      r := f !r (unsafe_get a i) (unsafe_get b i)
    done;
    !r

let fold_left3 f x a b c =
  let la = length a in
  if la <> length b || la <> length c then
    invalid_arg "ROArray.fold_left2: arrays must have the same length"
  else
    let r = ref x in
    for i = 0 to la - 1 do
      r := f !r (unsafe_get a i) (unsafe_get b i) (unsafe_get c i)
    done;
    !r

let for_all2 p a b =
  let n = length a in
  if n <> length b then
    invalid_arg "ROArray.for_all2: arrays must have the same length";
  let rec loop i =
    if i = n then true
    else if p (unsafe_get a i) (unsafe_get b i) then loop (succ i)
    else false in
  loop 0

