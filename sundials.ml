(* Aug 2010, Timothy Bourke (INRIA) *)

let extra_time_precision = ref false

let print_time (s1, s2) t =
  if !extra_time_precision
  then Printf.printf "%s%.15e%s" s1 t s2
  else Printf.printf "%s%e%s" s1 t s2

external format_float : string -> float -> string
    = "caml_format_float"

let floata = format_float "%a"

external get_big_real : unit -> float
    = "cvode_ml_big_real"

let big_real = get_big_real ()

external get_unit_roundoff : unit -> float
    = "cvode_ml_unit_roundoff"

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

    let length = Bigarray.Array1.dim

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

    let print_with_time' t v =
      print_time ("", "") t;
      app (Printf.printf "\t% .8f") v;
      print_newline ()

    let print_with_time t v =
      print_time ("", "") t;
      app (Printf.printf "\t% e") v;
      print_newline ()

    let clamp thres =
      let cf v = if abs_float v <= thres then 0.0 else v
      in
      if thres = 0.0 then (fun x -> ())
      else map cf
  end

(* root arrays *)

type int_array = (int32, Bigarray.int32_elt, Bigarray.c_layout) Bigarray.Array1.t
let make_int_array = Bigarray.Array1.create Bigarray.int32 Carray.layout

module Roots =
  struct
    type t = int_array
    type val_array = Carray.t

    let reset v = Bigarray.Array1.fill v 0l

    let create n =
      let a = make_int_array n in
      reset a;
      a

    let empty = create 0

    let length = Bigarray.Array1.dim

    let get roots i = roots.{i} <> 0l
    let get' roots i = Int32.to_int roots.{i}

    let set a i v = Bigarray.Array1.set a i (if v then 1l else 0l)

    let appi f v =
      for i = 0 to (length v - 1) do
        f i (v.{i} <> 0l)
      done

    let app f v =
      for i = 0 to (length v - 1) do
        f (v.{i} <> 0l)
      done

    let print vs =
      app (fun v -> print_string (if v then "\t1" else "\t0")) vs;
      print_newline ()

    let print' vs =
      Carray.appi (fun i v -> Printf.printf "\t% ld" v) vs;
      print_newline ()

    let fold_left f a vs =
      let len = Bigarray.Array1.dim vs - 1 in
      let a = ref a in
      for i = 0 to len do
        a := f !a (Int32.to_int vs.{i})
      done;
      !a

    let exists = fold_left (fun a x -> a || x <> 0) false
  end

