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

module Config = Sundials_Config

module Index = Sundials_Index

module Logfile = struct
  type t = Sundials_impl.Logfile.t

  let stderr   = Sundials_impl.Logfile.stderr
  let stdout   = Sundials_impl.Logfile.stdout
  let openfile = Sundials_impl.Logfile.openfile

  let output_string = Sundials_impl.Logfile.output_string
  let output_bytes  = Sundials_impl.Logfile.output_bytes

  let flush = Sundials_impl.Logfile.flush
  let close = Sundials_impl.Logfile.close
end

module Profiler = struct
  type t = Sundials_impl.Profiler.t

  let enabled = Sundials_configuration.profiling_enabled

  external make : string -> t
    = "sunml_profiler_make"

  external start : t -> string -> unit
    = "sunml_profiler_begin" [@@noalloc]

  external finish : t -> string -> unit
    = "sunml_profiler_end" [@@noalloc]

  external print : t -> Logfile.t -> unit
    = "sunml_profiler_print"
end

module Context = struct
  type t = Sundials_impl.Context.t
  exception ExternalProfilerInUse = Sundials_impl.Context.ExternalProfilerInUse
  let make    = Sundials_impl.Context.make
  let default = Sundials_impl.Context.default

  let get_profiler = Sundials_impl.Context.get_profiler
  let set_profiler = Sundials_impl.Context.set_profiler
end

exception RecoverableFailure
exception NonPositiveEwt

type 'f cfun = 'f Sundials_impl.Callback.cfun

let invoke = Sundials_impl.Callback.invoke

(* Note the type annotations are redundant because there's already a .mli, but
   explicit annotations improve performance for bigarrays.  *)
module RealArray = Sundials_RealArray

module RealArray2 = Sundials_RealArray2

(* Opaque arrays *)

module type ArrayBaseOps = sig (* {{{ *)
  type t
  type elt
  val error_name : string
  val create : int -> t
  val get : t -> int -> elt
  val set : t -> int -> elt -> unit
  val length : t -> int
end (* }}} *)

module ArrayLike (A : ArrayBaseOps) = struct (* {{{ *)
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

  let fill a x =
    for i = 0 to length a - 1 do
      set a i x
    done

  let blitn ~src ?(spos=0) ~dst ?(dpos=0) len =
    let na = length src
    and nb = length dst in
    if Sundials_configuration.safe &&
       (len < 0 || spos < 0 || na < spos+len || dpos < 0 || nb < dpos+len)
    then invalid_arg (Printf.sprintf "%s.blitn" error_name)
    else
      for i = 0 to len-1 do
        set dst (dpos + i) (get src (spos + i))
      done

  let blit ~src ~dst =
    for i = 0 to min (length src) (length dst) - 1 do
      set dst i (get src i)
    done

  let init n f =
    let a = create n in
    for i = 0 to n-1 do
      set a i (f i)
    done;
    a
end (* }}} *)

(* root arrays *)

module LintArray = Sundials_LintArray

module ROArray = Sundials_ROArray

module Roots = struct (* {{{ *)
  open Bigarray
  type t = (int32, int32_elt, c_layout) Array1.t

  type r =
    | NoRoot
    | Rising
    | Falling

  let root_of_int32 = function
    | 1l -> Rising
    | -1l -> Falling
    | 0l -> NoRoot
    | n ->
      failwith
        (Printf.sprintf
           "Roots.root_of_int32: invalid root event %ld" n)

  let int32_of_root x =
    match x with
    | NoRoot -> 0l
    | Rising -> 1l
    | Falling -> -1l

  let int_of_root x =
    match x with
    | NoRoot -> 0
    | Rising -> 1
    | Falling -> -1

  let reset v = Array1.fill v 0l

  let detected roots i = roots.{i} <> 0l
  let get roots i = root_of_int32 (roots.{i})
  let set a i v = a.{i} <- int32_of_root v

  let create n =
    let a = Array1.create int32 c_layout n in
    reset a;
    a

  let length = Array1.dim

  let string_of_root e =
    match e with
    | NoRoot  -> "_"
    | Rising  -> "R"
    | Falling -> "F"

  let pp fmt a =
    Format.pp_print_string fmt "[";
    Format.pp_open_hovbox fmt 0;
    for i = 0 to length a - 1 do
      if i > 0 then (
        Format.pp_print_string fmt ";";
        Format.pp_print_space fmt ();
      );
      Format.pp_print_string fmt (string_of_root (root_of_int32 a.{i}))
    done;
    Format.pp_close_box fmt ();
    Format.pp_print_string fmt "]"

  let ppi ?(start="[") ?(stop="]") ?(sep=";")
          ?(item=fun fmt i e ->
                    Format.pp_print_int fmt i;
                    Format.pp_print_string fmt "=";
                    Format.pp_print_string fmt (string_of_root e)) ()
          fmt a =
    Format.pp_print_string fmt start;
    Format.pp_open_hovbox fmt 0;
    for i = 0 to length a - 1 do
      if i > 0 then (
        Format.pp_print_string fmt sep;
        Format.pp_print_space fmt ();
      );
      item fmt i (root_of_int32 a.{i})
    done;
    Format.pp_close_box fmt ();
    Format.pp_print_string fmt stop

  module A = ArrayLike (struct (* {{{ *)
    type t = (int32, int32_elt, c_layout) Array1.t
    and elt = r
    let error_name = "Roots"
    let get = get
    let set = set
    let create = create
    let length = length
  end (* }}} *))

  let make n x = A.make n x
  let copy = A.copy
  let init = A.init

  let of_array = A.of_array
  let of_list = A.of_list
  let to_array = A.to_array
  let to_list = A.to_list

  let fill = A.fill

  let rising  roots i = roots.{i} = 1l
  let falling roots i = roots.{i} = -1l

  let set_noroot a i = set a i NoRoot
  let set_rising a i = set a i Rising
  let set_falling a i = set a i Falling

  let iteri = A.iteri
  let iter = A.iter

  let exists a =
    let n = length a in
    let rec go i = a.{i} <> 0l || (i < n-1 && go (i+1)) in
    go 0
end (* }}} *)

module RootDirs = struct (* {{{ *)
  open Bigarray
  type t = (int32, int32_elt, c_layout) Array1.t

  type d =
    | Increasing
    | Decreasing
    | IncreasingOrDecreasing

  let int32_of_rootdir x =
    match x with
    | Increasing -> 1l
    | Decreasing -> -1l
    | IncreasingOrDecreasing -> 0l

  let rootdir_of_int32 x =
    match x with
    | 1l -> Increasing
    | -1l -> Decreasing
    | 0l -> IncreasingOrDecreasing
    | n ->
      failwith
        (Printf.sprintf
           "Roots.rootdir_of_int32: invalid root direction %ld" n)

  let make n x =
    let a = Array1.create int32 c_layout n in
    Array1.fill a (int32_of_rootdir x);
    a

  let create n = make n IncreasingOrDecreasing

  let length a = Array1.dim a

  let string_of_rootdir e =
    match e with
    | Increasing -> "R"
    | Decreasing -> "F"
    | IncreasingOrDecreasing -> "E"

  let pp fmt a =
    Format.pp_print_string fmt "[";
    Format.pp_open_hovbox fmt 0;
    for i = 0 to length a - 1 do
      if i > 0 then (
        Format.pp_print_string fmt ";";
        Format.pp_print_space fmt ();
      );
      Format.pp_print_string fmt (string_of_rootdir (rootdir_of_int32 a.{i}))
    done;
    Format.pp_close_box fmt ();
    Format.pp_print_string fmt "]"

  let ppi ?(start="[") ?(stop="]") ?(sep=";")
          ?(item=fun fmt i e ->
                    Format.pp_print_int fmt i;
                    Format.pp_print_string fmt "=";
                    Format.pp_print_string fmt (string_of_rootdir e)) ()
          fmt a =
    Format.pp_print_string fmt start;
    Format.pp_open_hovbox fmt 0;
    for i = 0 to length a - 1 do
      if i > 0 then (
        Format.pp_print_string fmt sep;
        Format.pp_print_space fmt ();
      );
      item fmt i (rootdir_of_int32 a.{i})
    done;
    Format.pp_close_box fmt ();
    Format.pp_print_string fmt stop

  let copy n src =
    let nsrc = Array.length src in
    let a = Array1.create int32 c_layout n in
    if n > nsrc
    then Array1.fill a (int32_of_rootdir IncreasingOrDecreasing);
    for i = 0 to min n nsrc - 1 do
      a.{i} <- int32_of_rootdir src.(i)
    done;
    a

  let set a i v = a.{i} <- int32_of_rootdir v
  let get a i = rootdir_of_int32 a.{i}

  module A = ArrayLike (
    struct (* {{{ *)
      type t = (int32, int32_elt, c_layout) Array1.t
      and elt = d
      let error_name = "RootDirs"
      let create = create
      let set = set
      let get = get
      let length = length
    end (* }}} *))

  let of_array  = A.of_array
  let of_list   = A.of_list
  let to_array  = A.to_array
  let to_list   = A.to_list
  let fill      = A.fill
  let blitn     = A.blitn
  let blit      = A.blit

  let init = A.init

end (* }}} *)

module Constraint = struct (* {{{ *)
  let unconstrained = 0.0
  let geq_zero =  1.0
  let leq_zero = -1.0
  let gt_zero  =  2.0
  let lt_zero  = -2.0

  type t =
  | Unconstrained
  | GeqZero
  | LeqZero
  | GtZero
  | LtZero

  let of_float = function
    |  0.0 -> Unconstrained
    |  1.0 -> GeqZero
    | -1.0 -> LeqZero
    |  2.0 -> GtZero
    | -2.0 -> LtZero
    |    f -> raise (Invalid_argument
                    ("invalid constraint: " ^ string_of_float f))
  let to_float = function
    | Unconstrained -> 0.0
    | GeqZero       -> 1.0
    | LeqZero       -> -1.0
    | GtZero        -> 2.0
    | LtZero        -> -2.0
end (* }}} *)

module Matrix = Sundials_Matrix

module LinearSolver = Sundials_LinearSolver

module NonlinearSolver = Sundials_NonlinearSolver

module Util = struct (* {{{ *)

  type error_details = {
      error_code : int;
      module_name : string;
      function_name : string;
      error_message : string;
    }

  external format_float : string -> float -> string
      = "caml_format_float"

  let floata = format_float "%a"

  external c_compare_float : float -> float -> float -> int
      = "sunml_sundials_compare_tol"

  let compare_float ?(tol=10.0 *. Config.unit_roundoff) a b =
    c_compare_float a b tol = 0

end (* }}} *)

external c_init_module :
  (exn -> string -> unit)
  -> ('a Weak.t -> int -> 'a option)
  -> exn array
  -> unit
  = "sunml_sundials_init_module"

(* Used when user code throws an exception in a context where the
   underlying C library does not allow for clean handling.  *)
let warn_discarded_exn exn context =
  Printf.fprintf stderr
    ("WARNING: discarding exception %s raised in %s.\n%s")
    (Printexc.to_string exn) context
    (Printexc.get_backtrace ());
  flush stderr

let () =
  c_init_module warn_discarded_exn Weak.get
    (* Exceptions must be listed in the same order as
       sundials_exn_index.  *)
    [|RecoverableFailure;
      NonPositiveEwt;
      Config.NotImplementedBySundialsVersion|]

