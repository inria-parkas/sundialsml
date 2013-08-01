(* Sundials-specific pretty-printers.  *)
open Pprint
open Sundials

let pp_carray, dump_carray, show_carray, display_carray,
  print_carray, prerr_carray =
  printers_of_pp (fun fmt xs ->
    if !read_write_invariance
    then pp_array_like Carray.length Bigarray.Array1.get
           "(Carray.of_array [|" "|])" pp_float fmt xs
    else pp_array_like Carray.length Bigarray.Array1.get
           "[<" ">]" pp_float fmt xs)

let show_root_event x =
  (if !read_write_invariance then "Roots." else "")
  ^ Roots.string_of_root_event x

let pp_root_event, dump_root_event, show_root_event, display_root_event,
  print_root_event, prerr_root_event
    =
  printers_of_show show_root_event

let pp_root_info fmt xs =
  let get a i =
    (* If Roots is used improperly or if there's a bug in the binding, a root
       info array can contain garbage that doesn't correspond to any of NoRoot,
       Rising, or Falling.  Such values trigger a Failure in Roots.get.  *)
    try show_root_event (Roots.get' a i)
    with Failure _ -> "<garbage>"
  in
  if !read_write_invariance
  then pp_array_like Roots.length get "(Roots.of_array [|" "|])"
         pp_unquoted_string fmt xs
  else pp_array_like Roots.length get "[<" ">]" pp_unquoted_string fmt xs

let _, dump_root_info, show_root_info, display_root_info,
  print_root_info, prerr_root_info
    =
  printers_of_pp pp_root_info

let show_root_direction x =
  (if !read_write_invariance then "Roots." else "")
  ^ RootDirs.string_of_root_direction x

let pp_root_direction, dump_root_direction, show_root_direction,
  display_root_direction, print_root_direction, prerr_root_direction
    =
  printers_of_show show_root_direction

let pp_root_dirs fmt xs =
  let get a i =
    (* If Roots is used improperly or if there's a bug in the binding, a root
       info array can contain garbage that doesn't correspond to any of NoRoot,
       Rising, or Falling.  Such values trigger a Failure in Roots.get.  *)
    try RootDirs.string_of_root_direction (RootDirs.get a i)
    with Failure _ -> "<garbage>"
  in
  if !read_write_invariance
  then pp_array_like RootDirs.length get "(RootDirs.of_array [|" "|])"
         pp_unquoted_string fmt xs
  else pp_array_like RootDirs.length get "[<" ">]" pp_unquoted_string fmt xs

let _, dump_root_dirs, show_root_dirs, display_root_dirs,
  print_root_dirs, prerr_root_dirs
    =
  printers_of_pp pp_root_dirs
