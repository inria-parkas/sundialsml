(* expr_of functions for data types in sundials.mli.  *)
open Quickcheck
open Sundials
open Camlp4.PreCast

let _loc = Loc.ghost

let expr_of_carray v =
  let n = Carray.length v in
  if n = 0 then <:expr<Carray.create 0>>
  else <:expr<Carray.of_array
              [| $Ast.exSem_of_list (List.map (fun i -> <:expr<$`flo:v.{i}$>>)
                                      (enum 0 (n-1)))$ |]>>

let expr_of_root_event = function
  | Roots.Rising -> <:expr<Roots.Rising>>
  | Roots.Falling -> <:expr<Roots.Falling>>
  | Roots.NoRoot -> <:expr<Roots.NoRoot>>

let expr_of_root_info rs =
  <:expr<Roots.of_array [| $Ast.exSem_of_list (List.map expr_of_root_event
                                                 (Roots.to_list rs))$ |]>>

let expr_of_root_direction = function
  | RootDirs.Increasing -> <:expr<Ida.RootDirs.Increasing>>
  | RootDirs.Decreasing -> <:expr<Ida.RootDirs.Decreasing>>
  | RootDirs.IncreasingOrDecreasing ->
    <:expr<Ida.RootDirs.IncreasingOrDecreasing>>
