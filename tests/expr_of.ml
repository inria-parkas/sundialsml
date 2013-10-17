open Camlp4.PreCast.Syntax
open Ast
module Camlp4aux = Camlp4aux.Make (Camlp4.PreCast.Syntax)
open Camlp4aux
module Meta = Meta.Make (Meta.MetaLoc)

let _loc = Loc.ghost

let expr_of_int x = <:expr<$`int:x$>>
let expr_of_bool x = <:expr<$`bool:x$>>
let expr_of_int32 x = <:expr<$`int32:x$>>
let expr_of_int64 x = <:expr<$`int64:x$>>
let expr_of_naiveint x = <:expr<$`nativeint:x$>>
let expr_of_float x = <:expr<$`flo:x$>>
let expr_of_char x = <:expr<$`chr:x$>>
let expr_of_string x = <:expr<$`str:x$>>
let expr_of_unit () = <:expr<()>>
let expr_of_lazy_t e0 x = <:expr<lazy $e0 (Lazy.force x)$>>

let expr_of_array expr_of_elem x =
  let es = Array.fold_left (fun e x -> sem e (expr_of_elem x))
             (ExNil Loc.ghost) x
  in <:expr<[|$es$|]>>

let expr_of_list expr_of_elem =
  Meta.Expr.meta_list (fun _ -> expr_of_elem) _loc

let expr_of_pair expr_of_fst expr_of_snd (x,y) =
  <:expr<$expr_of_fst x$, $expr_of_snd y$>>

let expr_of_triple expr_of_fst expr_of_snd expr_of_thd (x,y,z) =
  <:expr<$expr_of_fst x$, $expr_of_snd y$, $expr_of_thd z$>>

let expr_of_option expr_of_contents = function
  | None -> <:expr<None>>
  | Some x -> <:expr<Some $expr_of_contents x$>>

let expr_of_exn_registry = ref []

let register_expr_of_exn f =
  expr_of_exn_registry := f :: !expr_of_exn_registry

let unregister_expr_of_exn f =
  expr_of_exn_registry := List.filter ((!=) f) !expr_of_exn_registry

let expr_of_exn exn =
  let exn_default = function
    | Not_found -> <:expr<Not_found>>
    | Out_of_memory -> <:expr<Out_of_memory>>
    | Stack_overflow -> <:expr<Stack_overflow>>
    | End_of_file -> <:expr<End_of_file>>
    | Division_by_zero -> <:expr<Division_by_zero>>
    | Sys_blocked_io -> <:expr<Sys_blocked_io>>
    | Match_failure (a,b,c) -> <:expr<Match_failure ($expr_of_string a$,
                                                     $expr_of_int b$,
                                                     $expr_of_int c$)>>
    | Assert_failure (a,b,c) -> <:expr<Assert_failure ($expr_of_string a$,
                                                       $expr_of_int b$,
                                                       $expr_of_int c$)>>
    | Invalid_argument s -> <:expr<Invalid_argument $`str:s$>>
    | Failure s -> <:expr<Failure $`str:s$>>
    | Sys_error s -> <:expr<Sys_Error $`str:s$>>
    | exn -> failwith ("expr_of_exn: no reifier found for "
                       ^ Printexc.to_string exn)
  in
  let rec go = function
    | [] -> exn_default
    | f::fs -> fun exn -> f (go fs) exn
  in go !expr_of_exn_registry exn
