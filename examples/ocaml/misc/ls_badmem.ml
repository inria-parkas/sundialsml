(* Compile with:
    ocamlc -o ls_badmem.byte -I +sundials -dllpath +sundials \
           bigarray.cma unix.cma sundials.cma ls_badmem.ml
 *)

let printf = Printf.printf

let f t y yd = yd.{0} <- 1.0

(* simulation *)

let y = Sundials.RealArray.of_array [| 0.0 |]
let y_nv = Nvector_serial.wrap y

let s1 = Cvode.(init Adams
                    (Newton
                      Spils.(solver (spgmr y_nv) prec_none))
                    default_tolerances
                    f 0.0 y_nv);;

let i = Cvode.Spils.get_num_lin_iters s1;;
printf "Spils.get_num_lin_iters s1 = %d\n" i;;

let s2 = Cvode.init Cvode.Adams
                    (Cvode.Newton Cvode.Diag.solver)
                    Cvode.default_tolerances
                    f 0.0 y_nv;;

let i = Cvode.Spils.get_num_lin_iters s2;;
printf "Spils.get_num_lin_iters s2 = %d\n" i;;

