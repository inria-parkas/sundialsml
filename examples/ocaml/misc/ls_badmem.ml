(* Compile with:
    ocamlc -o ls_badmem.byte -I +sundials -dllpath +sundials \
           bigarray.cma unix.cma sundials.cma ls_badmem.ml
 *)

let printf = Printf.printf

let f _ _ yd = yd.{0} <- 1.0

(* simulation *)

let y = Sundials.RealArray.of_array [| 0.0 |]
let y_nv = Nvector_serial.wrap y

let s1 = Cvode.(init Adams
                    default_tolerances
                    ~lsolver:Spils.(solver (spgmr y_nv) prec_none)
                    f 0.0 y_nv);;

let i = Cvode.Spils.get_num_lin_iters s1;;
printf "Spils.get_num_lin_iters s1 = %d\n" i;;

let s2 = Cvode.init Cvode.Adams
                    ~lsolver:Cvode.Diag.solver
                    Cvode.default_tolerances
                    f 0.0 y_nv;;

let i = Cvode.Diag.get_num_rhs_evals s2;;
printf "Diag.get_num_rhs_evals s2 = %d\n" i;;

