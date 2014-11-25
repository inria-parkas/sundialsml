#use "topfind";;
#require "sundialsml.mpi";;

let comm = Mpi.comm_world
let n = Mpi.comm_size comm
let my_id = Mpi.comm_rank comm
let pv = Nvector_parallel.make 1 n comm (float_of_int (my_id + 1));;

Printf.printf "%d: local=%f.\n" my_id (Nvector_parallel.local_array pv).{0};;
Printf.printf "Sum of abs. = %f\n" (Nvector_parallel.Ops.n_vl1norm pv);;
