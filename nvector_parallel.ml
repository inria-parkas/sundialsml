
type kind
type t = (Sundials.RealArray.t * int * MPI.Communicator, kind) nvector

exception IncorrectGlobalSize

let _ = Callback.register_exception
          "nvector_parallel_IncorrectGlobalSize" IncorrectGlobalSize

external wrap : (Sundials.RealArray.t * int * MPI.Communicator) -> t
  = "ml_nvec_wrap_parallel"

let make nl ng comm iv = wrap (Sundials.RealArray.init nl iv, ng, comm)

let unwrap nv =
  let data, _, _ = Sundials.unwrap nv in
  data

let global_length nv =
  let _, gl, _ = Sundials.unwrap nv in
  gl

let communicator nv =
  let _, _, comm = Sundials.unwrap nv in
  comm

