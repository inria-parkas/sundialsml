
type kind
type t = (Sundials.RealArray.t * int * Mpi.communicator, kind) Sundials.nvector

exception IncorrectGlobalSize

let _ = Callback.register_exception
          "nvector_parallel_IncorrectGlobalSize" IncorrectGlobalSize

external wrap : (Sundials.RealArray.t * int * Mpi.communicator) -> t
  = "ml_nvec_wrap_parallel"

let make nl ng comm iv = wrap (Sundials.RealArray.init nl iv, ng, comm)

let unwrap nv =
  let data, _, _ = Sundials.unvec nv in
  data

let global_length nv =
  let _, gl, _ = Sundials.unvec nv in
  gl

let communicator nv =
  let _, _, comm = Sundials.unvec nv in
  comm

