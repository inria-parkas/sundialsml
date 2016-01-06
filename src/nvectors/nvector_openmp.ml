
type data = Sundials.RealArray.t
type kind = Nvector_serial.kind
type t = (data, kind) Nvector.t

external c_wrap : int -> Sundials.RealArray.t
                    -> (Sundials.RealArray.t -> bool) -> t
  = "ml_nvec_wrap_openmp"

let wrap nthreads v =
  let len = Sundials.RealArray.length v in
  c_wrap nthreads v (fun v' -> len = Sundials.RealArray.length v')

let unwrap = Nvector.unwrap

let make nthreads n iv = wrap nthreads (Sundials.RealArray.make n iv)

external num_threads : t -> int
  = "ml_nvec_openmp_num_threads"

