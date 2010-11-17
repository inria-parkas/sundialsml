
val array_nvec_ops  : float array Nvector.Mutable.nvector_ops
val make            : int -> float -> float array Nvector.nvector
val wrap            : float array -> float array Nvector.nvector
val data            : 'a Nvector.nvector -> 'a

