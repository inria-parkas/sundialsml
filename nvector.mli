type 'a nvector

module Mutable :
  sig
    type 'a nvector_ops = {
      nvclone           : 'a -> 'a;
      nvdestroy         : ('a -> unit) option;
      nvspace           : ('a -> int * int) option;
      nvlinearsum       : float -> 'a -> float -> 'a -> 'a -> unit;
      nvconst           : float -> 'a -> unit;
      nvprod            : 'a -> 'a -> 'a -> unit;
      nvdiv             : 'a -> 'a -> 'a -> unit;
      nvscale           : float -> 'a -> 'a -> unit;
      nvabs             : 'a -> 'a -> unit;
      nvinv             : 'a -> 'a -> unit;
      nvaddconst        : 'a -> float -> 'a -> unit;
      nvmaxnorm         : 'a -> float;
      nvwrmsnorm        : 'a -> 'a -> float;
      nvmin             : 'a -> float;
      nvdotprod         : ('a -> 'a -> float) option;
      nvcompare         : (float -> 'a -> 'a -> unit) option;
      nvinvtest         : ('a -> 'a -> bool) option;
      nvwl2norm         : ('a -> 'a -> float) option;
      nvl1norm          : ('a -> float) option;
      nvwrmsnormmask    : ('a -> 'a -> 'a -> float) option;
      nvconstrmask      : ('a -> 'a -> 'a -> bool) option;
      nvminquotient     : ('a -> 'a -> float) option;
    }

    val make_nvector    : 'a nvector_ops -> 'a -> 'a nvector
    val nvector_data    : 'a nvector -> 'a
  end

module Immutable :
  sig
    type 'a nvector_ops = {
      nvclone           : 'a -> 'a;
      nvdestroy         : ('a -> unit) option;
      nvspace           : ('a -> int * int) option;
      nvlinearsum       : float -> 'a -> float -> 'a -> 'a;
      nvconst           : float -> 'a;
      nvprod            : 'a -> 'a -> 'a;
      nvdiv             : 'a -> 'a -> 'a;
      nvscale           : float -> 'a -> 'a;
      nvabs             : 'a -> 'a;
      nvinv             : 'a -> 'a;
      nvaddconst        : 'a -> float -> 'a;
      nvmaxnorm         : 'a -> float;
      nvwrmsnorm        : 'a -> 'a -> float;
      nvmin             : 'a -> float;
      nvdotprod         : ('a -> 'a -> float) option;
      nvcompare         : (float -> 'a -> 'a) option;
      nvinvtest         : ('a -> 'a -> bool) option;
      nvwl2norm         : ('a -> 'a -> float) option;
      nvl1norm          : ('a -> float) option;
      nvwrmsnormmask    : ('a -> 'a -> 'a -> float) option;
      nvconstrmask      : ('a -> 'a -> 'a -> bool) option;
      nvminquotient     : ('a -> 'a -> float) option;
    }

    val from_immutable  : 'a nvector_ops -> 'a ref Mutable.nvector_ops
    val make_nvector    : 'a nvector_ops -> 'a ref -> 'a ref nvector
    val nvector_data    : 'a ref nvector -> 'a
  end

