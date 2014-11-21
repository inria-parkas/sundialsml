Sundials/ML
===========

Sundials/ML is an interface to the Sundials suite of numerical solvers.
[Full documentation](http://inria-parkas.github.io/sundialsml/) is online.

[Sundials](http://computation.llnl.gov/casc/sundials/main.html) is a collection of five numerical solvers:
CVODE, CVODES, IDA, IDAS, and KINSOL.
This interface provides access to all features of the underlying library.

Installation
------------
1. [Download Sundials](http://computation.llnl.gov/casc/sundials/download/download.php), extract, and install it:
    1. `./configure --enable-examples --enable-shared`
    2. `make`
    3. `make install`
2. Either install from OPAM: `opam install sundialsml`, or
    1. [Download Sundials/ML](https://github.com/inria-parkas/sundialsml/releases), extract, and install it:
    2. `./configure`
    3. `make`
    4. `make install` or `make install-findlib`
3. [Start coding!](http://inria-parkas.github.io/sundialsml/#running)

Contact
-------
* [Support (public OCaml list)](mailto:caml-list@inria.fr?subject=Sundials/ML:)
* [Bug reports/Feature requests](https://github.com/inria-parkas/sundialsml/issues/new)
* [Pull requests](https://github.com/inria-parkas/sundialsml/compare)
