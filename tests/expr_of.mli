(** expr_of for primitive types.  This is to expr_of what pprint.ml is to
    pretty-printing.  All expressions will be tagged with the ghost location.
 *)
open Camlp4.PreCast.Ast

val expr_of_int : int -> expr
val expr_of_bool : bool -> expr
val expr_of_int32 : int32 -> expr
val expr_of_int64 : int64 -> expr
val expr_of_naiveint : nativeint -> expr
val expr_of_float : float -> expr
val expr_of_char : char -> expr
val expr_of_string : string -> expr
val expr_of_unit : unit -> expr
val expr_of_lazy_t : ('a -> expr) -> ('a lazy_t -> expr)

val expr_of_array : ('a -> expr) -> ('a array -> expr)
val expr_of_list : ('a -> expr) -> ('a list -> expr)

val expr_of_pair : ('a -> expr) -> ('b -> expr) -> ('a * 'b -> expr)
val expr_of_triple : ('a -> expr) -> ('b -> expr) -> ('c -> expr)
                   -> ('a * 'b * 'c -> expr)

val expr_of_option : ('a -> expr) -> ('a option -> expr)

(** Reifies an exception object.  Since exception objects are open-ended, this
    function can only handle exceptions that it knows about.  By default,
    [expr_of_exn exn] works only if [exn] is an exception defined in
    [Pervasives].  {!register_expr_of_exn} can be used to extend this function
    with new handlers.

    A handler [f : (exn -> expr) -> exn -> expr] is invoked as [f k exn] with a
    continuation [k] and should return either an [expr] if successful or
    otherwise call [k exn].  For example:

    {[
    exception Foo
    exception Bar
    let foo_bar_handler k = function
      | Foo -> <:expr<Foo>>
      | Bar -> <:expr<Bar>>
      | exn -> k exn
    ]}

    A handler is registered by calling [register_expr_of_exn f].  When
    [expr_of_exn exn] is evaluated, it calls the most recently registered
    handler with a continuation that invokes the second most recent handler.
    The second most recent handler receives a continuation that invokes the
    third most recent handler, and so on until all handlers have been tried.
    The oldest handler receives the default implementation as the continuation,
    which can reify any exception defined in [Pervasives].  If that also fails,
    then [expr_of_exn] raises [Failure].

    You can remove a handler by [unregister_expr_of_exn f] -- in this case [f]
    must be physically equal (i.e. pointer-equal) to a previously registered
    handler.  A handler registered more than once is removed completely.  *)
val expr_of_exn : exn -> expr

(** See {!expr_of_exn}.  *)
val register_expr_of_exn : ((exn -> expr) -> exn -> expr) -> unit

(** See {!expr_of_exn}.  *)
val unregister_expr_of_exn : ((exn -> expr) -> exn -> expr) -> unit
