(***********************************************************************)
(*                                                                     *)
(*               OCaml interface to (serial) Sundials                  *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a BSD 2-Clause License, refer to the file LICENSE.           *)
(*                                                                     *)
(***********************************************************************)

include Cvode

type 'a nvector = 'a Nvector.nvector

type root_array = Sundials.Roots.t
type root_val_array = Sundials.Roots.val_array

type 'a single_tmp = 'a
type 'a triple_tmp = 'a * 'a * 'a

type ('t, 'a) jacobian_arg =
  {
    jac_t   : float;
    jac_y   : 'a;
    jac_fy  : 'a;
    jac_tmp : 't
  }

type 'a prec_solve_arg =
  {
    rhs   : 'a;
    gamma : float;
    delta : float;
    left  : bool;
  }

(* Note this definition differs from the one in cvode_serial, so the tag
   values are different.  *)
type 'a linear_solver =
  | Diag
  | Spgmr of spils_params * 'a spils_callbacks
  | Spbcg of spils_params * 'a spils_callbacks
  | Sptfqmr of spils_params * 'a spils_callbacks

and spils_params = { maxl : int option;
                     prec_type : Spils.preconditioning_type; }

and 'a spils_callbacks =
  {
    prec_solve_fn : (('a single_tmp, 'a) jacobian_arg -> 'a prec_solve_arg
                     -> 'a nvector -> unit) option;
    prec_setup_fn : (('a triple_tmp, 'a) jacobian_arg -> bool -> float -> bool)
                    option;
    jac_times_vec_fn :
      (('a single_tmp, 'a) jacobian_arg
       -> 'a (* v *)
       -> 'a (* Jv *)
       -> unit) option;
  }

type 'a iter =
  | Newton of 'a linear_solver
  | Functional

type cvode_mem
type cvode_file
type c_weak_ref

type 'a session = {
        cvode      : cvode_mem;
        backref    : c_weak_ref;
        nroots     : int;
        err_file   : cvode_file;

        mutable exn_temp   : exn option;

        mutable rhsfn      : float -> 'a -> 'a -> unit;
        mutable rootsfn    : float -> 'a -> root_val_array -> unit;
        mutable errh       : Sundials.error_details -> unit;
        mutable errw       : 'a -> 'a -> unit;
        mutable presetupfn : ('a triple_tmp, 'a) jacobian_arg -> bool -> float -> bool;
        mutable presolvefn : ('a single_tmp, 'a) jacobian_arg -> 'a prec_solve_arg
                               -> 'a nvector -> unit;
        mutable jactimesfn : ('a single_tmp, 'a) jacobian_arg -> 'a -> 'a -> unit;

        mutable sensext    : Obj.t option (* Used by CVODES *)
      }

