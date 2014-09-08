/***********************************************************************
 *                                                                     *
 *                   OCaml interface to Sundials                       *
 *                                                                     *
 *  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *
 *                                                                     *
 *  Copyright 2014 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a BSD 2-Clause License, refer to the file LICENSE.           *
 *                                                                     *
 ***********************************************************************/

#ifndef _DLS_ML_H__
#define _DLS_ML_H__

#include <caml/mlvalues.h>

#define DLSMAT(v) (*(DlsMat *)Data_custom_val(v))

/* This enum must list exceptions in the same order as the call to
 * c_register_exns in dls.ml.  */
enum dls_exn_index {
    DLS_EXN_ZeroDiagonalElement = 0,
    DLS_EXN_SET_SIZE
};

#define DLS_EXN(name) (Field (Field (sundials_ml_exn_table,	\
				     DLS_EXN_SET),		\
			      DLS_EXN_ ## name))

#endif

