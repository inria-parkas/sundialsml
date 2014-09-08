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

#ifndef _SPILS_ML_H__
#define _SPILS_ML_H__

#include <caml/mlvalues.h>

enum gramschmidt_type_tag {
  VARIANT_SPILS_GRAMSCHMIDT_TYPE_MODIFIEDGS  = 0,
  VARIANT_SPILS_GRAMSCHMIDT_TYPE_CLASSICALGS,
};

enum preconditioning_type_tag {
  VARIANT_SPILS_PRECONDITIONING_TYPE_PRECNONE  = 0,
  VARIANT_SPILS_PRECONDITIONING_TYPE_PRECLEFT,
  VARIANT_SPILS_PRECONDITIONING_TYPE_PRECRIGHT,
  VARIANT_SPILS_PRECONDITIONING_TYPE_PRECBOTH,
};

int spils_precond_type(value);
int spils_gs_type(value);


/* This enum must list exceptions in the same order as the call to
 * c_register_exns in spils.ml.  */
enum spils_exn_index {
    SPILS_EXN_MemoryRequestFailure = 0,
    SPILS_EXN_ConvFailure,
    SPILS_EXN_QRfactFailure,
    SPILS_EXN_PSolveFailure,
    SPILS_EXN_ATimesFailure,
    SPILS_EXN_PSetFailure,
    SPILS_EXN_GSFailure,
    SPILS_EXN_QRSolFailure,
    SPILS_EXN_SET_SIZE
};

#define SPILS_EXN(name) (Field (Field (sundials_ml_exn_table,	\
				       SPILS_EXN_SET),		\
				SPILS_EXN_ ## name))

#endif

