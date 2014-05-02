/***********************************************************************
 *                                                                     *
 *               OCaml interface to (serial) Sundials                  *
 *                                                                     *
 *  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *
 *                                                                     *
 *  Copyright 2014 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a BSD 2-Clause License, refer to the file LICENSE.           *
 *                                                                     *
 ***********************************************************************/

/*
 * This module defines all constants and functions which do not depend on
 * the representation of continuous state vectors, i.e., those that are
 * shared between the bigarray and nvector versions of ida_ml_nvec.
 *
 */

#ifndef _KINSOL_ML_H__
#define _KINSOL_ML_H__

#include "sundials_ml.h"

/* Indices into the Kinsol_*.session type.  This enum must be in the same order as
 * the session type's member declaration.  The session data structure is shared
 * in four parts across the OCaml and C heaps.  See cvode_ml.h for details.
 */
enum kinsol_index {
    RECORD_KINSOL_SESSION_MEM = 0,
    RECORD_KINSOL_SESSION_BACKREF,
    // TODO
    RECORD_KINSOL_SESSION_SIZE	/* This has to come last. */
}

#endif
