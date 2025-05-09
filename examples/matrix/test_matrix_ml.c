/*
 * ----------------------------------------------------------------- 
 * Programmer(s): David J. Gardner @ LLNL
 *                Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * These test functions are designed to check a SUNMatrix module
 * implementation. 
 * -----------------------------------------------------------------
 */

#include <stdlib.h>

#if defined( SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
#include <time.h>
#include <unistd.h>
#endif

#include <caml/memory.h>
#include <caml/alloc.h>

/* ======================================================================
 * Private functions
 * ====================================================================*/

CAMLprim value ml_rand ()
{
    CAMLparam0 ();
    CAMLreturn (Val_int(rand()));
}

CAMLprim value ml_rand_max ()
{
    CAMLparam0 ();
    CAMLreturn (Val_int(RAND_MAX));
}

