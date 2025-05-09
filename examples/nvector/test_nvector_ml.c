/*
 * -----------------------------------------------------------------
 * $Revision: 4463 $
 * $Date: 2015-03-29 16:28:20 -0700 (Sun, 29 Mar 2015) $
 * ----------------------------------------------------------------- 
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------
 * Acknowledgements: These testing routines are based on an
 *                   NVECTOR testing routine by Daniel R. Reynolds
 *                   @ SMU.
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
 * These test functions are designed to check an NVECTOR module 
 * implementation. 
 *
 * NOTE: Many of these tests rely on the N_VGetArrayPointer routine 
 *       to get a pointer to the data component of an N_Vector. This 
 *       assumes the internal data is stored in a contiguous 
 *       realtype array.
 * -----------------------------------------------------------------
 */

#include <sundials/sundials_config.h>
#include <unistd.h>

#if defined( SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
#include <time.h>
#include <unistd.h>
#endif

#include <caml/memory.h>
#include <caml/alloc.h>

/* ======================================================================
 * Private functions
 * ====================================================================*/

/* Misc. */
CAMLprim value stdc_version ()
{
    CAMLparam0 ();

    /* C89 compilers don't define this, in which case it's replaced by
     * 0 in #if directives, which is what we're emulating here.  */
#ifndef __STDC_VERSION__
#define __STDC_VERSION__ 0
#endif
    CAMLreturn (Val_int (__STDC_VERSION__));
}
