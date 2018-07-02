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

#if defined( SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
#include <time.h>
#include <unistd.h>
#endif

#include <caml/memory.h>
#include <caml/alloc.h>

/* ======================================================================
 * Private functions
 * ====================================================================*/

#if defined( SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
time_t base_time_tv_sec = 0; /* Base time; makes time values returned
				by get_time easier to read when
				printed since they will be zero
				based.
			     */
#endif

CAMLprim value SetTiming(value onoff)
{
    CAMLparam1 (onoff);

#if defined( SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
    struct timespec spec;
    clock_gettime( CLOCK_MONOTONIC_RAW, &spec );
    base_time_tv_sec = spec.tv_sec;
#endif
    CAMLreturn (Val_unit);
}

/* ----------------------------------------------------------------------
 * Timer
 * --------------------------------------------------------------------*/
CAMLprim value get_time ()
{
    CAMLparam0 ();
#if defined( SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
    struct timespec spec;  
    clock_gettime( CLOCK_MONOTONIC_RAW, &spec );
    double time = (double)(spec.tv_sec - base_time_tv_sec) + ((double)(spec.tv_nsec) / 1E9);
#else
    double time = 0;
#endif
    CAMLreturn (caml_copy_double (time));
}

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

