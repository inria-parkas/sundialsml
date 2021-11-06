/***********************************************************************
 *                                                                     *
 *                   OCaml interface to Sundials                       *
 *                                                                     *
 *             Timothy Bourke, Jun Inoue, and Marc Pouzet              *
 *             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *
 *                                                                     *
 *  Copyright 2014 Institut National de Recherche en Informatique et   *
 *  en Automatique.  All rights reserved.  This file is distributed    *
 *  under a New BSD License, refer to the file LICENSE.                *
 *                                                                     *
 ***********************************************************************/

#ifndef __NVECTOR_ML_H__
#define __NVECTOR_ML_H__

#include <sundials/sundials_nvector.h>
#include <caml/mlvalues.h>
#include "../sundials/sundials_ml.h"

/* OCaml interface to Serial and Custom NVectors.

   Implementation principles.

   An NVector is represented abstractly on the OCaml side by the type:
	type ('data, 'kind) nvector

   where 'data is the underlying "data payload" and 'kind relates to the
   underlying representation (see below).

   The details of the underlying representation depend on two factors:

     I. How the nvector was created, either:

	created: the nvector was created from OCaml by calling a function
		 like Nvector_serial.wrap, Nvector_parallel.make,
		 Nvector_custom.make, etc.
	
	cloned:  the nvector was created by the Sundials C-code by cloning an
		 existing nvector.

    II. The 'kind of the nvector, one of:
	     Serial	    Nvector_serial.kind
	     Parallel	    Nvector_parallel.kind (nvector_par_ml.h/.c)
	     Custom	    Nvector_custom.kind

   The general schema is:

	     OCaml heap		:         C heap
    ----------------------------+---------------------------
				:
		 +------------------------------------+
		 |		:		      |
		 \   +-------+  :                     |
		  -->|payload|  :                     |
		 /   +-------+  :                     |
		 |		:                     |
      +--------+ |              :                     |
      | 'data  +-+              :                     |
      +--------+                :        +---------+  |
      |nvector +------------------------>|N_Vector |  |
      +--------+                :        +---------+  |
     ("caml-nvec")		:        |backlink |--+
				:        +---------+
				:         ("c-nvec")

   The c-nvec comprises two values placed contiguously in malloc-ed memory:
   a Sundials N_Vector and a backlink. The N_Vector itself contains two
   pointers to other malloc-ed memory: for the ops field and for the content
   field (the fields are not shown in the diagram). The N_Vector is passed
   to Sundials functions. The backlink is an OCaml Value that is registered
   as a root with the OCaml GC, it points to the OCaml representation of the
   nvector contents. The fine details of the c-nvec depend on the type of
   nvector, they are discussed below.

   The caml-nvec exists for a "created" nvector. It pairs a 'data field that
   references a payload box with a pointer to the NVector part of the
   associated c-nvec. The finalize function of the caml-nvec is responsible
   for cleaning up the associated c-nvec (deregistering the backlink GC root
   and free()-ing the memory). When an OCaml function is called with a
   caml-nvec, the interface passes the corresponding c-nvec to Sundials
   functions. Since the caml-nvec is held in a function argument, the c-nvec
   will not be GC-ed during such calls. Sundials never stocks such c-nvecs
   internally, it clones them, so the c-nvec will not be accessed after the
   call returns.

   For a "cloned" nvector, there is no caml-nvec. The backlink ensures that
   the payload is not GC-ed while the cloned nvector is still in use. When
   Sundials makes a callback with an N_Vector, the interface extracts the
   backlink to pass to the appropriate OCaml callback. Cloned nvectors are
   eventually freed by a call to their destroy operation. This operation is
   responsible for unregistering the global root at backlink and free-ing the
   memory.

   NB: It is possible to create an OCaml-C-loop, and thus prohibiting garbage
       collection of values, by referencing the caml-nvec from within the
       payload. We do not expect this to be done in normal use (since the
       payload is typically an array of datavalues of some sort).

   NB: None of these nvectors support the N_VCloneEmpty operation. This
       operation is used by some interface functions (from Fortran or Matlab)
       which replace the data field with a pointer to their own memory. This is
       not supported by the OCaml library!

   In summary, nvectors created from OCaml are:

     created by: ml_nvec_wrap_* (it allocates both caml-nvec and c-nvec).

     deleted by: finalizer during GC (when caml-nvec dies), no explicit
                 destruction allowed.

    and nvectors cloned from C (Sundials) are:

      created by: callml_vclone (it allocates just the c-nvec).

      deleted by: explicit call to nvdestroy field of N_Vector_Ops,
                  GC never initiates destruction of any part of the structure.


   Serial nvectors
   ---------------
   The payload is a Bigarray of floats and, as per usual, the underlying data
   is allocated in the C heap (it will not be moved by the GC).

   The N_Vector content->data field points to the data in the C heap which
   underlies the payload Bigarray.

   The N_Vector ops are identical to those of a standard serial N_Vector,
   except for nvclone, nvcloneempty, and nvdestroy which are functions,
   implemented in this file, to create the arrangement described here.

   Custom nvectors
   ---------------
   The payload is the value being wrapped. The content field is set to point
   to a Value containing the OCaml callback table, it is also registered as a
   global root.  The user must ensure the callbacks do not hold referenes to
   any particular nvector, for otherwise that nvector is never reclaimed.

   Parallel nvectors
   -----------------
   This is almost the same as serial nvectors, except the payload is a triple
   containing a Bigarray of floats (containing the local portion of the
   vector), global length, and an MPI communicator.

*/

/* must match the fields of Nvector_custom.nvector_ops */
enum nvector_ops_tag {
  NVECTOR_OPS_NVCHECK = 0,
  NVECTOR_OPS_NVCLONE,
  NVECTOR_OPS_NVSPACE,
  NVECTOR_OPS_NVGETLENGTH,
  NVECTOR_OPS_NVPRINTFILE,
  NVECTOR_OPS_NVLINEARSUM,
  NVECTOR_OPS_NVCONST,
  NVECTOR_OPS_NVPROD,
  NVECTOR_OPS_NVDIV,
  NVECTOR_OPS_NVSCALE,
  NVECTOR_OPS_NVABS,
  NVECTOR_OPS_NVINV,
  NVECTOR_OPS_NVADDCONST,
  NVECTOR_OPS_NVMAXNORM,
  NVECTOR_OPS_NVWRMSNORM,
  NVECTOR_OPS_NVMIN,

  NVECTOR_OPS_NVDOTPROD,
  NVECTOR_OPS_NVCOMPARE,
  NVECTOR_OPS_NVINVTEST,

  NVECTOR_OPS_NVWL2NORM,
  NVECTOR_OPS_NVL1NORM,
  NVECTOR_OPS_NVWRMSNORMMASK,
  NVECTOR_OPS_NVCONSTRMASK,
  NVECTOR_OPS_NVMINQUOTIENT,

  NVECTOR_OPS_NVGETCOMMUNICATOR,

  NVECTOR_OPS_NVLINEARCOMBINATION,
  NVECTOR_OPS_NVSCALEADDMULTI,
  NVECTOR_OPS_NVDOTPRODMULTI,

  NVECTOR_OPS_NVLINEARSUMVECTORARRAY,
  NVECTOR_OPS_NVSCALEVECTORARRAY,
  NVECTOR_OPS_NVCONSTVECTORARRAY,
  NVECTOR_OPS_NVWRMSNORMVECTORARRAY,
  NVECTOR_OPS_NVWRMSNORMMASKVECTORARRAY,
  NVECTOR_OPS_NVSCALEADDMULTIVECTORARRAY,
  NVECTOR_OPS_NVLINEARCOMBINATIONVECTORARRAY,

  NVECTOR_OPS_NVDOTPROD_LOCAL,
  NVECTOR_OPS_NVMAXNORM_LOCAL,
  NVECTOR_OPS_NVMIN_LOCAL,
  NVECTOR_OPS_NVL1NORM_LOCAL,
  NVECTOR_OPS_NVINVTEST_LOCAL,
  NVECTOR_OPS_NVCONSTRMASK_LOCAL,
  NVECTOR_OPS_NVMINQUOTIENT_LOCAL,
  NVECTOR_OPS_NVWSQRSUM_LOCAL,
  NVECTOR_OPS_NVWSQRSUMMASK_LOCAL,

  NVECTOR_OPS_SIZE
};

/* must match the declaration of Nvector.nvector_id */
enum nvector_id_tag {
  VARIANT_NVECTOR_ID_TAG_SERIAL	    = 0,
  VARIANT_NVECTOR_ID_TAG_PARALLEL,
  VARIANT_NVECTOR_ID_TAG_OPENMP,
  VARIANT_NVECTOR_ID_TAG_PTHREADS,
  VARIANT_NVECTOR_ID_TAG_PARHYP,
  VARIANT_NVECTOR_ID_TAG_PETSC,
  VARIANT_NVECTOR_ID_TAG_CUDA,
  VARIANT_NVECTOR_ID_TAG_RAJA,
  VARIANT_NVECTOR_ID_TAG_OPENMPDEV,
  VARIANT_NVECTOR_ID_TAG_TRILINOS,
  VARIANT_NVECTOR_ID_TAG_MANYVECTOR,
  VARIANT_NVECTOR_ID_TAG_MPIMANYVECTOR,
  VARIANT_NVECTOR_ID_TAG_MPIPLUSX,
  VARIANT_NVECTOR_ID_TAG_CUSTOM,
};

struct cnvec {
    struct _generic_N_Vector nvec;
    value backlink;
};

// Return the OCaml version of the nvector payload
#define NVEC_BACKLINK(nvec) (((struct cnvec *)nvec)->backlink)

enum nv_index {
  NVEC_PAYLOAD = 0,
  NVEC_CPTR,
  NVEC_CHECK,
  NVEC_CLONE,
  NVEC_SIZE, /* This has to come last. */
};

// NVEC_VAL turns an caml-nvec into a c-nvec
#define NVEC_CVAL(v) (*(N_Vector *)Data_custom_val(v))
#define NVEC_VAL(v) (NVEC_CVAL(Field(v, NVEC_CPTR)))

#define NVEC_TAG 0
#define NVEC_ALLOC() (caml_alloc(NVEC_SIZE, NVEC_TAG))

// Internal functions
N_Vector sunml_alloc_cnvec(size_t content_size, value backlink);
void sunml_clone_cnvec_ops(N_Vector dst, N_Vector src);
CAMLprim value sunml_alloc_caml_nvec(N_Vector nv, void (*finalizer)(value));
void sunml_free_cnvec(N_Vector nv);
CAMLprim void sunml_finalize_caml_nvec(value vnv);

void sunml_nvectors_into_array(int n, value vy, N_Vector *y);
value sunml_wrap_to_nvector_table(int n, N_Vector *y);
value sunml_wrap_to_nvector_tables(int n1, int n2, N_Vector **yy);

#if 400 <= SUNDIALS_LIB_VERSION
int sunml_arrays_of_nvectors(N_Vector *r[], int n, ...);
void sunml_arrays_of_nvectors2(int* nrows, int *ncols, N_Vector **vv[], int n, ...);
#endif

N_Vector *sunml_nvector_array_alloc(value vtable);
void sunml_nvector_array_free(N_Vector *nvarr);

// Creation functions
value ml_nvec_wrap_serial(value payload, value checkfn);
value ml_nvec_wrap_custom(value mlops, value payload, value checkfn);

#endif
