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

#ifndef __NVECTOR_ML_H__
#define __NVECTOR_ML_H__

#include <sundials/sundials_nvector.h>
#include <caml/mlvalues.h>

/* OCaml interface to Serial, Parallel, and Custom NVectors.

   The Parallel routines are only included when compiled with
	-D CVODE_ML_WITHMPI.

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
	     Parallel	    Nvector_parallel.kind
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
				:        |  ops    |
				:        +---------+
				:        |content  |
				:        +---------+
				:         ("c-nvec")

  The c-nvec comprises four values placed contiguously in malloc-ed memory:
  a Sundials N_Vector, a backlink, memory for the nvector ops field, memory
  for the nvector content field. The N_Vector is passed to Sundials
  functions (the ops and content fields are used from within the N_Vector).
  The backlink is an OCaml Value that is registered as a root with the OCaml
  GC, it points to the OCaml representation of the nvector contents. The
  find details of the c-nvec depend on the type of nvector, they are
  discussed below.

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

  In summary, for nvectors created from OCaml:

    on creation: ml_nvec_new is invoked (it allocates memory in the C heap).

    on deletion: finalize_nvec is invoked to free memory (it, in turn,
		 invokes callml_vdestroy).

   and for nvectors cloned from C (Sundials):

     on creation: callml_vclone is invoked (it allocates memory in the C
		  heap).

     on deletion: callml_vdestroy is invoked to free memory.

  Serial nvectors
  ---------------
  The payload is a Bigarray of floats and, as per usual, the underlying data
  is allocated in the C heap (it will not be moved by the GC).

  The N_Vector content->data field points into the data underlying the
  payload Bigarray.

  The N_Vector ops are identical to those of a standard serial N_Vector,
  except for nvclone, nvcloneempty, and nvdestroy which are functions,
  implemented in this file, to create the arrangement described here.

  Parallel nvectors
  -----------------
  The payload is a triple of Bigarray, int, and MPI communicator.

  The N_Vector content->data field points into the data underlying the
  payload Bigarray. The content->local_length field is set to the length of
  the Bigarray, the content->global_length field is set to the value of the
  int, and the C value corresponding to the MPI communicator is duplicated
  into content->comm.

  The N_Vector ops are identical to those of a standard serial N_Vector,
  except for nvclone, nvcloneempty, and nvdestroy which are functions,
  implemented in this file, to create the arrangement described here.

  Custom nvectors
  ---------------
  The payload is the value being wrapped. The content field is set to point
  to a Value containing the OCaml callback table, it is also registered as a
  global root.

*/

enum nvector_ops_tag {
  NVECTOR_OPS_NVCLONE = 0,
  NVECTOR_OPS_NVDESTROY,
  NVECTOR_OPS_NVSPACE,
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
  NVECTOR_OPS_SIZE
};

// Return the OCaml version of the nvector payload
#define NVEC_BACKLINK(nvec) (*(value *)(nvec + sizeof(struct _generic_N_Vector)))

// NVEC_VAL turns an caml-nvec into a c-nvec
#define NVEC_CVAL(v) (*(N_Vector *)Data_custom_val(v))
#define NVEC_VAL(v) (NVEC_CVAL(Field(v, 1)))

// Creation functions
value ml_nvec_wrap_serial(value payload);
value ml_nvec_wrap_parallel(value payload);
value ml_nvec_wrap_custom(value mlops, value payload);

// Custom operations
N_Vector callml_vclone(N_Vector w);
void callml_vdestroy(N_Vector v);
void callml_vspace(N_Vector v, long int *lrw, long int *liw);
void callml_vlinearsum(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);
void callml_vconst(realtype c, N_Vector z);
void callml_vprod(N_Vector x, N_Vector y, N_Vector z);
void callml_vdiv(N_Vector x, N_Vector y, N_Vector z);
void callml_vscale(realtype c, N_Vector x, N_Vector z);
void callml_vabs(N_Vector x, N_Vector z);
void callml_vinv(N_Vector x, N_Vector z);
void callml_vaddconst(N_Vector x, realtype b, N_Vector z);
realtype callml_vdotprod(N_Vector x, N_Vector y);
realtype callml_vmaxnorm(N_Vector x);
realtype callml_vwrmsnorm(N_Vector x, N_Vector w);
realtype callml_vwrmsnormmask(N_Vector x, N_Vector w, N_Vector id);
realtype callml_vmin(N_Vector x);
realtype callml_vwl2norm(N_Vector x, N_Vector w);
realtype callml_vl1norm(N_Vector x);
void callml_vcompare(realtype c, N_Vector x, N_Vector z);
booleantype callml_vinvtest(N_Vector x, N_Vector z);
booleantype callml_vconstrmask(N_Vector c, N_Vector x, N_Vector m);
realtype callml_vminquotient(N_Vector num, N_Vector denom);

#endif
