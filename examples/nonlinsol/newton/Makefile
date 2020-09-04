# This files is a wrapper around examples.mk that adds build rules for
# test_sunnonlinsol_*.ml.

FILES_TO_CLEAN=

USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=nonlinsol/newton
C_SUBDIR=sunnonlinsol/newton

EXAMPLES = $(if $(AT_LEAST_4_0_0), test_sunnonlinsol_newton.ml)

include ../nonlinsol.mk
