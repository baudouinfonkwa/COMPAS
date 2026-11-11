# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE
DEBUG      = FALSE

DIM        = 2
# DIM        = 3

COMP	   = gnu

USE_MPI    = TRUE
USE_OMP    = FALSE 
USE_CUDA   = FALSE
USE_HIP    = FALSE
USE_SYCL   = FALSE

# Optional phase count, isolated from the existing executable.
NPHASE ?= 2
ifneq ($(NPHASE),2)
USERSuffix := .N$(NPHASE)
endif

Bpack   := ./prob/Make.package 
Blocs   := . ./prob


include ../Make.exec


#############################################
# CODE SPECIFIC DEFINTIONS
#############################################

DEFINES += -DPHYSICS=SIXEQS_IE_NPHASE

DEFINES += -DNPHASE=$(NPHASE)

DEFINES += -DADVECTION=true
DIFFUSION ?= false
DEFINES += -DDIFFUSION=$(DIFFUSION)
DEFINES += -DNONCONSERVATIVE=true


