#*****************************************************************************
# Copyright 2018 Tristan Salles
#
# fillIT is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or any later version.
#
# fillIT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with fillIT.  If not, see <http://www.gnu.org/licenses/>.
#*****************************************************************************

## Makefile used to wrap Fortran functions used by fillIT

.PHONY: all clean

DIR=$(shell pwd)

F90 = gfortran
F2PY = f2py
FFLAGS = -shared -O2 -fPIC

# Fortran optimisation flags
FCFLAGS = -O3 -funroll-loops --param max-unroll-times=2 -cpp \
    -ffree-form -ftree-vectorize -ffast-math -lstdc++ -ffree-line-length-none

# F2PY optimisation flags
# F2PY_FLAGS = --f90exec=${F90} --opt='-O3 -ftree-vectorize '

# Depression filling algorithm using Planchon and Darboux
FMOD = queues

# LMOD = fillBarnes
# L2MOD = fillWei
# L3MOD = fillZhou

all:
	${F90} ${FFLAGS} -c ${FMOD}.f90

# ${F2PY} ${F2PY_FLAGS} -c ${LMOD}.pyf  ${FMOD}.f90 ${LMOD}.f90
# ${F2PY} ${F2PY_FLAGS} -c ${L2MOD}.pyf  ${FMOD}.f90 ${L2MOD}.f90
# ${F2PY} ${F2PY_FLAGS} -c ${L3MOD}.pyf  ${FMOD}.f90 ${L3MOD}.f90

clean:
	rm -fv *~ *.bak *.o *.mod *.original *.so
	rm -rfv *.dSYM
