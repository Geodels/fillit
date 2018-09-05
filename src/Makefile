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
FCFLAGS = -O3 -funroll-loops --param max-unroll-times=2 -cpp -shared \
    -ffree-form -ftree-vectorize -ffast-math -lstdc++ -ffree-line-length-none

# Depression filling algorithm using Planchon and Darboux
FMOD = queues

all:
	${F90} ${FFLAGS} -c ${FMOD}.f90

clean:
	rm -fv *~ *.bak *.o *.mod *.original *.so
	rm -rfv *.dSYM
