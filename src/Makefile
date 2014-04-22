program = ndpp
prefix = /usr/local

templates = $(wildcard templates/*.o)
xml_fort = xml-fortran/xmlparse.o \
           xml-fortran/read_xml_primitives.o \
           xml-fortran/write_xml_primitives.o

#===============================================================================
# Object Files
#===============================================================================

include OBJECTS

#===============================================================================
# User Options
#===============================================================================

COMPILER = gnu
DEBUG    = no
PROFILE  = no
OPTIMIZE = no
NATIVE   = no
OPENMP   = no
MPI      = no
HDF5     = no

#===============================================================================
# External Library Paths
#===============================================================================

# These are set up to work with Arch Linux; the user should revise these
# as necessary for their distribution
MPI_DIR   = /usr
HDF5_DIR  = /usr

#===============================================================================
# Add git SHA-1 hash
#===============================================================================

GIT_SHA1 = $(shell git log -1 2>/dev/null | head -n 1 | awk '{print $$2}')

#===============================================================================
# GNU Fortran compiler options
#===============================================================================

ifeq ($(COMPILER),gnu)
  F90 = gfortran
  F90FLAGS := -cpp -fbacktrace -DNO_F2008
  LDFLAGS =

  # Debugging
  ifeq ($(DEBUG),yes)
    F90FLAGS += -g -Wall -pedantic -std=f2008 -fbounds-check \
                -ffpe-trap=invalid,overflow,underflow
    LDFLAGS  += -g
  endif

  # Profiling
  ifeq ($(PROFILE),yes)
    F90FLAGS += -pg
    LDFLAGS  += -pg
  endif

  # Optimization
  ifeq ($(OPTIMIZE),yes)
    F90FLAGS += -O3
  endif

  # Native machine compilarion
  ifeq ($(NATIVE), yes)
  	F90FLAGS += -march=native -mtune=native
  endif
endif

#===============================================================================
# Intel Fortran compiler options
#===============================================================================

ifeq ($(COMPILER),intel)
  F90 = ifort
  F90FLAGS := -fpp -warn -assume byterecl -traceback
  LDFLAGS =

  # Debugging
  ifeq ($(DEBUG),yes)
    F90FLAGS += -g -ftrapuv -fp-stack-check -check all -fpe0
    LDFLAGS  += -g
  endif

  # Profiling
  ifeq ($(PROFILE),yes)
    F90FLAGS += -pg -prof-gen:srcpos
    LDFLAGS  += -pg -prof-gen:srcpos
  endif

  # Optimization
  ifeq ($(OPTIMIZE),yes)
    F90FLAGS += -O3
  endif

  # Native machine compilarion
  ifeq ($(NATIVE), yes)
  	F90FLAGS += -xHost -march=native -mtune=native -ipo
  endif
endif

#===============================================================================
# Setup External Libraries
#===============================================================================

# MPI for distributed-memory parallelism

ifeq ($(MPI),yes)
  F90 = $(MPI_DIR)/bin/mpif90
  F90FLAGS += -DMPI
endif

# OpenMP for shared-memory parallelism

ifeq ($(OPENMP),yes)
	ifeq ($(COMPILER),intel)
    F90FLAGS += -openmp -DOPENMP
    LDFLAGS += -openmp
	endif

	ifeq ($(COMPILER),gnu)
    F90FLAGS += -fopenmp -DOPENMP
    LDFLAGS += -fopenmp
	endif
endif

# HDF5 for I/O

ifeq ($(HDF5),yes)
	F90 = $(HDF5_DIR)/bin/h5pfc
  F90FLAGS += -DHDF5
endif

#===============================================================================
# Targets
#===============================================================================

all: xml-fortran $(program)
xml-fortran:
	cd xml-fortran; make MACHINE=$(MACHINE) F90=$(F90) F90FLAGS="$(F90FLAGS)" LDFLAGS="$(LDFLAGS)"
	cd templates; make F90=$(F90) F90FLAGS="$(F90FLAGS)" LDFLAGS="$(LDFLAGS)"
$(program): $(objects) $(omc_objects)
	$(F90) $(objects) $(omc_objects) $(templates) $(xml_fort) $(LDFLAGS) -o $@
tests:
	make all
	cd ../tests/test_scatt; make COMPILER=$(COMPILER) DEBUG=$(DEBUG) \
			OPTIMIZE=$(OPTIMIZE) PROFILE=$(PROFILE) MPI=$(MPI) HDF5=$(HDF5) \
			OPENMP=$(OPENMP)
install:
	@install -D $(program) $(DESTDIR)$(prefix)/bin/$(program)
	@install -D ../man/man1/ndpp.1 $(DESTDIR)$(prefix)/share/man/man1/ndpp.1
	@install -D ../LICENSE $(DESTDIR)$(prefix)/share/doc/$(program)/copyright
uninstall:
	@rm $(DESTDIR)$(prefix)/bin/$(program)
	@rm $(DESTDIR)$(prefix)/share/man/man1/ndpp.1
	@rm $(DESTDIR)$(prefix)/share/doc/$(program)/copyright
distclean: clean
	cd xml-fortran; make clean
	cd templates; make clean
	cd ../tests/test_scatt; make clean
clean:
	@rm -f *.o *.mod $(program)
neat:
	@rm -f *.o *.mod

#===============================================================================
# Rules
#===============================================================================

.SUFFIXES: .F90 .o
.PHONY: all xml-fortran tests install uninstall clean neat distclean

%.o: %.F90
	$(F90) $(F90FLAGS) -DGIT_SHA1="\"$(GIT_SHA1)\"" -Ixml-fortran -Itemplates -c $<

#===============================================================================
# Dependencies
#===============================================================================

include DEPENDENCIES
