templates = $(wildcard *.xml)
objects = $(templates:.xml=.o)

#===============================================================================
# Compiler Options
#===============================================================================

# Add include for subdirectory

override F90FLAGS += -I../xml-fortran

# Ignore unused variables

ifeq ($(F90),ifort)
  override F90FLAGS += -warn nounused
endif

#===============================================================================
# Targets
#===============================================================================

all: $(objects)
clean:
	@rm -f *.o *.mod *.out *.f90
	@rm -rf CodeCoverage *.dyn gmon.out pgopti* *.HTML

#===============================================================================
# Rules
#===============================================================================

.SUFFIXES: .f90 .o .xml
.PHONY: all clean
.PRECIOUS: %.f90

%.f90: %.xml
	../xml-fortran/xmlreader $(basename $@)

%.o: %.f90
	$(F90) $(F90FLAGS) -c $(basename $@).f90
