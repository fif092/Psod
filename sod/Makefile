include ../Pfof/trunk/Make.inc
#MPIFC = h5pfc

#OPTS = -O3 -g -cpp -DLONGINT -DREADHDF5 -traceback
#OPTS =   -O0 -g -traceback -fpe0 -ftrapuv -check bounds -cpp -DLONGINT -DREADHDF5 
#-warn all

#LDFLAGS = -g

EXE = sod_profile_rdfof_parallel_hdf5_zoom

VPATH = .:../Pfof/trunk/common/src
  
SRC =  	modconstant.f90     \
	modparameters.f90   \
        modvariables.f90    \
        modmpicom.f90       \
	modhdf5.f90         \
	modextract_zoom.f90      \
	modiocommons.f90    \
	modwritehalo_zoom.f90    \
	sod_profile_rdfof_parallel_hdf5_zoom.f90


FCFLAGS2=-DREADHDF5 -DWRITEHDF5
OBJ = $(SRC:%.f90=%.o)


%.o:%.f90
	@echo "--------------------------------"
	@echo building $*.o with parameters \
	FC=$(MPIFC), FCFLAGS=$(FCFLAGS) $(FCFLAGS2)
	$(MPIFC) -c $(FCFLAGS) $(FCFLAGS2) $<  

%.mod:%.f90
	@echo "--------------------------------"
	@echo building $*.mod with parameters \
	FC=$(MPIFC), FCFLAGS=$(FCFLAGS) $(FCFLAGS2)
	$(MPIFC) -c $(FCFLAGS) $(FCFLAGS2) $< 

all : $(EXE)

$(EXE) : $(OBJ)
	@echo "--------------------------------"
	@echo building $(EXE) 
	$(MPIFC) -o $(EXE) $(OBJ)  $(LDFLAGS)


clean :
	@echo "--------------------------------"
	@echo cleaning compiled objects, modules and executable 
	rm -f *.o *.mod *~ $(EXE)

