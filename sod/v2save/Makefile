MPIFC = h5pfc

#OPTS = -O3 -g -cpp -DLONGINT -DREADHDF5 -traceback
OPTS =   -O0 -g -traceback -fpe0 -ftrapuv -check bounds -cpp -DLONGINT -DREADHDF5 
#-warn all

LDFLAGS = -g

VPATH = .

EXE = sod_profile_rdfof_parallel_hdf5

SRC =  	modconstant.f90     \
	modparameters.f90   \
        modvariables.f90    \
        modmpicom.f90       \
	modhdf5.f90         \
	modextract.f90      \
	modwritehalo.f90    \
	sod_profile_rdfof_parallel_hdf5.f90

OBJ = $(SRC:%.f90=%.o)
FCFLAGS = $(OPTS)

%.o:%.f90
	@echo "--------------------------------"
	@echo building $*.o with parameters \
	FC=$(MPIFC), FCFLAGS=$(FCFLAGS)
	$(MPIFC) -c $(FCFLAGS) $<  

%.mod:%.f90
	@echo "--------------------------------"
	@echo building $*.mod with parameters \
	FC=$(MPIFC), FCFLAGS=$(FCFLAGS)
	$(MPIFC) -c $(FCFLAGS) $< 

all : $(EXE)

$(EXE) : $(OBJ)
	@echo "--------------------------------"
	@echo building $(EXE) 
	$(MPIFC) -o $(EXE) $(OBJ)  $(LDFLAGS)


clean :
	@echo "--------------------------------"
	@echo cleaning compiled objects, modules and executable 
	rm -f *.o *.mod *~ $(EXE)

