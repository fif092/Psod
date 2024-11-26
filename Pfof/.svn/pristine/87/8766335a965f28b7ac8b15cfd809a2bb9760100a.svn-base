include ../../../Make.inc
#MPIFC = h5pfc

#OPTS = -O3 -g -cpp -DLONGINT -Wall -fcheck=all 

#LDFLAGS = -g

VPATH = ../../../common/src:.

EXE = testextract

SRC = modhdf5.f90 		\
	modconstant.f90 	\
	modiocommons.f90	\
	modvariables.f90 	\
	modextract.f90 		\
	testextract.f90

OBJ = $(SRC:%.f90=%.o)
#FCFLAGS = $(OPTS)

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

