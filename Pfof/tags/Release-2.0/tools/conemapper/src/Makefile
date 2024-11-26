MPIFC = h5pfc

OPTS = -O2 -g  -Wall   -I/usr/local/include -cpp -DLONGINT -fbounds-check
LDFLAGS = -g

VPATH = ../../../common/src:.


SRC = modhdf5.f90 		\
	modparameters.f90	\
	modvariables.f90	\
	modio.f90		\
	modmap.f90		\
	conemapper.f90

OBJ = $(SRC:%.f90=%.o)
FCFLAGS = $(OPTS) -DWITHHDF5
EXE = conemapper


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

