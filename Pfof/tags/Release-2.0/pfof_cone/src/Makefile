MPIFC = h5pfc

OPTS = -O3 -g -Wall -cpp -DLONGINT
#OPTS = -O0 -g -Wall -cpp -DLONGINT -DDEBUG -fcheck=all -fbacktrace 
LDFLAGS = -g

VPATH = ../../common/src:.


SRC = modhdf5.f90 		\
	modconstant.f90 	\
	modparameters.f90	\
	modtiming.f90		\
	modvariables.f90	\
	modmpicom.f90		\
	modsortinterf.f90		\
	modxdmf.f90		\
	modwritehalo.f90	\
	modio.f90		\
	modhalo.f90		\
	modfofmpi.f90		\
	modfofpara.f90		\
	pfof_cone.f90

OBJ = $(SRC:%.f90=%.o)
FCFLAGS = $(OPTS) -DWITHHDF5
EXE = pfof_cone


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
	svnversion > pfof.version
clean :
	@echo "--------------------------------"
	@echo cleaning compiled objects, modules and executable 
	rm -f *.o *.mod *~ $(EXE) pfof.version

