MPIFC = h5pfc

OPTS = -O2 -g  -Wall   -I/usr/local/include -cpp -DLONGINT  -DWITHMPI
LDFLAGS = -g

VPATH = ../../../common/src:.

EXE = haloanalyzer
SRC = modhdf5.f90 \
	modvariable.f90 \
	modfunctions.f90 \
	modanalyzeengines.f90 \
	haloanalyzer.f90

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
	svnversion > haloanalyzer.version

clean :
	@echo "--------------------------------"
	@echo cleaning compiled objects, modules and executable 
	rm -f *.o *.mod *~ $(EXE)

