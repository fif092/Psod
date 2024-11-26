MPIFC = h5pfc

OPTS = -O2 -g  -I/usr/local/include -cpp -DLONGINT  -fcheck=all
LDFLAGS = -g

VPATH = ../../../common/src:.

EXE = conegravcreator 

SRC = modhdf5.f90 \
	modparam.f90 \
	modvariable.f90 \
	modio.f90 \
	modsortpart.f90 \
	conegravcreator.f90

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
	svnversion > conegravcreator.version


clean :
	@echo "--------------------------------"
	@echo cleaning compiled objects, modules and executable 
	rm -f *.o *.mod *~ $(EXE) conegravcreator.version

