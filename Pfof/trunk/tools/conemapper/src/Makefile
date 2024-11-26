include ../../../Make.inc

COMMONDIR = ../../../common/src/
VPATH = $(COMMONDIR):.

SRC = modhdf5.f90 		\
	modconstant.f90		\
	modiocommons.f90	\
	modreadparameters.f90	\
	modvariables.f90	\
	modio.f90		\
	modmap.f90		\
	conemapper.f90

OBJ = $(SRC:%.f90=%.o)

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

svnrev: 
	@SVNVERSION=$(shell svnversion -n . 2>/dev/null) ; \
	if [ $$? -eq 0 ] ; then cd $(COMMONDIR) ; \
	echo "#define SVNREV \"$(shell svnversion -n .)\"" > svnrev.h ; cd - ;\
	else cd $(COMMONDIR) ; echo "#define SVNREV \"236\"" > svnrev.h ; cd - ;\
	fi


.PHONY: svnrev clean all
