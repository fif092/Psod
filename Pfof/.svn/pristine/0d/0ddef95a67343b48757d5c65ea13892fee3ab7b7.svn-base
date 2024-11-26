include ../../../Make.inc

COMMONDIR = ../../../common/src/
VPATH= $(COMMONDIR):.

EXE = conegravcreator 

SRC = 	modconstant.f90 	\
	modmpicommons.f90 	\
	modhdf5.f90 		\
	modiocommons.f90 	\
	modreadinfo.f90 	\
	modvariables.f90	\
	modio.f90 		\
	modsortpart.f90 	\
	conegravcreator.f90

OBJ = $(SRC:%.f90=%.o)

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

$(EXE) : svnrev $(OBJ)
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
	else cd $(COMMONDIR) ; echo "#define SVNREV \"no_version_found\"" > svnrev.h ; cd - ;\
	fi


.PHONY: svnrev clean all
