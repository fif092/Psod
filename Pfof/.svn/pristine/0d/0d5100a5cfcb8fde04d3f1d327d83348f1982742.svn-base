include Make.inc

all: pfofsnap pfofcone conepartcreator conegravcreator conemapper

pfofsnap:
	cd pfof_snap/src ; make

pfofcone:
	cd pfof_cone/src ; make

conepartcreator:
	cd tools/conepartcreator/src ; make

conegravcreator:
	cd tools/conegravcreator/src ; make

conemapper:
	cd tools/conemapper/src ; make

clean:
	cd pfof_snap/src ; make clean
	cd pfof_cone/src ; make clean
	cd tools/conepartcreator/src ; make clean
	cd tools/conegravcreator/src ; make clean
	cd tools/conemapper/src ; make clean

install:
	cp pfof_snap/src/pfof_snap $(BINDIR)
	cp pfof_cone/src/pfof_cone $(BINDIR)
	cp tools/conegravcreator/src/conegravcreator $(BINDIR)
	cp tools/conepartcreator/src/conepartcreator $(BINDIR)
	cp tools/conemapper/src/conemapper $(BINDIR)
