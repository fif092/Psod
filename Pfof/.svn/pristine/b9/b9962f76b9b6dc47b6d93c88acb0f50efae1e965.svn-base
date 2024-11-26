include Make.inc

all: pfofhdf5 pfofcone conepartcreator conegravcreator conemapper

pfofhdf5:
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
