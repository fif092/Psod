MPIFC=h5pfc

# Intel release flags
FCFLAGS= -O3 -g -cpp -DLONGINT 
# Intel debug flags
#FCFLAGS= -O0 -g -warn all -C -traceback -cpp -DLONGINT -convert big_endian -DDEBUG

# GNU release flags
#FCFLAGS= -O3 -g -cpp -DLONGINT  #-fconvert=big-endian
# GNU debug flags
#FCFLAGS= -O0 -g -Wall -Wextra -std=f2003 -fcheck=all -fbacktrace -cpp -DDEBUG -DLONGINT -fconvert=big-endian

LDFLAGS=-g

BINDIR= 