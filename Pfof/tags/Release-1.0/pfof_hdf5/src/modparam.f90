!> @file
!! This file contains the definition of the input parameters.

!> This module contains the definition of the input parameters.
Module modparam

  Use modconst
  Character(len=3)   :: code_index     !< version of RAMSES used to produce the RAMSES files to analyze
  Character(len=200) :: pathinput      !< path to the directory containing the RAMSES files
  Character(len=200) :: namepart       !< base name of the files containing the particles information
  Character(len=200) :: nameinfo       !< name of the RAMSES info file
  Integer(kind=4)    :: grpsize        !< size of the I/O group used for the RAMSES simulation
  Logical            :: star           !< do the RAMSES files contain stars?
  Logical            :: metal          !< do the RAMSES files contain metalicities?
  Logical            :: readfromcube   !< should pFOF read particles information from previously created cube files?
  Real(kind=SP)      :: perco          !< FOF percolation length
  Integer(kind=4)    :: Mmin           !< minimum mass of the haloes
  Integer(kind=4)    :: Mmax           !< maximum mass of the haloes
  Logical            :: dofof          !< should pFOF perform FOF halo finding? (if not, it can still write cube files)
  Character(len=390) :: root           !< base name of the output files
  Logical            :: outcube        !< should pFOF write cube files?
  Logical            :: usehdf5        !< should pFOF use hdf5 file format for output (and also cube input)
  Integer(kind=4)    :: gatherwrite    !< gather parameter for parallel hdf5 output
  Integer(kind=4)    :: gatherread     !< gather parameter for parallel hdf5 input
  Logical            :: dotimings      !< should pFOF perform timings (imply extra synchronizations)

  ! Namelist for input file
  Namelist / ramsesinput / code_index, pathinput, namepart, nameinfo, grpsize, star, metal, readfromcube
  Namelist / fofparam / perco, Mmin, Mmax, dofof
  Namelist / outputparam / root, outcube, usehdf5, gatherwrite, gatherread, dotimings


End Module modparam
