!> @file
!! This file contains the definition of the input parameters.

!> This module contains the definition of the input parameters.
!> Authors: F. Roy, V. Bouillot
Module modparameters

  Implicit None
  Character(len=200) :: simulationName !< Simulation name to be written in output files
  Integer            :: outputNB       !< Ramses output number to be written in output files
  Character(len=3)   :: code_index     !< version of RAMSES used to produce the RAMSES files to analyze
  Character(len=200) :: pathinput      !< path to the directory containing the RAMSES files
  Character(len=200) :: namepart       !< base name of the files containing the particles information
  Character(len=200) :: nameinfo       !< name of the RAMSES info file
  Integer(kind=4)    :: grpsize        !< size of the I/O group used for the RAMSES simulation
  Logical(kind=4)    :: star           !< do the RAMSES files contain stars?
  Logical(kind=4)    :: metal          !< do the RAMSES files contain metalicities?
  Logical(kind=4)    :: potential      !< do the RAMSES files contain potential?
  Logical(kind=4)    :: force          !< do the RAMSES files contain force?
  Logical(kind=4)    :: readfromcube   !< should pFOF read particles information from previously created cube files?
  Logical(kind=4)    :: readfromsortedcube !< should pFOF read particles information from previously created sorted cube files?
  Real(kind=4)       :: perco          !< FOF percolation length
  Integer(kind=4)    :: Mmin           !< minimum mass of the haloes
  Integer(kind=4)    :: Mmax           !< maximum mass of the haloes
  Logical(kind=4)    :: dofof          !< should pFOF perform FOF halo finding? (if not, it can still write cube files)
  Logical(kind=4)    :: doUnbinding    !< Unbinding process? (not implemented yet)
  Logical(kind=4)    :: doSubHalo      !< SubHalo detection? (not implemented yet)
  Character(len=400) :: output_root    !< base name of the output files
  Logical(kind=4)    :: outcube        !< should pFOF write cube files?
  Integer(kind=4)    :: gatherwrite    !< gather parameter for parallel hdf5 output
  Integer(kind=4)    :: gatherread     !< gather parameter for parallel hdf5 input
  Logical(kind=4)    :: sortcube       !< sort the particles and write a 'sorted cube'
  Logical(kind=4)    :: dotimings      !< should pFOF perform timings (imply extra synchronizations)

  Character(len=50) :: origin          !< version number of the software (read from a .version text file)

  ! Namelist for input file
  Namelist / ramsesinput / simulationName, outputNB, code_index, pathinput, namepart, nameinfo, grpsize, &
       star, metal, potential, force, readfromcube, readfromsortedcube, gatherread
  Namelist / fofparam / perco, Mmin, Mmax, dofof, doUnbinding, doSubHalo
  Namelist / outputparam / output_root, outcube, gatherwrite,  sortcube, dotimings


End Module modparameters
