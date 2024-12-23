!> @file
!!This file contains the parameters read from pfof_cone.nml and the version number read from pfof.version. 

!> This module contains the parameters read from pfof_cone.nml and the version number read from pfof.version.
!>
!> Authors: F. Roy, V. Bouillot


Module modparameters

  Implicit None

  Character(len=200) :: simulationName !< Simulation name to be written in output files: contains the name of the cone (narrow or fullsky for instance)
  Integer            :: outputNB       !< Ramses output number to be written in output files: first output_ncoarse in the cone
  Character(len=200) :: shell_dir
  Character(len=200) :: shell_filename
  Integer(kind=4) :: shell_fid
  Integer(kind=4) :: shell_lid

  Real(kind=4) :: perco
  Integer(kind=4) :: Mmin
  Integer(kind=4) :: Mmax
  Integer(kind=4) :: nres
  Logical(kind=4) :: doUnbinding    !< Unbinding process? (not implemented yet)
  Logical(kind=4) :: doSubHalo      !< SubHalo detection? (not implemented yet)
  Logical(kind=4), parameter :: potential=.false.      !< Potential at the position of the particle? (here for compatibility reason with pfof
  Logical(kind=4), parameter :: force=.false.      !< Force at the position of the particle? (here for compatibility reason with pfof
  
  Character(len=400) :: output_root
  Logical(kind=4) :: dotimings

  Character(len=50) :: origin
  
  Namelist / shellparameters / simulationName, outputNB, shell_dir, shell_filename, shell_fid, shell_lid
  Namelist / fofparameters / perco, Mmin, Mmax, nres, doUnbinding, doSubHalo
  Namelist / outputparameters / output_root, dotimings

End Module modparameters
