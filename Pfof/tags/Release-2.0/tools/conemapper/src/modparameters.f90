!=======================================================================
! Author: Fabrice Roy (LUTH/CNRS/Observatoire de Paris)
! Fabrice.Roy@obspm.fr
!=======================================================================
Module modparameters

  Implicit None

  Character(len=200) :: simulationName
  Integer(kind=4) :: outputNB
  Character(len=200) :: shell_dir
  Character(len=196) :: shell_filename
  Integer(kind=4) :: shell_fid
  Integer(kind=4) :: shell_lid

  Real(kind=4) :: perco
  Integer(kind=4) :: Mmin
  Integer(kind=4) :: Mmax
  Integer(kind=4) :: nres
  Logical :: doUnbinding
  Logical :: doSubHalo
  
  Character(len=150) :: output_root
  Logical :: dotimings
  
  Namelist / shellparameters / simulationName, outputNB, shell_dir, shell_filename, shell_fid, shell_lid
  Namelist / fofparameters / perco, Mmin, Mmax, nres, doUnbinding, doSubHalo
  Namelist / outputparameters / output_root, dotimings

End Module modparameters
