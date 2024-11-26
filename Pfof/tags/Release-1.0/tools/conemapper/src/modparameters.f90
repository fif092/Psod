Module modparameters

  Implicit None

  Character(len=200) :: shell_dir
  Character(len=196) :: shell_filename
  Integer(kind=4) :: shell_fid
  Integer(kind=4) :: shell_lid

  Real(kind=4) :: perco
  Integer(kind=4) :: Mmin
  Integer(kind=4) :: Mmax
  Integer(kind=4) :: nres
  
  Character(len=150) :: output_root
  
  Namelist / shellparameters / shell_dir, shell_filename, shell_fid, shell_lid
  Namelist / fofparameters / perco, Mmin, Mmax, nres
  Namelist / outputparameters / output_root

End Module modparameters
