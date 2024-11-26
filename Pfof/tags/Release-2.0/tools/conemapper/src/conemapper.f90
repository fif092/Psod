!=======================================================================
! Author: Fabrice Roy (LUTH/CNRS/Observatoire de Paris)
! Fabrice.Roy@obspm.fr
!=======================================================================

Program conemapper

  Use modhdf5
  Use modio
  Use modmap

  Implicit None

  ! Initialization of HDF5
  Call hdf5_init()

  Call readparameters()

  Call h5readshellinfo()

  Call exploremap()

  ! Finalize HDF5
  Call hdf5_finalize()

End Program conemapper
