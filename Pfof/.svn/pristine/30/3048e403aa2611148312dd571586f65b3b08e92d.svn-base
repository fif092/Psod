!> pFOF_cone is a distributed implementation of the Friends-of-Friends halo finder algorithm designed for lightcones simulated with RAMSES. <br>
!! This parallel Friend of Friend implementation is based on a sequential
!! implementation written by Edouard Audit (CEA). <br>
!! It has been written by Fabrice Roy -- CNRS / LUTH (Observatoire de Paris) and  Vincent Bouillot. <br>
!! mail: fabrice.roy@obspm.fr <br>
!! mail: vincent.bouillot@obspm.fr <br>

Program pfof_cone

  Use modhdf5
  Use modio
  Use modvariables
  Use modparameters
  Use modmpicom
  Use modfofpara

  Implicit None

  ! Initialization of MPI
  Call Mpi_Init(mpierr)
  Call Mpi_Comm_Size(Mpi_Comm_World, procNB, mpierr)
  Call Mpi_Comm_Rank(Mpi_Comm_World, procID, mpierr)

  If(procID==0) Then
     Call title()
  End If

  ! Initialization of HDF5
  Call hdf5_init()

  Call readparameters()

  Call h5readmap()

  Call h5readshellinfo()

  Call h5readparticles()

  Call computeminmax()

  Call fofparacone()

  Call deallocall()

  ! Finalize HDF5
  Call hdf5_finalize()

  If(procID==0) Then
     Call theend()
  End If

  ! Finalize MPI
  Call Mpi_Finalize(mpierr)

End Program pfof_cone
