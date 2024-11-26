!=======================================================================
! Author: Fabrice Roy (LUTH/CNRS/Observatoire de Paris)
! Fabrice.Roy@obspm.fr
!=======================================================================

Program conecreator

  Use modparam
  Use modio
  Use modsortpart
  Use modvariable
  Use mpi
  Use modhdf5

  Implicit None


  Call Mpi_Init(mpierr)
  Call Mpi_Comm_Rank(Mpi_Comm_World, procID, mpierr)
  Call Mpi_Comm_Size(Mpi_Comm_World, procNB, mpierr)

  Call hdf5_init()


  Call readparameters

  Call readramsesfiles


  Call dividespace()

  Call tritas(idcube, pos, vel, npartloc)

#ifdef WITHMPI3
  Call Mpi_Wait(req_sumnpc, mpistat, mpierr)
#endif

  Call h5writecone()
     
  Deallocate(pos, vel, id)
  Deallocate(nparttab)

  Call hdf5_finalize()
  Call Mpi_Finalize(mpierr)

End Program conecreator
