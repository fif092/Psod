!=======================================================================
! Author: Fabrice Roy (LUTH/CNRS/Observatoire de Paris)
! Fabrice.Roy@obspm.fr
!=======================================================================

Program conecreator

!  Use modparam
  Use modio
  Use modsortpart
  Use modvariables
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

  If(param%do_read_ramses_part_id) Then
     If(.not.param%do_read_potential) Then
!        Call tritas(npartloc,idcube, pos, vel,ramsespartid)
        Call heapsort(npartloc, idcube, pos, vel, tid=ramsespartid)
     Else
        If(.not.param%do_read_gravitational_field) Then
           Call heapsort(npartloc, idcube, pos, vel, tp=pot, tid=ramsespartid)
        Else
           Call heapsort(npartloc, idcube, pos, vel, tf=field, tp=pot, tid=ramsespartid)
        End If
     End If
  Else
!     Call tritas(npartloc,idcube, pos, vel)
     If(.not.param%do_read_potential) Then
        Call heapsort(npartloc, idcube, pos, vel)
     Else
        If(.not.param%do_read_gravitational_field) Then
           Call heapsort(npartloc, idcube, pos, vel, tp=pot)
        Else
           Call heapsort(npartloc, idcube, pos, vel, tf=field, tp=pot)
        End If
     End If
     
  End If

#ifdef WITHMPI3
  Call Mpi_Wait(req_sumnpc, mpistat, mpierr)
#endif

  Call h5writecone()
     
  Deallocate(pos, vel, id)
  Deallocate(nparttab)

  Call hdf5_finalize()
  
  If(procID==0) Then
     Call theend()
  End If

  Call Mpi_Finalize(mpierr)

End Program conecreator
