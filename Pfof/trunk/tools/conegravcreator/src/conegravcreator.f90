!=======================================================================
! Author: Fabrice Roy (LUTH/CNRS/Observatoire de Paris)
! Fabrice.Roy@obspm.fr
!=======================================================================

Program conegravcreator

  Use modparam
  Use modio
  Use modsortpart
  Use modvariable
  Use mpi
  Use modhdf5

  Implicit None

  Integer :: ilvl, fc, lc

  Call Mpi_Init(mpierr)
  Call Mpi_Comm_Rank(Mpi_Comm_World, procID, mpierr)
  Call Mpi_Comm_Size(Mpi_Comm_World, procNB, mpierr)

  Call hdf5_init()

  Call readparameters

  Call readramsesfiles

  Call dividespace()

  fc = 1
  lc = 0
  ! we sort each level independantly
  Do ilvl = 1, param%nlevel
     lc = lc + ncellperlevel(ilvl)
     If(ncellperlevel(ilvl) > 0) Then
        Call tritas(ncellperlevel(ilvl),idcube(fc:lc),pos(:,fc:lc),acc(:,fc:lc),&
             phi(fc:lc),rho(fc:lc),refined(fc:lc))
        fc = fc + ncellperlevel(ilvl)
     End If
  End Do

#ifdef WITHMPI3
  Call Mpi_Wait(req_sumnpc, mpistat, mpierr)
#endif


  Call h5writeconegrav()
     
  Deallocate(pos, acc, rho, phi, refined)
  Deallocate(ncelltab)

  Call hdf5_finalize()

  If(procID==0) Then
     Call theend()
  End If

  Call Mpi_Finalize(mpierr)

End Program conegravcreator
