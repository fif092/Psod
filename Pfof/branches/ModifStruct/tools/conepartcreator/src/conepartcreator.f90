!==============================================================================
! Project: pFoF
! File: tools/conepartcreator/src/conepartcreator.f90
! Copyright Fabrice Roy (2015)
! Fabrice.Roy@obspm.fr
! 
! This software is a computer program whose purpose is to detect dark matter
! haloes in cosmological simulation with a parallel Friend of Friend algorithm. 
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can  use, 
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info". 
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability. 
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or 
! data to be ensured and, more generally, to use and operate it in the 
! same conditions as regards security. 
!
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL license and that you accept its terms.
!==============================================================================

!=======================================================================
! Author: Fabrice Roy (LUTH/CNRS/Observatoire de Paris)
! Fabrice.Roy@obspm.fr
!=======================================================================

Program conepartcreator

  Use modio
  Use modsortpart
  Use modvariables
  Use mpi
  Use modhdf5
  Use modmpicommons, only : procNB, procID

  Implicit None

  Integer(kind=4) :: mpierr

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
  Call Mpi_Finalize(mpierr)

End Program conepartcreator
