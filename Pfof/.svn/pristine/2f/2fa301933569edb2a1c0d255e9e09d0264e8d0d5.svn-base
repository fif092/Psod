!==============================================================================
! Project: pFoF
! File: pfof_cone/src/modmpicom.f90
! Copyright Fabrice Roy and Vincent Bouillot (2011)
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


!> @file
!! This file contains the subroutines and the variables used to prepare MPI communications 
!! and a routine used to broadcast the input parameters.

!> This module contains the subroutines and the variables used to prepare MPI communications 
!> and a routine used to broadcast the input parameters.
!>
!> Authors: F. Roy, V. Bouillot

Module modmpicom
  
  Use mpi
  Use modmpicommons

  Integer(kind=4) :: mpierr !< MPI error code

  Type(Type_info_process) :: info_proc

    
  Private
  Public :: info_proc, &
       settopology

Contains

  !=======================================================================
  !> This subroutine create a processes topology that maps the lightcone.
  Subroutine settopology()

    Implicit None

    Integer :: in

    Do in = 1, 6
       If(info_proc%global_comm%neighbours(in) == -1) Then
          info_proc%global_comm%neighbours(in) = MPI_PROC_NULL
       End If
    End Do

    ! global_comm def
    info_proc%global_comm%name = Mpi_Comm_World
    ! info_proc%global_comm%dims is read in modio.f90, in subroutine h5readmap
    info_proc%global_comm%pid = procID
    info_proc%global_comm%size = procNB
    info_proc%global_comm%periods = (/ .false., .false., .false. /)
    info_proc%global_comm%color = -1
    ! info_proc%global_comm%coords is read in modio.f90, in subroutine h5readmap
    ! info_proc%global_comm%neighbours is read in modio.f90, in subroutine h5readmap
    
    ! write_comm def
    info_proc%write_comm%color = info_proc%global_comm%coords(1) - 1
    
    Call Mpi_Comm_Split(info_proc%global_comm%name, info_proc%write_comm%color, procID, &
         info_proc%write_comm%name, mpierr)
    Call Mpi_Comm_Rank(info_proc%write_comm%name, info_proc%write_comm%pid, mpierr)
    Call Mpi_Comm_Size(info_proc%write_comm%name, info_proc%write_comm%size, mpierr) 
    
    info_proc%write_comm%dims = (/ 0, 0, 0 /)
    info_proc%write_comm%periods = (/ .false., .false., .false. /)
    info_proc%write_comm%coords = (/ 0, 0, 0 /)
    info_proc%write_comm%neighbours = (/ -1, -1, -1, -1, -1, -1 /)

#ifdef DEBUG
    Print *,'Proc',procID, ': cubecoord=',info_proc%global_comm%coords, &
         ' and ID in subcomm=',info_proc%write_comm%pid
#endif

    
    ! read_comm def: not used by pfof_cone for now
    info_proc%read_comm%name = -1
    info_proc%read_comm%dims = (/ -1, -1, -1 /)
    info_proc%read_comm%pid = -1
    info_proc%read_comm%size = -1
    info_proc%read_comm%periods = (/ .false., .false., .false. /)
    info_proc%read_comm%color = -1
    info_proc%read_comm%coords = (/ -1, -1, -1 /)
    info_proc%read_comm%neighbours = (/ -1, -1, -1, -1, -1, -1 /)
   
  End Subroutine settopology

End Module modmpicom
