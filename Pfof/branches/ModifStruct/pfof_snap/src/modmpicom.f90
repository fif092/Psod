!==============================================================================
! Project: pFoF
! File: pfof_snap/src/modmpicom.f90
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
!! This file contains the subroutines and the variables used for MPI communications.

!> This module contains the subroutines and the variables used for MPI communications.
!> Authors: F. Roy, V. Bouillot
Module modmpicom
  
  Use mpi
  Use modmpicommons
  Implicit none

  Type(Type_info_process) :: info_proc
  
  Integer(kind=4) :: dims(3) !< dimensions of the processes grid
  Integer(kind=4) :: MpiCube !< cartesian global communicator
  Integer(kind=4) :: MpiSubCubeWrite !< local communicator used for hdf5 parallel output
  Integer(kind=4) :: MpiSubCubeRead  !< local communicator used for hdf5 parallel input
  Integer(kind=4) :: CubeCoord(3)  !< coordinates of the process in the MpiCube communicator
  Integer(kind=4) :: commcolorWrite !< parameter used to create the MpiSubCubeWrite communicator
  Integer(kind=4) :: commcolorRead  !< parameter used to create the MpiSubCubeRead communicator
  Integer(kind=4) :: scprocIDWrite !< process id in the MpiSubCubeWrite communicator
  Integer(kind=4) :: scprocIDRead !< process id in the MpiSubCubeRead communicator
  Integer(kind=4) :: scprocNBWrite !< process number in the MpiSubCubeWrite communicator
  Integer(kind=4) :: scprocNBRead !< process number in the MpiSubCubeRead communicator
  Logical(kind=4) :: periods(3) !< periodicity of the global cartesian processes grid

  Private

  Public :: setcommunicators, &
       info_proc
!!$       dims, &
!!$       MpiCube, &
!!$       MpiSubCubeWrite, &
!!$       MpiSubCubeRead, &
!!$       scprocIDRead, &
!!$       scprocIDWrite, &
!!$       commcolorRead, &
!!$       commcolorWrite, &
!!$       CubeCoord

       

Contains

  !=======================================================================
  !> This subroutine creates the global cartesian grid of processes.
  !! It also creates the communicators used for the parallel hdf5 i/o if they are needed.
  Subroutine setcommunicators()

    Use modvariables, only : param
    Implicit None

    Integer(kind=4) :: mpierr

    ! Creation of a grid of processes.
    ! The grid is a cube, there are procNB**(1/3) processes in each dimension.
    ! The grid is periodic in each dimension.
    dims = int(procNB**(1./3.))
    periods = .true.
    
    Call Mpi_Cart_Create(Mpi_Comm_World,3,dims,periods,.true.,MPICube,mpierr)

!    Print *,'DIMS=',dims, '  MPI_CUBE=',MPICube, '   PROCID=',procID

    ! Process ID and number 
    Call Mpi_Comm_Rank  (MPICube, procID, mpierr)
    Call Mpi_Cart_Coords(MPICube, procID, 3, CubeCoord, mpierr)

    info_proc%global_comm%name = MPICube
    info_proc%global_comm%dims = dims
    info_proc%global_comm%pid = procID
    info_proc%global_comm%size = procNB
    info_proc%global_comm%periods = periods
    info_proc%global_comm%color = 0
    info_proc%global_comm%coords = CubeCoord
    
    
    ! Definition des neighbourss pour les 6 directions avec des conditions au bord periodiques
    ! Intialization of the 'neighboor processes id' 
    ! Neighboor : 1 = backward
    !             2 = forward
    !             3 = left
    !             4 = right
    !             5 = down
    !             6 = up
    
    Call Mpi_Cart_Shift(MpiCube,0,1,info_proc%global_comm%neighbours(1),&
         info_proc%global_comm%neighbours(2),mpierr)
    Call Mpi_Cart_Shift(MpiCube,1,1,info_proc%global_comm%neighbours(3),&
         info_proc%global_comm%neighbours(4),mpierr)
    Call Mpi_Cart_Shift(MpiCube,2,1,info_proc%global_comm%neighbours(5),&
         info_proc%global_comm%neighbours(6),mpierr)

!    info_proc%global_comm%neighbours = neighbours


    If(param%gatherwrite_factor > 1) Then
       If(procID==0) Then
          Print *, '*** Creating specific communicator for gathered cube output ***'  
       End If
       
       commcolorWrite = CubeCoord(1)/param%gatherwrite_factor + &
            (CubeCoord(2)/param%gatherwrite_factor) * (dims(1)/param%gatherwrite_factor) &
            + (CubeCoord(3)/param%gatherwrite_factor) * (dims(1)/param%gatherwrite_factor)*&
            (dims(1)/param%gatherwrite_factor)
       
       Call Mpi_Comm_Split(MpiCube, commcolorWrite, procID, MpiSubCubeWrite, mpierr)
       
       Call Mpi_Comm_Rank(MpiSubCubeWrite, scprocIDWrite, mpierr)
       Call Mpi_Comm_Size(MpiSubCubeWrite, scprocNBWrite, mpierr)

       info_proc%write_comm%name = MpiSubCubeWrite
       info_proc%write_comm%dims = (/ 0, 0, 0 /)
       info_proc%write_comm%pid = scprocIDWrite
       info_proc%write_comm%size = scprocNBWrite
       info_proc%write_comm%periods = (/ .false., .false., .false. /)
       info_proc%write_comm%color = commcolorWrite
       info_proc%write_comm%coords = (/ 0, 0, 0 /)
       info_proc%write_comm%neighbours = (/ 0, 0, 0, 0, 0, 0 /)

    Else
       info_proc%write_comm%name = -1
       info_proc%write_comm%dims = (/ -1, -1, -1 /)
       info_proc%write_comm%pid = -1
       info_proc%write_comm%size = -1
       info_proc%write_comm%periods = (/ .false., .false., .false. /)
       info_proc%write_comm%color = -1
       info_proc%write_comm%coords = (/ -1, -1, -1 /)
       info_proc%write_comm%neighbours = (/ -1, -1, -1, -1, -1, -1 /)

    End If

    If(param%gatherread_factor > 1) Then

       If(procID==0) Then
          Print *, '*** Creating specific communicator for gathered cube input ***'  
       End If

       commcolorRead = CubeCoord(1)/param%gatherread_factor +  &
            (CubeCoord(2)/param%gatherread_factor) * (dims(1)/param%gatherread_factor) &
            + (CubeCoord(3)/param%gatherread_factor) * (dims(1)/param%gatherread_factor)* &
            (dims(1)/param%gatherread_factor)
       
       Call Mpi_Comm_Split(MpiCube, commcolorRead, procID, MpiSubCubeRead, mpierr)
       
       Call Mpi_Comm_Rank(MpiSubCubeRead, scprocIDRead, mpierr)
       Call Mpi_Comm_Size(MpiSubCubeRead, scprocNBRead, mpierr)

       info_proc%read_comm%name = MpiSubCubeRead
       info_proc%read_comm%dims = (/ 0, 0, 0 /)
       info_proc%read_comm%pid = scprocIDRead
       info_proc%read_comm%size = scprocNBRead
       info_proc%read_comm%periods = (/ .false., .false., .false. /)
       info_proc%read_comm%color = commcolorRead
       info_proc%read_comm%coords = (/ 0, 0, 0 /)
       info_proc%read_comm%neighbours = (/ 0, 0, 0, 0, 0, 0 /)

    Else

       info_proc%read_comm%name = -1
       info_proc%read_comm%dims = (/ -1, -1, -1 /)
       info_proc%read_comm%pid = -1
       info_proc%read_comm%size = -1
       info_proc%read_comm%periods = (/ .false., .false., .false. /)
       info_proc%read_comm%color = -1
       info_proc%read_comm%coords = (/ -1, -1, -1 /)
       info_proc%read_comm%neighbours = (/ -1, -1, -1, -1, -1, -1 /)

    End If

  End Subroutine setcommunicators

End Module modmpicom
