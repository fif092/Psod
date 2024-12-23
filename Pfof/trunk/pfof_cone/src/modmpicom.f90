!> @file
!! This file contains the subroutines and the variables used to prepare MPI communications 
!! and a routine used to broadcast the input parameters.

!> This module contains the subroutines and the variables used to prepare MPI communications 
!> and a routine used to broadcast the input parameters.
!>
!> Authors: F. Roy, V. Bouillot

Module modmpicom
  
  Use mpi

  Integer(kind=4) :: mpierr !< MPI error code
  Integer(kind=4) :: mpireqs1,mpireqs2,mpireqs3,mpireqs4,mpireqs5,mpireqs6,mpireqs7
  Integer(kind=4) :: mpireqr1,mpireqr2,mpireqr3,mpireqr4,mpireqr5,mpireqr6,mpireqr7
  Integer(kind=4) :: procID !< process id in the global communicator
  Integer(kind=4) :: procNB !< process number in the global communicator
  Integer(kind=4) :: dims(3) !< dimensions of the processes grid
  Integer(kind=4) :: CubeCoord(3)  !< coordinates of the process in the MpiCube communicator
  Integer(kind=4), dimension(6) :: neighbours  !< array containing the id of the neighboor of the local process in the MpiCube communicator
  Integer(kind=4) :: commcolorWrite
  Integer(kind=4) :: scprocIDWrite
  Integer(kind=4) :: CommGatherWrite
  
Contains

  !=======================================================================
  !> This subroutine create a processes topology that maps the lightcone.
  Subroutine settopology()

    Implicit None

    Integer :: in

    Do in = 1, 6
       If(neighbours(in) == -1) Then
          neighbours(in) = MPI_PROC_NULL
       End If
    End Do

    commcolorWrite = CubeCoord(1) - 1
    Call Mpi_Comm_Split(Mpi_Comm_World, commcolorWrite, procID, CommGatherWrite, mpierr)       
    Call Mpi_Comm_Rank(CommGatherWrite, scprocIDWrite, mpierr)

#ifdef DEBUG
    Print *,'Proc',procID, ': cubecoord=',CubeCoord, ' and ID in subcomm=',scprocIDWrite
#endif
    

  End Subroutine settopology


  !=======================================================================
  !> The subroutine is used to abort the execution if something wrong is happening.
  Subroutine EmergencyStop(message,errcode)

    Implicit none

    Character(len=*), Intent(in) :: message
    Integer, Intent(in) :: errcode

    Write(*,1000) message,procID

1000 Format('*** ',A,' on process ',I5.5,'. ***')

    If(procID==0) Close(50)
    Close(54)

    Call Mpi_Abort(Mpi_Comm_World,errcode,mpierr)
    Stop

  End Subroutine EmergencyStop


End Module modmpicom
