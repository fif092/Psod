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
  Integer(kind=4) :: mpireqs1,mpireqs2,mpireqs3,mpireqs4,mpireqs5,mpireqs6
  Integer(kind=4) :: mpireqr1,mpireqr2,mpireqr3,mpireqr4,mpireqr5,mpireqr6
  Integer(kind=4) :: procID !< process id in the global communicator
  Integer(kind=4) :: procNB !< process number in the global communicator
  Integer(kind=4) :: dims(3) !< dimensions of the processes grid
  Integer(kind=4) :: CubeCoord(3)  !< coordinates of the process in the MpiCube communicator
  Integer(kind=4), dimension(6) :: neighbours  !< array containing the id of the neighboor of the local process in the MpiCube communicator

Contains


  !=======================================================================
  !> This subroutine broadcasts parameters from the process 0 to every other processes.
  Subroutine bcastparam()

    Use modparameters
    Implicit none

    Integer(kind=4) :: buffersize, nbbytes, b_pos
    Character, dimension(:),allocatable :: buffer

    ! Memory allocation for the input parameters broadcast buffer
    ! strings according to their length
    ! 6x Integer(kind=4)    outputNB, shell_fid, shell_lid, Mmin, Mmax, nres
    ! 1x Real(kind=4)       perco
    ! 3x Logical            doUnbinding, doSubHalo, dotimings

    ! buffer size in bits
    buffersize = len(shell_dir) + len(shell_filename) + len(output_root) + len(origin) + len(simulationName)
    
    ! size of integer(kind=4) + logical (same size)
    Call Mpi_Sizeof(shell_fid, nbbytes, mpierr)
    buffersize = buffersize + 9*nbbytes

    ! size of real(kind=SP)
    Call Mpi_Sizeof(perco, nbbytes, mpierr)
    buffersize = buffersize + nbbytes

    ! size of logical
!!$    Call Mpi_Sizeof(dotimings, nbbytes, mpierr)
!!$    buffersize = buffersize + 3*nbbytes

    allocate (buffer(0:buffersize-1))
    b_pos = 0
    
    If(procID==0) Then
       ! Pack input parameters in the buffer
       Call Mpi_Pack(simulationName, len(simulationName), Mpi_Character, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(     shell_dir,      len(shell_dir), Mpi_Character, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(shell_filename, len(shell_filename), Mpi_Character, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(   output_root,    len(output_root), Mpi_Character, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(        origin,         len(origin), Mpi_Character, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(      outputNB,                   1,   Mpi_Integer, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(     shell_fid,                   1,   Mpi_Integer, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(     shell_lid,                   1,   Mpi_Integer, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(          Mmin,                   1,   Mpi_Integer, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(          Mmax,                   1,   Mpi_Integer, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(         perco,                   1,      Mpi_Real, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(          nres,                   1,   Mpi_Integer, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(   doUnbinding,                   1,   Mpi_Logical, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(     doSubHalo,                   1,   Mpi_Logical, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(     dotimings,                   1,   Mpi_Logical, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
    End If

    ! Broadcast of input parameters
    Call Mpi_Bcast(buffer,buffersize,Mpi_Packed,0,Mpi_Comm_World,mpierr)

    ! Processes with ID != 0 unpack the input parameters
    If(procID /= 0) Then
       Call Mpi_Unpack(buffer, buffersize, b_pos, simulationName, len(simulationName), Mpi_Character, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,      shell_dir,      len(shell_dir), Mpi_Character, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos, shell_filename, len(shell_filename), Mpi_Character, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,    output_root,    len(output_root), Mpi_Character, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,         origin,         len(origin), Mpi_Character, Mpi_Comm_World, mpierr)
       Call Mpi_UnPack(buffer, buffersize, b_pos,       outputNB,                   1,   Mpi_Integer, Mpi_Comm_World, mpierr)
       Call Mpi_UnPack(buffer, buffersize, b_pos,      shell_fid,                   1,   Mpi_Integer, Mpi_Comm_World, mpierr)
       Call Mpi_UnPack(buffer, buffersize, b_pos,      shell_lid,                   1,   Mpi_Integer, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,           Mmin,                   1,   Mpi_Integer, Mpi_Comm_World, mpierr)
       Call Mpi_UnPack(buffer, buffersize, b_pos,           Mmax,                   1,   Mpi_Integer, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,          perco,                   1,      Mpi_Real, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,           nres,                   1,   Mpi_Integer, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,    doUnbinding,                   1,   Mpi_Logical, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,      doSubHalo,                   1,   Mpi_Logical, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,      dotimings,                   1,   Mpi_Logical, Mpi_Comm_World, mpierr)
    End If

    ! broadcast buffer deallocated
    Deallocate(buffer)

  End Subroutine bcastparam


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
