!> @file
!! This file contains the subroutines and the variables used for MPI communications.

!> This module contains the subroutines and the variables used for MPI communications.
!> Authors: F. Roy, V. Bouillot
Module modmpicom
  
  Use mpi

  Integer(kind=4) :: mpierr !< MPI error code
  Integer(kind=4) :: mpireqs1,mpireqs2,mpireqs3,mpireqs4,mpireqs5,mpireqs6, mpireqs7
  Integer(kind=4) :: mpireqr1,mpireqr2,mpireqr3,mpireqr4,mpireqr5,mpireqr6,mpireqr7
  Integer(kind=4) :: procID !< process id in the global communicator
  Integer(kind=4) :: procNB !< process number in the global communicator
  Integer(kind=4) :: dims(3) !< dimensions of the processes grid
  Integer(kind=4) :: MpiCube !< cartesian global communicator
  Integer(kind=4) :: MpiSubCubeWrite !< local communicator used for hdf5 parallel output
  Integer(kind=4) :: MpiSubCubeRead  !< local communicator used for hdf5 parallel input
  Integer(kind=4) :: CubeCoord(3)  !< coordinates of the process in the MpiCube communicator
  Integer(kind=4), dimension(6) :: neighbours  !< array containing the id of the neighboor of the local process in the MpiCube communicator
  Integer(kind=4) :: commcolorWrite !< parameter used to create the MpiSubCubeWrite communicator
  Integer(kind=4) :: commcolorRead  !< parameter used to create the MpiSubCubeRead communicator
  Integer(kind=4) :: scprocIDWrite !< process id in the MpiSubCubeWrite communicator
  Integer(kind=4) :: scprocIDRead !< process id in the MpiSubCubeRead communicator
  Integer(kind=4) :: scprocNBWrite !< process number in the MpiSubCubeWrite communicator
  Integer(kind=4) :: scprocNBRead !< process number in the MpiSubCubeRead communicator
  Logical(kind=4) :: periods(3) !< periodicity of the global cartesian processes grid



Contains

!!$  !=======================================================================
!!$  !> This subroutine broadcasts parameters from the process 0 to every other processes.
!!$  Subroutine bcastparam()
!!$
!!$    Use modparameters
!!$    Implicit none
!!$
!!$    Integer(kind=4) :: buffersize, nbbytes, b_pos
!!$    Character, dimension(:),allocatable :: buffer
!!$
!!$    ! Memory allocation for the input parameters broadcast buffer
!!$    ! length of the strings:  simulationName, pathinput, namepart, nameinfo, code_index
!!$    ! 6x Integer(kind=4)   grpsize, Mmin, Mmax, gatherwrite, gatherread, snapshot
!!$    ! 1x Real(kind=4)      perco
!!$    ! 11x Logical(kind=4)  star, metal, potential, force, outcube, dofof, readfromcube,
!!$    !                      dotimings, sortcube, 
!!$    !                      doUnbinding, doSubHalo
!!$
!!$    ! buffer size in bits
!!$    buffersize = len(simulation_name) + len(input_path) + &
!!$         len(part_inputfile) + len(info_inputfile) + len(code_index)
!!$
!!$    ! size of integer(kind=4) + size of logical(kind=4) (same size as integer(kind=4))
!!$    Call Mpi_Sizeof(grpsize, nbbytes, mpierr)
!!$    buffersize = buffersize + 17*nbbytes
!!$
!!$    ! size of real(kind=SP) 
!!$    Call Mpi_Sizeof(percolation_length, nbbytes, mpierr)
!!$    buffersize = buffersize + 1*nbbytes
!!$    
!!$    allocate (buffer(0:buffersize-1))
!!$    b_pos = 0
!!$
!!$    If(procID==0) Then
!!$       ! Pack input parameters in the buffer
!!$       Call Mpi_Pack(simulation_name, len(simulation_name), Mpi_Character, buffer, buffersize,&
!!$            b_pos, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Pack(snapshot, 1, Mpi_Integer, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Pack(code_index, len(code_index), Mpi_Character, buffer, buffersize, &
!!$            b_pos, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Pack(input_path, len(input_path), Mpi_Character, buffer, buffersize, &
!!$            b_pos, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Pack(part_inputfile, len(part_inputfile), Mpi_Character, buffer, buffersize, &
!!$            b_pos, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Pack(info_inputfile, len(info_inputfile), Mpi_Character, buffer, buffersize, &
!!$            b_pos, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Pack(percolation_length, 1, Mpi_Real, buffer, buffersize, &
!!$            b_pos, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Pack(grpsize, 1, Mpi_Integer, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Pack(mmin, 1, Mpi_Integer, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Pack(mmax, 1, Mpi_Integer, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Pack(gatherwrite_factor, 1, Mpi_Integer, buffer, buffersize, &
!!$            b_pos, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Pack(gatherread_factor, 1, Mpi_Integer, buffer, buffersize, &
!!$            b_pos, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Pack(do_skip_star, 1, Mpi_Logical, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Pack(do_skip_metal, 1, Mpi_Logical, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Pack(do_read_potential, 1, Mpi_Logical, buffer, buffersize, &
!!$            b_pos, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Pack(do_read_force, 1, Mpi_Logical, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Pack(do_write_cube, 1, Mpi_Logical, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Pack(do_sort_cube, 1, Mpi_Logical, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Pack(do_fof, 1, Mpi_Logical, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Pack(do_read_from_cube, 1, Mpi_Logical, buffer, buffersize, &
!!$            b_pos, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Pack(do_timings, 1, Mpi_Logical, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Pack(do_unbinding, 1, Mpi_Logical, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Pack(do_subhalo, 1, Mpi_Logical, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
!!$
!!$    End If
!!$
!!$    ! Broadcast of input parameters
!!$    Call Mpi_Bcast(buffer,buffersize,Mpi_Packed,0,Mpi_Comm_World,mpierr)
!!$
!!$    ! Processes with ID != 0 unpack the input parameters
!!$    If(procID /= 0) Then
!!$       Call Mpi_Unpack(buffer, buffersize, b_pos, simulation_name, len(simulation_name), &
!!$            Mpi_Character, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Unpack(buffer, buffersize, b_pos, snapshot, 1 , Mpi_Integer, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Unpack(buffer, buffersize, b_pos, code_index, len(code_index), &
!!$            Mpi_Character, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Unpack(buffer, buffersize, b_pos, input_path, len(input_path), &
!!$            Mpi_Character, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Unpack(buffer, buffersize, b_pos, part_inputfile, len(part_inputfile), &
!!$            Mpi_Character, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Unpack(buffer, buffersize, b_pos, info_inputfile, len(info_inputfile), &
!!$            Mpi_Character, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Unpack(buffer, buffersize, b_pos, percolation_length, 1, &
!!$            Mpi_Real, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Unpack(buffer, buffersize, b_pos, grpsize, 1, Mpi_Integer, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Unpack(buffer, buffersize, b_pos, mmin, 1, Mpi_Integer, Mpi_Comm_World, mpierr)
!!$       Call Mpi_UnPack(buffer, buffersize, b_pos, mmax, 1, Mpi_Integer, Mpi_Comm_World, mpierr)
!!$       Call Mpi_UnPack(buffer, buffersize, b_pos, gatherwrite_factor, 1, &
!!$            Mpi_Integer, Mpi_Comm_World, mpierr)
!!$       Call Mpi_UnPack(buffer, buffersize, b_pos, gatherread_factor, 1, &
!!$            Mpi_Integer, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Unpack(buffer, buffersize, b_pos, do_skip_star, 1, &
!!$            Mpi_Logical, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Unpack(buffer, buffersize, b_pos, do_skip_metal, 1, &
!!$            Mpi_Logical, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Unpack(buffer, buffersize, b_pos, do_read_potential, 1, &
!!$            Mpi_Logical, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Unpack(buffer, buffersize, b_pos, do_read_force, 1, &
!!$            Mpi_Logical, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Unpack(buffer, buffersize, b_pos, do_write_cube, 1, &
!!$            Mpi_Logical, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Unpack(buffer, buffersize, b_pos, do_sort_cube, 1, &
!!$            Mpi_Logical, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Unpack(buffer, buffersize, b_pos, do_fof, 1, Mpi_Logical, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Unpack(buffer, buffersize, b_pos, do_read_from_cube, 1, &
!!$            Mpi_Logical, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Unpack(buffer, buffersize, b_pos, do_timings, 1, Mpi_Logical, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Unpack(buffer, buffersize, b_pos, do_unbinding, 1, &
!!$            Mpi_Logical, Mpi_Comm_World, mpierr)
!!$       Call Mpi_Unpack(buffer, buffersize, b_pos, do_subhalo, 1, &
!!$            Mpi_Logical, Mpi_Comm_World, mpierr)
!!$    End If
!!$
!!$    ! broadcast buffer deallocated
!!$    Deallocate(buffer)
!!$
!!$  End Subroutine bcastparam

  !=======================================================================
  !> This subroutine creates the global cartesian grid of processes.
  !! It also creates the communicators used for the parallel hdf5 i/o if they are needed.
  Subroutine setcommunicators()

    Use modparameters
    Implicit None

    ! Creation of a grid of processes.
    ! The grid is a cube, there are procNB**(1/3) processes in each dimension.
    ! The grid is periodic in each dimension.
    dims = int(procNB**(1./3.))
    periods = .true. 
    
    Call Mpi_Cart_Create(Mpi_Comm_World,3,dims,periods,.true.,MPICube,mpierr)
    
    ! Process ID and number 
    Call Mpi_Comm_Rank  (MPICube, procID, mpierr)
    Call Mpi_Cart_Coords(MPICube, procID, 3, CubeCoord, mpierr)
    
    ! Definition des neighbourss pour les 6 directions avec des conditions au bord periodiques
    ! Intialization of the 'neighboor processes id' 
    ! Neighboor : 1 = backward
    !             2 = forward
    !             3 = left
    !             4 = right
    !             5 = down
    !             6 = up
    
    Call Mpi_Cart_Shift(MpiCube,0,1,neighbours(1),neighbours(2),mpierr)
    Call Mpi_Cart_Shift(MpiCube,1,1,neighbours(3),neighbours(4),mpierr)
    Call Mpi_Cart_Shift(MpiCube,2,1,neighbours(5),neighbours(6),mpierr)
    

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

    End If

  End Subroutine setcommunicators

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

    Call Mpi_Abort(MPICube,errcode,mpierr)
    Stop

  End Subroutine EmergencyStop

End Module modmpicom
