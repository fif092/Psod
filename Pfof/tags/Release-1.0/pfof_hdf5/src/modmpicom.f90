!> @file
!! This file contains the subroutines and the variables used for MPI communications.

!> This module contains the subroutines and the variables used for MPI communications.
Module modmpicom
  
  Use mpi
  Use modconst

  Integer(kind=4) :: mpierr !< MPI error code
  Integer(kind=4) :: mpireqs1,mpireqs2,mpireqs3,mpireqs4,mpireqr1,mpireqr2,mpireqr3,mpireqr4
  Integer(kind=4) :: procID !< process id in the global communicator
  Integer(kind=4) :: procNB !< process number in the global communicator
  Integer(kind=4) :: dims(3) !< dimensions of the processes grid
  Integer(kind=4) :: MpiCube !< cartesian global communicator
  Integer(kind=4) :: MpiSubCubeWrite !< local communicator used for hdf5 parallel output
  Integer(kind=4) :: MpiSubCubeRead  !< local communicator used for hdf5 parallel input
  Integer(kind=4) :: CubeCoord(3)  !< coordinates of the process in the MpiCube communicator
  Integer(kind=4) :: voisin(6)  !< array containing the id of the neighboor of the local process in the MpiCube communicator
  Integer(kind=4) :: commcolorWrite !< parameter used to create the MpiSubCubeWrite communicator
  Integer(kind=4) :: commcolorRead  !< parameter used to create the MpiSubCubeRead communicator
  Integer(kind=4) :: scprocIDWrite !< process id in the MpiSubCubeWrite communicator
  Integer(kind=4) :: scprocIDRead !< process id in the MpiSubCubeRead communicator
  Integer(kind=4) :: scprocNBWrite !< process number in the MpiSubCubeWrite communicator
  Integer(kind=4) :: scprocNBRead !< process number in the MpiSubCubeRead communicator
  Logical :: periods(3) !< periodicity of the global cartesian processes grid



Contains

  !=======================================================================
  !> This subroutine broadcasts parameters from the process 0 to every other processes.
  Subroutine bcastparam()

    Use modparam
    Implicit none

    Integer(kind=4) :: buffersize, nbbytes, b_pos
    Character, dimension(:),allocatable :: buffer

    ! Memory allocation for the input parameters broadcast buffer
    ! 80*4               : 4x Character(len=80) root, pathinput, namepart, nameinfo 
    ! 3                  : 1x Character(len=3)  code_index
    ! 3*bit_size(Mmin)/8 : 4x Integer(kind=4)   grpsize, Mmin, Mmax, gatherwrite, gatherread
    ! kind(perco)        : 1x Real(kind=4)      perco
    ! 7*4                : 7x Logical           star, metal, outcube, dofof, readfromcube, dotimings, usehdf5

    ! buffer size in bits
    buffersize = len(root) + len(pathinput) + len(namepart) + len(nameinfo) + len(code_index)

    ! size of integer(kind=4)
    Call Mpi_Sizeof(grpsize, nbbytes, mpierr)
    buffersize = buffersize + 5*nbbytes

    ! size of real(kind=SP)
    Call Mpi_Sizeof(perco, nbbytes, mpierr)
    buffersize = buffersize + nbbytes

    ! size of logical
    Call Mpi_Sizeof(star, nbbytes, mpierr)
    buffersize = buffersize + 7*nbbytes

    allocate (buffer(0:buffersize-1))
    b_pos = 0

    If(procID==0) Then
       ! Pack input parameters in the buffer
       Call Mpi_Pack(        root, len(root), Mpi_Character, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(  code_index, len(code_index), Mpi_Character, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(   pathinput, len(pathinput), Mpi_Character, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(    namepart, len(namepart), Mpi_Character, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(    nameinfo, len(nameinfo), Mpi_Character, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(       perco,  1,      Mpi_Real, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(     grpsize,  1,   Mpi_Integer, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(        Mmin,  1,   Mpi_Integer, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(        Mmax,  1,   Mpi_Integer, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack( gatherwrite,  1,   Mpi_Integer, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(  gatherread,  1,   Mpi_Integer, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(        star,  1,   Mpi_Logical, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(       metal,  1,   Mpi_Logical, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(     outcube,  1,   Mpi_Logical, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(       dofof,  1,   Mpi_Logical, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(readfromcube,  1,   Mpi_Logical, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(   dotimings,  1,   Mpi_Logical, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(     usehdf5,  1,   Mpi_Logical, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)

    End If

    ! Broadcast of input parameters
    Call Mpi_Bcast(buffer,buffersize,Mpi_Packed,0,Mpi_Comm_World,mpierr)

    ! Processes with ID != 0 unpack the input parameters
    If(procID /= 0) Then
       Call Mpi_Unpack(buffer, buffersize, b_pos,         root, len(root), Mpi_Character, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,   code_index, len(code_index), Mpi_Character, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,    pathinput, len(pathinput), Mpi_Character, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,     namepart, len(namepart), Mpi_Character, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,     nameinfo, len(nameinfo), Mpi_Character, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,        perco,  1,      Mpi_Real, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,      grpsize,  1,   Mpi_Integer, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,         Mmin,  1,   Mpi_Integer, Mpi_Comm_World, mpierr)
       Call Mpi_UnPack(buffer, buffersize, b_pos,         Mmax,  1,   Mpi_Integer, Mpi_Comm_World, mpierr)
       Call Mpi_UnPack(buffer, buffersize, b_pos,  gatherwrite,  1,   Mpi_Integer, Mpi_Comm_World, mpierr)
       Call Mpi_UnPack(buffer, buffersize, b_pos,   gatherread,  1,   Mpi_Integer, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,         star,  1,   Mpi_Logical, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,        metal,  1,   Mpi_Logical, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,      outcube,  1,   Mpi_Logical, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,        dofof,  1,   Mpi_Logical, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos, readfromcube,  1,   Mpi_Logical, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,    dotimings,  1,   Mpi_Logical, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,      usehdf5,  1,   Mpi_Logical, Mpi_Comm_World, mpierr)
    End If

    ! broadcast buffer deallocated
    Deallocate(buffer)

  End Subroutine bcastparam

  !=======================================================================
  !> This subroutine creates the global cartesian grid of processes.
  !! It also creates the communicators used for the parallel hdf5 i/o if they are needed.
  Subroutine setcommunicators()

    Use modparam
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
    
    ! Definition des voisins pour les 6 directions avec des conditions au bord periodiques
    ! Intialization of the 'neighboor processes id' 
    ! Neighboor : 1 = backward
    !             2 = forward
    !             3 = left
    !             4 = right
    !             5 = down
    !             6 = up
    
    Call Mpi_Cart_Shift(MpiCube,0,1,voisin(1),voisin(2),mpierr)
    Call Mpi_Cart_Shift(MpiCube,1,1,voisin(3),voisin(4),mpierr)
    Call Mpi_Cart_Shift(MpiCube,2,1,voisin(5),voisin(6),mpierr)
    

    If(gatherwrite > 1) Then
       If(procID==0) Then
          Print *, '*** Creating specific communicator for gathered cube output ***'  
       End If
       
       commcolorWrite = CubeCoord(1)/gatherwrite +  (CubeCoord(2)/gatherwrite) * (dims(1)/gatherwrite) &
            + (CubeCoord(3)/gatherwrite) * (dims(1)/gatherwrite)*(dims(1)/gatherwrite)
       
       Call Mpi_Comm_Split(MpiCube, commcolorWrite, procID, MpiSubCubeWrite, mpierr)
       
       Call Mpi_Comm_Rank(MpiSubCubeWrite, scprocIDWrite, mpierr)
       Call Mpi_Comm_Size(MpiSubCubeWrite, scprocNBWrite, mpierr)
    End If

    If(gatherread > 1) Then

       If(procID==0) Then
          Print *, '*** Creating specific communicator for gathered cube input ***'  
       End If

       commcolorRead = CubeCoord(1)/gatherread +  (CubeCoord(2)/gatherread) * (dims(1)/gatherread) &
            + (CubeCoord(3)/gatherread) * (dims(1)/gatherread)*(dims(1)/gatherread)
       
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
