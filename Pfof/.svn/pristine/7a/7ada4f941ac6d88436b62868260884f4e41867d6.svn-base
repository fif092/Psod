!=======================================================================
! Author: Fabrice Roy (LUTH/CNRS/Observatoire de Paris)
! Fabrice.Roy@obspm.fr
!=======================================================================

Module modvariable
  Use mpi

#ifdef LONGINT
  Integer, parameter :: PRI = 8
  Integer, parameter :: MPI_PRI = Mpi_Integer8
#else
  Integer, parameter :: PRI = 4
  Integer, parameter :: MPI_PRI = Mpi_Integer
#endif
  

  Integer(kind=PRI), dimension(:), allocatable :: id        ! id of the particles in the halo just read
  Real(kind=4), dimension(:,:), allocatable :: pos  ! position of the particles in the halo just read
  Real(kind=4), dimension(:,:), allocatable :: vel  ! velocities of the particles in the halo just read
  Integer(kind=PRI), dimension(:), allocatable :: ramsespartid

  Integer :: procID
  Integer :: procNB
  Integer :: mpierr

  Integer(kind=4), dimension(:), allocatable :: npartcubeloc, npartcube
  Integer(kind=4), dimension(:), allocatable :: idcube

  Integer(kind=4) :: ncx
  Integer(kind=4) :: ncy
  Integer(kind=4) :: ncz
  Integer(kind=4) :: ncube

  Real(kind=8) :: hy
  Real(kind=8) :: hz

  !! MPI Variables
  Integer :: req_sumnpc
  Integer(kind=4), dimension(Mpi_Status_Size) :: mpistat

End Module modvariable
