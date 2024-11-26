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
  

  Real(kind=4), dimension(:,:), allocatable :: pos  ! position of the cells in the halo just read
  Real(kind=4), dimension(:,:), allocatable :: acc  ! acceleration in the cells in the halo just read 
  Real(kind=4), dimension(:), allocatable :: phi  ! potential in the cells in the halo just read
  Real(kind=4), dimension(:), allocatable :: rho  ! density in the cells in the halo just read
  Integer(kind=4), dimension(:), allocatable :: refined ! is the cell refined or not

  Character(len=200), dimension(:), allocatable :: filelist    ! list of the filenames that we must read

  Integer :: procID
  Integer :: procNB
  Integer :: mpierr

  Integer(kind=4), dimension(:,:), allocatable :: ncellcubeloc, ncellcube
  Integer(kind=4), dimension(:), allocatable :: idcube

  Integer(kind=4) :: ncx
  Integer(kind=4) :: ncy
  Integer(kind=4) :: ncz
  Integer(kind=4) :: ncube

  Integer(kind=4), dimension(:), allocatable :: ncellperlevel

  Real(kind=8) :: hy
  Real(kind=8) :: hz

  !! MPI Variables
  Integer :: req_sumnpc
  Integer(kind=4), dimension(Mpi_Status_Size) :: mpistat

End Module modvariable
