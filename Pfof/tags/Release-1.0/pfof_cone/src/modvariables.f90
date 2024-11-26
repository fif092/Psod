Module modvariables
  Use mpi, only : Mpi_Integer, Mpi_Integer8

#ifdef LONGREAL
  Integer, parameter :: PR=8   !< Precision for real arrays (position/velocities)
#else
  Integer, parameter :: PR=4   !< Precision for real arrays (position/velocities)
#endif

#ifdef LONGINT
  Integer, parameter :: PRI = 8   !< Precision for integer arrays (id)
  Integer, parameter :: MPI_PRI = Mpi_Integer8  !< MPI precision for integer arrays
#else
  Integer, parameter :: PRI = 4   !< Precision for integer arrays (id)
  Integer, parameter :: MPI_PRI = Mpi_Integer   !< MPI precision for integer arrays
#endif



  Integer(kind=4) :: shell_nb
  Integer(kind=4), dimension(:,:), allocatable :: indexcube
!  Integer(kind=4), dimension(:), allocatable :: npartcube
  Integer(kind=4), dimension(:), allocatable :: ictable
  Integer(kind=4) :: shellcubes_size
  Integer(kind=4), dimension(6) :: neighbours
  Real(kind=8), dimension(6) :: boundaries 
  Integer(kind=4), dimension(:), allocatable :: cubepershell
  Integer(kind=4) :: mynpart
  Integer(kind=PRI) :: npart

  Real(kind=PR), dimension(:,:), allocatable :: pos
  Real(kind=PR), dimension(:,:), allocatable :: vel
  Integer(kind=PRI), dimension(:), allocatable :: id

  Real(kind=PR), dimension(3) :: partmin, partmax

Contains

  Subroutine deallocall()

    Implicit None

!    If(Allocated(npartcube)) Deallocate(npartcube)
    If(Allocated(indexcube)) Deallocate(indexcube)
    If(Allocated(ictable)) Deallocate(ictable)
    If(Allocated(cubepershell)) Deallocate(cubepershell)
    If(Allocated(pos)) Deallocate(pos)
    If(Allocated(vel)) Deallocate(vel)
    If(Allocated(id)) Deallocate(id)

  End Subroutine deallocall

End Module modvariables
