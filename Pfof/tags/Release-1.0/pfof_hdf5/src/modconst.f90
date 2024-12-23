!> @file
!! This file contains some constant parameters definition

!> Constant parameters definition
Module modconst

!  Include 'mpif.h'
  Use mpi

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

  Integer, parameter :: DP = kind(1.d0)   !< Definition of double precision
  Integer, parameter :: SP = kind(1.e0)   !< Definition of simple precision

  ! Output Units 
  Integer, parameter :: Ulog=50   !< I/O unit for text log file
  Integer, parameter :: Ucub=51   !< I/O unit for binary cube file
  Integer, parameter :: Umas=52   !< I/O unit for binary mass file
  Integer, parameter :: Ustr=53   !< I/O unit for binary haloes file
  Integer, parameter :: Uopa=54   !< I/O unit for text input parameter


End Module modconst
