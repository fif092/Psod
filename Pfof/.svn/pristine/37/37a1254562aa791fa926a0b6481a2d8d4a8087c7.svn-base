!=======================================================================
! Author: Fabrice Roy (LUTH/CNRS/Observatoire de Paris)
! Fabrice.Roy@obspm.fr
!=======================================================================
Module modvariable

#ifdef WITHMPI
  Use mpi
#endif
  
  ! If you do not know the type of integer used for the ID, use h5dump on the first halo file.
  ! h5dump -d haloID test_halo_00000.h5
  ! This will print the dataset haloID containing the ID of each halo in the file. 
  ! You should see something like
  ! HDF5 "test_halo_00000.h5" {
  ! DATASET "haloID" {
  !    DATATYPE  H5T_STD_I64LE
  ! The datatype should contain the number of bits used, i.e. 32 for kind=4 or 64 for kind=8.
#ifdef LONGINT
  Integer, parameter :: PRI = 8
#ifdef WITHMPI
  Integer, parameter :: MPI_PRI = Mpi_Integer8
#endif
#else
  Integer, parameter :: PRI = 4
#ifdef WITHMPI
  Integer, parameter :: MPI_PRI = Mpi_Integer
#endif
#endif

  !! Variables read from the HDF5 file
  Integer(kind=PRI), dimension(:), allocatable :: id        ! id of the particles in the halo just read
  Real(kind=4), dimension(:,:), allocatable :: pos  ! position of the particles in the halo just read
  Real(kind=4), dimension(:,:), allocatable :: vel  ! velocities of the particles in the halo just read
  Real(kind=4), dimension(:), allocatable :: mp 


  Integer(kind=PRI) :: currenthaloID ! haloID of the halo analyzed in your function
  Integer(kind=4) :: currenthalo     ! index of the halo which varies from 1 to N where N is the total number of haloes you are analyzing
  Integer(kind=4) :: nbhaloanalyzed  ! total number of haloes that you want to analyze 

#ifdef WITHMPI
  Integer :: procID
  Integer :: procNB
  Integer :: mpierr
#endif

End Module modvariable
