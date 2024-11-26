!=======================================================================
! Author: Fabrice Roy (LUTH/CNRS/Observatoire de Paris)
! Fabrice.Roy@obspm.fr
!=======================================================================

!=======================================================================
! Generic program reading a set of hdf5 halo files produced by pFoF
! and analyzing them
! You can analyze either a subset of haloes by providing thier ID (analyzelist)
! or each halo in a subset of files by providing the name of the files (analyzefile)
!=======================================================================
Program haloanalyzer

  Use modanalyzeengines
  Use modfunctions
  Implicit none

  Type(ListOfVar) :: var         ! list of the datasets needed for your analyzis
  Character(len=400) :: filename  ! name of the parameter file: and example of the 2 different kind of files can be found with the source files

  ! You can add your own variables here
  ! for example loop indices 
  Integer :: i
  ! or variables required to open/write an output file
  Integer(hid_t) :: file_id
  Character(len=16) :: name
  Character(len=50) :: origin
  Character(len=15) :: codeversion
  Integer :: ierr


#ifdef WITHMPI
  Call Mpi_Init(mpierr)
#endif

  Call hdf5_init()
  

  ! You can define the datasets you need to read to perform your analyzis
  var%pos = .true.    ! default is .true.
  var%vel = .false.   ! default is .true.
  var%id = .false.    ! default is .true.

  ! The first argument is the name of the file containing the list of the ID of the haloes that you want to analyze
  ! The second argument is the name of your analyzis subroutine (defined in modfunctions.f90) 
  ! compos is an example that you can find in modfunctions.f90
  ! The third argument (optional, used in the following example) is the list of variables (datasets) that you need to read from the HDF5 files
  filename = 'halolist.nml'
  Call analyzelist(filename,compos,var)

  ! Output of the results of your analyzes here
  Print *,'Position of the center of mass (list) :'
  Do i=1,size(cmpos,2)
     Print *,'POS: ',cmpos(:,i)
  End Do

  ! Deallocate output arrays before you perform other analyzes
  If(Allocated(cmpos)) Deallocate(cmpos)

  
  ! The first argument is the name of the file containing the list of the files that you want to analyze
  ! The second argument is the name of your analyzis subroutine (defined in modfunctions.f90) 
  ! comvel is an example that you can find in modfunctions.f90
  ! The third argument (optional, not used in the following example) is the list of variables (datasets) that you need to read from the HDF5 files
  filename = 'filelist.nml'
  Call analyzefile(filename, comvel)

  ! Output of the results of your analyzes here
  Print *,'Velocity of the center of mass (all) :'
  Do i=1,size(cmvel,2)
     Print *, 'VEL: ',cmvel(:,i)
  End Do


  ! Example of parallel hdf5 output
  ! Use parallel HDF5 only if you compile with -DWITHHDF5
  name='cmvel'
  filename = 'cmvel.h5'
  Open(Unit=10, File='haloanalyzer.version', status='old', Iostat=ierr)
  If(ierr/=0) Then
     Print *, 'version of haloanalyzer not defined: file not found'
     codeversion = 'undefined'
  Else
     Read(10,*,Iostat=ierr) codeversion
     If(ierr/=0) Then
        Print *, 'version of haloanalyzer not defined: file empty'
        codeversion='undefined'
     End If
     Close(10)
  End If
  origin = 'Created with haloanalyzer version '//codeversion
  Call hdf5_create_mpi_file(filename, Mpi_Comm_World, file_id, origin)
  Call hdf5_write_mpi_data(file_id, name, 3, size(cmvel,2),cmvel,Mpi_Comm_World) 
  Call hdf5_close_mpi_file(file_id)

  If(Allocated(cmvel)) Deallocate(cmvel)

  Call hdf5_finalize()

#ifdef WITHMPI
  Call Mpi_Finalize(mpierr)
#endif

End Program haloanalyzer
