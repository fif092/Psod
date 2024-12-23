!=======================================================================
! Author: Fabrice Roy (LUTH/CNRS/Observatoire de Paris)
! Fabrice.Roy@obspm.fr
!=======================================================================

!> @file
!! This file contains wrappers for HDF5 serial and MPI I/O functions.

!> This module contains wrappers for HDF5 serial and MPI I/O functions.
Module modhdf5

  Use mpi
  Use hdf5
  Use iso_c_binding

  Implicit None

  Private
  Integer :: h5err !< hdf5 error code
  Integer :: mpierr
    
  Public :: hdf5_open_file, &
       hdf5_open_mpi_file, &
       hdf5_create_mpi_file, &
       hdf5_close_mpi_file, &
       hdf5_create_file, &
       hdf5_close_file, &
       hdf5_open_group, &
       hdf5_create_group, &
       hdf5_close_group, &
       hdf5_write_data, &
       hdf5_read_data, &
       hdf5_write_attr, &
       hdf5_read_attr, &
       hdf5_write_mpi_data, &
       hdf5_read_mpi_data, &
       hdf5_init, &
       hdf5_finalize

  Public :: hid_t, hsize_t
  
  !> Generic inteface used to write attributes in a hdf5 file.
  !> @param[in] id id of the file/group where the attribute will be written
  !> @param[in] name name of the attribute
  !> @param[in] n1 first dimension of the attribute array (optional)
  !> @param[in] n2 second dimension of the attribute array (optional)
  !> @param[in] data 0-D, 1-D or 2-D attribute of type Integer(4), Real(4), Real(8) or Character
  Interface hdf5_write_attr
     Module Procedure hdf5_write_int4_attr0D   ! (id, name, data)
     Module Procedure hdf5_write_int4_attr1D   ! (id, name, n1, data)
     Module Procedure hdf5_write_int4_attr2D   ! (id, name, n1, n2, data)
     Module Procedure hdf5_write_real4_attr0D  ! (id, name, data)
     Module Procedure hdf5_write_real4_attr1D  ! (id, name, n1, data)
     Module Procedure hdf5_write_real4_attr2D  ! (id, name, n1, n2, data)
     Module Procedure hdf5_write_real8_attr0D  ! (id, name, data)
     Module Procedure hdf5_write_real8_attr1D  ! (id, name, n1, data)
     Module Procedure hdf5_write_real8_attr2D  ! (id, name, n1, n2, data)
     Module Procedure hdf5_write_char_attr     ! (id, name, data)
  End Interface hdf5_write_attr

  !> Generic inteface used to read attributes from a hdf5 file.
  !> @param[in] id id of the file/group where the attribute will be read
  !> @param[in] name name of the attribute
  !> @param[in] n1 first dimension of the attribute array (optional)
  !> @param[in] n2 second dimension of the attribute array (optional)
  !> @param[in,out] data 0-D, 1-D or 2-D attribute of type Integer(4), Real(4), Real(8) or Character
  Interface hdf5_read_attr
     Module Procedure hdf5_read_int4_attr0D   ! (id, name, data)
     Module Procedure hdf5_read_int4_attr1D   ! (id, name, n1, data)
     Module Procedure hdf5_read_int4_attr2D   ! (id, name, n1, n2, data)
     Module Procedure hdf5_read_real4_attr0D  ! (id, name, data)
     Module Procedure hdf5_read_real4_attr1D  ! (id, name, n1, data)
     Module Procedure hdf5_read_real4_attr2D  ! (id, name, n1, n2, data)
     Module Procedure hdf5_read_real8_attr0D  ! (id, name, data)
     Module Procedure hdf5_read_real8_attr1D  ! (id, name, n1, data)
     Module Procedure hdf5_read_real8_attr2D  ! (id, name, n1, n2, data)
     Module Procedure hdf5_read_char_attr     ! (id, name, data)
  End Interface hdf5_read_attr
  
  !> Generic inteface used to write data in a hdf5 file.
  !> @param[in] id id of the file/group where the dataset will be written
  !> @param[in] name name of the dataset
  !> @param[in] n1 first dimension of the data array
  !> @param[in] n2 second dimension of the data array (optional)
  !> @param[in] data 1-D or 2-D data array of type Integer(4), Integer(8), Real(4) or Real(8)
  Interface hdf5_write_data
     Module Procedure hdf5_write_int4_1D       ! (id, name, n1, data)
     Module Procedure hdf5_write_int4_2D       ! (id, name, n1, n2, data)
     Module Procedure hdf5_write_int8_1D       ! (id, name, n1, data)
     Module Procedure hdf5_write_int8_2D       ! (id, name, n1, n2, data)
     Module Procedure hdf5_write_real4_1D      ! (id, name, n1, data)
     Module Procedure hdf5_write_real4_2D      ! (id, name, n1, n2, data)
     Module Procedure hdf5_write_real8_1D      ! (id, name, n1, data)
     Module Procedure hdf5_write_real8_2D      ! (id, name, n1, n2, data)
  End Interface hdf5_write_data

  !> Generic inteface used to read data from a hdf5 file.
  !> @param[in] id id of the file/group where the dataset will be read
  !> @param[in] name name of the dataset
  !> @param[in] n1 first dimension of the data array
  !> @param[in] n2 second dimension of the data array (optional)
  !> @param[in,out] data 1-D or 2-D data array of type Integer(4), Integer(8), Real(4) or Real(8)
  Interface hdf5_read_data
     Module Procedure hdf5_read_int4_1D        ! (id, name, n1, data)
     Module Procedure hdf5_read_int4_2D        ! (id, name, n1, n2, data)
     Module Procedure hdf5_read_int8_1D        ! (id, name, n1, data)
     Module Procedure hdf5_read_int8_2D        ! (id, name, n1, n2, data)
     Module Procedure hdf5_read_real4_1D       ! (id, name, n1, data)
     Module Procedure hdf5_read_real4_2D       ! (id, name, n1, n2, data)
     Module Procedure hdf5_read_real8_1D       ! (id, name, n1, data)
     Module Procedure hdf5_read_real8_2D       ! (id, name, n1, n2, data)
  End Interface hdf5_read_data

  !> Generic inteface used to write data in a parallel hdf5 file.
  !> @param[in] id id of the file/group where the dataset will be written
  !> @param[in] name name of the dataset
  !> @param[in] n1 first dimension of the data array
  !> @param[in] n2 second dimension of the data array (optional)
  !> @param[in] data 1-D or 2-D data array of type Integer(4), Integer(8), Real(4) or Real(8)
  !> @param[in] comm MPI communicator used for the output (the same communicator should have been used to create the file)
  Interface hdf5_write_mpi_data
     Module Procedure hdf5_write_mpi_int4_1D       ! (id, name, n1, data)
     Module Procedure hdf5_write_mpi_int4_2D       ! (id, name, n1, n2, data)
     Module Procedure hdf5_write_mpi_int8_1D       ! (id, name, n1, data)
     Module Procedure hdf5_write_mpi_int8_2D       ! (id, name, n1, n2, data)
     Module Procedure hdf5_write_mpi_real4_1D      ! (id, name, n1, data)
     Module Procedure hdf5_write_mpi_real4_2D      ! (id, name, n1, n2, data)
     Module Procedure hdf5_write_mpi_real8_1D      ! (id, name, n1, data)
     Module Procedure hdf5_write_mpi_real8_2D      ! (id, name, n1, n2, data)
  End Interface hdf5_write_mpi_data

  !> Generic inteface used to read data from a parallel hdf5 file.
  !> @param[in] id id of the file/group where the dataset will be read
  !> @param[in] name name of the dataset
  !> @param[in] n1 first dimension of the data array
  !> @param[in] n2 second dimension of the data array (optional)
  !> @param[in,out] data 1-D or 2-D data array of type Integer(4), Integer(8), Real(4) or Real(8)
  !> @param[in] comm MPI communicator used for the input (the same communicator should have been used to open the file)
  Interface hdf5_read_mpi_data
     Module Procedure hdf5_read_mpi_int4_1D
     Module Procedure hdf5_read_mpi_int4_2D
     Module Procedure hdf5_read_mpi_int8_1D
     Module Procedure hdf5_read_mpi_int8_2D
     Module Procedure hdf5_read_mpi_real4_1D
     Module Procedure hdf5_read_mpi_real4_2D
     Module Procedure hdf5_read_mpi_real8_1D
     Module Procedure hdf5_read_mpi_real8_2D
  End Interface hdf5_read_mpi_data
  !

Contains

  !===============================================!
  !         INIT/CLOSE ROUTINES
  !===============================================!

  !=======================================================================
  !> Open the hdf5 interface
  Subroutine hdf5_init()

    Implicit None

    Call h5open_f(h5err)

#ifdef DEBUG
    If(h5err /= 0) Then
       Print *,'HDF5 error while initializing the HDF5 interface'
       Call Mpi_Abort(Mpi_Comm_World, h5err, mpierr)
    End If
#endif

  End Subroutine hdf5_init


  !=======================================================================
  !> Close the hdf5 interface
  Subroutine hdf5_finalize()

    Implicit None

    Call h5close_f(h5err)

#ifdef DEBUG
    If(h5err /= 0) Then
       Print *,'HDF5 error while closing the HDF5 interface'
       Call Mpi_Abort(Mpi_Comm_World, h5err, mpierr)
    End If
#endif

  End Subroutine hdf5_finalize


  !===============================================!
  !         FILE MANAGEMENT ROUTINES
  !===============================================!

  !=======================================================================
  !> Open an existing hdf5 file
  Subroutine hdf5_open_file(filename, file_id)

    Implicit None

    Character(len=400), intent(in) :: filename           !< name of the file to open
    Integer(kind=hid_t), intent(out) :: file_id         !< hdf5 id of the opened file

    ! open the file 
    Call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, h5err)

#ifdef DEBUG
    If(h5err /= 0) Then
       Print *,'HDF5 error while opening file ',trim(filename)
       Call Mpi_Abort(Mpi_Comm_World, h5err, mpierr)
    End If
#endif

  End Subroutine hdf5_open_file
  

  !=======================================================================
  !> Open an existing 'parallel' hdf5 file shared by the communicator comm
  Subroutine hdf5_open_mpi_file(filename, comm, file_id)

    Implicit None

    Character(len=400), intent(in) :: filename           !< name of the parallel file to open
    Integer(kind=4), intent(in) :: comm                 !< MPI communicator used for the file access
    Integer(kind=hid_t), intent(out) :: file_id         !< hdf5 id of the opened file

    Integer(hid_t) :: plist_id                          ! Property list identifier for parallel file
    
    ! Properties of the file
    Call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, h5err)

#ifdef DEBUG
    If(h5err /= 0) Then
       Print *,'HDF5 error while creating the property list for parallel file ',trim(filename)
       Call Mpi_Abort(Mpi_Comm_World, h5err, mpierr)
    End If
#endif

    Call h5pset_fapl_mpio_f(plist_id, comm, Mpi_Info_Null, h5err)

#ifdef DEBUG
    If(h5err /= 0) Then
       Print *,'HDF5 error while setting the property list for parallel file ',trim(filename)
       Call Mpi_Abort(Mpi_Comm_World, h5err, mpierr)
    End If
#endif

    ! open the file 
    Call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, h5err, access_prp = plist_id)

#ifdef DEBUG
    If(h5err /= 0) Then
       Print *,'HDF5 error while opening the parallel file ',trim(filename)
       Call Mpi_Abort(Mpi_Comm_World, h5err, mpierr)
    End If
#endif

    Call h5pclose_f(plist_id, h5err)

  End Subroutine hdf5_open_mpi_file


  !=======================================================================
  !> Create and open a 'parallel' hdf5 file shared by the communicator comm
  Subroutine hdf5_create_mpi_file(filename, comm, file_id, origin)

    Implicit None

    Character(len=400), intent(in) :: filename           !< name of the parallel file to create
    Integer, intent(in) :: comm                         !< MPI communicator used for the file access
    Integer(hid_t), intent(out) :: file_id              !< hdf5 id of the created file
    Character(len=50), optional :: origin   !< origin of the file (code+version) : optional argument

    Character(len=50) :: deforigin                      ! default value for origin 
    Integer :: rank                                     ! Nb of dimensions for the attributes
    Integer(size_t) :: alen                             ! Attribute length in bytes
    Integer(hsize_t), dimension(1) :: adims             ! Dimensions of the attributes
    Integer(hid_t) :: plist_id                          ! Property list identifier for parallel file
    Integer(hid_t) :: aspace_id                         ! Attribute dataspace identifier
    Integer(hid_t) :: attr_id                           ! Attribute identifier
    Integer(hid_t) :: atype_id                          ! Attribute datatype identifier
    Integer(hid_t) :: gr_root_id                        ! Root group identifier
    Character(len=20) :: aname                          ! Attribute name
    Character(len=8) :: simdate
    Character(len=10) :: simtime
    Character(len=10) :: attrdate
    Character(len=12) :: attrtime



    ! Properties of the file
    Call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, h5err)
       
    Call h5pset_fapl_mpio_f(plist_id, comm, Mpi_Info_Null, h5err)

    ! Creation of the file 
    Call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5err, access_prp = plist_id)
    Call h5pclose_f(plist_id, h5err)

    ! Open the '/' group to add some info as attributes
    Call h5gopen_f(file_id, "/", gr_root_id, h5err)

    rank = 1
    adims(1) = 1

    ! Get date and time
    Call date_and_time(date=simdate,time=simtime)

    attrdate(1:4) = simdate(1:4)
    attrdate(5:5) = '/'
    attrdate(6:7) = simdate(5:6)
    attrdate(8:8) = '/'
    attrdate(9:10) = simdate(7:8)

    attrtime(1:2) = simtime(1:2)
    attrtime(3:3) = 'h'
    attrtime(4:5) = simtime(3:4)
    attrtime(6:6) = 'm'
    attrtime(7:8) = simtime(5:6)
    attrtime(9:9) = 's'
    attrtime(10:12) = simtime(8:10)

    alen = len(attrdate)
    

    ! Write date as attribute
    Call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5err)
    Call h5tset_size_f(atype_id, alen, h5err)
    Call h5screate_f(H5S_SCALAR_F, aspace_id, h5err)
    aname = 'date'
    Call h5acreate_f(gr_root_id, aname, atype_id, aspace_id, attr_id, h5err)
    Call h5awrite_f(attr_id, atype_id, attrdate, adims, h5err)
    Call h5aclose_f(attr_id, h5err)
    Call h5sclose_f(aspace_id, h5err)
    Call h5tclose_f(atype_id, h5err)

    alen = len(attrtime)

    ! Write time as attribute
    Call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5err)
    Call h5tset_size_f(atype_id, alen, h5err)
    Call h5screate_f(H5S_SCALAR_F, aspace_id, h5err)
    aname = 'time'
    Call h5acreate_f(gr_root_id, aname, atype_id, aspace_id, attr_id, h5err)
    Call h5awrite_f(attr_id, atype_id, attrtime, adims, h5err)
    Call h5aclose_f(attr_id, h5err)      
    Call h5sclose_f(aspace_id, h5err)
    Call h5tclose_f(atype_id, h5err)

    ! Write origin as attribute    
    alen = len(deforigin)
    deforigin = 'Origin not defined '

    If(Present(origin) ) Then
       deforigin = origin
    End If

    Call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5err)
    Call h5tset_size_f(atype_id,alen, h5err)
    Call h5screate_f(H5S_SCALAR_F, aspace_id, h5err)
    aname = 'Origin'
    Call h5acreate_f(gr_root_id, aname, atype_id, aspace_id, attr_id, h5err)
    Call h5awrite_f(attr_id, atype_id, deforigin, adims, h5err)
    Call h5aclose_f(attr_id, h5err)
    Call h5sclose_f(aspace_id, h5err)
    Call h5tclose_f(atype_id, h5err)

    ! Close the group.
    Call h5gclose_f(gr_root_id, h5err)

  End Subroutine hdf5_create_mpi_file


  !=======================================================================
  !> Close a 'parallel' hdf5 file
  Subroutine hdf5_close_mpi_file(file_id)
    
    Implicit None

    Integer(hid_t), intent(in) :: file_id              !< hdf5 identifier of the parallel file to close 

    ! Terminate access to the file.
    Call h5fclose_f(file_id, h5err)

#ifdef DEBUG
    If(h5err /= 0) Then
       Print *,'HDF5 error while closing the parallel file'
       Call Mpi_Abort(Mpi_Comm_World, h5err, mpierr)
    End If
#endif

  End Subroutine hdf5_close_mpi_file


  !=======================================================================
  !> Create and open a hdf5 file
  Subroutine hdf5_create_file(filename, file_id, origin)

    Implicit None

    Character(len=400), intent(in) :: filename           !< name of the parallel file to create
    Integer(hid_t), intent(out) :: file_id              !< hdf5 identifier of the parallel file 
    Character(len=50), optional :: origin   !< origin of the file (code+version) : optional argument
    
    Integer :: rank                                     ! Nb of dimensions for the attributes
    Integer(size_t) :: alen                             ! Attribute length in bytes
    Integer(hsize_t), dimension(1) :: adims             ! Dimensions of the attributes
    Integer(hid_t) :: aspace_id                         ! Attribute dataspace identifier
    Integer(hid_t) :: attr_id                           ! Attribute identifier
    Integer(hid_t) :: atype_id                          ! Attribute datatype identifier
    Integer(hid_t) :: gr_root_id                        ! Root group identifier
    Character(len=20) :: aname                          ! Attribute name
    Character(len=8) :: simdate                         ! date of creation of the file
    Character(len=10) :: simtime                        ! time of creation of the file
    Character(len=10) :: attrdate
    Character(len=12) :: attrtime
    Character(len=50) :: deforigin

    ! Create a new file using default properties.
    Call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5err)

    ! Open the '/' group to add some info as attributes
    Call h5gopen_f(file_id, "/", gr_root_id, h5err)

    rank = 1
    adims(1) = 1

    ! Get date and time
    Call date_and_time(date=simdate,time=simtime)

    attrdate(1:4) = simdate(1:4)
    attrdate(5:5) = '/'
    attrdate(6:7) = simdate(5:6)
    attrdate(8:8) = '/'
    attrdate(9:10) = simdate(7:8)

    attrtime(1:2) = simtime(1:2)
    attrtime(3:3) = 'h'
    attrtime(4:5) = simtime(3:4)
    attrtime(6:6) = 'm'
    attrtime(7:8) = simtime(5:6)
    attrtime(9:9) = 's'
    attrtime(10:12) = simtime(8:10)


    alen = len(attrdate)

    ! Write date as attribute
    Call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5err)
    Call h5tset_size_f(atype_id, alen, h5err)
    Call h5screate_f(H5S_SCALAR_F, aspace_id, h5err)
    aname = 'date'
    Call h5acreate_f(gr_root_id, aname, atype_id, aspace_id, attr_id, h5err)
    Call h5awrite_f(attr_id, atype_id, attrdate, adims, h5err)
    Call h5aclose_f(attr_id, h5err)
    Call h5sclose_f(aspace_id, h5err)
    Call h5tclose_f(atype_id, h5err)

    alen = len(attrtime)

    ! Write time as attribute
    Call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5err)
    Call h5tset_size_f(atype_id, alen, h5err)
    Call h5screate_f(H5S_SCALAR_F, aspace_id, h5err)
    aname = 'time'
    Call h5acreate_f(gr_root_id, aname, atype_id, aspace_id, attr_id, h5err)
    Call h5awrite_f(attr_id, atype_id, attrtime, adims, h5err)
    Call h5aclose_f(attr_id, h5err)      
    Call h5sclose_f(aspace_id, h5err)
    Call h5tclose_f(atype_id, h5err)

    ! Write origin as attribute    
    alen = len(deforigin)
    deforigin = 'Origin not defined '

    If(Present(origin) ) Then
       deforigin = origin
    End If
     
    Call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5err)
    Call h5tset_size_f(atype_id,alen, h5err)
    Call h5screate_f(H5S_SCALAR_F, aspace_id, h5err)
    aname = 'Origin'
    Call h5acreate_f(gr_root_id, aname, atype_id, aspace_id, attr_id, h5err)
    Call h5awrite_f(attr_id, atype_id, deforigin, adims, h5err)
    Call h5aclose_f(attr_id, h5err)
    Call h5sclose_f(aspace_id, h5err)
    Call h5tclose_f(atype_id, h5err)


    ! Close the group.
    Call h5gclose_f(gr_root_id, h5err)

  End Subroutine hdf5_create_file


  !=======================================================================
  !> Close a h5 file
  Subroutine hdf5_close_file(file_id)

    Implicit None

    Integer(hid_t), intent(in) :: file_id              !< hdf5 identifier of the file to close

    ! Terminate access to the file.
    Call h5fclose_f(file_id, h5err)

#ifdef DEBUG
    If(h5err /= 0) Then
       Print *,'HDF5 error while closing the file'
       Call Mpi_Abort(Mpi_Comm_World, h5err, mpierr)
    End If
#endif

  End Subroutine hdf5_close_file


  !===============================================!
  !         GROUP MANAGEMENT ROUTINES
  !===============================================!

  !=======================================================================
  !> Open an existing hdf5 group 
  Subroutine hdf5_open_group(id, name, group_id)

    Implicit None

    Integer(hid_t), intent(in) :: id                        !< "parent" identifier
    Character(len=16), intent(in) :: name                   !< group name
    Integer(hid_t), intent(out) :: group_id                 !< group identifier

    ! Open the group
    Call h5gopen_f(id, name, group_id, h5err)

#ifdef DEBUG
    If(h5err /= 0) Then
       Print *,'HDF5 error while opening the group ',trim(name)
       Call Mpi_Abort(Mpi_Comm_World, h5err, mpierr)
    End If
#endif

  End Subroutine hdf5_open_group


  !=======================================================================
  !> Create and open a hdf5 group 
  Subroutine hdf5_create_group(id, name, group_id)

    Implicit None

    Integer(hid_t), intent(in) :: id                        !< "parent" identifier
    Character(len=16), intent(in) :: name                   !< group name
    Integer(hid_t), intent(out) :: group_id                 !< group identifier

    ! Open the '/' group to add some info as attributes
    Call h5gcreate_f(id, name, group_id, h5err)    

#ifdef DEBUG
    If(h5err /= 0) Then
       Print *,'HDF5 error while creating the group ',trim(name)
       Call Mpi_Abort(Mpi_Comm_World, h5err, mpierr)
    End If
#endif

  End Subroutine hdf5_create_group


  !=======================================================================
  !> Close a hdf5 group 
  Subroutine hdf5_close_group(group_id)

    Implicit None

    Integer(hid_t), intent(in) :: group_id                 !< group identifier

    ! close group
    Call h5gclose_f(group_id, h5err)

  End Subroutine hdf5_close_group

  
  !===============================================!
  !         SERIAL DATA OUTPUT ROUTINES
  !===============================================!


  !=======================================================================
  !> Write a 1-D integer(kind=4) array in a serial file
  Subroutine hdf5_write_int4_1D(id, name, n1, data)

    Implicit None

    Integer(hid_t), intent(in) :: id                            !< id of the file/group where the dataset will be written 
    Character(len=16), intent(in) :: name                       !< name of the dataset                                    
    Integer, intent(in) :: n1                                   !< dimension of the data                                  
    Integer(kind=4), dimension(n1), intent(in), target :: data  !< data array                                             

    Integer, parameter :: rank = 1                   ! dimension of the array = 1
    Integer(hsize_t), dimension(1) :: dim1D          ! number of elements = size
    Integer(hid_t) :: dset_id                        ! id of the dataset
    Integer(hid_t) :: dspace_id                      ! id of the dataspace
    Integer(hid_t) :: h5_kind                        ! HDF5 real type
    Type(c_ptr) :: ptr1D                             ! pointer to the array


    h5_kind = h5kind_to_type(4,H5_INTEGER_KIND)
    ptr1D = C_LOC(data(1))
    dim1D(1) = n1
    ! Create the dataspace.
    Call h5screate_simple_f(rank, dim1D, dspace_id, h5err)

    ! Create the dataset with default properties.
    Call h5dcreate_f(id, name, h5_kind, dspace_id, &
         dset_id, h5err)

    CALL h5dwrite_f(dset_id, h5_kind, ptr1D, h5err)

    ! End access to the dataset and release resources used by it.
    Call h5dclose_f(dset_id, h5err)

    ! Terminate access to the data space.
    Call h5sclose_f(dspace_id, h5err)

  End Subroutine hdf5_write_int4_1D


  !=======================================================================
  !> Write a 1-D integer(kind=8) array in a serial file
  Subroutine hdf5_write_int8_1D(id, name, n1, data)

    Implicit None

    Integer(hid_t), intent(in) :: id                           !< id of the file/group where the dataset will be written 
    Character(len=16), intent(in) :: name                      !< name of the dataset                                    
    Integer, intent(in) :: n1                                  !< dimension of the data                                  
    Integer(kind=8), dimension(n1), intent(in), target :: data !< data array                                             

    Integer, parameter :: rank = 1                   ! dimension of the array = 1
    Integer(hsize_t), dimension(1) :: dim1D          ! number of elements = size
    Integer(hid_t) :: dset_id                        ! id of the dataset
    Integer(hid_t) :: dspace_id                      ! id of the dataspace
    Integer(hid_t) :: h5_kind                        ! HDF5 real type
    Type(c_ptr) :: ptr1D                             ! pointer to the array


    h5_kind = h5kind_to_type(8,H5_INTEGER_KIND)
    ptr1D = C_LOC(data(1))
    dim1D(1) = n1
    ! Create the dataspace.
    Call h5screate_simple_f(rank, dim1D, dspace_id, h5err)

    ! Create the dataset with default properties.
    Call h5dcreate_f(id, name, h5_kind, dspace_id, &
         dset_id, h5err)

    CALL h5dwrite_f(dset_id, h5_kind, ptr1D, h5err)

    ! End access to the dataset and release resources used by it.
    Call h5dclose_f(dset_id, h5err)

    ! Terminate access to the data space.
    Call h5sclose_f(dspace_id, h5err)

  End Subroutine hdf5_write_int8_1D


  !=======================================================================
  !> Write a 2-D integer(kind=4) array in a serial file
  Subroutine hdf5_write_int4_2D(id, name, n1, n2, data)

    Implicit None

    Integer(hid_t), intent(in) :: id                              !< id of the file/group where the dataset will be written 
    Character(len=16), intent(in) :: name                         !< name of the dataset                                    
    Integer, intent(in) :: n1                                     !< first dimension of the data                           
    Integer, intent(in) :: n2                                     !< second dimension of the data                          
    Integer(kind=4), dimension(n1,n2), intent(in), target :: data !< data array                                             

    Integer, parameter :: rank = 2                   ! dimension of the array = 1
    Integer(hsize_t), dimension(2) :: dim2D          ! number of elements = size
    Integer(hid_t) :: dset_id                        ! id of the dataset
    Integer(hid_t) :: dspace_id                      ! id of the dataspace
    Integer(hid_t) :: h5_kind                        ! HDF5 real type
    Type(c_ptr) :: ptr2D                             ! pointer to the array


    h5_kind = h5kind_to_type(4,H5_INTEGER_KIND)
    ptr2D = C_LOC(data(1,1))
    dim2D(1) = n1
    dim2D(2) = n2
    ! Create the dataspace.
    Call h5screate_simple_f(rank, dim2D, dspace_id, h5err)

    ! Create the dataset with default properties.
    Call h5dcreate_f(id, name, h5_kind, dspace_id, &
         dset_id, h5err)

    CALL h5dwrite_f(dset_id, h5_kind, ptr2D, h5err)

    ! End access to the dataset and release resources used by it.
    Call h5dclose_f(dset_id, h5err)

    ! Terminate access to the data space.
    Call h5sclose_f(dspace_id, h5err)

  End Subroutine hdf5_write_int4_2D


  !=======================================================================
  !> Write a 2-D integer(kind=8) array in a serial file
  Subroutine hdf5_write_int8_2D(id, name, n1, n2, data)

    Implicit None

    Integer(hid_t),    intent(in) :: id                            !< id of the file/group where the dataset will be written
    Character(len=16), intent(in) :: name                          !< name of the dataset                                   
    Integer,           intent(in) :: n1                            !< first dimension of the data                           
    Integer,           intent(in) :: n2                            !< second dimension of the data                          
    Integer(kind=8), dimension(n1,n2), intent(in), target :: data  !< data array                                            

    Integer, parameter :: rank = 2                   ! dimension of the array = 1
    Integer(hsize_t), dimension(2) :: dim2D          ! number of elements = size
    Integer(hid_t) :: dset_id                        ! id of the dataset
    Integer(hid_t) :: dspace_id                      ! id of the dataspace
    Integer(hid_t) :: h5_kind                        ! HDF5 real type
    Type(c_ptr) :: ptr2D                             ! pointer to the array


    h5_kind = h5kind_to_type(8,H5_INTEGER_KIND)
    ptr2D = C_LOC(data(1,1))
    dim2D(1) = n1
    dim2D(2) = n2
    ! Create the dataspace.
    Call h5screate_simple_f(rank, dim2D, dspace_id, h5err)

    ! Create the dataset with default properties.
    Call h5dcreate_f(id, name, h5_kind, dspace_id, &
         dset_id, h5err)

    CALL h5dwrite_f(dset_id, h5_kind, ptr2D, h5err)

    ! End access to the dataset and release resources used by it.
    Call h5dclose_f(dset_id, h5err)

    ! Terminate access to the data space.
    Call h5sclose_f(dspace_id, h5err)

  End Subroutine hdf5_write_int8_2D


  !=======================================================================
  !> Write a 1-D real(kind=4) array in a serial file
  Subroutine hdf5_write_real4_1D(id, name, n1, data)

    Implicit None

    Integer(hid_t), intent(in) :: id                        !< id of the file/group where the dataset will be written 
    Character(len=16), intent(in) :: name                   !< name of the dataset                                    
    Integer, intent(in) :: n1                               !< dimension of the data                                  
    Real(kind=4), dimension(n1), intent(in), target :: data !< data array                                             

    Integer, parameter :: rank = 1                   ! dimension of the array = 1
    Integer(hsize_t), dimension(1) :: dim1D          ! number of elements = size
    Integer(hid_t) :: dset_id                        ! id of the dataset
    Integer(hid_t) :: dspace_id                      ! id of the dataspace
    Integer(hid_t) :: h5_kind                        ! HDF5 real type
    Type(c_ptr) :: ptr1D                             ! pointer to the array


    h5_kind = h5kind_to_type(4,H5_REAL_KIND)
    ptr1D = C_LOC(data(1))
    dim1D(1) = n1
    ! Create the dataspace.
    Call h5screate_simple_f(rank, dim1D, dspace_id, h5err)

    ! Create the dataset with default properties.
    Call h5dcreate_f(id, name, h5_kind, dspace_id, &
         dset_id, h5err)

    CALL h5dwrite_f(dset_id, h5_kind, ptr1D, h5err)

    ! End access to the dataset and release resources used by it.
    Call h5dclose_f(dset_id, h5err)

    ! Terminate access to the data space.
    Call h5sclose_f(dspace_id, h5err)

  End Subroutine hdf5_write_real4_1D


  !=======================================================================
  !> Write a 1-D real(kind=8) array in a serial file
  Subroutine hdf5_write_real8_1D(id, name, n1, data)

    Implicit None

    Integer(hid_t), intent(in) :: id                         !< id of the file/group where the dataset will be written 
    Character(len=16), intent(in) :: name                    !< name of the dataset                                    
    Integer, intent(in) :: n1                                !< dimension of the data                                  
    Real(kind=8), dimension(n1), intent(in), target :: data  !< data array                                             

    Integer, parameter :: rank = 1                   ! dimension of the array = 1
    Integer(hsize_t), dimension(1) :: dim1D          ! number of elements = size
    Integer(hid_t) :: dset_id                        ! id of the dataset
    Integer(hid_t) :: dspace_id                      ! id of the dataspace
    Integer(hid_t) :: h5_kind                        ! HDF5 real type
    Type(c_ptr) :: ptr1D                             ! pointer to the array


    h5_kind = h5kind_to_type(8,H5_REAL_KIND)
    ptr1D = C_LOC(data(1))
    dim1D(1) = n1
    ! Create the dataspace.
    Call h5screate_simple_f(rank, dim1D, dspace_id, h5err)

    ! Create the dataset with default properties.
    Call h5dcreate_f(id, name, h5_kind, dspace_id, &
         dset_id, h5err)

    CALL h5dwrite_f(dset_id, h5_kind, ptr1D, h5err)

    ! End access to the dataset and release resources used by it.
    Call h5dclose_f(dset_id, h5err)

    ! Terminate access to the data space.
    Call h5sclose_f(dspace_id, h5err)

  End Subroutine hdf5_write_real8_1D


  !=======================================================================
  !> Write a 2-D real(kind=4) array in a serial file
  Subroutine hdf5_write_real4_2D(id, name, n1, n2, data)

    Implicit None

    Integer(hid_t), intent(in) :: id                           !< id of the file/group where the dataset will be written
    Character(len=16), intent(in) :: name                      !< name of the dataset                                   
    Integer, intent(in) :: n1                                  !< first dimension of the data                           
    Integer, intent(in) :: n2                                  !< second dimension of the data                          
    Real(kind=4), dimension(n1,n2), intent(in), target :: data !< data array                                            

    Integer, parameter :: rank = 2                   ! dimension of the array = 1
    Integer(hsize_t), dimension(2) :: dim2D          ! number of elements = size
    Integer(hid_t) :: dset_id                        ! id of the dataset
    Integer(hid_t) :: dspace_id                      ! id of the dataspace
    Integer(hid_t) :: h5_kind                        ! HDF5 real type
    Type(c_ptr) :: ptr2D                             ! pointer to the array


    h5_kind = h5kind_to_type(4,H5_REAL_KIND)
    ptr2D = C_LOC(data(1,1))
    dim2D(1) = n1
    dim2D(2) = n2
    ! Create the dataspace.
    Call h5screate_simple_f(rank, dim2D, dspace_id, h5err)

    ! Create the dataset with default properties.
    Call h5dcreate_f(id, name, h5_kind, dspace_id, &
         dset_id, h5err)

    CALL h5dwrite_f(dset_id, h5_kind, ptr2D, h5err)

    ! End access to the dataset and release resources used by it.
    Call h5dclose_f(dset_id, h5err)

    ! Terminate access to the data space.
    Call h5sclose_f(dspace_id, h5err)

  End Subroutine hdf5_write_real4_2D


  !=======================================================================
  !> Write a 2-D real(kind=8) array in a serial file
  Subroutine hdf5_write_real8_2D(id, name, n1, n2, data)

    Implicit None

    Integer(hid_t), intent(in) :: id                            !< id of the file/group where the dataset will be written
    Character(len=16), intent(in) :: name                       !< name of the dataset                                   
    Integer, intent(in) :: n1                                   !< first dimension of the data                           
    Integer, intent(in) :: n2                                   !< second dimension of the data                          
    Real(kind=8), dimension(n1,n2), intent(in), target :: data  !< data array                                            

    Integer, parameter :: rank = 2                   ! dimension of the array = 1
    Integer(hsize_t), dimension(2) :: dim2D          ! number of elements = size
    Integer(hid_t) :: dset_id                        ! id of the dataset
    Integer(hid_t) :: dspace_id                      ! id of the dataspace
    Integer(hid_t) :: h5_kind                        ! HDF5 real type
    Type(c_ptr) :: ptr2D                             ! pointer to the array


    h5_kind = h5kind_to_type(8,H5_REAL_KIND)
    ptr2D = C_LOC(data(1,1))
    dim2D(1) = n1
    dim2D(2) = n2
    ! Create the dataspace.
    Call h5screate_simple_f(rank, dim2D, dspace_id, h5err)

    ! Create the dataset with default properties.
    Call h5dcreate_f(id, name, h5_kind, dspace_id, &
         dset_id, h5err)

    CALL h5dwrite_f(dset_id, h5_kind, ptr2D, h5err)

    ! End access to the dataset and release resources used by it.
    Call h5dclose_f(dset_id, h5err)

    ! Terminate access to the data space.
    Call h5sclose_f(dspace_id, h5err)

  End Subroutine hdf5_write_real8_2D


  !===============================================!
  !           MPI DATA OUTPUT ROUTINES
  !===============================================!


  !=======================================================================
  !> Write a 1-D integer(kind=4) array in a parallel file
  Subroutine hdf5_write_mpi_int4_1D(id, name, n1, data, comm)

    Implicit None

    Integer(hid_t), intent(in) :: id                           !< id of the file/group where the dataset will be written
    Character(len=16), intent(in) :: name                      !< name of the dataset                                   
    Integer(kind=4), intent(in) :: n1                          !< dimension of the local data                           
    Integer(kind=4), dimension(n1), intent(in), target :: data !< local data array                                      
    Integer(kind=4), intent(in) :: comm                        !< MPI communicator used

    Integer(hsize_t), dimension(1) :: local_dims               ! Local dataset dimensions
    Integer(hsize_t), dimension(1) :: global_dims              ! Global dataspace dimensions
    Integer(hsize_t), dimension(1) :: offset                   ! Offset
    Type(c_ptr) :: ptr                                         ! Pointer to the array
    Integer, parameter :: rank = 1                             ! Number of dimensions = 1
    Integer(hid_t) :: dset_id                                  ! Dataset identifier
    Integer(hid_t) :: plist_id                                 ! Property list identifier for parallel file
    Integer(hid_t) :: dspace_id                                ! Dataspace identifier
    Integer(hid_t) :: mem_id                                   ! Local memory dataspace identifier for parallel I/O
    Integer(hid_t) :: h5_kind                                  ! HDF5 integer type
    Integer(kind=4), dimension(:), allocatable :: ntab         ! tmp tab used for global dims and offset calculation
    Integer(kind=4) :: i
    Integer(kind=4) :: nproc, idproc

    Call Mpi_Comm_Size(comm, nproc, mpierr)
    Call Mpi_Comm_Rank(comm, idproc, mpierr)

    ptr = C_LOC(data(1))
    
    Allocate(ntab(nproc))

    ! Creation du type hdf5 correspondant au type des donnees a ecrire
    h5_kind = h5kind_to_type(4,H5_INTEGER_KIND)

    ! Offset and global dimensions have to be computed
    Call Mpi_Allgather(n1, 1, Mpi_Integer, ntab, 1, Mpi_Integer, comm, mpierr)

    global_dims(1) = 0

    Do i = 1, nproc
       global_dims(1) = global_dims(1) + ntab(i)
    End Do
    
    local_dims(1) = n1
    
    offset(1) = 0
    
    Do i = 1, idproc
       offset(1) = offset(1) + ntab(i)
    End Do

    ! Create the dataspace.
    Call h5screate_simple_f(rank, global_dims, dspace_id, h5err)      

    ! Create the dataset with default properties.
    Call h5dcreate_f(id, name, h5_kind, dspace_id, &
         dset_id, h5err)

    ! Create the local memory dataspace
    Call h5screate_simple_f(rank, local_dims, mem_id, h5err)

    Call h5dget_space_f(dset_id, dspace_id, h5err)
    Call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, local_dims, h5err)

    ! Create property list for collective dataset write
    Call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, h5err) 
    Call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, h5err)

    Call h5dwrite_f(dset_id, h5_kind, ptr, h5err, &
         file_space_id = dspace_id, mem_space_id = mem_id, xfer_prp = plist_id)

    ! End access to the dataset and release resources used by it.
    Call h5dclose_f(dset_id, h5err)

    ! Terminate access to the data space.
    Call h5sclose_f(dspace_id, h5err)

    ! Terminate access to the property list.
    Call h5pclose_f(plist_id, h5err)

    ! End access to the local memory dataspace
    Call h5sclose_f(mem_id, h5err)
    
    Deallocate(ntab)

  End Subroutine hdf5_write_mpi_int4_1D


  !=======================================================================
  !> Write a 2-D integer(kind=4) array in a parallel file
  !! The array is distributed along the 2nd dimension
  Subroutine hdf5_write_mpi_int4_2D(id, name, n1, n2, data, comm)

    Implicit None

    Integer(hid_t), intent(in) :: id                              !< id of the file/group where the dataset will be written 
    Character(len=16), intent(in) :: name                         !< name of the dataset                                    
    Integer(kind=4), intent(in) :: n1                             !< first dimension of the local data                      
    Integer(kind=4), intent(in) :: n2                             !< second dimension of the local data                     
    Integer(kind=4), dimension(n1,n2), intent(in), target :: data !< local data array                                       
    Integer(kind=4), intent(in) :: comm                           !< MPI communicator used

    Integer(hsize_t), dimension(2) :: local_dims               ! Local dataset dimensions
    Integer(hsize_t), dimension(2) :: global_dims              ! Global dataspace dimensions
    Integer(hsize_t), dimension(2) :: offset                   ! Offset
    Type(c_ptr) :: ptr                                         ! Pointer to the array
    Integer, parameter :: rank = 2                             ! Number of dimensions = 1
    Integer(hid_t) :: dset_id                                  ! Dataset identifier
    Integer(hid_t) :: plist_id                                 ! Property list identifier for parallel file
    Integer(hid_t) :: dspace_id                                ! Dataspace identifier
    Integer(hid_t) :: mem_id                                   ! Local memory dataspace identifier for parallel I/O
    Integer(hid_t) :: h5_kind                                  ! HDF5 integer type
    Integer(kind=4), dimension(:), allocatable :: ntab         ! tmp tab used for global dims and offset calculation
    Integer(kind=4) :: i
    Integer(kind=4) :: nproc, idproc

    Call Mpi_Comm_Size(comm, nproc, mpierr)
    Call Mpi_Comm_Rank(comm, idproc, mpierr)

    ptr = C_LOC(data(1,1))
    
    Allocate(ntab(nproc))

    ! Creation du type hdf5 correspondant au type des donnees a ecrire
    h5_kind = h5kind_to_type(4,H5_INTEGER_KIND)

    ! Offset and global dimensions have to be computed
    Call Mpi_Allgather(n2, 1, Mpi_Integer, ntab, 1, Mpi_Integer, comm, mpierr)

    global_dims(1) = n1
    global_dims(2) = 0

    Do i = 1, nproc
       global_dims(2) = global_dims(2) + ntab(i)
    End Do
    
    local_dims(1) = n1
    local_dims(2) = n2

    offset(1) = 0
    offset(2) = 0
    
    Do i = 1, idproc
       offset(2) = offset(2) + ntab(i)
    End Do

    ! Create the dataspace.
    Call h5screate_simple_f(rank, global_dims, dspace_id, h5err)      

    ! Create the dataset with default properties.
    Call h5dcreate_f(id, name, h5_kind, dspace_id, &
         dset_id, h5err)

    ! Create the local memory dataspace
    Call h5screate_simple_f(rank, local_dims, mem_id, h5err)

    Call h5dget_space_f(dset_id, dspace_id, h5err)
    Call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, local_dims, h5err)

    ! Create property list for collective dataset write
    Call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, h5err) 
    Call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, h5err)

    Call h5dwrite_f(dset_id, h5_kind, ptr, h5err, &
         file_space_id = dspace_id, mem_space_id = mem_id, xfer_prp = plist_id)

    ! End access to the dataset and release resources used by it.
    Call h5dclose_f(dset_id, h5err)

    ! Terminate access to the data space.
    Call h5sclose_f(dspace_id, h5err)

    ! Terminate access to the property list.
    Call h5pclose_f(plist_id, h5err)

    ! End access to the local memory dataspace
    Call h5sclose_f(mem_id, h5err)
    
    Deallocate(ntab)

  End Subroutine hdf5_write_mpi_int4_2D


  !=======================================================================
  !> Write a 1-D integer(kind=8) array in a parallel file
  Subroutine hdf5_write_mpi_int8_1D(id, name, n1, data, comm)

    Implicit None

    Integer(hid_t), intent(in) :: id                           !< id of the file/group where the dataset will be written
    Character(len=16), intent(in) :: name                      !< name of the dataset                                   
    Integer(kind=4), intent(in) :: n1                          !< dimension of the local data                           
    Integer(kind=8), dimension(n1), intent(in), target :: data !< local data array                                      
    Integer(kind=4), intent(in) :: comm                        !< MPI communicator used                                 

    Integer(hsize_t), dimension(1) :: local_dims               ! Local dataset dimensions
    Integer(hsize_t), dimension(1) :: global_dims              ! Global dataspace dimensions
    Integer(hsize_t), dimension(1) :: offset                   ! Offset
    Type(c_ptr) :: ptr                                         ! Pointer to the array
    Integer, parameter :: rank = 1                             ! Number of dimensions = 1
    Integer(hid_t) :: dset_id                                  ! Dataset identifier
    Integer(hid_t) :: plist_id                                 ! Property list identifier for parallel file
    Integer(hid_t) :: dspace_id                                ! Dataspace identifier
    Integer(hid_t) :: mem_id                                   ! Local memory dataspace identifier for parallel I/O
    Integer(hid_t) :: h5_kind                                  ! HDF5 integer type
    Integer(kind=4), dimension(:), allocatable :: ntab         ! tmp tab used for global dims and offset calculation
    Integer(kind=4) :: i
    Integer(kind=4) :: nproc, idproc

    Call Mpi_Comm_Size(comm, nproc, mpierr)
    Call Mpi_Comm_Rank(comm, idproc, mpierr)

    ptr = C_LOC(data(1))
    
    Allocate(ntab(nproc))

    ! Creation du type hdf5 correspondant au type des donnees a ecrire
    h5_kind = h5kind_to_type(8,H5_INTEGER_KIND)

    ! Offset and global dimensions have to be computed
    Call Mpi_Allgather(n1, 1, Mpi_Integer, ntab, 1, Mpi_Integer, comm, mpierr)

    global_dims(1) = 0

    Do i = 1, nproc
       global_dims(1) = global_dims(1) + ntab(i)
    End Do
    
    local_dims(1) = n1
    
    offset(1) = 0
    
    Do i = 1, idproc
       offset(1) = offset(1) + ntab(i)
    End Do

    ! Create the dataspace.
    Call h5screate_simple_f(rank, global_dims, dspace_id, h5err)      

    ! Create the dataset with default properties.
    Call h5dcreate_f(id, name, h5_kind, dspace_id, &
         dset_id, h5err)

    ! Create the local memory dataspace
    Call h5screate_simple_f(rank, local_dims, mem_id, h5err)

    Call h5dget_space_f(dset_id, dspace_id, h5err)
    Call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, local_dims, h5err)

    ! Create property list for collective dataset write
    Call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, h5err) 
    Call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, h5err)

    Call h5dwrite_f(dset_id, h5_kind, ptr, h5err, &
         file_space_id = dspace_id, mem_space_id = mem_id, xfer_prp = plist_id)

    ! End access to the dataset and release resources used by it.
    Call h5dclose_f(dset_id, h5err)

    ! Terminate access to the data space.
    Call h5sclose_f(dspace_id, h5err)

    ! Terminate access to the property list.
    Call h5pclose_f(plist_id, h5err)

    ! End access to the local memory dataspace
    Call h5sclose_f(mem_id, h5err)
    
    Deallocate(ntab)

  End Subroutine hdf5_write_mpi_int8_1D


  !=======================================================================
  !> Write a 2-D integer(kind=8) array in a parallel file
  !> The array is distributed along the 2nd dimension
  Subroutine hdf5_write_mpi_int8_2D(id, name, n1, n2, data, comm)

    Implicit None

    Integer(hid_t), intent(in) :: id                              !< id of the file/group where the dataset will be written
    Character(len=16), intent(in) :: name                         !< name of the dataset                                   
    Integer(kind=4), intent(in) :: n1                             !< first dimension of the local data                           
    Integer(kind=4), intent(in) :: n2                             !< second dimension of the local data                           
    Integer(kind=8), dimension(n1,n2), intent(in), target :: data !< local data array                                      
    Integer(kind=4), intent(in) :: comm                           !< MPI communicator used                                 

    Integer(hsize_t), dimension(2) :: local_dims               ! Local dataset dimensions
    Integer(hsize_t), dimension(2) :: global_dims              ! Global dataspace dimensions
    Integer(hsize_t), dimension(2) :: offset                   ! Offset
    Type(c_ptr) :: ptr                                         ! Pointer to the array
    Integer, parameter :: rank = 2                             ! Number of dimensions = 1
    Integer(hid_t) :: dset_id                                  ! Dataset identifier
    Integer(hid_t) :: plist_id                                 ! Property list identifier for parallel file
    Integer(hid_t) :: dspace_id                                ! Dataspace identifier
    Integer(hid_t) :: mem_id                                   ! Local memory dataspace identifier for parallel I/O
    Integer(hid_t) :: h5_kind                                  ! HDF5 integer type
    Integer(kind=4), dimension(:), allocatable :: ntab         ! tmp tab used for global dims and offset calculation
    Integer(kind=4) :: i
    Integer(kind=4) :: nproc, idproc

    Call Mpi_Comm_Size(comm, nproc, mpierr)
    Call Mpi_Comm_Rank(comm, idproc, mpierr)

    ptr = C_LOC(data(1,1))
    
    Allocate(ntab(nproc))

    ! Creation du type hdf5 correspondant au type des donnees a ecrire
    h5_kind = h5kind_to_type(8,H5_INTEGER_KIND)

    ! Offset and global dimensions have to be computed
    Call Mpi_Allgather(n2, 1, Mpi_Integer, ntab, 1, Mpi_Integer, comm, mpierr)

    global_dims(1) = n1
    global_dims(2) = 0

    Do i = 1, nproc
       global_dims(2) = global_dims(2) + ntab(i)
    End Do
    
    local_dims(1) = n1
    local_dims(2) = n2
    
    offset(1) = 0
    offset(2) = 0
    
    Do i = 1, idproc
       offset(2) = offset(2) + ntab(i)
    End Do

    ! Create the dataspace.
    Call h5screate_simple_f(rank, global_dims, dspace_id, h5err)      

    ! Create the dataset with default properties.
    Call h5dcreate_f(id, name, h5_kind, dspace_id, &
         dset_id, h5err)

    ! Create the local memory dataspace
    Call h5screate_simple_f(rank, local_dims, mem_id, h5err)

    Call h5dget_space_f(dset_id, dspace_id, h5err)
    Call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, local_dims, h5err)

    ! Create property list for collective dataset write
    Call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, h5err) 
    Call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, h5err)

    Call h5dwrite_f(dset_id, h5_kind, ptr, h5err, &
         file_space_id = dspace_id, mem_space_id = mem_id, xfer_prp = plist_id)

    ! End access to the dataset and release resources used by it.
    Call h5dclose_f(dset_id, h5err)

    ! Terminate access to the data space.
    Call h5sclose_f(dspace_id, h5err)

    ! Terminate access to the property list.
    Call h5pclose_f(plist_id, h5err)

    ! End access to the local memory dataspace
    Call h5sclose_f(mem_id, h5err)
    
    Deallocate(ntab)

  End Subroutine hdf5_write_mpi_int8_2D


  !=======================================================================
  !> Write a 1-D real(kind=4) array in a parallel file
  Subroutine hdf5_write_mpi_real4_1D(id, name, n1, data, comm)

    Implicit None

    Integer(hid_t), intent(in) :: id                           !< id of the file/group where the dataset will be written
    Character(len=16), intent(in) :: name                      !< name of the dataset                                   
    Integer(kind=4), intent(in) :: n1                          !< dimension of the local data                           
    Real(kind=4), dimension(n1), intent(in), target :: data    !< local data array                                      
    Integer(kind=4), intent(in) :: comm                        !< MPI communicator used                                 

    Integer(hsize_t), dimension(1) :: local_dims               ! Local dataset dimensions
    Integer(hsize_t), dimension(1) :: global_dims              ! Global dataspace dimensions
    Integer(hsize_t), dimension(1) :: offset                   ! Offset
    Type(c_ptr) :: ptr                                         ! Pointer to the array
    Integer, parameter :: rank = 1                             ! Number of dimensions = 1
    Integer(hid_t) :: dset_id                                  ! Dataset identifier
    Integer(hid_t) :: plist_id                                 ! Property list identifier for parallel file
    Integer(hid_t) :: dspace_id                                ! Dataspace identifier
    Integer(hid_t) :: mem_id                                   ! Local memory dataspace identifier for parallel I/O
    Integer(hid_t) :: h5_kind                                  ! HDF5 integer type
    Integer(kind=4), dimension(:), allocatable :: ntab         ! tmp tab used for global dims and offset calculation
    Integer(kind=4) :: i
    Integer(kind=4) :: nproc, idproc

    Call Mpi_Comm_Size(comm, nproc, mpierr)
    Call Mpi_Comm_Rank(comm, idproc, mpierr)

    ptr = C_LOC(data(1))
    
    Allocate(ntab(nproc))

    ! Creation du type hdf5 correspondant au type des donnees a ecrire
    h5_kind = h5kind_to_type(4,H5_REAL_KIND)

    ! Offset and global dimensions have to be computed
    Call Mpi_Allgather(n1, 1, Mpi_Integer, ntab, 1, Mpi_Integer, comm, mpierr)

    global_dims(1) = 0

    Do i = 1, nproc
       global_dims(1) = global_dims(1) + ntab(i)
    End Do
    
    local_dims(1) = n1
    
    offset(1) = 0
    
    Do i = 1, idproc
       offset(1) = offset(1) + ntab(i)
    End Do

    ! Create the dataspace.
    Call h5screate_simple_f(rank, global_dims, dspace_id, h5err)      

    ! Create the dataset with default properties.
    Call h5dcreate_f(id, name, h5_kind, dspace_id, &
         dset_id, h5err)

    ! Create the local memory dataspace
    Call h5screate_simple_f(rank, local_dims, mem_id, h5err)

    Call h5dget_space_f(dset_id, dspace_id, h5err)
    Call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, local_dims, h5err)

    ! Create property list for collective dataset write
    Call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, h5err) 
    Call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, h5err)

    Call h5dwrite_f(dset_id, h5_kind, ptr, h5err, &
         file_space_id = dspace_id, mem_space_id = mem_id, xfer_prp = plist_id)

    ! End access to the dataset and release resources used by it.
    Call h5dclose_f(dset_id, h5err)

    ! Terminate access to the data space.
    Call h5sclose_f(dspace_id, h5err)

    ! Terminate access to the property list.
    Call h5pclose_f(plist_id, h5err)

    ! End access to the local memory dataspace
    Call h5sclose_f(mem_id, h5err)
    
    Deallocate(ntab)

  End Subroutine hdf5_write_mpi_real4_1D


  !=======================================================================
  !> Write a 2-D real(kind=4) array in a parallel file
  !> The array is distributed along the 2nd dimension
  Subroutine hdf5_write_mpi_real4_2D(id, name, n1, n2, data, comm)

    Implicit None

    Integer(hid_t), intent(in) :: id                            !< id of the file/group where the dataset will be written
    Character(len=16), intent(in) :: name                       !< name of the dataset                                   
    Integer(kind=4), intent(in) :: n1                           !< first dimension of the local data                     
    Integer(kind=4), intent(in) :: n2                           !< second dimension of the local data                    
    Real(kind=4), dimension(n1,n2), intent(in), target :: data  !< local data array                                      
    Integer(kind=4), intent(in) :: comm                         !< MPI communicator used                                 

    Integer(hsize_t), dimension(2) :: local_dims               ! Local dataset dimensions
    Integer(hsize_t), dimension(2) :: global_dims              ! Global dataspace dimensions
    Integer(hsize_t), dimension(2) :: offset                   ! Offset
    Type(c_ptr) :: ptr                                         ! Pointer to the array
    Integer, parameter :: rank = 2                             ! Number of dimensions = 1
    Integer(hid_t) :: dset_id                                  ! Dataset identifier
    Integer(hid_t) :: plist_id                                 ! Property list identifier for parallel file
    Integer(hid_t) :: dspace_id                                ! Dataspace identifier
    Integer(hid_t) :: mem_id                                   ! Local memory dataspace identifier for parallel I/O
    Integer(hid_t) :: h5_kind                                  ! HDF5 integer type
    Integer(kind=4), dimension(:), allocatable :: ntab         ! tmp tab used for global dims and offset calculation
    Integer(kind=4) :: i
    Integer(kind=4) :: nproc, idproc

    Call Mpi_Comm_Size(comm, nproc, mpierr)
    Call Mpi_Comm_Rank(comm, idproc, mpierr)

    ptr = C_LOC(data(1,1))
    
    Allocate(ntab(nproc))

    ! Creation du type hdf5 correspondant au type des donnees a ecrire
    h5_kind = h5kind_to_type(4,H5_REAL_KIND)

    ! Offset and global dimensions have to be computed
    Call Mpi_Allgather(n2, 1, Mpi_Integer, ntab, 1, Mpi_Integer, comm, mpierr)

    global_dims(1) = n1
    global_dims(2) = 0

    Do i = 1, nproc
       global_dims(2) = global_dims(2) + ntab(i)
    End Do
    
    local_dims(1) = n1
    local_dims(2) = n2

    offset(1) = 0
    offset(2) = 0
    
    Do i = 1, idproc
       offset(2) = offset(2) + ntab(i)
    End Do

    ! Create the dataspace.
    Call h5screate_simple_f(rank, global_dims, dspace_id, h5err)      

    ! Create the dataset with default properties.
    Call h5dcreate_f(id, name, h5_kind, dspace_id, &
         dset_id, h5err)

    ! Create the local memory dataspace
    Call h5screate_simple_f(rank, local_dims, mem_id, h5err)

    Call h5dget_space_f(dset_id, dspace_id, h5err)
    Call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, local_dims, h5err)

    ! Create property list for collective dataset write
    Call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, h5err) 
    Call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, h5err)

    Call h5dwrite_f(dset_id, h5_kind, ptr, h5err, &
         file_space_id = dspace_id, mem_space_id = mem_id, xfer_prp = plist_id)

    ! End access to the dataset and release resources used by it.
    Call h5dclose_f(dset_id, h5err)

    ! Terminate access to the data space.
    Call h5sclose_f(dspace_id, h5err)

    ! Terminate access to the property list.
    Call h5pclose_f(plist_id, h5err)

    ! End access to the local memory dataspace
    Call h5sclose_f(mem_id, h5err)
    
    Deallocate(ntab)

  End Subroutine hdf5_write_mpi_real4_2D


  !=======================================================================
  !> Write a 1-D real(kind=8) array in a parallel file
  Subroutine hdf5_write_mpi_real8_1D(id, name, n1, data, comm)

    Implicit None

    Integer(hid_t), intent(in) :: id                           !< id of the file/group where the dataset will be written
    Character(len=16), intent(in) :: name                      !< name of the dataset                                   
    Integer(kind=4), intent(in) :: n1                          !< dimension of the local data                           
    Real(kind=8), dimension(n1), intent(in), target :: data    !< local data array                                      
    Integer(kind=4), intent(in) :: comm                        !< MPI communicator used                                 

    Integer(hsize_t), dimension(1) :: local_dims               ! Local dataset dimensions
    Integer(hsize_t), dimension(1) :: global_dims              ! Global dataspace dimensions
    Integer(hsize_t), dimension(1) :: offset                   ! Offset
    Type(c_ptr) :: ptr                                         ! Pointer to the array
    Integer, parameter :: rank = 1                             ! Number of dimensions = 1
    Integer(hid_t) :: dset_id                                  ! Dataset identifier
    Integer(hid_t) :: plist_id                                 ! Property list identifier for parallel file
    Integer(hid_t) :: dspace_id                                ! Dataspace identifier
    Integer(hid_t) :: mem_id                                   ! Local memory dataspace identifier for parallel I/O
    Integer(hid_t) :: h5_kind                                  ! HDF5 integer type
    Integer(kind=4), dimension(:), allocatable :: ntab         ! tmp tab used for global dims and offset calculation
    Integer(kind=4) :: i
    Integer(kind=4) :: nproc, idproc

    Call Mpi_Comm_Size(comm, nproc, mpierr)
    Call Mpi_Comm_Rank(comm, idproc, mpierr)

    ptr = C_LOC(data(1))
    
    Allocate(ntab(nproc))

    ! Creation du type hdf5 correspondant au type des donnees a ecrire
    h5_kind = h5kind_to_type(8,H5_REAL_KIND)

    ! Offset and global dimensions have to be computed
    Call Mpi_Allgather(n1, 1, Mpi_Integer, ntab, 1, Mpi_Integer, comm, mpierr)

    global_dims(1) = 0

    Do i = 1, nproc
       global_dims(1) = global_dims(1) + ntab(i)
    End Do
    
    local_dims(1) = n1
    
    offset(1) = 0
    
    Do i = 1, idproc
       offset(1) = offset(1) + ntab(i)
    End Do

    ! Create the dataspace.
    Call h5screate_simple_f(rank, global_dims, dspace_id, h5err)      

    ! Create the dataset with default properties.
    Call h5dcreate_f(id, name, h5_kind, dspace_id, &
         dset_id, h5err)

    ! Create the local memory dataspace
    Call h5screate_simple_f(rank, local_dims, mem_id, h5err)

    Call h5dget_space_f(dset_id, dspace_id, h5err)
    Call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, local_dims, h5err)

    ! Create property list for collective dataset write
    Call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, h5err) 
    Call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, h5err)

    Call h5dwrite_f(dset_id, h5_kind, ptr, h5err, &
         file_space_id = dspace_id, mem_space_id = mem_id, xfer_prp = plist_id)

    ! End access to the dataset and release resources used by it.
    Call h5dclose_f(dset_id, h5err)

    ! Terminate access to the data space.
    Call h5sclose_f(dspace_id, h5err)

    ! Terminate access to the property list.
    Call h5pclose_f(plist_id, h5err)

    ! End access to the local memory dataspace
    Call h5sclose_f(mem_id, h5err)
    
    Deallocate(ntab)

  End Subroutine hdf5_write_mpi_real8_1D


  !=======================================================================
  !> Write a 2-D real(kind=8) array in a parallel file
  !> The array is distributed along the 2nd dimension
  Subroutine hdf5_write_mpi_real8_2D(id, name, n1, n2, data, comm)

    Implicit None

    Integer(hid_t), intent(in) :: id                           !< id of the file/group where the dataset will be written
    Character(len=16), intent(in) :: name                      !< name of the dataset                                   
    Integer(kind=4), intent(in) :: n1                          !< first dimension of the local data                     
    Integer(kind=4), intent(in) :: n2                          !< second dimension of the local data                    
    Real(kind=8), dimension(n1,n2), intent(in), target :: data !< local data array                                      
    Integer(kind=4), intent(in) :: comm                        !< MPI communicator used                                 

    Integer(hsize_t), dimension(2) :: local_dims               ! Local dataset dimensions
    Integer(hsize_t), dimension(2) :: global_dims              ! Global dataspace dimensions
    Integer(hsize_t), dimension(2) :: offset                   ! Offset
    Type(c_ptr) :: ptr                                         ! Pointer to the array
    Integer, parameter :: rank = 2                             ! Number of dimensions = 1
    Integer(hid_t) :: dset_id                                  ! Dataset identifier
    Integer(hid_t) :: plist_id                                 ! Property list identifier for parallel file
    Integer(hid_t) :: dspace_id                                ! Dataspace identifier
    Integer(hid_t) :: mem_id                                   ! Local memory dataspace identifier for parallel I/O
    Integer(hid_t) :: h5_kind                                  ! HDF5 integer type
    Integer(kind=4), dimension(:), allocatable :: ntab         ! tmp tab used for global dims and offset calculation
    Integer(kind=4) :: i
    Integer(kind=4) :: nproc, idproc

    Call Mpi_Comm_Size(comm, nproc, mpierr)
    Call Mpi_Comm_Rank(comm, idproc, mpierr)

    ptr = C_LOC(data(1,1))
    
    Allocate(ntab(nproc))

    ! Creation du type hdf5 correspondant au type des donnees a ecrire
    h5_kind = h5kind_to_type(8,H5_REAL_KIND)

    ! Offset and global dimensions have to be computed
    Call Mpi_Allgather(n2, 1, Mpi_Integer, ntab, 1, Mpi_Integer, comm, mpierr)

    global_dims(1) = n1
    global_dims(2) = 0

    Do i = 1, nproc
       global_dims(2) = global_dims(2) + ntab(i)
    End Do
    
    local_dims(1) = n1
    local_dims(2) = n2
    
    offset(1) = 0
    offset(2) = 0
    
    Do i = 1, idproc
       offset(2) = offset(2) + ntab(i)
    End Do

    ! Create the dataspace.
    Call h5screate_simple_f(rank, global_dims, dspace_id, h5err)      

    ! Create the dataset with default properties.
    Call h5dcreate_f(id, name, h5_kind, dspace_id, &
         dset_id, h5err)

    ! Create the local memory dataspace
    Call h5screate_simple_f(rank, local_dims, mem_id, h5err)

    Call h5dget_space_f(dset_id, dspace_id, h5err)
    Call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, local_dims, h5err)

    ! Create property list for collective dataset write
    Call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, h5err) 
    Call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, h5err)

    Call h5dwrite_f(dset_id, h5_kind, ptr, h5err, &
         file_space_id = dspace_id, mem_space_id = mem_id, xfer_prp = plist_id)

    ! End access to the dataset and release resources used by it.
    Call h5dclose_f(dset_id, h5err)

    ! Terminate access to the data space.
    Call h5sclose_f(dspace_id, h5err)

    ! Terminate access to the property list.
    Call h5pclose_f(plist_id, h5err)

    ! End access to the local memory dataspace
    Call h5sclose_f(mem_id, h5err)
    
    Deallocate(ntab)

  End Subroutine hdf5_write_mpi_real8_2D


  !===============================================!
  !         SERIAL DATA INPUT ROUTINES
  !===============================================!


  !=======================================================================
  !> Read a 1-D integer(kind=4) array from a serial file
  Subroutine hdf5_read_int4_1D(id, name, n1, data)

    Implicit None

    Integer(hid_t), intent(in) :: id                               !< id of the file/group where the dataset will be read
    Character(len=16), intent(in) :: name                          !< name of the dataset                                
    Integer(kind=4), intent(in) :: n1                              !< dimension of the array to read                     
    Integer(kind=4), dimension(n1), intent(inout), target :: data  !< array                                               

    Type(c_ptr) :: ptr
    Integer(hid_t) :: dset_id                               ! id of the dataset
    Integer(hid_t) :: h5_kind                               ! HDF5 integer type

    ptr = C_LOC(data(1))
    ! hdf5 type corresponding to the integer type to read
    h5_kind = h5kind_to_type(4,H5_INTEGER_KIND)

    ! open the dataset
    Call h5dopen_f(id, name, dset_id, h5err)
    
    ! read the dataset
    Call h5dread_f(dset_id, h5_kind, ptr, h5err)

    ! close the dataset
    Call h5dclose_f(dset_id, h5err)

  End Subroutine hdf5_read_int4_1D


  !=======================================================================
  !> Read a 1-D integer(kind=8) array from a serial file
  Subroutine hdf5_read_int8_1D(id, name, n1, data)

    Implicit None

    Integer(hid_t), intent(in) :: id                               !< id of the file/group where the dataset will be read 
    Character(len=16), intent(in) :: name                          !< name of the dataset                                 
    Integer(kind=4), intent(in) :: n1                              !< dimension of the array to read                
    Integer(kind=8), dimension(n1), intent(inout), target :: data  !< array                                              
                                                                                                                  
    Type(c_ptr) :: ptr
    Integer(hid_t) :: dset_id                               ! id of the dataset
    Integer(hid_t) :: h5_kind                               ! HDF5 integer type

    ptr = C_LOC(data(1))
    ! hdf5 type corresponding to the integer type to read
    h5_kind = h5kind_to_type(8,H5_INTEGER_KIND)

    ! open the dataset
    Call h5dopen_f(id, name, dset_id, h5err)
    
    ! read the dataset
    Call h5dread_f(dset_id, h5_kind, ptr, h5err)

    ! close the dataset
    Call h5dclose_f(dset_id, h5err)

  End Subroutine hdf5_read_int8_1D


  !=======================================================================
  !> Read a 2-D integer(kind=4) array from a serial file
  Subroutine hdf5_read_int4_2D(id, name, n1, n2, data)

    Implicit None

    Integer(hid_t), intent(in) :: id                                   !< id of the file/group where the dataset will be read 
    Character(len=16), intent(in) :: name                              !< name of the dataset                                 
    Integer(kind=4), intent(in) :: n1                                  !< first dimension of the array to read                
    Integer(kind=4), intent(in) :: n2                                  !< second dimension of the array to read               
    Integer(kind=4), dimension(n1,n2), intent(inout), target :: data   !< array                                               

    Type(c_ptr) :: ptr
    Integer(hid_t) :: dset_id                               ! id of the dataset
    Integer(hid_t) :: h5_kind                               ! HDF5 integer type

    ptr = C_LOC(data(1,1))
    ! hdf5 type corresponding to the integer type to read
    h5_kind = h5kind_to_type(4,H5_INTEGER_KIND)

    ! open the dataset
    Call h5dopen_f(id, name, dset_id, h5err)
    
    ! read the dataset
    Call h5dread_f(dset_id, h5_kind, ptr, h5err)

    ! close the dataset
    Call h5dclose_f(dset_id, h5err)

  End Subroutine hdf5_read_int4_2D


  !=======================================================================
  !> Read a 2-D integer(kind=8) array from a serial file
  Subroutine hdf5_read_int8_2D(id, name, n1, n2, data)

    Implicit None

    Integer(hid_t), intent(in) :: id                                  !< id of the file/group where the dataset will be read 
    Character(len=16), intent(in) :: name                             !< name of the dataset                                 
    Integer(kind=4), intent(in) :: n1                                 !< first dimension of the array to read                
    Integer(kind=4), intent(in) :: n2                                 !< second dimension of the array to read               
    Integer(kind=8), dimension(n1,n2), intent(inout), target :: data  !< array                                                

    Type(c_ptr) :: ptr
    Integer(hid_t) :: dset_id                               ! id of the dataset
    Integer(hid_t) :: h5_kind                               ! HDF5 integer type

    ptr = C_LOC(data(1,1))
    ! hdf5 type corresponding to the integer type to read
    h5_kind = h5kind_to_type(8,H5_INTEGER_KIND)

    ! open the dataset
    Call h5dopen_f(id, name, dset_id, h5err)
    
    ! read the dataset
    Call h5dread_f(dset_id, h5_kind, ptr, h5err)

    ! close the dataset
    Call h5dclose_f(dset_id, h5err)

  End Subroutine hdf5_read_int8_2D


  !=======================================================================
  !> Read a 1-D real(kind=4) array from a serial file
  Subroutine hdf5_read_real4_1D(id, name, n1, data)

    Implicit None

    Integer(hid_t), intent(in) :: id                            !< id of the file/group where the dataset will be read
    Character(len=16), intent(in) :: name                       !< name of the dataset                                
    Integer(kind=4), intent(in) :: n1                           !< dimension of the array to read                     
    Real(kind=4), dimension(n1), intent(inout), target :: data  !< array                                               

    Type(c_ptr) :: ptr
    Integer(hid_t) :: dset_id                               ! id of the dataset
    Integer(hid_t) :: h5_kind                               ! HDF5 integer type

    ptr = C_LOC(data(1))
    ! hdf5 type corresponding to the integer type to read
    h5_kind = h5kind_to_type(4,H5_REAL_KIND)

    ! open the dataset
    Call h5dopen_f(id, name, dset_id, h5err)
    
    ! read the dataset
    Call h5dread_f(dset_id, h5_kind, ptr, h5err)

    ! close the dataset
    Call h5dclose_f(dset_id, h5err)

  End Subroutine hdf5_read_real4_1D


  !=======================================================================
  !> Read a 1-D real(kind=8) array from a serial file
  Subroutine hdf5_read_real8_1D(id, name, n1, data)

    Implicit None

    Integer(hid_t), intent(in) :: id                            !< id of the file/group where the dataset will be read
    Character(len=16), intent(in) :: name                       !< name of the dataset                                
    Integer(kind=4), intent(in) :: n1                           !< dimension of the array to read                     
    Real(kind=8), dimension(n1), intent(inout), target :: data  !< array                                               

    Type(c_ptr) :: ptr
    Integer(hid_t) :: dset_id                               ! id of the dataset
    Integer(hid_t) :: h5_kind                               ! HDF5 integer type

    ptr = C_LOC(data(1))
    ! hdf5 type corresponding to the integer type to read
    h5_kind = h5kind_to_type(8,H5_REAL_KIND)

    ! open the dataset
    Call h5dopen_f(id, name, dset_id, h5err)
    
    ! read the dataset
    Call h5dread_f(dset_id, h5_kind, ptr, h5err)

    ! close the dataset
    Call h5dclose_f(dset_id, h5err)

  End Subroutine hdf5_read_real8_1D


  !=======================================================================
  !> Read a 2-D real(kind=4) array from a serial file
  Subroutine hdf5_read_real4_2D(id, name, n1, n2, data)

    Implicit None

    Integer(hid_t), intent(in) :: id                               !< id of the file/group where the dataset will be read 
    Character(len=16), intent(in) :: name                          !< name of the dataset                                 
    Integer(kind=4), intent(in) :: n1                              !< first dimension of the array to read                
    Integer(kind=4), intent(in) :: n2                              !< second dimension of the array to read               
    Real(kind=4), dimension(n1,n2), intent(inout), target :: data  !< array                                                

    Type(c_ptr) :: ptr
    Integer(hid_t) :: dset_id                               ! id of the dataset
    Integer(hid_t) :: h5_kind                               ! HDF5 integer type

    ptr = C_LOC(data(1,1))
    ! hdf5 type corresponding to the integer type to read
    h5_kind = h5kind_to_type(4,H5_REAL_KIND)

    ! open the dataset
    Call h5dopen_f(id, name, dset_id, h5err)
    
    ! read the dataset
    Call h5dread_f(dset_id, h5_kind, ptr, h5err)

    ! close the dataset
    Call h5dclose_f(dset_id, h5err)

  End Subroutine hdf5_read_real4_2D


  !=======================================================================
  !> Read a 2-D real(kind=8) array from a serial file
  Subroutine hdf5_read_real8_2D(id, name, n1, n2, data)

    Implicit None

    Integer(hid_t), intent(in) :: id                               !< id of the file/group where the dataset will be read 
    Character(len=16), intent(in) :: name                          !< name of the dataset                                 
    Integer(kind=4), intent(in) :: n1                              !< first dimension of the array to read                
    Integer(kind=4), intent(in) :: n2                              !< second dimension of the array to read               
    Real(kind=8), dimension(n1,n2), intent(inout), target :: data  !< array                                                

    Type(c_ptr) :: ptr
    Integer(hid_t) :: dset_id                               ! id of the dataset
    Integer(hid_t) :: h5_kind                               ! HDF5 integer type

    ptr = C_LOC(data(1,1))
    ! hdf5 type corresponding to the integer type to read
    h5_kind = h5kind_to_type(8,H5_REAL_KIND)

    ! open the dataset
    Call h5dopen_f(id, name, dset_id, h5err)
    
    ! read the dataset
    Call h5dread_f(dset_id, h5_kind, ptr, h5err)

    ! close the dataset
    Call h5dclose_f(dset_id, h5err)

  End Subroutine hdf5_read_real8_2D


  !===============================================!
  !         MPI DATA INPUT ROUTINES
  !===============================================!


  !=======================================================================
  !> Read a 1-D integer4 array from a parallel file
  !! The array is distributed along the 2nd dimension
  Subroutine hdf5_read_mpi_int4_1D(id, name, n1, data, comm)

    Implicit None

    Integer(hid_t), intent(in) :: id                              !< id of the file/group where the dataset will be read 
    Character(len=16), intent(in) :: name                         !< name of the dataset                                 
    Integer(kind=4), intent(in) :: n1                             !< dimension of the local array to read                      
    Integer(kind=4), dimension(n1), intent(inout), target :: data !< array                                               
    Integer(kind=4), intent(in) :: comm                           !< MPI communicator used                               

    Type(c_ptr) :: ptr                                         ! Pointer to the array
    Integer(hid_t) :: dset_id                                  ! id of the dataset
    Integer(hid_t) :: plist_id                                 ! id of the dataset
    Integer(hid_t) :: dspace_id                                ! id of the dataspace
    Integer(hid_t) :: mem_id                                   ! id of the dataspace
    Integer(hid_t) :: h5_kind                                  ! HDF5 real type
    Integer(hsize_t), dimension(1) :: local_dims               ! local array dimensions
    Integer(hsize_t), dimension(1) :: offset                   ! offset used for the parallel read
    Integer(kind=4), dimension(:), allocatable :: ntab         ! tmp array used for offset computation
    Integer(kind=4) :: nproc, idproc                           ! number of process and id of local process
    Integer(kind=4) :: i

    ptr = C_LOC(data(1))
    
    ! hdf5 type corresponding to the integer type to read
    h5_kind = h5kind_to_type(4,H5_INTEGER_KIND)

    local_dims(1) = n1

    Call Mpi_Comm_Size(comm, nproc, mpierr)
    Call Mpi_Comm_Rank(comm, idproc, mpierr)
    
    Allocate(ntab(nproc))

    ! Offset has to be computed
    Call Mpi_Allgather(n1, 1, Mpi_Integer, ntab, 1, Mpi_Integer, comm, mpierr)
    
    offset(1) = 0
    
    Do i = 1, idproc
       offset(1) = offset(1) + ntab(i)
    End Do

    ! create memory dataspace
    Call h5screate_simple_f(1, local_dims, mem_id, h5err)
    ! open dataset
    Call h5dopen_f(id, name, dset_id, h5err)
    ! read the file dataspace
    Call h5dget_space_f(dset_id, dspace_id, h5err)

    Call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, local_dims, h5err)

    ! Create property list for collective dataset write
    Call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, h5err) 
    Call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, h5err)

    ! read the dataset
    Call h5dread_f(dset_id,h5_kind,ptr,h5err,&
         mem_space_id=mem_id,file_space_id=dspace_id,xfer_prp=plist_id)
    
    Call h5pclose_f(plist_id, h5err)
    ! close the file dataspace
    Call h5sclose_f(dspace_id, h5err)
    ! close the dataset
    Call h5sclose_f(mem_id, h5err)
    ! close the dataset
    Call h5dclose_f(dset_id, h5err)

    Deallocate(ntab)

  End Subroutine hdf5_read_mpi_int4_1D


  !=======================================================================
  !> Read a 2-D integer4 array from a parallel file
  !! The array is distributed along the 2nd dimension
  Subroutine hdf5_read_mpi_int4_2D(id, name, n1, n2, data, comm)

    Implicit None

    Integer(hid_t), intent(in) :: id                                  !< id of the file/group where the dataset will be read 
    Character(len=16), intent(in) :: name                             !< name of the dataset                                 
    Integer(kind=4), intent(in) :: n1                                 !< first dimension of the local array to read          
    Integer(kind=4), intent(in) :: n2                                 !< second dimension of the local array to read         
    Integer(kind=4), dimension(n1,n2), intent(inout), target :: data  !< array                                               
    Integer(kind=4), intent(in) :: comm                               !< MPI communicator used                               

    Type(c_ptr) :: ptr                                         ! Pointer to the array
    Integer(hid_t) :: dset_id                                  ! id of the dataset
    Integer(hid_t) :: plist_id                                 ! id of the dataset
    Integer(hid_t) :: dspace_id                                ! id of the dataspace
    Integer(hid_t) :: mem_id                                   ! id of the dataspace
    Integer(hid_t) :: h5_kind                                  ! HDF5 real type
    Integer(hsize_t), dimension(2) :: local_dims               ! local array dimensions
    Integer(hsize_t), dimension(2) :: offset                   ! offset used for the parallel read
    Integer(kind=4), dimension(:), allocatable :: ntab         ! tmp array used for offset computation
    Integer(kind=4) :: nproc, idproc                           ! number of process and id of local process
    Integer(kind=4) :: i


    ptr = C_LOC(data(1,1))
    
    ! hdf5 type corresponding to the integer type to read
    h5_kind = h5kind_to_type(4,H5_INTEGER_KIND)

    local_dims(1) = n1
    local_dims(2) = n2

    Call Mpi_Comm_Size(comm, nproc, mpierr)
    Call Mpi_Comm_Rank(comm, idproc, mpierr)
    
    Allocate(ntab(nproc))

    ! Offset has to be computed
    Call Mpi_Allgather(n2, 1, Mpi_Integer, ntab, 1, Mpi_Integer, comm, mpierr)
    
    offset(1) = 0
    offset(2) = 0
    
    Do i = 1, idproc
       offset(2) = offset(2) + ntab(i)
    End Do

    ! create memory dataspace
    Call h5screate_simple_f(2, local_dims, mem_id, h5err)
    ! open dataset
    Call h5dopen_f(id, name, dset_id, h5err)
    ! read the file dataspace
    Call h5dget_space_f(dset_id, dspace_id, h5err)

    Call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, local_dims, h5err)

    ! Create property list for collective dataset write
    Call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, h5err) 
    Call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, h5err)

    ! read the dataset
    Call h5dread_f(dset_id,h5_kind,ptr,h5err,&
         mem_space_id=mem_id,file_space_id=dspace_id,xfer_prp=plist_id)
    
    Call h5pclose_f(plist_id, h5err)
    ! close the file dataspace
    Call h5sclose_f(dspace_id, h5err)
    ! close the dataset
    Call h5sclose_f(mem_id, h5err)
    ! close the dataset
    Call h5dclose_f(dset_id, h5err)

    Deallocate(ntab)

  End Subroutine hdf5_read_mpi_int4_2D


  !=======================================================================
  !> Read a 1-D integer8 array from a parallel file 
  !! The array is distributed along the 2nd dimension
  Subroutine hdf5_read_mpi_int8_1D(id, name, n1, data, comm)

    Implicit None

    Integer(hid_t), intent(in) :: id                               !< id of the file/group where the dataset will be read 
    Character(len=16), intent(in) :: name                          !< name of the dataset                                 
    Integer(kind=4), intent(in) :: n1                              !< dimension of the local array to read                      
    Integer(kind=8), dimension(n1), intent(inout), target :: data  !< array                                               
    Integer(kind=4), intent(in) :: comm                            !< MPI communicator used                               

    Type(c_ptr) :: ptr                                         ! Pointer to the array
    Integer(hid_t) :: dset_id                                  ! id of the dataset
    Integer(hid_t) :: plist_id                                 ! id of the dataset
    Integer(hid_t) :: dspace_id                                ! id of the dataspace
    Integer(hid_t) :: mem_id                                   ! id of the dataspace
    Integer(hid_t) :: h5_kind                                  ! HDF5 real type
    Integer(hsize_t), dimension(1) :: local_dims               ! local array dimensions
    Integer(hsize_t), dimension(1) :: offset                   ! offset used for the parallel read
    Integer(kind=4), dimension(:), allocatable :: ntab         ! tmp array used for offset computation
    Integer(kind=4) :: nproc, idproc                           ! number of process and id of local process
    Integer(kind=4) :: i


    ptr = C_LOC(data(1))
    
    ! hdf5 type corresponding to the integer type to read
    h5_kind = h5kind_to_type(8,H5_INTEGER_KIND)

    local_dims(1) = n1

    Call Mpi_Comm_Size(comm, nproc, mpierr)
    Call Mpi_Comm_Rank(comm, idproc, mpierr)
    
    Allocate(ntab(nproc))

    ! Offset has to be computed
    Call Mpi_Allgather(n1, 1, Mpi_Integer, ntab, 1, Mpi_Integer, comm, mpierr)
    
    offset(1) = 0
    
    Do i = 1, idproc
       offset(1) = offset(1) + ntab(i)
    End Do

    ! create memory dataspace
    Call h5screate_simple_f(1, local_dims, mem_id, h5err)
    ! open dataset
    Call h5dopen_f(id, name, dset_id, h5err)
    ! read the file dataspace
    Call h5dget_space_f(dset_id, dspace_id, h5err)

    Call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, local_dims, h5err)

    ! Create property list for collective dataset write
    Call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, h5err) 
    Call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, h5err)

    ! read the dataset
    Call h5dread_f(dset_id,h5_kind,ptr,h5err,&
         mem_space_id=mem_id,file_space_id=dspace_id,xfer_prp=plist_id)
    
    Call h5pclose_f(plist_id, h5err)
    ! close the file dataspace
    Call h5sclose_f(dspace_id, h5err)
    ! close the dataset
    Call h5sclose_f(mem_id, h5err)
    ! close the dataset
    Call h5dclose_f(dset_id, h5err)

    Deallocate(ntab)

  End Subroutine hdf5_read_mpi_int8_1D


  !=======================================================================
  !> Read a 2-D integer8 array from a parallel file 
  !! The array is distributed along the 2nd dimension
  Subroutine hdf5_read_mpi_int8_2D(id, name, n1, n2, data, comm)

    Implicit None

    Integer(hid_t), intent(in) :: id                                  !< id of the file/group where the dataset will be read 
    Character(len=16), intent(in) :: name                             !< name of the dataset                                 
    Integer(kind=4), intent(in) :: n1                                 !< first dimension of the local array to read          
    Integer(kind=4), intent(in) :: n2                                 !< second dimension of the local array to read              
    Integer(kind=8), dimension(n1,n2), intent(inout), target :: data  !< array                                               
    Integer(kind=4), intent(in) :: comm                               !< MPI communicator used                               

    Type(c_ptr) :: ptr                                         ! Pointer to the array
    Integer(hid_t) :: dset_id                                  ! id of the dataset
    Integer(hid_t) :: plist_id                                 ! id of the dataset
    Integer(hid_t) :: dspace_id                                ! id of the dataspace
    Integer(hid_t) :: mem_id                                   ! id of the dataspace
    Integer(hid_t) :: h5_kind                                  ! HDF5 real type
    Integer(hsize_t), dimension(2) :: local_dims               ! local array dimensions
    Integer(hsize_t), dimension(2) :: offset                   ! offset used for the parallel read
    Integer(kind=4), dimension(:), allocatable :: ntab         ! tmp array used for offset computation
    Integer(kind=4) :: nproc, idproc                           ! number of process and id of local process
    Integer(kind=4) :: i


    ptr = C_LOC(data(1,1))
    
    ! hdf5 type corresponding to the integer type to read
    h5_kind = h5kind_to_type(8,H5_INTEGER_KIND)

    local_dims(1) = n1
    local_dims(2) = n2

    Call Mpi_Comm_Size(comm, nproc, mpierr)
    Call Mpi_Comm_Rank(comm, idproc, mpierr)
    
    Allocate(ntab(nproc))

    ! Offset has to be computed
    Call Mpi_Allgather(n2, 1, Mpi_Integer, ntab, 1, Mpi_Integer, comm, mpierr)
    
    offset(1) = 0
    offset(2) = 0
    
    Do i = 1, idproc
       offset(2) = offset(2) + ntab(i)
    End Do

    ! create memory dataspace
    Call h5screate_simple_f(2, local_dims, mem_id, h5err)
    ! open dataset
    Call h5dopen_f(id, name, dset_id, h5err)
    ! read the file dataspace
    Call h5dget_space_f(dset_id, dspace_id, h5err)

    Call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, local_dims, h5err)

    ! Create property list for collective dataset write
    Call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, h5err) 
    Call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, h5err)

    ! read the dataset
    Call h5dread_f(dset_id,h5_kind,ptr,h5err,&
         mem_space_id=mem_id,file_space_id=dspace_id,xfer_prp=plist_id)
    
    Call h5pclose_f(plist_id, h5err)
    ! close the file dataspace
    Call h5sclose_f(dspace_id, h5err)
    ! close the dataset
    Call h5sclose_f(mem_id, h5err)
    ! close the dataset
    Call h5dclose_f(dset_id, h5err)

    Deallocate(ntab)

  End Subroutine hdf5_read_mpi_int8_2D


  !=======================================================================
  !> Read a 1-D real4 array from a parallel file
  !! The array is distributed along the 2nd dimension
  Subroutine hdf5_read_mpi_real4_1D(id, name, n1, data, comm)

    Implicit None

    Integer(hid_t), intent(in) :: id                            !< id of the file/group where the dataset will be read 
    Character(len=16), intent(in) :: name                       !< name of the dataset                                 
    Integer(kind=4), intent(in) :: n1                           !< dimension of the local array to read                      
    Real(kind=4), dimension(n1), intent(inout), target :: data  !< array                                               
    Integer(kind=4), intent(in) :: comm                         !< MPI communicator used                               

    Type(c_ptr) :: ptr                                         ! Pointer to the array
    Integer(hid_t) :: dset_id                                  ! id of the dataset
    Integer(hid_t) :: plist_id                                 ! id of the dataset
    Integer(hid_t) :: dspace_id                                ! id of the dataspace
    Integer(hid_t) :: mem_id                                   ! id of the dataspace
    Integer(hid_t) :: h5_kind                                  ! HDF5 real type
    Integer(hsize_t), dimension(1) :: local_dims               ! local array dimensions
    Integer(hsize_t), dimension(1) :: offset                   ! offset used for the parallel read
    Integer(kind=4), dimension(:), allocatable :: ntab         ! tmp array used for offset computation
    Integer(kind=4) :: nproc, idproc                           ! number of process and id of local process
    Integer(kind=4) :: i

    ptr = C_LOC(data(1))
    
    ! hdf5 type corresponding to the integer type to read
    h5_kind = h5kind_to_type(4,H5_REAL_KIND)

    local_dims(1) = n1

    Call Mpi_Comm_Size(comm, nproc, mpierr)
    Call Mpi_Comm_Rank(comm, idproc, mpierr)
    
    Allocate(ntab(nproc))

    ! Offset has to be computed
    Call Mpi_Allgather(n1, 1, Mpi_Integer, ntab, 1, Mpi_Integer, comm, mpierr)
    
    offset(1) = 0
    
    Do i = 1, idproc
       offset(1) = offset(1) + ntab(i)
    End Do

    ! create memory dataspace
    Call h5screate_simple_f(1, local_dims, mem_id, h5err)
    ! open dataset
    Call h5dopen_f(id, name, dset_id, h5err)
    ! read the file dataspace
    Call h5dget_space_f(dset_id, dspace_id, h5err)

    Call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, local_dims, h5err)

    ! Create property list for collective dataset write
    Call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, h5err) 
    Call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, h5err)

    ! read the dataset
    Call h5dread_f(dset_id,h5_kind,ptr,h5err,&
         mem_space_id=mem_id,file_space_id=dspace_id,xfer_prp=plist_id)
    
    Call h5pclose_f(plist_id, h5err)
    ! close the file dataspace
    Call h5sclose_f(dspace_id, h5err)
    ! close the dataset
    Call h5sclose_f(mem_id, h5err)
    ! close the dataset
    Call h5dclose_f(dset_id, h5err)

    Deallocate(ntab)

  End Subroutine hdf5_read_mpi_real4_1D


  !=======================================================================
  !> Read a 2-D real4 array from a parallel file
  !! The array is distributed along the 2nd dimension
  Subroutine hdf5_read_mpi_real4_2D(id, name, n1, n2, data, comm)

    Implicit None

    Integer(hid_t), intent(in) :: id                               !< id of the file/group where the dataset will be read 
    Character(len=16), intent(in) :: name                          !< name of the dataset                                 
    Integer(kind=4), intent(in) :: n1                              !< first dimension of the local array to read          
    Integer(kind=4), intent(in) :: n2                              !< second dimension of the local array to read         
    Real(kind=4), dimension(n1,n2), intent(inout), target :: data  !< array                                               
    Integer(kind=4), intent(in) :: comm                            !< MPI communicator used                               

    Type(c_ptr) :: ptr                                         ! Pointer to the array
    Integer(hid_t) :: dset_id                                  ! id of the dataset
    Integer(hid_t) :: plist_id                                 ! id of the dataset
    Integer(hid_t) :: dspace_id                                ! id of the dataspace
    Integer(hid_t) :: mem_id                                   ! id of the dataspace
    Integer(hid_t) :: h5_kind                                  ! HDF5 real type
    Integer(hsize_t), dimension(2) :: local_dims               ! local array dimensions
    Integer(hsize_t), dimension(2) :: offset                   ! offset used for the parallel read
    Integer(kind=4), dimension(:), allocatable :: ntab         ! tmp array used for offset computation
    Integer(kind=4) :: nproc, idproc                           ! number of process and id of local process
    Integer(kind=4) :: i


    ptr = C_LOC(data(1,1))
    
    ! hdf5 type corresponding to the integer type to read
    h5_kind = h5kind_to_type(4,H5_REAL_KIND)

    local_dims(1) = n1
    local_dims(2) = n2

    Call Mpi_Comm_Size(comm, nproc, mpierr)
    Call Mpi_Comm_Rank(comm, idproc, mpierr)
    
    Allocate(ntab(nproc))

    ! Offset has to be computed
    Call Mpi_Allgather(n2, 1, Mpi_Integer, ntab, 1, Mpi_Integer, comm, mpierr)
    
    offset(1) = 0
    offset(2) = 0
    
    Do i = 1, idproc
       offset(2) = offset(2) + ntab(i)
    End Do

    ! create memory dataspace
    Call h5screate_simple_f(2, local_dims, mem_id, h5err)
    ! open dataset
    Call h5dopen_f(id, name, dset_id, h5err)
    ! read the file dataspace
    Call h5dget_space_f(dset_id, dspace_id, h5err)

    Call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, local_dims, h5err)

    ! Create property list for collective dataset write
    Call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, h5err) 
    Call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, h5err)

    ! read the dataset
    Call h5dread_f(dset_id,h5_kind,ptr,h5err,&
         mem_space_id=mem_id,file_space_id=dspace_id,xfer_prp=plist_id)
    
    Call h5pclose_f(plist_id, h5err)
    ! close the file dataspace
    Call h5sclose_f(dspace_id, h5err)
    ! close the dataset
    Call h5sclose_f(mem_id, h5err)
    ! close the dataset
    Call h5dclose_f(dset_id, h5err)

    Deallocate(ntab)

  End Subroutine hdf5_read_mpi_real4_2D


  !=======================================================================
  !> Read a 1-D real8 array from a parallel file
  !! The array is distributed along the 2nd dimension
  Subroutine hdf5_read_mpi_real8_1D(id, name, n1, data, comm)

    Implicit None

    Integer(hid_t), intent(in) :: id                            !< id of the file/group where the dataset will be read 
    Character(len=16), intent(in) :: name                       !< name of the dataset                                 
    Integer(kind=4), intent(in) :: n1                           !< dimension of the local array to read                      
    Real(kind=8), dimension(n1), intent(inout), target :: data  !< array                                               
    Integer(kind=4), intent(in) :: comm                         !< MPI communicator used                               

    Type(c_ptr) :: ptr                                         ! Pointer to the array
    Integer(hid_t) :: dset_id                                  ! id of the dataset
    Integer(hid_t) :: plist_id                                 ! id of the dataset
    Integer(hid_t) :: dspace_id                                ! id of the dataspace
    Integer(hid_t) :: mem_id                                   ! id of the dataspace
    Integer(hid_t) :: h5_kind                                  ! HDF5 real type
    Integer(hsize_t), dimension(1) :: local_dims               ! local array dimensions
    Integer(hsize_t), dimension(1) :: offset                   ! offset used for the parallel read
    Integer(kind=4), dimension(:), allocatable :: ntab         ! tmp array used for offset computation
    Integer(kind=4) :: nproc, idproc                           ! number of process and id of local process
    Integer(kind=4) :: i


    ptr = C_LOC(data(1))
    
    ! hdf5 type corresponding to the integer type to read
    h5_kind = h5kind_to_type(8,H5_REAL_KIND)

    local_dims(1) = n1

    Call Mpi_Comm_Size(comm, nproc, mpierr)
    Call Mpi_Comm_Rank(comm, idproc, mpierr)
    
    Allocate(ntab(nproc))

    ! Offset has to be computed
    Call Mpi_Allgather(n1, 1, Mpi_Integer, ntab, 1, Mpi_Integer, comm, mpierr)
    
    offset(1) = 0
    
    Do i = 1, idproc
       offset(1) = offset(1) + ntab(i)
    End Do

    ! create memory dataspace
    Call h5screate_simple_f(1, local_dims, mem_id, h5err)
    ! open dataset
    Call h5dopen_f(id, name, dset_id, h5err)
    ! read the file dataspace
    Call h5dget_space_f(dset_id, dspace_id, h5err)

    Call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, local_dims, h5err)

    ! Create property list for collective dataset write
    Call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, h5err) 
    Call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, h5err)

    ! read the dataset
    Call h5dread_f(dset_id,h5_kind,ptr,h5err,&
         mem_space_id=mem_id,file_space_id=dspace_id,xfer_prp=plist_id)
    
    Call h5pclose_f(plist_id, h5err)
    ! close the file dataspace
    Call h5sclose_f(dspace_id, h5err)
    ! close the dataset
    Call h5sclose_f(mem_id, h5err)
    ! close the dataset
    Call h5dclose_f(dset_id, h5err)

    Deallocate(ntab)

  End Subroutine hdf5_read_mpi_real8_1D


  !=======================================================================
  !> Read a 2-D real8 array from a parallel file
  !! The array is distributed along the 2nd dimension
  Subroutine hdf5_read_mpi_real8_2D(id, name, n1, n2, data, comm)

    Implicit None

    Integer(hid_t), intent(in) :: id                                !< id of the file/group where the dataset will be read 
    Character(len=16), intent(in) :: name                           !< name of the dataset                                 
    Integer(kind=4), intent(in) :: n1                               !< first dimension of the local array to read          
    Integer(kind=4), intent(in) :: n2                               !< second dimension of the local array to read         
    Real(kind=8), dimension(n1,n2), intent(inout), target :: data   !< array                                               
    Integer(kind=4), intent(in) :: comm                             !< MPI communicator used                               

    Type(c_ptr) :: ptr                                         ! Pointer to the array
    Integer(hid_t) :: dset_id                                  ! id of the dataset
    Integer(hid_t) :: plist_id                                 ! id of the dataset
    Integer(hid_t) :: dspace_id                                ! id of the dataspace
    Integer(hid_t) :: mem_id                                   ! id of the dataspace
    Integer(hid_t) :: h5_kind                                  ! HDF5 real type
    Integer(hsize_t), dimension(2) :: local_dims               ! local array dimensions
    Integer(hsize_t), dimension(2) :: offset                   ! offset used for the parallel read
    Integer(kind=4), dimension(:), allocatable :: ntab         ! tmp array used for offset computation
    Integer(kind=4) :: nproc, idproc                           ! number of process and id of local process
    Integer(kind=4) :: i


    ptr = C_LOC(data(1,1))
    
    ! hdf5 type corresponding to the real type to read
    h5_kind = h5kind_to_type(8,H5_REAL_KIND)

    local_dims(1) = n1
    local_dims(2) = n2

    Call Mpi_Comm_Size(comm, nproc, mpierr)
    Call Mpi_Comm_Rank(comm, idproc, mpierr)
    
    Allocate(ntab(nproc))

    ! Offset has to be computed
    Call Mpi_Allgather(n2, 1, Mpi_Integer, ntab, 1, Mpi_Integer, comm, mpierr)
    
    offset(1) = 0
    offset(2) = 0
    
    Do i = 1, idproc
       offset(2) = offset(2) + ntab(i)
    End Do

    ! create memory dataspace
    Call h5screate_simple_f(2, local_dims, mem_id, h5err)
    ! open dataset
    Call h5dopen_f(id, name, dset_id, h5err)
    ! read the file dataspace
    Call h5dget_space_f(dset_id, dspace_id, h5err)

    Call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, local_dims, h5err)

    ! Create property list for collective dataset write
    Call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, h5err) 
    Call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, h5err)

    ! read the dataset
    Call h5dread_f(dset_id,h5_kind,ptr,h5err,&
         mem_space_id=mem_id,file_space_id=dspace_id,xfer_prp=plist_id)
    
    Call h5pclose_f(plist_id, h5err)
    ! close the file dataspace
    Call h5sclose_f(dspace_id, h5err)
    ! close the dataset
    Call h5sclose_f(mem_id, h5err)
    ! close the dataset
    Call h5dclose_f(dset_id, h5err)

    Deallocate(ntab)

  End Subroutine hdf5_read_mpi_real8_2D
  

  !===============================================!
  !       SERIAL ATTRIBUTE OUTPUT ROUTINES
  !===============================================!


  !=======================================================================
  !> Write an integer4 as attribute
  Subroutine hdf5_write_int4_attr0D(id, name, data)
    
    Implicit None

    Integer(hid_t), intent(in) :: id                           !< id of the file/group where the attribute will be written
    Character(len=16), intent(in) :: name                      !< name of the attribute                                   
    Integer(kind=4), intent(in) :: data                        !< attribute value                                         

    Integer(hsize_t), dimension(1) :: dim                      ! Local dataset dimensions
    Integer :: rank 
    Integer(hid_t) :: aspace_id                                ! id of the dataspace
    Integer(hid_t) :: atype_id                                 ! id of the dataspace
    Integer(hid_t) :: attr_id                                  ! id of the dataspace
    Integer(size_t) :: alen                                    ! Attribute length in bytes

    
    dim(1) = 1
    rank = 1
    ! Create scalar data space for the attribute.
    Call h5screate_simple_f(rank, dim, aspace_id, h5err)

    ! Create datatype for the attribute.
    alen = 4
    Call h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, h5err)
    Call h5tset_size_f(atype_id, alen, h5err)
    ! Create dataset attribute.
    Call h5acreate_f(id, name, atype_id, aspace_id, attr_id, h5err)  
    ! Write the attribute data.
    Call h5awrite_f(attr_id, atype_id, data, dim, h5err)
    ! Close the attribute.
    Call h5aclose_f(attr_id, h5err)
    ! Terminate access to the data space.
    Call h5sclose_f(aspace_id, h5err)
    Call h5tclose_f(atype_id, h5err)

  End Subroutine hdf5_write_int4_attr0D


  !=======================================================================
  !> Write a 1-D array of integer4 as attribute
  Subroutine hdf5_write_int4_attr1D(id, name, n1, data)
    
    Implicit None

    Integer(hid_t), intent(in) :: id                           !< id of the file/group where the attribute will be written
    Character(len=16), intent(in) :: name                      !< name of the attribute                                   
    Integer(kind=4), intent(in) :: n1                          !< dimension of the attribute array                        
    Integer(kind=4), dimension(n1), intent(in) :: data         !< attribute value                                         

    Integer(hsize_t), dimension(1) :: dim                      ! Local dataset dimensions
    Integer :: rank 
    Integer(hid_t) :: aspace_id                                ! id of the dataspace
    Integer(hid_t) :: atype_id                                 ! id of the dataspace
    Integer(hid_t) :: attr_id                                  ! id of the dataspace
    Integer(size_t) :: alen                                    ! Attribute length in bytes

    
    dim(1) = n1
    rank = 1
    Call h5screate_simple_f(rank, dim, aspace_id, h5err)

    alen = 4
    Call h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, h5err)
    Call h5tset_size_f(atype_id, alen, h5err)
    Call h5acreate_f(id, name, atype_id, aspace_id, attr_id, h5err)  
    Call h5awrite_f(attr_id, atype_id, data, dim, h5err)
    Call h5aclose_f(attr_id, h5err)
    Call h5sclose_f(aspace_id, h5err)
    Call h5tclose_f(atype_id, h5err)

  End Subroutine hdf5_write_int4_attr1D


  !=======================================================================
  ! Write a 2-D array of integer4 as attribute
  Subroutine hdf5_write_int4_attr2D(id, name, n1, n2, data)
    
    Implicit None

    Integer(hid_t), intent(in) :: id                           !< id of the file/group where the attribute will be written
    Character(len=16), intent(in) :: name                      !< name of the attribute                                   
    Integer(kind=4), intent(in) :: n1                          !< first dimension of the attribute array                        
    Integer(kind=4), intent(in) :: n2                          !< second dimension of the attribute array                 
    Integer(kind=4), dimension(n1,n2), intent(in) :: data      !< attribute value                                         

    Integer(hsize_t), dimension(2) :: dim                      ! Local dataset dimensions
    Integer :: rank 
    Integer(hid_t) :: aspace_id                                ! id of the dataspace
    Integer(hid_t) :: atype_id                                 ! id of the dataspace
    Integer(hid_t) :: attr_id                                  ! id of the dataspace
    Integer(size_t) :: alen                                    ! Attribute length in bytes

    
    dim(1) = n1
    dim(2) = n2
    rank = 2
    Call h5screate_simple_f(rank, dim, aspace_id, h5err)

    alen = 4
    Call h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, h5err)
    Call h5tset_size_f(atype_id, alen, h5err)
    Call h5acreate_f(id, name, atype_id, aspace_id, attr_id, h5err)  
    Call h5awrite_f(attr_id, atype_id, data, dim, h5err)
    Call h5aclose_f(attr_id, h5err)
    Call h5sclose_f(aspace_id, h5err)
    Call h5tclose_f(atype_id, h5err)

  End Subroutine hdf5_write_int4_attr2D

  !=======================================================================
  !> Write a real4 as attribute
  Subroutine hdf5_write_real4_attr0D(id, name, data)
    
    Implicit None

    Integer(hid_t), intent(in) :: id                           !< id of the file/group where the attribute will be written
    Character(len=16), intent(in) :: name                      !< name of the attribute                                   
    Real(kind=4), intent(in) :: data                           !< attribute value                                         

    Integer(hsize_t), dimension(1) :: dim                      ! Local dataset dimensions
    Integer :: rank 
    Integer(hid_t) :: aspace_id                                ! id of the dataspace
    Integer(hid_t) :: atype_id                                 ! id of the dataspace
    Integer(hid_t) :: attr_id                                  ! id of the dataspace
    Integer(size_t) :: alen                                    ! Attribute length in bytes

    
    dim(1) = 1
    rank = 1
    Call h5screate_simple_f(rank, dim, aspace_id, h5err)

    alen = 4
    Call h5tcopy_f(H5T_NATIVE_REAL, atype_id, h5err)
    Call h5tset_size_f(atype_id, alen, h5err)
    Call h5acreate_f(id, name, atype_id, aspace_id, attr_id, h5err)  
    Call h5awrite_f(attr_id, atype_id, data, dim, h5err)
    Call h5aclose_f(attr_id, h5err)
    Call h5sclose_f(aspace_id, h5err)
    Call h5tclose_f(atype_id, h5err)

  End Subroutine hdf5_write_real4_attr0D


  !=======================================================================
  !> Write a 1-D array of real4 as attribute
  Subroutine hdf5_write_real4_attr1D(id, name, n1, data)
    
    Implicit None

    Integer(hid_t), intent(in) :: id                           !< id of the file/group where the attribute will be written
    Character(len=16), intent(in) :: name                      !< name of the attribute                                   
    Integer(kind=4), intent(in) :: n1                          !< dimension of the attribute array                        
    Real(kind=4), dimension(n1), intent(in) :: data            !< attribute value                                         

    Integer(hsize_t), dimension(1) :: dim                      ! Local dataset dimensions
    Integer :: rank 
    Integer(hid_t) :: aspace_id                                ! id of the dataspace
    Integer(hid_t) :: atype_id                                 ! id of the dataspace
    Integer(hid_t) :: attr_id                                  ! id of the dataspace
    Integer(size_t) :: alen                                    ! Attribute length in bytes

    
    dim(1) = n1
    rank = 1
    Call h5screate_simple_f(rank, dim, aspace_id, h5err)

    alen = 4
    Call h5tcopy_f(H5T_NATIVE_REAL, atype_id, h5err)
    Call h5tset_size_f(atype_id, alen, h5err)
    Call h5acreate_f(id, name, atype_id, aspace_id, attr_id, h5err)  
    Call h5awrite_f(attr_id, atype_id, data, dim, h5err)
    Call h5aclose_f(attr_id, h5err)
    Call h5sclose_f(aspace_id, h5err)
    Call h5tclose_f(atype_id, h5err)

  End Subroutine hdf5_write_real4_attr1D


  !=======================================================================
  ! Write a 2-D array of real4 as attribute
  Subroutine hdf5_write_real4_attr2D(id, name, n1, n2, data)
    
    Implicit None

    Integer(hid_t), intent(in) :: id                           !< id of the file/group where the attribute will be written
    Character(len=16), intent(in) :: name                      !< name of the attribute                                   
    Integer(kind=4), intent(in) :: n1                          !< first dimension of the attribute array                  
    Integer(kind=4), intent(in) :: n2                          !< second dimension of the attribute array                 
    Real(kind=4), dimension(n1,n2), intent(in) :: data         !< attribute value                                         

    Integer(hsize_t), dimension(2) :: dim                      ! Local dataset dimensions
    Integer :: rank 
    Integer(hid_t) :: aspace_id                                ! id of the dataspace
    Integer(hid_t) :: atype_id                                 ! id of the dataspace
    Integer(hid_t) :: attr_id                                  ! id of the dataspace
    Integer(size_t) :: alen                                    ! Attribute length in bytes

    
    dim(1) = n1
    dim(2) = n2
    rank = 2
    Call h5screate_simple_f(rank, dim, aspace_id, h5err)

    alen = 4
    Call h5tcopy_f(H5T_NATIVE_REAL, atype_id, h5err)
    Call h5tset_size_f(atype_id, alen, h5err)
    Call h5acreate_f(id, name, atype_id, aspace_id, attr_id, h5err)  
    Call h5awrite_f(attr_id, atype_id, data, dim, h5err)
    Call h5aclose_f(attr_id, h5err)
    Call h5sclose_f(aspace_id, h5err)
    Call h5tclose_f(atype_id, h5err)

  End Subroutine hdf5_write_real4_attr2D


  !=======================================================================
  !> Write a real8 as attribute
  Subroutine hdf5_write_real8_attr0D(id, name, data)
    
    Implicit None

    Integer(hid_t), intent(in) :: id                           !< id of the file/group where the attribute will be written
    Character(len=16), intent(in) :: name                      !< name of the attribute                                   
    Real(kind=8), intent(in) :: data                           !< attribute value                                         

    Integer(hsize_t), dimension(1) :: dim                      ! Local dataset dimensions
    Integer :: rank 
    Integer(hid_t) :: aspace_id                                ! id of the dataspace
    Integer(hid_t) :: atype_id                                 ! id of the dataspace
    Integer(hid_t) :: attr_id                                  ! id of the dataspace
    Integer(size_t) :: alen                                    ! Attribute length in bytes

    
    dim(1) = 1
    rank = 1
    Call h5screate_simple_f(rank, dim, aspace_id, h5err)

    alen = 8
    Call h5tcopy_f(H5T_NATIVE_DOUBLE, atype_id, h5err)
    Call h5tset_size_f(atype_id, alen, h5err)
    Call h5acreate_f(id, name, atype_id, aspace_id, attr_id, h5err)  
    Call h5awrite_f(attr_id, atype_id, data, dim, h5err)
    Call h5aclose_f(attr_id, h5err)
    Call h5sclose_f(aspace_id, h5err)
    Call h5tclose_f(atype_id, h5err)

  End Subroutine hdf5_write_real8_attr0D


  !=======================================================================
  !> Write a 1-D array of real8 as attribute
  Subroutine hdf5_write_real8_attr1D(id, name, n1, data)
    
    Implicit None

    Integer(hid_t), intent(in) :: id                           !< id of the file/group where the attribute will be written
    Character(len=16), intent(in) :: name                      !< name of the attribute                                   
    Integer(kind=4), intent(in) :: n1                          !< dimension of the attribute array                        
    Real(kind=8), dimension(n1), intent(in) :: data            !< attribute value                                         

    Integer(hsize_t), dimension(1) :: dim                      ! Local dataset dimensions
    Integer :: rank 
    Integer(hid_t) :: aspace_id                                ! id of the dataspace
    Integer(hid_t) :: atype_id                                 ! id of the dataspace
    Integer(hid_t) :: attr_id                                  ! id of the dataspace
    Integer(size_t) :: alen                                    ! Attribute length in bytes

    
    dim(1) = n1
    rank = 1
    Call h5screate_simple_f(rank, dim, aspace_id, h5err)

    alen = 8
    Call h5tcopy_f(H5T_NATIVE_DOUBLE, atype_id, h5err)
    Call h5tset_size_f(atype_id, alen, h5err)
    Call h5acreate_f(id, name, atype_id, aspace_id, attr_id, h5err)  
    Call h5awrite_f(attr_id, atype_id, data, dim, h5err)
    Call h5aclose_f(attr_id, h5err)
    Call h5sclose_f(aspace_id, h5err)
    Call h5tclose_f(atype_id, h5err)

  End Subroutine hdf5_write_real8_attr1D


  !=======================================================================
  ! Write a 2-D array of real8 as attribute
  Subroutine hdf5_write_real8_attr2D(id, name, n1, n2, data)
    
    Implicit None

    Integer(hid_t), intent(in) :: id                           !< id of the file/group where the attribute will be written
    Character(len=16), intent(in) :: name                      !< name of the attribute                                   
    Integer(kind=4), intent(in) :: n1                          !< first dimension of the attribute array                  
    Integer(kind=4), intent(in) :: n2                          !< second dimension of the attribute array                 
    Real(kind=8), dimension(n1,n2), intent(in) :: data         !< attribute value                                         

    Integer(hsize_t), dimension(2) :: dim                      ! Local dataset dimensions
    Integer :: rank 
    Integer(hid_t) :: aspace_id                                ! id of the dataspace
    Integer(hid_t) :: atype_id                                 ! id of the dataspace
    Integer(hid_t) :: attr_id                                  ! id of the dataspace
    Integer(size_t) :: alen                                    ! Attribute length in bytes

    
    dim(1) = n1
    dim(2) = n2
    rank = 2
    Call h5screate_simple_f(rank, dim, aspace_id, h5err)

    alen = 8
    Call h5tcopy_f(H5T_NATIVE_DOUBLE, atype_id, h5err)
    Call h5tset_size_f(atype_id, alen, h5err)
    Call h5acreate_f(id, name, atype_id, aspace_id, attr_id, h5err)  
    Call h5awrite_f(attr_id, atype_id, data, dim, h5err)
    Call h5aclose_f(attr_id, h5err)
    Call h5sclose_f(aspace_id, h5err)
    Call h5tclose_f(atype_id, h5err)

  End Subroutine hdf5_write_real8_attr2D


  !=======================================================================
  !> Write a string as attribute
  Subroutine hdf5_write_char_attr(id, name, data)
    
    Implicit None

    Integer(hid_t), intent(in) :: id                           !< id of the file/group where the attribute will be written
    Character(len=16), intent(in) :: name                      !< name of the attribute                                   
    Character(*), intent(in) :: data                           !< attribute value                                         

    Integer(hsize_t), dimension(1) :: dim                      ! Local dataset dimensions
    Integer(hid_t) :: aspace_id                                ! id of the dataspace
    Integer(hid_t) :: atype_id                                 ! id of the dataspace
    Integer(hid_t) :: attr_id                                  ! id of the dataspace
    Integer(size_t) :: alen                                    ! Attribute length in bytes

    
    dim(1) = 1
    alen = len(data)
    Call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5err)
    Call h5tset_size_f(atype_id, alen, h5err)
    Call h5screate_f(H5S_SCALAR_F, aspace_id, h5err)
    Call h5acreate_f(id, name, atype_id, aspace_id, attr_id, h5err)
    Call h5awrite_f(attr_id, atype_id, data, dim, h5err)
    Call h5aclose_f(attr_id, h5err)
    Call h5sclose_f(aspace_id, h5err)
    Call h5tclose_f(atype_id, h5err)

  End Subroutine hdf5_write_char_attr


  !===============================================!
  !         ATTRIBUTE INPUT ROUTINES
  !===============================================!

  !=======================================================================
  !> Read an integer4 attribute in a hdf5 file
  Subroutine hdf5_read_int4_attr0D(id, name, attr)

    Implicit None

    Integer(kind=hid_t), intent(in) :: id              !< id of the file/group where the attribute will be read
    Character(len=16), intent(in) :: name              !< name of the attribute                              
    Integer(kind=4), intent(inout) :: attr             !< attribute value                                    

    Integer(kind=hsize_t), dimension(1) :: dim1D       !< dimension of the attribute: here = 1
    Integer(kind=hid_t) :: attr_id                     !< attribute id
    Integer(kind=hid_t) :: atype_id                    !< attribute type id

    dim1D(1) = 1

    Call h5aopen_name_f(id, name, attr_id, h5err) 
    Call h5aget_type_f(attr_id, atype_id, h5err) 
    Call h5aread_f(attr_id, atype_id, attr, dim1D, h5err)
    Call h5aclose_f(attr_id, h5err)

  End Subroutine hdf5_read_int4_attr0D


  !=======================================================================
  !> Read a 1-D integer4 attribute in a hdf5 file
  Subroutine hdf5_read_int4_attr1D(id, name, n1, attr)

    Implicit None

    Integer(kind=hid_t), intent(in) :: id                   !< id of the file/group where the attribute will be read
    Character(len=16), intent(in) :: name                   !< name of the attribute
    Integer(kind=4), intent(in) :: n1                       !< dimension of the attribute array
    Integer(kind=4), dimension(n1), intent(inout) :: attr   !< attribute value

    Integer(kind=hsize_t), dimension(1) :: dim1D            !< dimension of the attribute
    Integer(kind=hid_t) :: attr_id                          !< attribute id
    Integer(kind=hid_t) :: atype_id                         !< attribute type id

    dim1D(1) = n1

    Call h5aopen_name_f(id, name, attr_id, h5err) 
    Call h5aget_type_f(attr_id, atype_id, h5err) 
    Call h5aread_f(attr_id, atype_id, attr, dim1D, h5err)
    Call h5aclose_f(attr_id, h5err)

  End Subroutine hdf5_read_int4_attr1D


  !=======================================================================
  !> Read a 2-D integer4 attribute in a hdf5 file
  Subroutine hdf5_read_int4_attr2D(id, name, n1, n2, attr)

    Implicit None

    Integer(kind=hid_t), intent(in) :: id                      !< id of the file/group where the attribute will be read 
    Character(len=16), intent(in) :: name                      !< name of the attribute
    Integer(kind=4), intent(in) :: n1                          !< first dimension of the attribute array
    Integer(kind=4), intent(in) :: n2                          !< second dimension of the attribute array
    Integer(kind=4), dimension(n1,n2), intent(inout) :: attr   !< attribute value

    Integer(kind=hsize_t), dimension(2) :: dim2D               !< dimensions of the attribute
    Integer(kind=hid_t) :: attr_id                             !< attribute id
    Integer(kind=hid_t) :: atype_id                            !< attribute type id

    dim2D(1) = n1
    dim2D(2) = n2

    Call h5aopen_name_f(id, name, attr_id, h5err) 
    Call h5aget_type_f(attr_id, atype_id, h5err) 
    Call h5aread_f(attr_id, atype_id, attr, dim2D, h5err)
    Call h5aclose_f(attr_id, h5err)

  End Subroutine hdf5_read_int4_attr2D


  !=======================================================================
  !> Read a real4 attribute in a hdf5 file
  Subroutine hdf5_read_real4_attr0D(id, name, attr)

    Implicit None

    Integer(kind=hid_t), intent(in) :: id            !< id of the file/group where the attribute will be read
    Character(len=16), intent(in) :: name            !< name of the attribute
    Real(kind=4), intent(inout) :: attr              !< attribute value

    Integer(kind=hsize_t), dimension(1) :: dim1D     !< dimension of the attribute: here = 1
    Integer(kind=hid_t) :: attr_id                   !< attribute id
    Integer(kind=hid_t) :: atype_id                  !< attribute type id

    dim1D(1) = 1

    Call h5aopen_name_f(id, name, attr_id, h5err) 
    Call h5aget_type_f(attr_id, atype_id, h5err) 
    Call h5aread_f(attr_id, atype_id, attr, dim1D, h5err)
    Call h5aclose_f(attr_id, h5err)

  End Subroutine hdf5_read_real4_attr0D


  !=======================================================================
  !> Read a 1-D real4 attribute in a hdf5 file
  Subroutine hdf5_read_real4_attr1D(id, name, n1, attr)

    Implicit None

    Integer(kind=hid_t), intent(in) :: id                 !< id of the file/group where the attribute will be read
    Character(len=16), intent(in) :: name                 !< name of the attribute
    Integer(kind=4), intent(in) :: n1                     !< dimension of the attribute array
    Real(kind=4), dimension(n1), intent(inout) :: attr    !< attribute value

    Integer(kind=hsize_t), dimension(1) :: dim1D          !< dimension of the attribute
    Integer(kind=hid_t) :: attr_id                        !< attribute id
    Integer(kind=hid_t) :: atype_id                       !< attribute type id

    dim1D(1) = n1

    Call h5aopen_name_f(id, name, attr_id, h5err) 
    Call h5aget_type_f(attr_id, atype_id, h5err) 
    Call h5aread_f(attr_id, atype_id, attr, dim1D, h5err)
    Call h5aclose_f(attr_id, h5err)

  End Subroutine hdf5_read_real4_attr1D


  !=======================================================================
  !> Read a 2-D real4 attribute in a hdf5 file
  Subroutine hdf5_read_real4_attr2D(id, name, n1, n2, attr)

    Implicit None

    Integer(kind=hid_t), intent(in) :: id                   !< id of the file/group where the attribute will be read
    Character(len=16), intent(in) :: name                   !< name of the attribute
    Integer(kind=4), intent(in) :: n1                       !< first dimension of the attribute array
    Integer(kind=4), intent(in) :: n2                       !< second dimension of the attribute array
    Real(kind=4), dimension(n1,n2), intent(inout) :: attr   !< attribute value

    Integer(kind=hsize_t), dimension(2) :: dim2D            !< dimensions of the attribute
    Integer(kind=hid_t) :: attr_id                          !< attribute id
    Integer(kind=hid_t) :: atype_id                         !< attribute type id

    dim2D(1) = n1
    dim2D(2) = n2

    Call h5aopen_name_f(id, name, attr_id, h5err) 
    Call h5aget_type_f(attr_id, atype_id, h5err) 
    Call h5aread_f(attr_id, atype_id, attr, dim2D, h5err)
    Call h5aclose_f(attr_id, h5err)

  End Subroutine hdf5_read_real4_attr2D


  !=======================================================================
  !> Read a real8 attribute in a hdf5 file
  Subroutine hdf5_read_real8_attr0D(id, name, attr)

    Implicit None

    Integer(kind=hid_t), intent(in) :: id            !< id of the file/group where the attribute will be read
    Character(len=16), intent(in) :: name            !< name of the attribute
    Real(kind=8), intent(inout) :: attr              !< attribute value

    Integer(kind=hsize_t), dimension(1) :: dim1D     !< dimension of the attribute: here = 1
    Integer(kind=hid_t) :: attr_id                   !< attribute id
    Integer(kind=hid_t) :: atype_id                  !< attribute type id

    dim1D(1) = 1

    Call h5aopen_name_f(id, name, attr_id, h5err) 
    Call h5aget_type_f(attr_id, atype_id, h5err) 
    Call h5aread_f(attr_id, atype_id, attr, dim1D, h5err)
    Call h5aclose_f(attr_id, h5err)

  End Subroutine hdf5_read_real8_attr0D


  !=======================================================================
  !> Read a 1-D real8 attribute in a hdf5 file
  Subroutine hdf5_read_real8_attr1D(id, name, n1, attr)

    Implicit None

    Integer(kind=hid_t), intent(in) :: id                 !< id of the file/group where the attribute will be read
    Character(len=16), intent(in) :: name                 !< name of the attribute
    Integer(kind=4), intent(in) :: n1                     !< dimension of the attribute array
    Real(kind=8), dimension(n1), intent(inout) :: attr    !< attribute value

    Integer(kind=hsize_t), dimension(1) :: dim1D          !< dimension of the attribute
    Integer(kind=hid_t) :: attr_id                        !< attribute id
    Integer(kind=hid_t) :: atype_id                       !< attribute type id

    dim1D(1) = n1

    Call h5aopen_name_f(id, name, attr_id, h5err) 
    Call h5aget_type_f(attr_id, atype_id, h5err) 
    Call h5aread_f(attr_id, atype_id, attr, dim1D, h5err)
    Call h5aclose_f(attr_id, h5err)

  End Subroutine hdf5_read_real8_attr1D


  !=======================================================================
  !> Read a 2-D real4 attribute in a hdf5 file
  Subroutine hdf5_read_real8_attr2D(id, name, n1, n2, attr)

    Implicit None

    Integer(kind=hid_t), intent(in) :: id                   !< id of the file/group where the attribute will be read
    Character(len=16), intent(in) :: name                   !< name of the attribute
    Integer(kind=4), intent(in) :: n1                       !< first dimension of the attribute array
    Integer(kind=4), intent(in) :: n2                       !< second dimension of the attribute array
    Real(kind=8), dimension(n1,n2), intent(inout) :: attr   !< attribute value

    Integer(kind=hsize_t), dimension(2) :: dim2D            !< dimensions of the attribute
    Integer(kind=hid_t) :: attr_id                          !< attribute id
    Integer(kind=hid_t) :: atype_id                         !< attribute type id

    dim2D(1) = n1
    dim2D(2) = n2

    Call h5aopen_name_f(id, name, attr_id, h5err) 
    Call h5aget_type_f(attr_id, atype_id, h5err) 
    Call h5aread_f(attr_id, atype_id, attr, dim2D, h5err)
    Call h5aclose_f(attr_id, h5err)

  End Subroutine hdf5_read_real8_attr2D


  !=======================================================================
  !> Read string attribute in a hdf5 file
  Subroutine hdf5_read_char_attr(id, name, n1, attr)

    Implicit None

    Integer(kind=hid_t), intent(in) :: id                 !< id of the file/group where the attribute will be read
    Character(len=16), intent(in) :: name                 !< name of the attribute
    Integer(kind=4), intent(in) :: n1                     !< dimension of the attribute string
    Character(len=n1), intent(inout) :: attr              !< attribute value

    Integer(kind=hsize_t), dimension(1) :: dim1D          !< dimension of the attribute
    Integer(kind=hid_t) :: attr_id                        !< attribute id
    Integer(kind=hid_t) :: atype_id                       !< attribute type id

    dim1D(1) = n1

    Call h5aopen_name_f(id, name, attr_id, h5err) 
    Call h5aget_type_f(attr_id, atype_id, h5err) 
    Call h5aread_f(attr_id, atype_id, attr, dim1D, h5err)
    Call h5aclose_f(attr_id, h5err)

  End Subroutine hdf5_read_char_attr

End Module modhdf5

