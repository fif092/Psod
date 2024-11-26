!> @file
!! This file contains routines to read HDF5 halo files created with pfof and pfof_cone.

!> This module contains routines to read HDF5 halo files created with pfof and pfof_cone.
!>
!> Authors: F. Roy

Module modreadhalo
  Integer, parameter :: PRIID = 8 !< Precision used for ID; must match the value used to create the HDF5 file you are reading

  Real(kind=8), dimension(:,:), allocatable :: hf_halocomPos, hf_halocomVel
  Integer(kind=PRIID), dimension(:), allocatable :: hf_haloID 
  Integer(kind=4), dimension(:), allocatable :: hf_haloMass
  

  
Contains
  !> Reads a HDF5 halomass file created with pfof
  Subroutine h5readhalomass(filename)
    
    Use modhdf5

    Implicit None
    
    Character(len=400), intent(in) :: filename !< Name of the file
    Character(len=H5STRLEN) :: groupname         ! name of the group
    Integer(hid_t) :: file_id,gr_id
    Character(len=H5STRLEN) :: dname
    Integer(kind=4) :: haloNB

    !! You will probably have to define the following variables in another module
    !! and use this module here as well as in the subroutine from which you want to call h5readhalomass.
!    Real(kind=8), dimension(:,:), allocatable :: hf_halocomPos, hf_halocomVel
!    Integer(kind=PRIID), dimension(:), allocatable :: hf_haloID 
!    Integer(kind=4), dimension(:), allocatable :: hf_haloMass
    !!

    ! Open the file
    Call hdf5_open_file(filename, file_id)
    
    groupname = 'metadata'
    Call hdf5_open_group(file_id,groupname, gr_id)
  
    ! Read the number of halos
    dname = 'nhalo_file'
    Call hdf5_read_attr(gr_id, dname, haloNB)

    Call hdf5_close_group(gr_id)
    
    
    ! Allocate arrays for the position, velocity, ID and mass of the halos
    Allocate(hf_halocomPos(3,haloNB), hf_halocomVel(3,haloNB), hf_haloID(haloNB), hf_haloMass(haloNB))

    groupname = 'data'
    Call hdf5_open_group(file_id,groupname, gr_id)

    ! Read position of the halos
    dname = 'position_halo'
    Call hdf5_read_data(gr_id, dname, 3, haloNB, hf_halocomPos)

    ! Read velocity of the halos
    dname = 'velocity_halo'
    Call hdf5_read_data(gr_id, dname, 3, haloNB, hf_halocomVel)

    ! Read ID of the halos
    dname = 'identity_halo'
    Call hdf5_read_data(gr_id, dname, haloNB, hf_haloID)

    ! Read mass of the halos
    dname = 'npart_halo'
    Call hdf5_read_data(gr_id, dname, haloNB, hf_haloMass)
    
    Call hdf5_close_group(gr_id)

    ! Close file
    Call hdf5_close_file(file_id)

  End Subroutine h5readhalomass

End Module modreadhalo
