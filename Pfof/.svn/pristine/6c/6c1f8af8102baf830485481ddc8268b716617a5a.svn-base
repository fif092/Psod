!> @file
!! This file contains routines to read HDF5 halo files created with pfof and pfof_cone.

!> This module contains routines to read HDF5 halo files created with pfof and pfof_cone.
!>
!> Authors: F. Roy

Module modreadhalo

  Integer, parameter :: PRI = 8 !< Precision used for ID; must match the value used to create the HDF5 file you are reading

Contains
  !> Reads a HDF5 halomass file created with pfof
  Subroutine h5readhalomass(filename)
    
    Use modhdf5

    Implicit None
    
    Character(len=400), intent(in) :: filename !< Name of the file
    Integer(hid_t) :: file_id
    Character(len=16) :: dname
    Integer(kind=4) :: haloNB

    !! You will probably have to define the following variables in another module
    !! and use this module here as well as in the subroutine from which you want to call h5readhalomass.
    Real(kind=8), dimension(:,:), allocatable :: halocomPos, halocomVel
    Integer(kind=PRI), dimension(:), allocatable :: haloID 
    Integer(kind=4), dimension(:), allocatable :: haloMass
    !!

    ! Open the file
    Call hdf5_open_file(filename, file_id)
      
    ! Read the number of halos
    dname = 'haloNB'
    Call hdf5_read_attr(file_id, dname, haloNB)

    ! Allocate arrays for the position, velocity, ID and mass of the halos
    Allocate(halocomPos(3,haloNB), halocomVel(3,haloNB), haloID(haloNB), haloMass(haloNB))

    ! Read position of the halos
    dname = 'halocomPos'
    Call hdf5_read_data(file_id, dname, 3, haloNB, halocomPos)

    ! Read velocity of the halos
    dname = 'halocomVel'
    Call hdf5_read_data(file_id, dname, 3, haloNB, halocomVel)

    ! Read ID of the halos
    dname = 'haloID'
    Call hdf5_read_data(file_id, dname, haloNB, haloID)

    ! Read mass of the halos
    dname = 'haloMass'
    Call hdf5_read_data(file_id, dname, haloNB, haloMass)

    ! Close file
    Call hdf5_close_file(file_id)

  End Subroutine h5readhalomass

End Module modreadhalo
