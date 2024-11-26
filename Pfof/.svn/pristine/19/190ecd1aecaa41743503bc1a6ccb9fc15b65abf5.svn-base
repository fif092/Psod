!==============================================================================
! Project: pFoF
! File: common/src/modreadhalo.f90
! Copyright Fabrice Roy (2015)
! Fabrice.Roy@obspm.fr
! 
! This software is a computer program whose purpose is to detect dark matter
! haloes in cosmological simulation with a parallel Friend of Friend algorithm.
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can  use, 
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info". 
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability. 
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or 
! data to be ensured and, more generally, to use and operate it in the 
! same conditions as regards security. 
!
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL license and that you accept its terms.
!==============================================================================


!> @file
!! This file contains routines to read HDF5 halo files created with pfof and pfof_cone.

!> This module contains routines to read HDF5 halo files created with pfof and pfof_cone.
!>
!> Authors: F. Roy

Module modreadhalo

  Use modconstant

  Private
  Public :: h5read_halo_hfprop

Contains
  !> Reads a HDF5 halo hfprop file created with pfof
  Subroutine h5read_halo_hfprop(filename, common_metadata, parameter_pfof, info_ramses, info_cone, &
       position_halo, velocity_halo, rmax_halo, identity_halo, npart_halo)
    
    Use modhdf5
    Use modiocommons

    Implicit None
    
    ! input/output variables
    Character(len=400), intent(in) :: filename !< Name of the file
    Type(Type_common_metadata), intent(out) :: common_metadata
    Real(kind=8), dimension(:,:), allocatable, intent(out), optional :: position_halo, velocity_halo
    Real(kind=8), dimension(:), allocatable, intent(out), optional :: rmax_halo
    Integer(kind=PRI), dimension(:), allocatable, intent(out), optional :: identity_halo 
    Integer(kind=4), dimension(:), allocatable, intent(out), optional :: npart_halo
    Class(Type_parameter_pfof), intent(out), optional :: parameter_pfof
    Type(Type_info_ramses), intent(out), optional :: info_ramses
    Type(Type_info_cone_part), intent(out), optional :: info_cone
    
    ! local variables
    Integer(hid_t) :: file_id, data_id, meta_id
    Character(len=H5STRLEN) :: dname, groupname
    Integer(kind=4) :: nhalo

    
    ! Open the file
    Call hdf5_open_file(filename, file_id)

    ! open metadata group
    groupname='metadata'
    Call hdf5_open_group(file_id, groupname, meta_id)
    
    ! read metadata
    ! Read the number of halos
    dname = 'nhalo_file'
    Call hdf5_read_attr(meta_id, dname, nhalo)
#ifdef DEBUG
    Print *,'DEBUG: nhalo_file=',nhalo
#endif

    ! close group 'metadata'
    Call hdf5_close_group(meta_id)

    Call h5read_meta_common_metadata(file_id, common_metadata)

    If(present(info_ramses)) Then
       Call h5read_meta_info_ramses(file_id, info_ramses)
    End If

    If(present(info_cone)) Then
       Call h5read_meta_info_cone(file_id, info_cone)
    End If

    If(present(parameter_pfof)) Then
       Call h5read_meta_pfof_parameter(file_id, parameter_pfof)
    End If


    ! open group 'data'
    groupname='data'
    Call hdf5_open_group(file_id, groupname, data_id)

    ! read data
    If(present(position_halo)) Then
       If(.not.Allocated(position_halo)) Allocate(position_halo(3,nhalo))
       ! Read position of the halos
       dname = 'position_halo'
       Call hdf5_read_data(data_id, dname, 3, nhalo, position_halo)
    End If

    If(present(velocity_halo)) Then
       If(.not.Allocated(velocity_halo)) Allocate(velocity_halo(3,nhalo))
       ! Read velocity of the halos
       dname = 'velocity_halo'
       Call hdf5_read_data(data_id, dname, 3, nhalo, velocity_halo)
    End If

    If(present(identity_halo)) Then
       If(.not.Allocated(identity_halo)) Allocate(identity_halo(nhalo))
       ! Read identity of the halos
       dname = 'identity_halo'
       Call hdf5_read_data(data_id, dname, nhalo, identity_halo)
    End If
 
    If(present(npart_halo)) Then
       If(.not.Allocated(npart_halo)) Allocate(npart_halo(nhalo))
       ! Read npart in each halo
       dname = 'npart_halo'
       Call hdf5_read_data(data_id, dname, nhalo, npart_halo)
    End If

    If(present(rmax_halo)) Then
       If(.not.Allocated(rmax_halo)) Allocate(rmax_halo(nhalo))
       ! Read max. radius of the halos
       dname = 'rmax_halo'
       Call hdf5_read_data(data_id, dname, nhalo, rmax_halo)
    End If

    ! close group 'data'
    Call hdf5_close_group(data_id)

    ! Close file
    Call hdf5_close_file(file_id)

  End Subroutine h5read_halo_hfprop

End Module modreadhalo
