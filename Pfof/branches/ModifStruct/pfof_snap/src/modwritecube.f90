!==============================================================================
! Project: pFoF
! File: pfof_snap/src/modwritecube.f90
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
!! This file contains subroutines used to write HDF5 cubes.

!> This module contains subroutines used to write HDF5 cubes.
!> Authors: F. Roy
Module modwritecube

  Use mpi
  Use modconstant
  Use modvarcommons
  Use modhdf5
  Use modmpicom
  Use modvariables
  Use modiocommons
  Use modmpicommons

  Implicit None

  Procedure(), Pointer :: writecube

  Private
  Public :: selectwritecubetype, &
       writecube


Contains

  !=======================================================================
  Subroutine selectwritecubetype()

    Implicit none
    
    If(param%gatherwrite_factor==1) Then
       If(param%do_sort_cube) Then
          writecube => h5writesortedcube
       Else
          writecube => h5writecube
       End If
    Else
       If(param%do_sort_cube) Then
          writecube => mpih5writesortedcube
       Else
          writecube => mpih5writecube
       End If
    End If
    

  End Subroutine selectwritecubetype

  !=======================================================================
  !> This subroutine writes the position, velocity and id of each particle on the process in a hdf5 file.
  !! One file is written per MPI process
  Subroutine h5writecube()

    Implicit none

    Character(len=400) :: filecube
    Character(len=5)  :: pid_char

    Character(len=H5STRLEN) :: dsetname                           ! Dataset name
    Character(len=H5STRLEN) :: aname                              ! Attribute name
    Character(len=H5STRLEN) :: adata
    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: codename

    Integer(hid_t) :: file_id                               ! File identifier
    Integer(hid_t) :: gr_id                                 ! Group identifier

    Integer(kind=4) :: mpierr
    Real(kind=4), dimension(6) :: boundaries
    Integer(kind=4), dimension(procNB) :: nparttab
    Integer(kind=8) :: npart8
    Logical(kind=4) :: islast
    
#ifdef DEBUG
    Print *,"Enter h5writecube on process ",procID
#endif    
    
    Write(pid_char(1:5),'(I5.5)') procID
    filecube = 'pfof_cube_snap_part_data_'//trim(param%simulation_name)//'_'//pid_char//'.h5'
    
    Call Mpi_Gather(local_npart, 1, Mpi_Integer, nparttab, 1, Mpi_Integer, 0, &
         info_proc%global_comm%name, mpierr)

    ! create the hdf5 file
    Call hdf5_create_file(filecube, file_id)

    codename='pfof_snap'
    npart8 = 2**(inforamses%levelmin)
    npart8 = npart8**3
    Call h5write_meta_common(file_id, codename, npart8)
    Call h5write_meta_pfof_parameter(file_id, param)
    islast = .false.
    Call h5write_meta_info_ramses(file_id, inforamses, islast)


    ! open the root group
    groupname = 'metadata'
    Call hdf5_open_group(file_id,groupname, gr_id)

    ! Write type as attribute    
    aname = 'file_type'
    adata = 'cube'
    Call hdf5_write_attr(gr_id, aname, adata)
        
    ! Write the number of particles in this file as an integer(kind=8) dataset
    dsetname = 'npart_file'
    npart8 = int(local_npart, kind=8)
    Call hdf5_write_data(gr_id, dsetname, npart8)

    ! Write the process ID as an attribute
    aname = 'procID'
    Call hdf5_write_attr(gr_id, aname, procID)

    aname = 'nfile'
    Call hdf5_write_attr(gr_id, aname, procNB)
    
    ! Write the boundaries of the cube as an attribute
    boundaries(1) = xmin
    boundaries(2) = xmax
    boundaries(3) = ymin
    boundaries(4) = ymax
    boundaries(5) = zmin
    boundaries(6) = zmax
    aname = 'boundaries'

    Call hdf5_write_attr(gr_id, aname, 6, boundaries)

    Call hdf5_close_group(gr_id)

    groupname = 'data'
    Call hdf5_create_group(file_id,groupname,gr_id)

    ! Write the position of the particles
    dsetname='position_part'
    Call hdf5_write_data(gr_id, dsetname, 3, local_npart, position)

    ! Write the velocity of the particles
    dsetname='velocity_part'
    Call hdf5_write_data(gr_id, dsetname, 3, local_npart, velocity) 

    ! If we use potential, write potential of the particles
    If(param%do_read_gravitational_field) Then
       dsetname = 'gravitational_field_part'
       Call hdf5_write_data(gr_id, dsetname, 3, local_npart, field)
    End If   

    ! If we use potential, write potential of the particles
    If(param%do_read_potential) Then
       dsetname = 'potential_part'
       Call hdf5_write_data(gr_id, dsetname, local_npart, potential)
    End If

    ! Write the ID of the particles
    dsetname='identity_part'
    Call hdf5_write_data(gr_id, dsetname, local_npart, pfof_id) 

    Call hdf5_close_group(gr_id)

    Call hdf5_close_file(file_id)

#ifdef DEBUG
    Print *,'Exit h5writecube on process ',procID
#endif    

  End Subroutine h5writecube


  !=======================================================================
  !> This subroutine writes the position, the velocity and the id of each particle 
  !! on the process in a hdf5 file. 
  !! The particles are gathered 
  Subroutine mpih5writecube()

    Implicit none

    Character(len=400) :: filecube
    Character(len=5)  :: pid_char

    Character(len=H5STRLEN) :: dsetname                           ! Dataset name
    Character(len=H5STRLEN) :: aname                              ! Attribute name
    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: adata
    Character(len=H5STRLEN) :: codename

    Integer(hid_t) :: file_id                               ! File identifier
    Integer(hid_t) :: gr_id                                 ! Group identifier

    Integer(kind=4) :: procperfile
    Integer(kind=4) :: npart
    Integer(kind=4) :: nfile
    Integer(kind=4), dimension(:), allocatable :: partnb_tab

    Real(kind=4), dimension(6) :: boundaries
    Real(kind=4), dimension(:,:), allocatable :: bound_tab
    Integer(kind=8) :: npart8
    Integer(kind=4) :: mpierr
    Logical(kind=4) :: islast
    
#ifdef DEBUG
    Print *,"Enter mpih5writecube on process ",procID
#endif    

    ! number of processes writing in the same file
    procperfile = param%gatherwrite_factor**3
        
    Allocate(partnb_tab(procperfile))
    Allocate(bound_tab(6,procperfile))

    Write(pid_char(1:5),'(I5.5)') info_proc%write_comm%color
    filecube = 'pfof_cube_snap_part_data_'//trim(param%simulation_name)//'_'//pid_char//'.h5'

    ! create the hdf5 file
    Call hdf5_create_mpi_file(filecube, info_proc%write_comm%name, file_id)

    codename='pfof_snap'
    npart8 = 2**(inforamses%levelmin)
    npart8 = npart8**3
    Call h5write_meta_common(file_id, codename, npart8)
    Call h5write_meta_pfof_parameter(file_id, param)
    islast = .false.
    Call h5write_meta_info_ramses(file_id, inforamses, islast)

    ! open the root group
    groupname = 'metadata'
    Call hdf5_open_group(file_id, groupname, gr_id)

    ! Offset and global dimensions have to be computed
    Call Mpi_Allgather(local_npart, 1, Mpi_Integer, partnb_tab, 1, Mpi_Integer, info_proc%write_comm%name, mpierr)

    aname='npart_cube_array'
    Call hdf5_write_attr(gr_id, aname, procperfile, partnb_tab)

    ! Write type as attribute    
    aname = 'file_type'
    adata = 'mpicube'
    Call hdf5_write_attr(gr_id, aname, adata)

    aname = 'gatherwrite_factor'
    Call hdf5_write_attr(gr_id, aname, param%gatherwrite_factor)

    npart = sum(partnb_tab)
    npart8 = int(npart, kind=8)
    dsetname = 'npart_file'
    Call hdf5_write_data(gr_id, dsetname, npart8)

    aname = 'nfile'
    nfile = procNB / procperfile
    Call hdf5_write_attr(gr_id, aname, nfile)

    aname = 'procID'
    Call hdf5_write_attr(gr_id, aname, info_proc%write_comm%color) 

    ! Write the boundaries of the cube as an attribute
    boundaries(1) = xmin
    boundaries(2) = xmax
    boundaries(3) = ymin
    boundaries(4) = ymax
    boundaries(5) = zmin
    boundaries(6) = zmax

    Call Mpi_Allgather(boundaries,6,Mpi_Real, bound_tab, 6, Mpi_Real, info_proc%write_comm%name, mpierr)

    aname = 'boundaries_array'
    Call hdf5_write_attr(gr_id, aname, 6, procperfile, bound_tab)

    boundaries(1) = minval(bound_tab(1,:))
    boundaries(2) = maxval(bound_tab(2,:))
    boundaries(3) = minval(bound_tab(3,:))
    boundaries(4) = maxval(bound_tab(4,:))
    boundaries(5) = minval(bound_tab(5,:))
    boundaries(6) = maxval(bound_tab(6,:))
    aname = 'boundaries'
    Call hdf5_write_attr(gr_id, aname, 6, boundaries)

    Call hdf5_close_group(gr_id)

    groupname = 'data'
    Call hdf5_create_group(file_id, groupname, gr_id)

    dsetname = 'position_part'
    Call hdf5_write_mpi_data(gr_id, dsetname, 3, local_npart, position, info_proc%write_comm%name)

    dsetname = 'velocity_part'
    Call hdf5_write_mpi_data(gr_id, dsetname, 3, local_npart, velocity, info_proc%write_comm%name)

    dsetname = 'identity_part'
    Call hdf5_write_mpi_data(gr_id, dsetname, local_npart, pfof_id, info_proc%write_comm%name)

    If(param%do_read_potential) Then
       dsetname = 'potential_part'
       Call hdf5_write_mpi_data(gr_id, dsetname, local_npart, potential, info_proc%write_comm%name)
    End If

    If(param%do_read_gravitational_field) Then
       dsetname = 'gravitational_field_part'
       Call hdf5_write_mpi_data(gr_id, dsetname, 3, local_npart, field, info_proc%write_comm%name)
    End If

    ! Close the root group.
    Call hdf5_close_group(gr_id)

    Deallocate(partnb_tab, bound_tab)

    ! Close h5 file
    Call hdf5_close_mpi_file(file_id)

#ifdef DEBUG
    Print *,"Exit mpih5writecube on process ",procID
#endif    

  End Subroutine mpih5writecube


  !=======================================================================

  Subroutine h5writesortedcube()

    Use modsort

    Implicit none

    Character(len=400) :: filecube
    Character(len=5)  :: pid_char
    Character(len=8) :: charic

    Character(len=H5STRLEN) :: dsetname                           ! Dataset name
    Character(len=H5STRLEN) :: aname                              ! Attribute name
    Character(len=H5STRLEN) :: adata
    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: codename

    Integer(hid_t) :: file_id                               ! File identifier
    Integer(hid_t) :: gr_id
    Integer(hid_t) :: gr_data_id

    Real(kind=4), dimension(6) :: boundaries

    Integer(kind=4) :: ic, ix, iy, iz, ip, deb, fin, nc, ncdim, deltam1
    Integer(kind=4), dimension(:), allocatable :: npartcube
    Integer(kind=8), dimension(:), allocatable :: ictable
    Logical(kind=4) :: fileexist, fileopened

    Integer(kind=8) :: npart8
    Logical(kind=4) :: islast

#ifdef DEBUG
    Print *,"Enter h5writesortedcube on process ",procID
#endif    

    ! each FoF cube is divided into nc=512=8x8x8 groups, with ncdim=8
    ncdim = 8
    nc = ncdim*ncdim*ncdim
    Allocate(ictable(local_npart))
    Allocate(npartcube(nc))
    npartcube = 0
    ! Ramses coarse grid is composed of nres^3 cells 
    ! on a process: nres^3 / procNB cells 
    ! => there is (nres/(ncdim*dims(1))^3 coarse cells in each group
    ! the size of a group is nres/(ncdim*dims(1))*(1/nres) where 1/nres is the size of 1 Ramses coarse cell
    ! => delta = 1/(ncdim*dims(1)) => 1/delta = ncdim*dims(1)
    deltam1 = ncdim*info_proc%global_comm%dims(1)

    ! We compute the "group" index of each particle and the number of particle in each group
    Do ip=1, local_npart
       ix = int((position(1,ip) - xmin)*deltam1 + 1)
       iy = int((position(2,ip) - ymin)*deltam1 + 1)
       iz = int((position(3,ip) - zmin)*deltam1 + 1)
       ! rounding issue
       If(ix>ncdim) ix=ncdim
       If(iy>ncdim) iy=ncdim
       If(iz>ncdim) iz=ncdim
       ic = ix + (iy-1)*ncdim + (iz-1)*ncdim*ncdim
       ictable(ip) = ic
       npartcube(ic) = npartcube(ic) + 1
    End Do

    ! We sort the particles along their group id
    If(param%do_read_potential .and. param%do_read_gravitational_field) Then
       Call heapsort(local_npart,ictable,position,velocity,field,potential,pfof_id)
    Else If (param%do_read_potential .and. .not. param%do_read_gravitational_field) Then
       Call heapsort(local_npart,ictable,position,velocity,potential,pfof_id)
    Else If (param%do_read_gravitational_field .and. .not. param%do_read_potential) Then
       Call heapsort(local_npart,ictable,position,velocity,field,pfof_id)
    Else
       Call heapsort(local_npart,ictable,position,velocity,pfof_id)
    End If

    ! We open a sortedcube file and write the groups into it
    Write(pid_char(1:5),'(I5.5)') procID
    filecube = 'pfof_cube_snap_part_data_'//trim(param%simulation_name)//'_'//pid_char//'.h5'

    Inquire(File=filecube,exist=fileexist, opened=fileopened)
    If(fileopened) Then
       Print *,'Error on process ',procID,' : File ',filecube,' already opened'
    End If

    ! create the hdf5 file
    Call hdf5_create_file(filecube, file_id)

    codename='pfof_snap'
    npart8 = 2**(inforamses%levelmin)
    npart8 = npart8**3
    Call h5write_meta_common(file_id, codename, npart8)
    Call h5write_meta_pfof_parameter(file_id, param)
    islast = .false.
    Call h5write_meta_info_ramses(file_id, inforamses, islast)

    ! open the root group
    groupname = 'metadata'
    Call hdf5_open_group(file_id, groupname, gr_id)

    ! Write type as attribute    
    aname = 'file_type'
    adata = 'sortedcube'
    Call hdf5_write_attr(gr_id, aname, adata)

    ! Write the number of particles as an attribute
    dsetname = 'npart_file'
    npart8 = int(local_npart, kind=8)
    Call hdf5_write_data(gr_id, dsetname, npart8)

    ! Write the process ID as an attribute
    aname = 'procID'
    Call hdf5_write_attr(gr_id, aname, procID)

    aname = 'nfile'
    Call hdf5_write_attr(gr_id, aname, procNB)
    
    ! Write the boundaries of the cube as an attribute
    boundaries(1) = xmin
    boundaries(2) = xmax
    boundaries(3) = ymin
    boundaries(4) = ymax
    boundaries(5) = zmin
    boundaries(6) = zmax
    aname = 'boundaries'
    Call hdf5_write_attr(gr_id, aname, 6, boundaries)

    aname='ngroup'
    Call hdf5_write_attr(gr_id, aname, nc)

    aname='1/groupsize'
    Call hdf5_write_attr(gr_id, aname, deltam1)

    dsetname = 'npart_grp_array'
    Call hdf5_write_data(gr_id, dsetname, nc, npartcube)

    Call hdf5_close_group(gr_id)

    groupname = 'data'
    Call hdf5_create_group(file_id, groupname, gr_data_id)

    deb = 1
    ! For each non empty group we create an HDF5 group and write dataset into it
    Do ic = 1, nc
       If(npartcube(ic) /= 0) Then
          Write(charic(1:8),'(I8.8)') ic
          groupname = 'group'//charic
          
          Call hdf5_create_group(gr_data_id,groupname,gr_id)
       
          fin = deb + npartcube(ic) - 1

          ! Write the position of the particles
          dsetname='position_part'
          Call hdf5_write_data(gr_id, dsetname, 3, npartcube(ic), position(:,deb:fin))
          
          ! Write the velocity of the particles
          dsetname='velocity_part'
          Call hdf5_write_data(gr_id, dsetname, 3, npartcube(ic), velocity(:,deb:fin)) 
          
          ! Write the ID of the particles
          dsetname='identity_part'
          Call hdf5_write_data(gr_id, dsetname, npartcube(ic), pfof_id(deb:fin)) 
          
          ! Write potential if it is used
          If(param%do_read_potential) Then
             dsetname = 'potential_part'
             Call hdf5_write_data(gr_id, dsetname, npartcube(ic), potential(deb:fin))
          End If

          ! Write force if it is used
          If(param%do_read_gravitational_field) Then
             dsetname = 'gravitational_field_part'
             Call hdf5_write_data(gr_id, dsetname, 3, npartcube(ic), field(:,deb:fin))
          End If
          
          Call hdf5_close_group(gr_id)
          deb = fin + 1

       End If
    End Do
    
    Call hdf5_close_group(gr_data_id)

    Call hdf5_close_file(file_id)

#ifdef DEBUG
    Print *,'Exit h5writesortedcube on process ',procID
#endif    

  End Subroutine h5writesortedcube


  !=======================================================================
  !> This subroutine writes the position, the velocity and the id of each particle 
  !! on the process in a hdf5 file. 
  !! The particles are gathered 
  Subroutine mpih5writesortedcube()

    Use modsort

    Implicit none

    Character(len=400) :: filecube
    Character(len=5)  :: pid_char
    Character(len=8) :: charic

    Character(len=H5STRLEN) :: dsetname                           ! Dataset name
    Character(len=H5STRLEN) :: aname                              ! Attribute name
    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: adata
    Character(len=H5STRLEN) :: codename

    Integer(hid_t) :: file_id                               ! File identifier
    Integer(hid_t) :: gr_id
    Integer(hid_t) :: gr_data_id

    Integer(kind=4) :: procperfile
    Integer(kind=4) :: npart
    Integer(kind=4) :: nfile
    Integer(kind=4) :: ic, ix, iy, iz, ip, deb, fin, nc, ncdim, deltam1, fic
    Integer(kind=4), dimension(:), allocatable :: npartcube
    Integer(kind=8), dimension(:), allocatable :: ictable

    Integer(kind=4), dimension(:), allocatable :: partnb_tab

    Real(kind=4), dimension(6) :: boundaries
    Real(kind=4), dimension(:,:), allocatable :: bound_tab

    Logical :: empty

    Integer(kind=4) :: mpierr
    Integer(kind=8) :: npart8
    Logical(kind=4) :: islast

#ifdef DEBUG
    Print *,"Enter mpih5writesortedcube on process ",procID
#endif    


    ! each FoF cube is divided into nc=512=8x8x8 groups, with ncdim=8
    ncdim = 8
    nc = ncdim*ncdim*ncdim
    Allocate(ictable(local_npart))
    Allocate(npartcube(nc))
    npartcube = 0
    ! Ramses coarse grid is composed of nres^3 cells 
    ! on a process: nres^3 / procNB cells 
    ! => there is (nres/(ncdim*dims(1))^3 coarse cells in each group
    ! the size of a group is nres/(ncdim*dims(1))*(1/nres) where 1/nres is the size of 1 Ramses coarse cell
    ! => delta = 1/(ncdim*dims(1)) => 1/delta = ncdim*dims(1)
    deltam1 = ncdim*info_proc%global_comm%dims(1)

    ! We compute the "group" index of each particle and the number of particle in each group
    Do ip=1, local_npart
       ix = int((position(1,ip) - xmin)*deltam1 + 1)
       iy = int((position(2,ip) - ymin)*deltam1 + 1)
       iz = int((position(3,ip) - zmin)*deltam1 + 1)
       ! rounding issue
       If(ix>ncdim) ix=ncdim
       If(iy>ncdim) iy=ncdim
       If(iz>ncdim) iz=ncdim
       ic = ix + (iy-1)*ncdim + (iz-1)*ncdim*ncdim
       ictable(ip) = ic
       npartcube(ic) = npartcube(ic) + 1
    End Do

    ! We sort the particles along their group id
    If(param%do_read_potential .and. param%do_read_gravitational_field) Then
       Call heapsort(local_npart,ictable,position,velocity,field,potential,pfof_id)
    Else If (param%do_read_potential .and. .not. param%do_read_gravitational_field) Then
       Call heapsort(local_npart,ictable,position,velocity,potential,pfof_id)
    Else If (param%do_read_gravitational_field .and. .not. param%do_read_potential) Then
       Call heapsort(local_npart,ictable,position,velocity,field,pfof_id)
    Else
       Call heapsort(local_npart,ictable,position,velocity,pfof_id)
    End If

#ifdef DEBUG
    Print *,'Particles sorted'
#endif

    ! number of processes writing in the same file
    procperfile = param%gatherwrite_factor**3
        
    Allocate(partnb_tab(procperfile))
    Allocate(bound_tab(6,procperfile))

    Write(pid_char(1:5),'(I5.5)') info_proc%write_comm%color
    filecube = 'pfof_cube_snap_part_data_'//trim(param%simulation_name)//'_'//pid_char//'.h5'

    ! create the hdf5 file
    Call hdf5_create_mpi_file(filecube, info_proc%write_comm%name, file_id)

    codename = 'pfof_snap'
    npart8 = 2**(inforamses%levelmin)
    npart8 = npart8**3
    Call h5write_meta_common(file_id, codename, npart8)
    Call h5write_meta_pfof_parameter(file_id, param)
    islast = .false.
    Call h5write_meta_info_ramses(file_id, inforamses, islast)

    ! open the root group
    groupname = 'metadata'
    Call hdf5_open_group(file_id, groupname, gr_id)

    ! Offset and global dimensions have to be computed
    Call Mpi_Allgather(local_npart, 1, Mpi_Integer, partnb_tab, 1, Mpi_Integer, info_proc%write_comm%name, mpierr)

    aname='npart_cube_array'
    Call hdf5_write_attr(gr_id, aname, procperfile, partnb_tab)

    aname='ngroup'
    Call hdf5_write_attr(gr_id, aname, nc*procperfile)

    aname='1/groupsize'
    Call hdf5_write_attr(gr_id, aname, deltam1)

    dsetname = 'npart_grp_array'
    Call hdf5_write_mpi_data(gr_id, dsetname, nc, npartcube, info_proc%write_comm%name)

    aname = 'gatherwrite_factor'
    Call hdf5_write_attr(gr_id, aname, param%gatherwrite_factor)

    ! Write type as attribute    
    aname = 'file_type'
    adata = 'mpisortedcube'
    Call hdf5_write_attr(gr_id, aname, adata)

    npart = sum(partnb_tab)
    dsetname = 'npart_file'
    npart8 = int(npart, kind=8)
    Call hdf5_write_data(gr_id, dsetname, npart8)
    
    aname = 'nfile'
    nfile = procNB / procperfile
    Call hdf5_write_attr(gr_id, aname, nfile)

    aname = 'procID'
    Call hdf5_write_attr(gr_id, aname, info_proc%write_comm%color) 

    ! Write the boundaries of the cube as an attribute
    boundaries(1) = xmin
    boundaries(2) = xmax
    boundaries(3) = ymin
    boundaries(4) = ymax
    boundaries(5) = zmin
    boundaries(6) = zmax

    Call Mpi_Allgather(boundaries,6,Mpi_Real, bound_tab, 6, Mpi_Real, info_proc%write_comm%name, mpierr)

    aname = 'boundaries_array'
    Call hdf5_write_attr(gr_id, aname, 6, procperfile, bound_tab)

    boundaries(1) = minval(bound_tab(1,:))
    boundaries(2) = maxval(bound_tab(2,:))
    boundaries(3) = minval(bound_tab(3,:))
    boundaries(4) = maxval(bound_tab(4,:))
    boundaries(5) = minval(bound_tab(5,:))
    boundaries(6) = maxval(bound_tab(6,:))
    aname = 'boundaries'
    Call hdf5_write_attr(gr_id, aname, 6, boundaries)

    Call hdf5_close_group(gr_id)

    deb = 1
    fin = 1
    fic = nc*info_proc%write_comm%pid

    groupname = 'data'
    Call hdf5_create_group(file_id, groupname, gr_data_id)

    ! For each non empty group we create an HDF5 group and write dataset into it
    Do ic = 1, nc*procperfile

       Write(charic(1:8),'(I8.8)') ic
       groupname = 'group'//charic
       Call hdf5_create_group(gr_data_id, groupname, gr_id)

       fic = ic - nc*info_proc%write_comm%pid
       If(fic>=1 .and. fic<=nc) Then
          If(npartcube(fic) == 0) Then
             empty = .true.
             If(deb>local_npart) Then
                deb=1
                fin=1
             End If
          Else
             empty = .false.
             fin = deb + npartcube(fic) - 1
#ifdef DEBUG
             Print *,'group ',groupname,' opened:', procID, fic, deb, fin
#endif
          End If
       Else
          empty = .true.
          If(deb>local_npart .or. fin > local_npart) Then
             deb=1
             fin=1
          End If
       End If

       ! Write the position of the particles
       dsetname='position_part'
       Call hdf5_write_mpi_data(gr_id, dsetname, 3, fin-deb+1, position(:,deb:fin), &
            info_proc%write_comm%name, empty)
       
       ! Write the velocity of the particles
       dsetname='velocity_part'
       Call hdf5_write_mpi_data(gr_id, dsetname, 3, fin-deb+1, velocity(:,deb:fin),&
            info_proc%write_comm%name, empty) 
       
       ! Write the ID of the particles
       dsetname='identity_part'
       Call hdf5_write_mpi_data(gr_id, dsetname, fin-deb+1, pfof_id(deb:fin),&
            info_proc%write_comm%name,empty) 

       ! Write potential if it is used
       If(param%do_read_potential) Then
          dsetname = 'potential_part'
          Call hdf5_write_mpi_data(gr_id, dsetname, fin-deb+1, potential(deb:fin),&
               info_proc%write_comm%name,empty)
       End If

       ! Write force if it is used
       If(param%do_read_gravitational_field) Then
          dsetname = 'gravitational_field_part'
          Call hdf5_write_mpi_data(gr_id, dsetname, 3, fin-deb+1, field(:,deb:fin),&
               info_proc%write_comm%name,empty)
       End If
       
       Call hdf5_close_group(gr_id)

       If(fic>=1 .and. fic<=nc) Then
          If(npartcube(fic) /= 0) Then
             deb = fin + 1
#ifdef DEBUG
             Print *,'group ',groupname,' closed'
#endif
          End If
       End If
          
    End Do

    Call hdf5_close_group(gr_data_id)
    Deallocate(partnb_tab, bound_tab)

    ! Close h5 file
    Call hdf5_close_mpi_file(file_id)

#ifdef DEBUG
    Print *,"Exit mpih5writesortedcube on process ",procID
#endif    

  End Subroutine mpih5writesortedcube


End Module modwritecube
