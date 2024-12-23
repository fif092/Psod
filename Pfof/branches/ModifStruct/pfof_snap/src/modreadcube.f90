!==============================================================================
! Project: pFoF
! File: pfof_snap/src/modreadcube.f90
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
!! This file contains functions to read HDF5 cube files.

!> This module contains functions to read HDF5 cube files.
!> Authors: F. Roy

Module modreadcube

  Use modconstant
  Use modhdf5
  Use modmpicom
  Use modmpicommons
  Use modvariables
  Use modiocommons
  Use modvarcommons
  Use mpi
  

  Implicit none

  Procedure(), Pointer :: readcube

  Private
  Public :: selectreadcubetype, &
       readcube

Contains

  !=======================================================================
  Subroutine selectreadcubetype()

    Implicit none

    Integer(kind=hid_t) :: file_id
    Integer(kind=hid_t) :: gr_id
    Character(len=400) :: filename
    Character(len=H5STRLEN) :: h5name
    Character(len=H5STRLEN) :: filetype
    Integer(kind=4) :: mpierr
    
    If(procID==0) Then
       filename = trim(param%input_path)//trim(param%part_input_file)
       Call hdf5_open_file(filename, file_id)
       h5name = 'metadata'
       Call hdf5_open_group(file_id, h5name, gr_id)
       h5name = 'file_type'
       Call hdf5_read_attr(gr_id, h5name, len(filetype), filetype)
       Call hdf5_close_group(gr_id)
       Call hdf5_close_file(file_id)
    End If

    Call Mpi_Bcast(filetype, len(filetype), Mpi_Char, 0, Mpi_Comm_World, mpierr)

    Select case (filetype)
    Case ("cube")
       readcube => h5readcube
    Case ("mpicube")
       readcube => mpih5readcube
    Case ("sortedcube")
       readcube => h5readsortedcube
    Case ("mpisortedcube")
       readcube => mpih5readsortedcube
    Case default 
       Call EmergencyStop("Error while reading cube type",31)
    End Select

  End Subroutine selectreadcubetype

  !=======================================================================
  !> This subroutine reads hdf5 cube files created by pFOF.
  Subroutine h5readcube()

    Implicit none

    Character(len=400) :: filename                          ! File name
    Integer(hid_t) :: file_id                               ! File identifier
    Integer(hid_t) :: gr_id                                 ! Group identifier
    Character(len=H5STRLEN) :: aname                              ! Attribute name
    Character(len=H5STRLEN) :: dsetname                           ! Dataset name
    Character(len=H5STRLEN) :: groupname
    Real(kind=4), dimension(6) :: boundaries
    Integer(kind=8) :: npart8
    Integer(kind=4) :: length
    Logical(kind=4) :: islast

#ifdef DEBUG
    Print *,"Enter h5readcube on process ",procID
#endif

    filename = trim(param%input_path)//trim(param%part_input_file)
    length=len_trim(filename)
    Write(filename(length-7:length-3),'(I5.5)') procID

    ! open the file
    Call hdf5_open_file(filename, file_id)

    Call h5read_meta_common(file_id, common_meta)
    islast = .false.
    Call h5read_meta_info_ramses(file_id, inforamses, islast)

    nres = 2**inforamses%levelmin
    ngrid = int(nres,kind=8)**3
    global_npart = ngrid

#ifdef DEBUG
    Print *,'levelmin=',inforamses%levelmin, ' on process ',procID
    Print *,'nres=',nres, ' on process ',procID
    Print *,'ngrid=',ngrid, ' on process ',procID
#endif

    ! open the root group
    groupname = 'metadata'
    Call hdf5_open_group(file_id,groupname, gr_id)

    ! read attribute partNB
    dsetname = 'npart_file'
    Call hdf5_read_data(gr_id, dsetname, npart8)
    local_npart = int(npart8, kind=4)

    ! read attribute boundaries
    aname = 'boundaries'
    Call hdf5_read_attr(gr_id, aname, 6, boundaries)

    xmin=boundaries(1)
    xmax=boundaries(2)
    ymin=boundaries(3)
    ymax=boundaries(4)
    zmin=boundaries(5)
    zmax=boundaries(6)

    Call hdf5_close_group(gr_id)

    If(.not.Allocated(pfof_id)) Allocate(pfof_id(local_npart))
    If(.not.Allocated(position))  Allocate(position(3,local_npart))
    If(.not.Allocated(velocity))  Allocate(velocity(3,local_npart))

    groupname = 'data'
    Call hdf5_open_group(file_id, groupname, gr_id)

    ! read position of the particles
    dsetname = 'position_part'
    Call hdf5_read_data(gr_id, dsetname, 3, local_npart, position)

    ! read velocity of the particles
    dsetname = 'velocity_part'
    Call hdf5_read_data(gr_id, dsetname, 3, local_npart, velocity)

    ! read id of the particles
    dsetname = 'identity_part'
    Call hdf5_read_data(gr_id, dsetname, local_npart, pfof_id)
    
    If(param%do_read_potential) Then
       If(.not.Allocated(potential)) Allocate(potential(local_npart))
       dsetname = 'potential_part'
       Call hdf5_read_data(gr_id, dsetname, local_npart, potential)
    End If

    If(param%do_read_gravitational_field) Then
       If(.not.Allocated(field)) Allocate(field(3,local_npart))
       dsetname = 'gravitational_field_part'
       Call hdf5_read_data(gr_id, dsetname, 3, local_npart, field)
    End If

    ! Close the root group.
    Call hdf5_close_group(gr_id)

    Call hdf5_close_file(file_id)


#ifdef DEBUG
    Print *,"Exit h5readcube on process ",procID
#endif

  End Subroutine h5readcube



  !=======================================================================
  !> This subroutine reads hdf5 sorted cube files created by pFOF.
  Subroutine h5readsortedcube()

    Implicit none

    Character(len=400) :: filename                          ! File name
    Integer(hid_t) :: file_id                               ! File identifier
    Integer(hid_t) :: gr_data_id                            ! Group identifier
    Integer(hid_t) :: gr_id                                 ! Group identifier
    Character(len=8)  :: gid_char
    Character(len=H5STRLEN) :: aname                              ! Attribute name
    Character(len=H5STRLEN) :: dsetname                           ! Dataset name
    Character(len=H5STRLEN) :: groupname
    Real(kind=4), dimension(6) :: boundaries
    Integer(kind=4) :: ngroup
    Integer(kind=4), dimension(:), allocatable :: npartpergroup
    Integer(kind=4) :: igroup
    Integer(kind=4) :: indbeg, indend
    Integer(kind=4) :: length
    Integer(kind=8) :: npart8
    Logical(kind=4) :: islast
    
#ifdef DEBUG
    Print *,"Enter h5readsortedcube on process ",procID
#endif

    filename = trim(param%input_path)//trim(param%part_input_file)
    length=len_trim(filename)
    Write(filename(length-7:length-3),'(I5.5)') procID

    ! open the file
    Call hdf5_open_file(filename, file_id)

    Call h5read_meta_common(file_id, common_meta)
    islast = .false.
    Call h5read_meta_info_ramses(file_id, inforamses, islast)

    nres = 2**inforamses%levelmin
    ngrid = int(nres,kind=8)**3
    global_npart = ngrid

    ! open the root group
    groupname = 'metadata'
    Call hdf5_open_group(file_id,groupname, gr_id)

    ! read attribute partNB
    dsetname = 'npart_file'
    Call hdf5_read_data(gr_id, dsetname, npart8)
    local_npart = int(npart8,kind=4)

    ! read attribute boundaries
    aname = 'boundaries'
    Call hdf5_read_attr(gr_id, aname, 6, boundaries)

    xmin=boundaries(1)
    xmax=boundaries(2)
    ymin=boundaries(3)
    ymax=boundaries(4)
    zmin=boundaries(5)
    zmax=boundaries(6)

    ! read attribute ngroup
    aname = 'ngroup'
    Call hdf5_read_attr(gr_id, aname, ngroup)

    Allocate(npartpergroup(ngroup))
    aname = 'npart_grp_array'
    Call hdf5_read_data(gr_id, aname, ngroup, npartpergroup)

    Call hdf5_close_group(gr_id)

    groupname='data'
    Call hdf5_open_group(file_id, groupname, gr_data_id)

    If(.not.Allocated(pfof_id)) Allocate(pfof_id(local_npart))
    If(.not.Allocated(position))  Allocate(position(3,local_npart))
    If(.not.Allocated(velocity))  Allocate(velocity(3,local_npart))

    indbeg = 1
    indend = 0
    ! loop over the groups
    Do igroup=1, ngroup
       If(npartpergroup(igroup) /= 0) Then ! there is at least one part. to read in this group
          Write(gid_char(1:8),'(I8.8)') igroup
          groupname='group'//gid_char
          Call hdf5_open_group(gr_data_id,groupname,gr_id)
          
          indend = indbeg + npartpergroup(igroup) - 1
          ! read position of the particles
          dsetname = 'position_part'
          Call hdf5_read_data(gr_id, dsetname, 3, local_npart, position(:,indbeg:indend))
          
          ! read velocity of the particles
          dsetname = 'velocity_part'
          Call hdf5_read_data(gr_id, dsetname, 3, local_npart, velocity(:,indbeg:indend))
          
          ! read id of the particles
          dsetname = 'identity_part'
          Call hdf5_read_data(gr_id, dsetname, local_npart, pfof_id(indbeg:indend))

          ! read potential if requested
          If(param%do_read_potential) Then
             If(.not.Allocated(potential)) Allocate(potential(local_npart))
             dsetname = 'potential_part'
             Call hdf5_read_data(gr_id, dsetname, local_npart, potential(indbeg:indend))
          End If

          ! read force if requested
          If(param%do_read_gravitational_field) Then
             If(.not.Allocated(field)) Allocate(field(3,local_npart))
             dsetname = 'gravitational_field_part'
             Call hdf5_read_data(gr_id, dsetname, 3, local_npart, field(:,indbeg:indend))
          End If
          
          indbeg = indend + 1

          Call hdf5_close_group(gr_id)
       End If
    End Do

    If(indend /= local_npart) Then
       Print *,'Error while reading particles from file ',filename
       Call EmergencyStop('Error in h5readsortedcube',32)
    End If

    ! Close the root group.
    Call hdf5_close_group(gr_data_id)

    Call hdf5_close_file(file_id)

    Deallocate(npartpergroup)

#ifdef DEBUG
    Print *,"Exit h5readsortedcube on process ",procID
#endif

  End Subroutine h5readsortedcube


  !=======================================================================
 
  Subroutine mpih5readcube()
    
#ifdef DEBUG
    Use modmpicommons, only : procID
#endif

    Implicit none

    Character(len=400) :: filename

    Character(len=H5STRLEN) :: dsetname                           ! Dataset name
    Character(len=H5STRLEN) :: aname                              ! Attribute name
    Character(len=H5STRLEN) :: groupname

    Integer(hid_t) :: file_id                               ! File identifier
    Integer(hid_t) :: gr_id                                 ! Group identifier
    
    Integer(kind=4) :: procperfile
    Integer(kind=4), dimension(:), allocatable :: partnb_tab

    Real(kind=4), dimension(:,:), allocatable :: bound_tab
    Integer(kind=8) :: npart8
    Integer(kind=4) :: length
    Logical(kind=4) :: islast
    
#ifdef DEBUG
    Integer(kind=4) :: mpierr

    Print *,"Enter mpih5readcube on process ",procID
#endif    
    
    ! number of processes writing in the same file
    procperfile = info_proc%read_comm%size !param%gatherread_factor**3
    
    Allocate(partnb_tab(procperfile))
    Allocate(bound_tab(6,procperfile))
        
    filename = trim(param%input_path)//trim(param%part_input_file)
    length=len_trim(filename)
    Write(filename(length-7:length-3),'(I5.5)') info_proc%read_comm%color

    Call hdf5_open_mpi_file(filename, info_proc%read_comm%name, file_id)

    Call h5read_meta_common(file_id, common_meta)
    islast = .false.
    Call h5read_meta_info_ramses(file_id, inforamses, islast)

    nres = 2**inforamses%levelmin
    ngrid = int(nres,kind=8)**3
    global_npart = ngrid

    ! open the root group
    groupname = 'metadata'
    Call hdf5_open_group(file_id,groupname, gr_id)

    ! read attribute partNB
    dsetname = 'npart_file'
    Call hdf5_read_data(gr_id, dsetname, npart8)
    local_npart = int(npart8,kind=4)

    ! read attribute partNB
    aname = 'npart_cube_array'
    Call hdf5_read_attr(gr_id, aname, 6, partnb_tab)

#ifdef DEBUG
    Print *,"partNB_tab read on process ",procID
    Call Mpi_Barrier(Mpi_Comm_World, mpierr)
#endif    

    local_npart = partnb_tab(info_proc%read_comm%pid+1)

    ! read attribute boundaries
    aname = 'boundaries_array'
    Call hdf5_read_attr(gr_id, aname, 6, procperfile, bound_tab)

#ifdef DEBUG
    Print *,"boundaries_tab read on process ",procID
    Call Mpi_Barrier(Mpi_Comm_World, mpierr)
#endif    

    xmin=bound_tab(1,info_proc%read_comm%pid+1)
    xmax=bound_tab(2,info_proc%read_comm%pid+1)
    ymin=bound_tab(3,info_proc%read_comm%pid+1)
    ymax=bound_tab(4,info_proc%read_comm%pid+1)
    zmin=bound_tab(5,info_proc%read_comm%pid+1)
    zmax=bound_tab(6,info_proc%read_comm%pid+1)

    Call hdf5_close_group(gr_id)

    If(.not.Allocated(pfof_id)) Allocate(pfof_id(local_npart))
    If(.not.Allocated(position))  Allocate(position(3,local_npart))
    If(.not.Allocated(velocity))  Allocate(velocity(3,local_npart))

    groupname = 'data'
    Call hdf5_open_group(file_id, groupname, gr_id)

    dsetname = 'position_part'
    Call hdf5_read_mpi_data(gr_id, dsetname, 3, local_npart, position, info_proc%read_comm%name)

#ifdef DEBUG
    Print *,"position read on process ",procID
    Call Mpi_Barrier(Mpi_Comm_World, mpierr)
#endif    

    dsetname = 'velocity_part'
    Call hdf5_read_mpi_data(gr_id, dsetname, 3, local_npart, velocity, info_proc%read_comm%name)

#ifdef DEBUG
    Print *,"velocity read on process ",procID
    Call Mpi_Barrier(Mpi_Comm_World, mpierr)
#endif    

    dsetname = 'identity_part'
    Call hdf5_read_mpi_data(gr_id, dsetname, local_npart, pfof_id, info_proc%read_comm%name)

#ifdef DEBUG
    Print *,"ID read on process ",procID
    Call Mpi_Barrier(Mpi_Comm_World, mpierr)
#endif    

    If(param%do_read_potential) Then
       If(.not.Allocated(potential)) Allocate(potential(local_npart))
       dsetname = 'potential_part'
       Call hdf5_read_mpi_data(gr_id, dsetname, local_npart, potential, info_proc%read_comm%name)
    End If

    ! if requested, try to read the force attribute
    If(param%do_read_gravitational_field) Then
       If(.not.Allocated(field)) Allocate(field(3,local_npart))
       dsetname = 'gravitational_field_part'
       Call hdf5_read_mpi_data(gr_id, dsetname, 3, local_npart, field, info_proc%read_comm%name)
    End If

    ! Close the root group.
    Call hdf5_close_group(gr_id)

    Call hdf5_close_file(file_id)
    
    Deallocate(partnb_tab)
    Deallocate(bound_tab)


#ifdef DEBUG
    Print *,"Exit mpih5readcube on process ",procID
#endif    
    
    
  End Subroutine mpih5readcube


  !=======================================================================
  Subroutine mpih5readsortedcube()
    
    Implicit none

    Character(len=400) :: filename
    Character(len=8)  :: gid_char

    Character(len=H5STRLEN) :: dsetname                           ! Dataset name
    Character(len=H5STRLEN) :: aname                              ! Attribute name
    Character(len=H5STRLEN) :: groupname

    Integer(hid_t) :: file_id                               ! File identifier
    Integer(hid_t) :: gr_data_id                            ! Group identifier
    Integer(hid_t) :: gr_id                                 ! Group identifier
    
    Integer(kind=4) :: procperfile
    Integer(kind=4), dimension(:), allocatable :: partnb_tab
    Integer(kind=4) :: igroup, firstgroup, lastgroup
    Integer(kind=4) :: indbeg, indend
    Integer(kind=4) :: nc, ngroup
    
    Integer(kind=4), dimension(:), allocatable :: npartpergroup

    Real(kind=4), dimension(:,:), allocatable :: bound_tab
    Integer(kind=8) :: npart8
    Integer(kind=4) :: length
    Logical(kind=4) :: islast
    
#ifdef DEBUG
    Integer(kind=4) :: mpierr

    Print *,"Enter mpih5readsortedcube on process ",procID
#endif    
    
    ! number of processes reading in the same file
    procperfile = param%gatherread_factor**3
    
    Allocate(partnb_tab(procperfile))
    Allocate(bound_tab(6,procperfile))
        
    filename = trim(param%input_path)//trim(param%part_input_file)
    length=len_trim(filename)
    Write(filename(length-7:length-3),'(I5.5)') info_proc%read_comm%color

    ! we open the file for a serial read
    Call hdf5_open_file(filename, file_id)

    Call h5read_meta_common(file_id, common_meta)
    islast = .false.
    Call h5read_meta_info_ramses(file_id, inforamses, islast)

    nres = 2**inforamses%levelmin
    ngrid = int(nres,kind=8)**3
    global_npart = ngrid

    ! open the root group
    groupname = 'metadata'
    Call hdf5_open_group(file_id,groupname, gr_id)

    ! read attribute partNB
    dsetname = 'npart_file'
    Call hdf5_read_data(gr_id, dsetname, npart8)
    local_npart = int(npart8,kind=4)

#ifdef DEBUG
    Print *,"Value of nres read = ",nres, " on process ",procID
    Call Mpi_Barrier(Mpi_Comm_World, mpierr)
#endif

    ! read attribute partNB
    aname = 'npart_cube_array'
    Call hdf5_read_attr(gr_id, aname, procperfile, partnb_tab)

#ifdef DEBUG
    Print *,"partNB_tab read on process ",procID
    Call Mpi_Barrier(Mpi_Comm_World, mpierr)
#endif    

    local_npart = partnb_tab(info_proc%read_comm%pid+1)

    ! read attribute boundaries
    aname = 'boundaries_array'
    Call hdf5_read_attr(gr_id, aname, 6, procperfile, bound_tab)

#ifdef DEBUG
    Print *,"boundaries_tab read on process ",procID
    Call Mpi_Barrier(Mpi_Comm_World, mpierr)
#endif    

    xmin=bound_tab(1,info_proc%read_comm%pid+1)
    xmax=bound_tab(2,info_proc%read_comm%pid+1)
    ymin=bound_tab(3,info_proc%read_comm%pid+1)
    ymax=bound_tab(4,info_proc%read_comm%pid+1)
    zmin=bound_tab(5,info_proc%read_comm%pid+1)
    zmax=bound_tab(6,info_proc%read_comm%pid+1)

    ! read attribute ngroup
    aname = 'ngroup'
    Call hdf5_read_attr(gr_id, aname, ngroup)

    Allocate(npartpergroup(ngroup))
    aname = 'npart_grp_array'
    Call hdf5_read_data(gr_id, aname, ngroup, npartpergroup)

    nc = ngroup / procperfile
    firstgroup = nc * info_proc%read_comm%pid + 1
    lastgroup  = nc *(info_proc%read_comm%pid+1)

    Call hdf5_close_group(gr_id)

    If(.not.Allocated(pfof_id)) Allocate(pfof_id(local_npart))
    If(.not.Allocated(position))  Allocate(position(3,local_npart))
    If(.not.Allocated(velocity))  Allocate(velocity(3,local_npart))

    groupname ='data'
    Call hdf5_open_group(file_id, groupname, gr_data_id)

    indbeg = 1
    indend = 0
    ! loop over the groups
    Do igroup=firstgroup, lastgroup
       If(npartpergroup(igroup) /= 0) Then ! there is at least one part. to read in this group
          Write(gid_char(1:8),'(I8.8)') igroup
          groupname='group'//gid_char
          Call hdf5_open_group(gr_data_id,groupname,gr_id)
          
          indend = indbeg + npartpergroup(igroup) - 1
          ! read position of the particles
          dsetname = 'position_part'
          Call hdf5_read_data(gr_id, dsetname, 3, local_npart, position(:,indbeg:indend))
          
          ! read velocity of the particles
          dsetname = 'velocity_part'
          Call hdf5_read_data(gr_id, dsetname, 3, local_npart, velocity(:,indbeg:indend))
          
          ! read id of the particles
          dsetname = 'identity_part'
          Call hdf5_read_data(gr_id, dsetname, local_npart, pfof_id(indbeg:indend))
          
          ! read potential if requested
          If(param%do_read_potential) Then
             If(.not.Allocated(potential)) Allocate(potential(local_npart))
             dsetname = 'potential_part'
             Call hdf5_read_data(gr_id, dsetname, local_npart, potential(indbeg:indend))
          End If

          ! read force if requested
          If(param%do_read_gravitational_field) Then
             If(.not.Allocated(field)) Allocate(field(3,local_npart))
             dsetname = 'gravitational_field_part'
             Call hdf5_read_data(gr_id, dsetname, 3, local_npart, field(:,indbeg:indend))
          End If
          
          indbeg = indend + 1

          Call hdf5_close_group(gr_id)
       End If
    End Do


    ! Close the root group.
    Call hdf5_close_group(gr_data_id)

    Call hdf5_close_file(file_id)
    
    Deallocate(partnb_tab)
    Deallocate(bound_tab)

#ifdef DEBUG
    Print *,"Exit mpih5readsortedcube on process ",procID
#endif    
        
  End Subroutine mpih5readsortedcube

End Module modreadcube
