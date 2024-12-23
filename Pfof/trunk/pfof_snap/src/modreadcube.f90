! Copyright (c) 2007-2016 CNRS, Fabrice Roy, Vincent Bouillot
! Author: Fabrice Roy (LUTH/CNRS/PSL), fabrice.roy@obspm.fr
! Vincent Bouillot (LUTH/CNRS/PSL)
! Distributed under the CeCILL 2.0 license. For full terms see the file LICENSE.

!> @file
!!
!! This file contains the subroutines used to read "cube" files produced by pfof_snap.
!!
!! There are 4 differents cube types, there is one reading subroutine per cube type.
!!
!! There is a subroutine used to select the right reading function according to the file_type
!! metadata read from the HDF5 cube file.
!!
!! @author Fabrice Roy


Module modreadcube

  Use modconstant
  Use modhdf5
  Use modmpicom
  Use modparameters
  Use modvariables
  Use modiocommons

  Implicit none

  Procedure(), Pointer :: readcube !< Procedure pointer pointing to the cube reading subroutine selected from the "file_type" metadata in the HDF5 cube file.

Contains

  !=======================================================================
  !> @function
  !! Opens a cube file, reads the file_type metadata, and selects the corresponding reading subroutine.
  !! Sets the readcube procedure pointer to the right subroutine.
  Subroutine selectcubetype()

    Implicit none

    Integer(kind=hid_t) :: file_id
    Integer(kind=hid_t) :: gr_id
    Character(len=400) :: filename
    Character(len=H5STRLEN) :: h5name
    Character(len=H5STRLEN) :: filetype
    
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
    Print *,'FILETYPE=',filetype

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

  End Subroutine selectcubetype

  !=======================================================================
  !> Reads a standard hdf5 cube file created by pfof_snap.
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

    islast = .false.
    Call h5readcommonmetadata(file_id, islast, inforamses)

    nres = 2**inforamses%levelmin
    ngrid = int(nres,kind=8)**3
    npart = ngrid

#ifdef DEBUG
    Print *,'NRES=',nres, ' on process ',procID
#endif
    
    ! open the root group
    groupname = 'metadata'
    Call hdf5_open_group(file_id,groupname, gr_id)

    ! read attribute partNB
    dsetname = 'npart_file'
    Call hdf5_read_data(gr_id, dsetname, npart8)
    mynpart = int(npart8, kind=4)

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

    If(.not.Allocated(id)) Allocate(id(mynpart))
    If(.not.Allocated(pos))  Allocate(pos(3,mynpart))
    If(.not.Allocated(vel))  Allocate(vel(3,mynpart))

    groupname = 'data'
    Call hdf5_open_group(file_id, groupname, gr_id)

    ! read position of the particles
    dsetname = 'position_part'
    Call hdf5_read_data(gr_id, dsetname, 3, mynpart, pos)

    ! read velocity of the particles
    dsetname = 'velocity_part'
    Call hdf5_read_data(gr_id, dsetname, 3, mynpart, vel)

    ! read id of the particles
    dsetname = 'identity_part'
    Call hdf5_read_data(gr_id, dsetname, mynpart, id)
    
    If(param%do_read_potential) Then
       If(.not.Allocated(pot)) Allocate(pot(mynpart))
       dsetname = 'potential_part'
       Call hdf5_read_data(gr_id, dsetname, mynpart, pot)
    End If

    If(param%do_read_gravitational_field) Then
       If(.not.Allocated(for)) Allocate(for(3,mynpart))
       dsetname = 'gravitational_field_part'
       Call hdf5_read_data(gr_id, dsetname, 3, mynpart, for)
    End If

    ! Close the root group.
    Call hdf5_close_group(gr_id)

    Call hdf5_close_file(file_id)


#ifdef DEBUG
    Print *,"Exit h5readcube on process ",procID
#endif

  End Subroutine h5readcube



  !=======================================================================
  !> Reads a hdf5 sorted cube file created by pfof_snap
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

    islast = .false.
    Call h5readcommonmetadata(file_id, islast, inforamses)

    nres = 2**inforamses%levelmin
    ngrid = int(nres,kind=8)**3
    npart = ngrid

    ! open the root group
    groupname = 'metadata'
    Call hdf5_open_group(file_id,groupname, gr_id)

    ! read attribute partNB
    dsetname = 'npart_file'
    Call hdf5_read_data(gr_id, dsetname, npart8)
    mynpart = int(npart8,kind=4)

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

    If(.not.Allocated(id)) Allocate(id(mynpart))
    If(.not.Allocated(pos))  Allocate(pos(3,mynpart))
    If(.not.Allocated(vel))  Allocate(vel(3,mynpart))

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
          Call hdf5_read_data(gr_id, dsetname, 3, mynpart, pos(:,indbeg:indend))
          
          ! read velocity of the particles
          dsetname = 'velocity_part'
          Call hdf5_read_data(gr_id, dsetname, 3, mynpart, vel(:,indbeg:indend))
          
          ! read id of the particles
          dsetname = 'identity_part'
          Call hdf5_read_data(gr_id, dsetname, mynpart, id(indbeg:indend))

          ! read potential if requested
          If(param%do_read_potential) Then
             If(.not.Allocated(pot)) Allocate(pot(mynpart))
             dsetname = 'potential_part'
             Call hdf5_read_data(gr_id, dsetname, mynpart, pot(indbeg:indend))
          End If

          ! read force if requested
          If(param%do_read_gravitational_field) Then
             If(.not.Allocated(for)) Allocate(for(3,mynpart))
             dsetname = 'gravitational_field_part'
             Call hdf5_read_data(gr_id, dsetname, 3, mynpart, for(:,indbeg:indend))
          End If
          
          indbeg = indend + 1

          Call hdf5_close_group(gr_id)
       End If
    End Do

    If(indend /= mynpart) Then
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
  !> Reads a parallel HDF5 cube file created by pfof_snap. 
  Subroutine mpih5readcube()
    
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
    Print *,"Enter mpih5readcube on process ",procID
#endif    
    
    ! number of processes writing in the same file
    procperfile = param%gatherread_factor**3
    
    Allocate(partnb_tab(procperfile))
    Allocate(bound_tab(6,procperfile))
        
    filename = trim(param%input_path)//trim(param%part_input_file)
    length=len_trim(filename)
    Write(filename(length-7:length-3),'(I5.5)') commcolorRead

    Call hdf5_open_mpi_file(filename, MpiSubCubeRead, file_id)

    islast = .false.
    Call h5readcommonmetadata(file_id, islast, inforamses)

    nres = 2**inforamses%levelmin
    ngrid = int(nres,kind=8)**3
    npart = ngrid

    ! open the root group
    groupname = 'metadata'
    Call hdf5_open_group(file_id,groupname, gr_id)

    ! read attribute partNB
    dsetname = 'npart_file'
    Call hdf5_read_data(gr_id, dsetname, npart8)
    mynpart = int(npart8,kind=4)

    ! read attribute partNB
    aname = 'npart_cube_array'
    Call hdf5_read_attr(gr_id, aname, 6, partnb_tab)

#ifdef DEBUG
    Print *,"partNB_tab read on process ",procID
    Call Mpi_Barrier(Mpi_Comm_World, mpierr)
#endif    

    mynpart = partnb_tab(scprocIDRead+1)

    ! read attribute boundaries
    aname = 'boundaries_array'
    Call hdf5_read_attr(gr_id, aname, 6, procperfile, bound_tab)

#ifdef DEBUG
    Print *,"boundaries_tab read on process ",procID
    Call Mpi_Barrier(Mpi_Comm_World, mpierr)
#endif    

    xmin=bound_tab(1,scprocIDRead+1)
    xmax=bound_tab(2,scprocIDRead+1)
    ymin=bound_tab(3,scprocIDRead+1)
    ymax=bound_tab(4,scprocIDRead+1)
    zmin=bound_tab(5,scprocIDRead+1)
    zmax=bound_tab(6,scprocIDRead+1)

    Call hdf5_close_group(gr_id)

    If(.not.Allocated(id)) Allocate(id(mynpart))
    If(.not.Allocated(pos))  Allocate(pos(3,mynpart))
    If(.not.Allocated(vel))  Allocate(vel(3,mynpart))

    groupname = 'data'
    Call hdf5_open_group(file_id, groupname, gr_id)

    dsetname = 'position_part'
    Call hdf5_read_mpi_data(gr_id, dsetname, 3, mynpart, pos, MpiSubCubeRead)

#ifdef DEBUG
    Print *,"pos read on process ",procID
    Call Mpi_Barrier(Mpi_Comm_World, mpierr)
#endif    

    dsetname = 'velocity_part'
    Call hdf5_read_mpi_data(gr_id, dsetname, 3, mynpart, vel, MpiSubCubeRead)

#ifdef DEBUG
    Print *,"vel read on process ",procID
    Call Mpi_Barrier(Mpi_Comm_World, mpierr)
#endif    

    dsetname = 'identity_part'
    Call hdf5_read_mpi_data(gr_id, dsetname, mynpart, id, MpiSubCubeRead)

#ifdef DEBUG
    Print *,"ID read on process ",procID
    Call Mpi_Barrier(Mpi_Comm_World, mpierr)
#endif    

    If(param%do_read_potential) Then
       If(.not.Allocated(pot)) Allocate(pot(mynpart))
       dsetname = 'potential_part'
       Call hdf5_read_mpi_data(gr_id, dsetname, mynpart, pot, MpiSubCubeRead)
    End If

    ! if requested, try to read the force attribute
    If(param%do_read_gravitational_field) Then
       If(.not.Allocated(for)) Allocate(for(3,mynpart))
       dsetname = 'gravitational_field_part'
       Call hdf5_read_mpi_data(gr_id, dsetname, 3, mynpart, for, MpiSubCubeRead)
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
  !> Reads a parallel and sorted HDF5 cube file created by pfof_snap.   
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
    integer::i

#ifdef DEBUG
    Print *,"Enter mpih5readsortedcube on process ",procID
#endif    
    
    ! number of processes reading in the same file
    procperfile = param%gatherread_factor**3
    
    Allocate(partnb_tab(procperfile))
    Allocate(bound_tab(6,procperfile))
        
    filename = trim(param%input_path)//trim(param%part_input_file)
    length=len_trim(filename)
    Write(filename(length-7:length-3),'(I5.5)') commcolorRead

    ! we open the file for a serial read
    Call hdf5_open_file(filename, file_id)

    islast = .false.
    Call h5readcommonmetadata(file_id, islast, inforamses)

    nres = 2**inforamses%levelmin
    ngrid = int(nres,kind=8)**3
    npart = ngrid

    ! open the root group
    groupname = 'metadata'
    Call hdf5_open_group(file_id,groupname, gr_id)

    ! read attribute partNB
    dsetname = 'npart_file'
    Call hdf5_read_data(gr_id, dsetname, npart8)
    mynpart = int(npart8,kind=4)

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

    mynpart = partnb_tab(scprocIDRead+1)

    ! read attribute boundaries
    aname = 'boundaries_array'
    Call hdf5_read_attr(gr_id, aname, 6, procperfile, bound_tab)

#ifdef DEBUG
    Print *,"boundaries_tab read on process ",procID
    Call Mpi_Barrier(Mpi_Comm_World, mpierr)
#endif    

    xmin=bound_tab(1,scprocIDRead+1)
    xmax=bound_tab(2,scprocIDRead+1)
    ymin=bound_tab(3,scprocIDRead+1)
    ymax=bound_tab(4,scprocIDRead+1)
    zmin=bound_tab(5,scprocIDRead+1)
    zmax=bound_tab(6,scprocIDRead+1)

    ! read attribute ngroup
    aname = 'ngroup'
    Call hdf5_read_attr(gr_id, aname, ngroup)

    Allocate(npartpergroup(ngroup))
    aname = 'npart_grp_array'
    Call hdf5_read_data(gr_id, aname, ngroup, npartpergroup)

    nc = ngroup / procperfile
    firstgroup = nc * scprocIDRead + 1
    lastgroup  = nc *(scprocIDRead+1)

    Call hdf5_close_group(gr_id)

    If(.not.Allocated(id)) Allocate(id(mynpart))
    If(.not.Allocated(pos))  Allocate(pos(3,mynpart))
    If(.not.Allocated(vel))  Allocate(vel(3,mynpart))

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
          Call hdf5_read_data(gr_id, dsetname, 3, mynpart, pos(:,indbeg:indend))
          
          ! read velocity of the particles
          dsetname = 'velocity_part'
          Call hdf5_read_data(gr_id, dsetname, 3, mynpart, vel(:,indbeg:indend))
          
          ! read id of the particles
          dsetname = 'identity_part'
          Call hdf5_read_data(gr_id, dsetname, mynpart, id(indbeg:indend))
          
          ! read potential if requested
          If(param%do_read_potential) Then
             If(.not.Allocated(pot)) Allocate(pot(mynpart))
             dsetname = 'potential_part'
             Call hdf5_read_data(gr_id, dsetname, mynpart, pot(indbeg:indend))
          End If

          ! read force if requested
          If(param%do_read_gravitational_field) Then
             If(.not.Allocated(for)) Allocate(for(3,mynpart))
             dsetname = 'gravitational_field_part'
             Call hdf5_read_data(gr_id, dsetname, 3, mynpart, for(:,indbeg:indend))
          End If
          
          indbeg = indend + 1

          Call hdf5_close_group(gr_id)
       End If
    End Do


    ! Close the root group.
    Call hdf5_close_group(gr_data_id)

    Call hdf5_close_file(file_id)


        ! WARNING  TEST
    !do i=1,mynpart
    !   if(pos(3,i) < 1.0e-8) then
    !      print *,procID,' FICHIER CUBE :',filename,mynpart
    !      print *,'Z = 0 : ',i, pos(3,i),scprocIDRead,partnb_tab
    !      print*,pos(1,i),pos(2,i),pos(3,i)
    !   End if
    !end do


    
    Deallocate(partnb_tab)
    Deallocate(bound_tab)

#ifdef DEBUG
    Print *,"Exit mpih5readsortedcube on process ",procID
#endif    
        
  End Subroutine mpih5readsortedcube

End Module modreadcube
