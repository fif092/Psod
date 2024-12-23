!==============================================================================
! Project: pFoF
! File: pfof_cone/src/modio.f90
! Copyright Fabrice Roy and Vincent Bouillot (2011)
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
!!This file contains the subroutines and common variables used for I/O in pFoF_cone.

!> This module contains the subroutines and common variables used for I/O in pFoF_cone.
!>
!> Authors: F. Roy, V. Bouillot


Module modio
  
  Use modmpicommons
  Use modvarcommons
  Use modconstant
  Use modhdf5
  Use modvariables
  Use mpi

  Integer(kind=4) :: ioerr

  Private
  Public :: title, &
       readparameters, &
       h5readparticles, &
       h5readmap, &
       h5readshellinfo
  
Contains
  Subroutine title()
    
    Character(len=12) :: charpnb

    Write(charpnb(1:12),*) procNB
    
    Print *,' '
    Print *,'        /\ \       /\ \       /\ \         /\ \    '
    Print *,'       /  \ \     /  \ \     /  \ \       /  \ \   '
    Print *,'      / /\ \ \   / /\ \ \   / /\ \ \     / /\ \ \  '
    Print *,'     / / /\ \_\ / / /\ \_\ / / /\ \ \   / / /\ \_\ '
    Print *,'    / / /_/ / // /_/_ \/_// / /  \ \_\ / /_/_ \/_/ '
    Print *,'   / / /__\/ // /____/\  / / /   / / // /____/\    '
    Print *,'  / / /_____// /\____\/ / / /   / / // /\____\/    '
    Print *,' / / /      / / /      / / /___/ / // / /          '
    Print *,'/ / /      / / /      / / /____\/ // / /           '
    Print *,'\/_/       \/_/       \/_________/ \/_/            '
    Print *,' '
    Print *,'Code written by F.Roy and V.Bouillot'
    Print *,'based on a serial implementation written by E.Audit'
    Print *,'(see A&A 564, A13 (2014))'
    Print *,' '
    Print *,'Cone version run'
    Print *,'Number of processes: '//trim(adjustl(charpnb))
    Print *,' '

  End Subroutine title

  !=======================================================================

  Subroutine theend()

    Implicit none
    
    Print *,' '
    Print *,'Run Completed!'
    Print *,' '

  End Subroutine theend

  !=======================================================================
  !> This subroutine reads the parameters from the file pfof_cone.nml
  Subroutine readparameters()
    
    Use modvariables, only : param
    Use modmpicommons
    Use modreadparameters
    
    Implicit none

    Character(len=400) :: filename
    Integer :: ioerr
    Character(len=84) :: filelog 
    Integer(kind=4) :: mpierr
    Character(len=500) :: errormessage
    

    ! The process 0 read the input parameters and pack them for the broadcast.
    If(procID==0) Then
       filename='pfof_cone.nml'
       Call read_pfof_cone_parameters(filename, param, ioerr, errormessage, .true.)
       
       If(ioerr>0) Then
          Print *,errormessage
          Call EmergencyStop(errormessage,ioerr)
       End If

       Print *,'Parallel FoF lightcone version'
       Print *,procNB,' processes:'

       ! Open log file
       filelog = 'pfof_log_'//trim(param%simulation_name)//'.log'
       Open(Unit=Ulog,file=filelog)
       Write(Ulog,*) 'Parallel FoF'
       Write(Ulog,*) procNB,' processes:'

    End If

    Call create_mpi_type_param_pfof_cone()
    
    Call Mpi_Bcast(param, 1, Mpi_Type_parameter_pfof_cone, 0, Mpi_Comm_World, mpierr)

    ! nb of shell files 
    shell_nb = param%shell_last_id - param%shell_first_id + 1
    
  End Subroutine readparameters


  !=======================================================================
  !> This subroutine reads the cone size and the number of particles in each
  !> cubic group from each shell file, then writes an index of the cubes found
  !> in each shell file
  Subroutine h5readshellinfo()

    Use modiocommons
    Implicit None
    
    Character(len=H5STRLEN) :: dataname
    Character(len=H5STRLEN) :: groupname
    Character(len=5) :: charish
    Character(len=400) :: shellname
    Integer(kind=4) :: ic
    Integer(kind=4), dimension(3) :: nctab
    Integer(kind=4) :: ish, ishell
    Integer(kind=4), dimension(:), allocatable :: nptemp
    Integer(kind=4), dimension(:,:), allocatable :: indextemp
    Integer(kind=4) :: my_shell_fid, my_shell_lid, my_shell_nb
    Integer(kind=4), dimension(:,:), allocatable :: indexfile
    Integer(kind=4) :: icube
    Integer(kind=hid_t) :: file_id
    Integer(kind=hid_t) :: gr_meta_id
    Integer(kind=4) :: rootlast, roottmp
    Logical(kind=4) :: islast
    Integer(kind=4) :: mpierr


    ! initialization
    rootlast = 0
    islast = .false.

    Call create_mpi_type_info_ramses()
    Call create_mpi_type_info_cone_part()

    
    If(procID == 0) Then
       Print *,'Building index of dataset to read from shell files'
       Print *,' '
    End If


    ! 3 cases:
    If(shell_nb == procNB) Then   ! case 1: exactly one shell file per process
       my_shell_fid = param%shell_first_id + procID
       my_shell_lid = my_shell_fid
       my_shell_nb = 1
    Else If(shell_nb > procNB) Then ! case 2: one or more shell file per process
       my_shell_nb = shell_nb / procNB
       If(procID < mod(shell_nb, procNB)) Then
          my_shell_nb = my_shell_nb + 1
          my_shell_fid = param%shell_first_id + procID*my_shell_nb
          my_shell_lid = my_shell_fid + my_shell_nb - 1
       Else
          my_shell_fid = param%shell_first_id + mod(shell_nb,procNB) + procID*my_shell_nb
          my_shell_lid = my_shell_fid + my_shell_nb - 1
       End If
    Else ! case 3: some processes don't read any shell file here
       If(procID < shell_nb) Then
          my_shell_fid = param%shell_first_id + procID
          my_shell_lid = my_shell_fid
          my_shell_nb = 1
       Else
          my_shell_fid = 0
          my_shell_lid = -1
          my_shell_nb = 0
       End If
    End If
    
    Do ish = my_shell_fid, my_shell_lid
       
       Write(charish(1:5),'(I5.5)') ish
       shellname = trim(param%input_path)//trim(param%part_input_file)//'_'//charish//'.h5'
       Call hdf5_open_file(shellname, file_id)
       

       !Call hdf5_read_info_ramses
       !  Subroutine h5read_meta_info_ramses(file_id, inforamses)
       Call h5read_meta_common(file_id, common_meta)
       Call h5read_meta_info_cone(file_id, infocone, islast)
       Call h5read_meta_info_ramses(file_id, inforamses, islast)
       nres = 2**inforamses%levelmin

       ! this is the last shell: the metadata read from this shell will be written in the output files
       If(ish==param%shell_last_id) Then
          inforamseslast = inforamses
          infoconelast = infocone
          rootlast = procID
       End If
              
       groupname = 'metadata'
       Call hdf5_open_group(file_id, groupname, gr_meta_id)

       dataname = 'ncube_array'
       Call hdf5_read_attr(gr_meta_id, dataname, 3, nctab)
    
       ncube = nctab(1)*nctab(2)*nctab(3)
       If(.not.Allocated(indextemp)) Then
          Allocate(indextemp(ncube,shell_nb))
          indextemp = 0
       End If

       
       If(.not.Allocated(nptemp)) Allocate(nptemp(ncube))
       dataname = 'npart_cube_array'
       Call hdf5_read_data(gr_meta_id, dataname, ncube, nptemp)

       Do ic = 1, ncube
          If(nptemp(ic) /= 0) Then
             ishell = ish - param%shell_first_id + 1
             indextemp(ic,ishell) = ish
          End If
       End Do

       Call hdf5_close_group(gr_meta_id)
       Call hdf5_close_file(file_id)
    End Do
    
    If(procNB > shell_nb) Then
       Call Mpi_Bcast(ncube, 1, Mpi_Integer, 0, Mpi_Comm_World, mpierr)
       If(.not.Allocated(indextemp)) Then
          Allocate(indextemp(ncube,shell_nb))
          indextemp = 0
       End If
    End If
    
    Call Mpi_Allreduce(rootlast, roottmp, 1, Mpi_Integer, Mpi_Sum, Mpi_Comm_World, mpierr)
    rootlast = roottmp
    Call Mpi_Bcast(inforamseslast, 1, Mpi_Type_info_ramses, rootlast, Mpi_Comm_World, mpierr)
    Call Mpi_Bcast(infoconelast, 1, Mpi_Type_info_cone_part, rootlast, Mpi_Comm_World, mpierr)

    If(procID >= shell_nb) Then
       inforamses = inforamseslast
       infocone = infoconelast
       nres = 2**inforamses%levelmin       
    End If

    
    Allocate(indexfile(ncube,shell_nb))
    Call Mpi_Allreduce(indextemp, indexfile, ncube*shell_nb, Mpi_Integer, Mpi_Sum, Mpi_Comm_World, mpierr)

    Deallocate(indextemp)

    Allocate(indexcube(shell_nb, shellcubes_size))
    Allocate(cubepershell(shell_nb))
    cubepershell = 0
    Do ic = 1, shellcubes_size
       icube = ictable(ic)
       Do ish = 1, shell_nb
          If(indexfile(icube,ish) /= 0) Then
             cubepershell(ish) = cubepershell(ish)+1
             indexcube(ish,cubepershell(ish)) = icube
          End If
       End Do
    End Do


    Deallocate(indexfile)


  End Subroutine h5readshellinfo


  !=======================================================================
  !> This subroutine reads the map file 
  Subroutine h5readmap()

    Use modmpicom

    Implicit None
    
    Character(len=400) :: filename
    Character(len=H5STRLEN) :: dataname
    Character(len=H5STRLEN) :: groupname
    Character(len=8) :: charpid
    Integer(kind=hid_t) :: file_id
    Integer(kind=hid_t) :: gr_id
    Integer(kind=4), dimension(:), allocatable :: ictmp
    Integer(kind=4) :: mpierr

    If(procID==0) Then
       Print *,'Reading process mapping generated by conemapper from procmap.h5 file'
       Print *,' '
    End If

    filename = 'procmap.h5'
    Call hdf5_open_file(filename, file_id)

    dataname = 'ncube_shell'
    Call hdf5_read_attr(file_id,dataname,shellcubes_size)

    dataname = 'commdims'
    Call hdf5_read_attr(file_id,dataname,3, info_proc%global_comm%dims)

    Allocate(ictmp(shellcubes_size))

    Write(charpid(1:8),'(I8.8)') procID
    groupname = 'process_'//charpid
    Call hdf5_open_group(file_id, groupname, gr_id)

    dataname = 'npart_process'
    Call hdf5_read_attr(gr_id, dataname, local_npart)

    dataname = 'boundaries'
    Call hdf5_read_attr(gr_id, dataname, 6, boundaries)

    dataname = 'neighbours'
    Call hdf5_read_data(gr_id, dataname, 6, info_proc%global_comm%neighbours)

    dataname = 'coords'
    Call hdf5_read_data(gr_id, dataname, 3, info_proc%global_comm%coords)

    dataname = 'idc_array'
    Call hdf5_read_data(gr_id, dataname, shellcubes_size, ictmp)

    Call hdf5_close_group(gr_id)

    Call hdf5_close_file(file_id)

    ! on enleve les 0
    shellcubes_size = count(ictmp/=0)
    Allocate(ictable(shellcubes_size))
    ictable(:) = ictmp(1:shellcubes_size)

    Allocate(npart_tab(procNB))
    Call Mpi_Allgather(local_npart, 1, Mpi_Integer, npart_tab, 1, Mpi_Integer, Mpi_Comm_World, mpierr) 
    Deallocate(ictmp)

    Call settopology() 

  End Subroutine h5readmap

  !=======================================================================
  !> This subroutine reads the particles from the shell files 
  Subroutine h5readparticles()

    Implicit None
    
    Character(len=8) :: charic
    Character(len=5) :: charish
    Character(len=400) :: shellname
    Character(len=H5STRLEN) :: dataname
    Character(len=H5STRLEN) :: groupname
    
    Integer(kind=hid_t) :: file_id
    Integer(kind=hid_t) :: gr_id
    Integer(kind=hid_t) :: gr_meta_id
    Integer(kind=hid_t) :: gr_data_id
    Integer(kind=4) :: ic, ish
    Integer(kind=4) :: deb, fin, n
    Integer(kind=4) :: ip
    Integer(kind=PRI) :: fid

    Integer(kind=4), dimension(:), allocatable :: npartloc

    Allocate(position(3, local_npart))
    Allocate(velocity(3, local_npart))
    Allocate(pfof_id(local_npart))
    If(param%do_read_ramses_part_id) Allocate(ramses_id(local_npart))
    If(param%do_read_potential) Allocate(potential(local_npart))
    If(param%do_read_gravitational_field) Allocate(field(3,local_npart))

    Allocate(npartloc(ncube))

    If(procID==0) Then
       Print *,'Reading particles position and velocity from shell files'
       Print *,' '
    End If

#ifdef DEBUG
    n = 0
    ! check
    Do ish=1, shell_nb
       Write(charish(1:5),'(I5.5)') ish + param%shell_first_id - 1
       shellname = trim(param%input_path)//trim(param%part_input_file)//'_'//charish//'.h5'
       Call hdf5_open_file(shellname, file_id)

       groupname='metadata'
       Call hdf5_open_group(file_id, groupname, gr_meta_id)
       
       !Read nb of part in each cube for this file
       dataname = 'npart_cube_array'
       Call hdf5_read_data(gr_meta_id, dataname, ncube, npartloc)
       
       Do ic = 1, shellcubes_size
          n = n + npartloc(ictable(ic))
       End Do
       
       Call hdf5_close_group(gr_meta_id)
       
       Call hdf5_close_file(file_id)
    End Do
    Print *,'Check particles number before reading:',n, ' on process',procID
#endif
    

    deb = 1
    ! Loop over the files found in index
    Do ish = 1, shell_nb
       If(cubepershell(ish) /= 0) Then
          Write(charish(1:5),'(I5.5)') ish + param%shell_first_id - 1
          shellname = trim(param%input_path)//trim(param%part_input_file)//'_'//charish//'.h5'
          Call hdf5_open_file(shellname, file_id)

          groupname='metadata'
          Call hdf5_open_group(file_id, groupname, gr_meta_id)
       
          !Read nb of part in each cube for this file
          dataname = 'npart_cube_array'
          Call hdf5_read_data(gr_meta_id, dataname, ncube, npartloc)

          Call hdf5_close_group(gr_meta_id)
          groupname='data'
          Call hdf5_open_group(file_id, groupname, gr_data_id)

          Do ic = 1, cubepershell(ish)
             Write(charic(1:8),'(I8.8)') indexcube(ish,ic)
             groupname = 'cube'//charic
             Call hdf5_open_group(gr_data_id, groupname, gr_id)
          
             dataname = 'position_part'
             n = npartloc(indexcube(ish,ic))
             fin = deb + n - 1          
             Call hdf5_read_data(gr_id, dataname, 3, n, position(:,deb:fin))
             dataname = 'velocity_part'
             Call hdf5_read_data(gr_id, dataname, 3, n, velocity(:,deb:fin))
          
             If(param%do_read_ramses_part_id) Then
                dataname='identity_part_ramses'
                Call hdf5_read_data(gr_id, dataname, n, ramses_id(deb:fin))
             End If

             If(param%do_read_potential) Then
                dataname='potential_part'
                Call hdf5_read_data(gr_id, dataname, n, potential(deb:fin))
             End If

             If(param%do_read_gravitational_field) Then
                dataname='gravitational_field_part'
                Call hdf5_read_data(gr_id, dataname, 3, n, field(:,deb:fin))
             End If

             deb = fin + 1
             Call hdf5_close_group(gr_id)

          End Do

          Call hdf5_close_group(gr_data_id)

          Call hdf5_close_file(file_id)
       End If
    End Do  

#ifdef DEBUG
    Print '(I8,A,I8,A,I6)', fin, ' / ', local_npart, ' particles read from shell files on process ',procID 
#endif

    ! Attribution d'un ID unique allant de 1 a npartglobal
    global_npart = 0
    fid = 0
    Do ip = 1, procNB
       global_npart = global_npart + npart_tab(ip)
    End Do
    Do ip = 1, procID
       fid = fid + npart_tab(ip)
    End Do

    Do ip = 1, local_npart
       pfof_id(ip) = fid + ip
    End Do

    If(procID==0) Then
       Print *, 'Total number of particles in the cone: ', global_npart
       Print *, ' '
    End If

    Deallocate(npartloc)    
    
  End Subroutine h5readparticles


  
End Module modio
