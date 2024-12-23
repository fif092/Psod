!> @file
!!This file contains the subroutines and common variables used for I/O in pFoF_cone.

!> This module contains the subroutines and common variables used for I/O in pFoF_cone.
!>
!> Authors: F. Roy, V. Bouillot


Module modio
  
  Use modparameters
  Use modhdf5
  Use modvariables
  Use modmpicom

  Integer(kind=4) :: ncube
  Integer(kind=4) :: ioerr
  Integer(kind=4), dimension(:), allocatable :: npart_tab

Contains
  Subroutine title()

    Implicit none
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
    
    Use modparameters
    Use modmpicommons
    Use modmpicom
    
    Implicit none

    ! Input parameters
    Character(len=200) :: input_path
    Character(len=200) :: part_input_file
    Logical(kind=4)    :: do_read_ramses_part_id
    Logical(kind=4)    :: do_read_potential
    Logical(kind=4)    :: do_read_gravitational_field
    Integer(kind=4)    :: shell_first_id
    Integer(kind=4)    :: shell_last_id

    ! Friend of friend parameters
    Real(kind=4)       :: percolation_length
    Integer(kind=4)    :: mmin
    Integer(kind=4)    :: mmax
    Logical(kind=4)    :: do_unbinding
    Logical(kind=4)    :: do_subHalo
    
    ! Output parameters
    Character(len=200) :: simulation_name
    Logical(kind=4)    :: do_timings
    Logical(kind=4)    :: do_gather_halo
    
    Integer :: ioerr
    Character(len=84) :: fileopa 
    Character(len=84) :: filelog 
    
    ! Namelist for input file
    Namelist / input_parameters  / input_path, part_input_file, do_read_ramses_part_id, &
         do_read_potential, do_read_gravitational_field, shell_first_id, shell_last_id
    Namelist / fof_parameters    / percolation_length, mmin, mmax, do_unbinding, do_subhalo
    Namelist / output_parameters / simulation_name, do_gather_halo, do_timings
    

    ! The process 0 read the input parameters and pack them for the broadcast.
    If(procID==0) Then
       
       ! Read input parameters in file 'pfof_cone.nml'
       Open(10, file='pfof_cone.nml', iostat=ioerr) 
       
       If(ioerr>0) Then
          Print *,'** Error opening input file pfof_cone.nml. Please check this file. **'
          Call Mpi_Abort(Mpi_Comm_World,21,mpierr)
       End If
       
       Read(10, nml=input_parameters)
       Read(10, nml=fof_parameters)
       Read(10, nml=output_parameters)
       Close(10)
       
       Print *,'Parallel FoF'
       Print *,procNB,' processes:'

       param%input_path = input_path 
       param%part_input_file = part_input_file 
       param%do_read_ramses_part_id = do_read_ramses_part_id 
       param%do_read_potential = do_read_potential 
       param%do_read_gravitational_field = do_read_gravitational_field
       param%shell_first_id = shell_first_id
       param%shell_last_id = shell_last_id
       param%percolation_length = percolation_length 
       param%mmin = mmin 
       param%mmax = mmax 
       param%do_unbinding = do_unbinding 
       param%do_subhalo = do_subhalo
       param%simulation_name = simulation_name
       param%do_gather_halo = do_gather_halo
       param%do_timings = do_timings


       ! Open file .nml and write input parameters
       fileopa = 'pfof_parameters_'//trim(simulation_name)//'.nml'
       Open(Unit=Uopa,file=fileopa)
       Write(Uopa, nml=input_parameters)
       Write(Uopa, nml=fof_parameters)
       Write(Uopa, nml=output_parameters)
       Close(Uopa)    

       ! Open log file
       filelog = 'pfof_log_'//trim(param%simulation_name)//'.log'
       Open(Unit=Ulog,file=filelog)
       Write(Ulog,*) 'Parallel FoF'
       Write(Ulog,*) procNB,' processes:'
       
       ! Write input parameters in .log file
       Write(Ulog,*) 'Input parameters'
       Write(Ulog, nml=input_parameters)
       Write(Ulog, nml=fof_parameters)
       Write(Ulog, nml=output_parameters)
       
    End If

    Call create_mpi_type_info()
    
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
    Use modmpicommons
    Implicit None
    
    Character(len=H5STRLEN) :: dataname
    Character(len=H5STRLEN) :: groupname
    Character(len=5) :: charish
    Character(len=400) :: shellname
    Integer :: ic
    Integer, dimension(3) :: nctab
    Integer :: ish, ishell
    Integer, dimension(:), allocatable :: nptemp
    Integer, dimension(:,:), allocatable :: indextemp
    Integer(kind=4) :: my_shell_fid, my_shell_lid, my_shell_nb
    Integer(kind=4), dimension(:,:), allocatable :: indexfile
    Integer(kind=4) :: icube
    Integer(kind=hid_t) :: file_id
    Integer(kind=hid_t) :: gr_meta_id
    Integer(kind=4) :: rootlast, roottmp
    Logical(kind=4) :: islast

    ! initialization
    rootlast = 0
    islast = .false.

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
       Call h5readcommonmetadata(file_id, islast, inforamses, infocone)

       ! this is the last shell: the metadata read from this shell will be written in the output files
       If(ish==param%shell_last_id) Then
          inforamseslast = inforamses
          infoconelast = infocone
          rootlast = procID
       End If

       nres = 2**inforamses%levelmin

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
    Call Mpi_Bcast(inforamseslast, 1, Mpi_Type_inforamses, rootlast, Mpi_Comm_World, mpierr)
    Call Mpi_Bcast(infoconelast, 1, Mpi_Type_infocone_part, rootlast, Mpi_Comm_World, mpierr)

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

    Implicit None
    
    Character(len=400) :: filename
    Character(len=H5STRLEN) :: dataname
    Character(len=H5STRLEN) :: groupname
    Character(len=8) :: charpid
    Integer(kind=hid_t) :: file_id
    Integer(kind=hid_t) :: gr_id
    Integer(kind=4), dimension(:), allocatable :: ictmp
    Integer(kind=4) :: process_number

    If(procID==0) Then
       Print *,'Reading process mapping generated by conemapper from procmap.h5 file'
       Print *,' '
    End If

    filename = 'procmap.h5'
    Call hdf5_open_file(filename, file_id)

    dataname = 'process_number'
    Call hdf5_read_attr(file_id, dataname, process_number)
    
    If(process_number /= procNB) Then
       Print *,'Error: the number of processes used by pfof_cone does not match the number of processes in the procmap.h5 file'
       Call Mpi_Abort(Mpi_Comm_World, 1001, mpierr)
    End If

    dataname = 'ncube_shell'
    Call hdf5_read_attr(file_id,dataname,shellcubes_size)

    dataname = 'commdims'
    Call hdf5_read_attr(file_id,dataname,3, dims)

    Allocate(ictmp(shellcubes_size))

    Write(charpid(1:8),'(I8.8)') procID
    groupname = 'process_'//charpid
    Call hdf5_open_group(file_id, groupname, gr_id)

    dataname = 'npart_process'
    Call hdf5_read_attr(gr_id, dataname, mynpart)

    dataname = 'boundaries'
    Call hdf5_read_attr(gr_id, dataname, 6, boundaries)

    dataname = 'neighbours'
    Call hdf5_read_data(gr_id, dataname, 6, neighbours)

    dataname = 'coords'
    Call hdf5_read_data(gr_id, dataname, 3, CubeCoord)

    dataname = 'idc_array'
    Call hdf5_read_data(gr_id, dataname, shellcubes_size, ictmp)

    Call hdf5_close_group(gr_id)

    Call hdf5_close_file(file_id)

    ! on enleve les 0
    shellcubes_size = count(ictmp/=0)
    Allocate(ictable(shellcubes_size))
    ictable(:) = ictmp(1:shellcubes_size)

    Allocate(npart_tab(procNB))
    Call Mpi_Allgather(mynpart, 1, Mpi_Integer, npart_tab, 1, Mpi_Integer, Mpi_Comm_World, mpierr) 
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

    Allocate(pos(3, mynpart))
    Allocate(vel(3, mynpart))
    Allocate(id(mynpart))
    If(param%do_read_ramses_part_id) Allocate(ramsesid(mynpart))
    If(param%do_read_potential) Allocate(pot(mynpart))
    If(param%do_read_gravitational_field) Allocate(for(3,mynpart))

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
             Call hdf5_read_data(gr_id, dataname, 3, n, pos(:,deb:fin))
             dataname = 'velocity_part'
             Call hdf5_read_data(gr_id, dataname, 3, n, vel(:,deb:fin))
          
             If(param%do_read_ramses_part_id) Then
                dataname='identity_part_ramses'
                Call hdf5_read_data(gr_id, dataname, n, ramsesid(deb:fin))
             End If

             If(param%do_read_potential) Then
                dataname='potential_part'
                Call hdf5_read_data(gr_id, dataname, n, pot(deb:fin))
             End If

             If(param%do_read_gravitational_field) Then
                dataname='gravitational_field_part'
                Call hdf5_read_data(gr_id, dataname, 3, n, for(:,deb:fin))
             End If

             deb = fin + 1
             Call hdf5_close_group(gr_id)

          End Do

          Call hdf5_close_group(gr_data_id)

          Call hdf5_close_file(file_id)
       End If
    End Do  

#ifdef DEBUG
    Print '(I8,A,I8,A,I6)', fin, ' / ', mynpart, ' particles read from shell files on process ',procID 
#endif

    ! Attribution d'un ID unique allant de 1 a npartglobal
    npart = 0
    fid = 0
    Do ip = 1, procNB
       npart = npart + npart_tab(ip)
    End Do
    Do ip = 1, procID
       fid = fid + npart_tab(ip)
    End Do

    Do ip = 1, mynpart
       id(ip) = fid + ip
    End Do

    If(procID==0) Then
       Print *, 'Total number of particles in the cone: ', npart
       Print *, ' '
    End If

    Deallocate(npartloc)    
    
  End Subroutine h5readparticles


  
End Module modio
