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
  !> This subroutine reads the parameters from the file pfof_cone.nml
  Subroutine readparameters()

    Implicit None
    
    Character(len=15) :: codeversion


    ! process 0 reads the parameters file
    If(procID == 0) Then
       Print *, 'Reading parameters in pfof_cone.nml file'
       Print *,' '
       Open(Unit=10, file='pfof_cone.nml', iostat=ioerr) 
       If(ioerr>0) Then
          Print *,'** Error opening input file pfof_cone.nml. Please check this file. **'
          Call Mpi_Abort(Mpi_Comm_World,ioerr,mpierr)
       End If
       Read(10, nml=shellparameters)
       Read(10, nml=fofparameters)
       Read(10, nml=outputparameters)
       Close(10)

       ! Read code version
       Open(Unit=10, File='pfof.version', status='old', Iostat=ioerr)
       If(ioerr/=0) Then
          Print *, 'Version of pfof not defined: pfof.version file not found'
          Print *,' '
          codeversion = 'undefined'
       Else
          Read(10,*,Iostat=ioerr) codeversion
          If(ioerr/=0) Then
             Print *, 'Version of pfof not defined: pfof.version file empty'
             Print *,' ' 
             codeversion='undefined'
          Else
             Print *,'Version of pfof used: '//codeversion
             Print *,' '
          End If
          Close(10)
       End If
       origin = 'Created with pfof version '//codeversion
       
       
    End If

    ! process 0 broadcasts the parameters
    Call bcastparam()

    ! nb of shell files 
    shell_nb = shell_lid - shell_fid + 1
    
  End Subroutine readparameters


  !=======================================================================
  !> This subroutine reads the cone size and the number of particles in each
  !> cubic group from each shell file, then writes an index of the cubes found
  !> in each shell file
  Subroutine h5readshellinfo()
    
    Implicit None
    
    Character(len=16) :: dataname
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

    If(procID == 0) Then
       Print *,'Building index of dataset to read from shell files'
       Print *,' '
    End If


    ! 3 cases:
    If(shell_nb == procNB) Then   ! case 1: exactly one shell file per process
       my_shell_fid = shell_fid + procID
       my_shell_lid = my_shell_fid
       my_shell_nb = 1
    Else If(shell_nb > procNB) Then ! case 2: one or more shell file per process
       my_shell_nb = shell_nb / procNB
       If(procID < mod(shell_nb, procNB)) Then
          my_shell_nb = my_shell_nb + 1
          my_shell_fid = shell_fid + procID*my_shell_nb
          my_shell_lid = my_shell_fid + my_shell_nb - 1
       Else
          my_shell_fid = shell_fid + mod(shell_nb,procNB) + procID*my_shell_nb
          my_shell_lid = my_shell_fid + my_shell_nb - 1
       End If
    Else ! case 3: some processes don't read any shell file here
       If(procID < shell_nb) Then
          my_shell_fid = shell_fid + procID
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
       shellname = trim(shell_dir)//trim(shell_filename)//'_'//charish//'.h5'
       Call hdf5_open_file(shellname, file_id)
       
       dataname = 'nctab'
       Call hdf5_read_attr(file_id, dataname, 3, nctab)
    
       ncube = nctab(1)*nctab(2)*nctab(3)
       If(.not.Allocated(indextemp)) Then
          Allocate(indextemp(ncube,shell_nb))
          indextemp = 0
       End If

       
       If(.not.Allocated(nptemp)) Allocate(nptemp(ncube))
       dataname = 'npartcube'
       Call hdf5_read_data(file_id, dataname, ncube, nptemp)

       Do ic = 1, ncube
          If(nptemp(ic) /= 0) Then
             ishell = ish - shell_fid + 1
             indextemp(ic,ishell) = ish
          End If
       End Do
       
       Call hdf5_close_file(file_id)
    End Do
    
    If(procNB > shell_nb) Then
       Call Mpi_Bcast(ncube, 1, Mpi_Integer, Mpi_Comm_World, 0, mpierr)
       If(.not.Allocated(indextemp)) Then
          Allocate(indextemp(ncube,shell_nb))
          indextemp = 0
       End If
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
    Character(len=16) :: dataname
    Character(len=16) :: groupname
    Character(len=8) :: charpid
    Integer(kind=hid_t) :: file_id
    Integer(kind=hid_t) :: gr_id
    Integer(kind=4), dimension(:), allocatable :: ictmp

    If(procID==0) Then
       Print *,'Reading process mapping generated by conemapper from procmap.h5 file'
       Print *,' '
    End If

    filename = 'procmap.h5'
    Call hdf5_open_file(filename, file_id)

    dataname = 'shellcubes_size'
    Call hdf5_read_attr(file_id,dataname,shellcubes_size)

    dataname = 'commdims'
    Call hdf5_read_attr(file_id,dataname,3, dims)

    Allocate(ictmp(shellcubes_size))

    Write(charpid(1:8),'(I8.8)') procID
    groupname = 'process_'//charpid
    Call hdf5_open_group(file_id, groupname, gr_id)

    dataname = 'my_npart'
    Call hdf5_read_attr(gr_id, dataname, mynpart)

    dataname = 'boundaries'
    Call hdf5_read_attr(gr_id, dataname, 6, boundaries)

    dataname = 'neighbours'
    Call hdf5_read_data(gr_id, dataname, 6, neighbours)

    dataname = 'shellcubes'
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
    Character(len=16) :: dataname
    Character(len=16) :: groupname
    
    Integer(kind=hid_t) :: file_id
    Integer(kind=hid_t) :: gr_id
    Integer(kind=4) :: ic, ish
    Integer(kind=4) :: deb, fin, n
    Integer(kind=4) :: ip
    Integer(kind=PRI) :: fid

    Integer(kind=4), dimension(:), allocatable :: npartloc

    Allocate(pos(3, mynpart))
    Allocate(vel(3, mynpart))
    Allocate(id(mynpart))

    Allocate(npartloc(ncube))

    If(procID==0) Then
       Print *,'Reading particles position and velocity from shell files'
       Print *,' '
    End If

#ifdef DEBUG
    n = 0
    ! check
    Do ish=1, shell_nb
       Write(charish(1:5),'(I5.5)') ish + shell_fid - 1
       shellname = trim(shell_dir)//trim(shell_filename)//'_'//charish//'.h5'
       Call hdf5_open_file(shellname, file_id)

       !Read nb of part in each cube for this file
       dataname = 'npartcube'
       Call hdf5_read_data(file_id, dataname, ncube, npartloc)
       
       Do ic = 1, shellcubes_size
          n = n + npartloc(ictable(ic))
       End Do
       
       Call hdf5_close_file(file_id)
    End Do
    Print *,'Check particles number before reading:',n, ' on process',procID
#endif
    

    deb = 1
    ! Loop over the files found in index
    Do ish = 1, shell_nb
       If(cubepershell(ish) /= 0) Then
          Write(charish(1:5),'(I5.5)') ish + shell_fid - 1
          shellname = trim(shell_dir)//trim(shell_filename)//'_'//charish//'.h5'
          Call hdf5_open_file(shellname, file_id)

          !Read nb of part in each cube for this file
          dataname = 'npartcube'
          Call hdf5_read_data(file_id, dataname, ncube, npartloc)

          Do ic = 1, cubepershell(ish)
             Write(charic(1:8),'(I8.8)') indexcube(ish,ic)
             groupname = 'cube'//charic
             Call hdf5_open_group(file_id, groupname, gr_id)
          
             dataname = 'pos'
             n = npartloc(indexcube(ish,ic))
             fin = deb + n - 1          
             Call hdf5_read_data(gr_id, dataname, 3, n, pos(:,deb:fin))
             dataname = 'vel'
             Call hdf5_read_data(gr_id, dataname, 3, n, vel(:,deb:fin))
             deb = fin + 1
          
             Call hdf5_close_group(gr_id)

          End Do

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



!!$  !=======================================================================
!!$  !> This subroutine writes for each halo the position, the velocity and the id of each particle in this halo
!!$  !! in a hdf5 file.
!!$  !! One file is written per MPI process.
!!$  !! There is one group per halo, the name of the group is the ID of the halo.
!!$  Subroutine h5writehalopart(haloNB, halopartNB, haloMass, haloID, halopartPos, halopartVel, halopartID)
!!$
!!$    Use modhdf5
!!$    Use modparameters
!!$    Use modmpicom
!!$    Use modvariables
!!$    Use modxdmf
!!$    Implicit none
!!$    
!!$    Integer(kind=4),                            intent(in) :: haloNB, halopartNB
!!$    Integer(kind=PRI), dimension(haloNB),       intent(in) :: haloID
!!$    Integer(kind=PRI), dimension(halopartNB),   intent(in), target :: halopartID
!!$    Integer(kind=4),   dimension(haloNB),       intent(in) :: haloMass
!!$    Real   (kind=4),  dimension(3,halopartNB), intent(in), target :: halopartPos, halopartVel
!!$
!!$    Integer(kind=4) :: ih, hptr
!!$    Character(len=400) :: filestrct
!!$    Character(len=390) :: filebase
!!$    Character(len=5)  :: pid_char
!!$    Character(len=16) :: groupname
!!$    Character(len=16) :: dsetname                           ! Dataset name
!!$    Character(len=16) :: aname                              ! Attribute name
!!$
!!$    Integer(hid_t) :: file_id                               ! File identifier
!!$    Integer(hid_t) :: gr_root_id                            ! Group identifier
!!$    Integer(hid_t) :: gr_halo_id                            ! Group identifier
!!$
!!$    Integer(kind=4), dimension(procNB) :: haloNBtab, dspl
!!$    Integer(kind=4), dimension(:), allocatable :: haloMasstab
!!$    Integer(kind=PRI), dimension(:), allocatable :: haloIDtab
!!$    Integer(kind=PRI), dimension(2) :: IDminmax
!!$    Integer(kind=4) :: haloNBall, p
!!$
!!$    
!!$#ifdef DEBUG
!!$    Print *,"Enter h5writehalopart on process ",procID
!!$#endif    
!!$
!!$    Write(pid_char(1:5),'(I5.5)') procID
!!$    filestrct = trim(output_root)//'_halo_'//pid_char//'.h5'
!!$    filebase = trim(output_root)//'_halo'
!!$
!!$    Call Mpi_Gather(haloNB, 1, Mpi_Integer, haloNBtab, 1, Mpi_Integer, 0, Mpi_Comm_World, mpierr)
!!$    If(procID==0) Then
!!$       haloNBall = haloNBtab(1)
!!$       dspl(1) = 0
!!$       Do p=2, procNB 
!!$          dspl(p) = haloNBall
!!$          haloNBall= haloNBall + haloNBtab(p)
!!$       End Do
!!$       Allocate(haloMasstab(haloNBall), haloIDtab(haloNBall))
!!$    Else
!!$       Allocate(haloMasstab(0), haloIDtab(0))
!!$    End If
!!$    Call Mpi_Gatherv(haloMass, haloNB, Mpi_Integer, haloMasstab, haloNBtab, dspl, Mpi_Integer, 0, Mpi_Comm_World, mpierr)
!!$    Call Mpi_Gatherv(haloID, haloNB, Mpi_PRI, haloIDtab, haloNBtab, dspl, Mpi_PRI, 0, Mpi_Comm_World, mpierr)
!!$
!!$    If(procID==0) Then
!!$       Call writesinglehaloxdmf(procNB,filebase,haloNBtab,haloMasstab,haloIDtab)
!!$    End If
!!$
!!$    ! create the hdf5 file
!!$    Call hdf5_create_file(filestrct, file_id, origin)
!!$
!!$    ! open the root group
!!$    groupname = '/'
!!$    Call hdf5_open_group(file_id,groupname, gr_root_id)
!!$
!!$    aname = 'nfile'
!!$    Call hdf5_write_attr(gr_root_id, aname, procNB)
!!$    
!!$    ! write the number of haloes as an attribute
!!$    aname = "haloNB"
!!$    Call hdf5_write_attr(gr_root_id, aname, haloNB)
!!$
!!$    !! write the halo ID as data and not attribute: it seems that we cannot write integer(kind=8) attribute 
!!$    aname = 'haloID'
!!$    Call hdf5_write_data(gr_root_id, aname, haloNB, haloID)
!!$
!!$    aname = 'haloIDminmax'
!!$    IDminmax(1) = haloID(1)
!!$    IDminmax(2) = haloID(haloNB)
!!$    Call hdf5_write_data(gr_root_id, aname, 2, IDminmax)
!!$
!!$    ! pointer to the current halo
!!$    hptr = 1
!!$    
!!$    Do ih = 1, haloNB
!!$       ! create a group for each halo
!!$       groupname = "halo_00000000000"
!!$       Write(groupname(6:16),'(I11.11)') haloID(ih)
!!$       Call hdf5_create_group(gr_root_id, groupname, gr_halo_id)
!!$       ! create an attribute containing the number of particles in the halo
!!$       aname = "halopartNB"
!!$       Call hdf5_write_attr(gr_halo_id, aname, haloMass(ih))
!!$       
!!$       dsetname="pos"
!!$       Call hdf5_write_data(gr_halo_id, dsetname, 3, haloMass(ih), halopartPos(:,hptr:hptr+haloMass(ih)-1)) 
!!$       
!!$       dsetname="vel"
!!$       Call hdf5_write_data(gr_halo_id, dsetname, 3, haloMass(ih), halopartVel(:,hptr:hptr+haloMass(ih)-1)) 
!!$       
!!$       dsetname = "ID"
!!$       Call hdf5_write_data(gr_halo_id, dsetname, haloMass(ih), halopartID(hptr:hptr+haloMass(ih)-1)) 
!!$       
!!$       ! Close the halo group.
!!$       Call hdf5_close_group(gr_halo_id)
!!$       
!!$       ! move the pointer to the next halo
!!$       hptr = hptr + haloMass(ih)
!!$    End Do
!!$    
!!$    ! Close the root group.
!!$    Call hdf5_close_group(gr_root_id)
!!$       
!!$    Call hdf5_close_file(file_id)
!!$
!!$#ifdef DEBUG
!!$    Print *,"Exit h5writehalopart on process ",procID
!!$#endif
!!$
!!$  End Subroutine h5writehalopart


!!$  !=======================================================================
!!$  !> This subroutine writes, for each halo, its mass (number of particles), 
!!$  !! the position of its center of mass and its ID
!!$  !! in only one hdf5 file using parallel HDF5.
!!$  Subroutine mpih5writehalomass(haloNB_all, haloNB, haloMass, halocomPos, halocomVel, haloID)
!!$
!!$    Use modparameters
!!$    Use modvariables
!!$    Use modmpicom
!!$    Use modhdf5
!!$
!!$    Implicit none
!!$    
!!$    Integer(kind=4),                        intent(in) :: haloNB_all
!!$    Integer(kind=4),                        intent(in) :: haloNB
!!$    Integer(kind=4),   dimension(haloNB),   intent(in), target :: haloMass
!!$    Real(kind=8),     dimension(3,haloNB), intent(in), target :: halocomPos
!!$    Real(kind=8),     dimension(3,haloNB), intent(in), target :: halocomVel
!!$    Integer(kind=PRI), dimension(haloNB),   intent(in), target :: haloID
!!$
!!$    Character(len=400) :: filename
!!$    Character(len=16) :: aname                           ! Attribute name
!!$    Character(len=16) :: dsetname                        ! Dataset name
!!$    Integer(hid_t) :: file_id                            ! File identifier
!!$    
!!$
!!$    filename = trim(output_root)//'_halomass.h5'
!!$
!!$    ! Create h5 parallel file
!!$    Call hdf5_create_mpi_file(filename, Mpi_Comm_World, file_id, origin)
!!$
!!$    ! Write number of halos as attribute
!!$    aname = 'haloNB'
!!$    Call hdf5_write_attr(file_id, aname, haloNB_all)
!!$
!!$    dsetname = 'halocomPos'
!!$    Call hdf5_write_mpi_data(file_id, dsetname, 3, haloNB, halocomPos, Mpi_Comm_World)
!!$
!!$    dsetname = 'halocomVel'
!!$    Call hdf5_write_mpi_data(file_id, dsetname, 3, haloNB, halocomVel, Mpi_Comm_World)
!!$
!!$    dsetname = 'haloID'
!!$    Call hdf5_write_mpi_data(file_id, dsetname, haloNB, haloID, Mpi_Comm_World)
!!$
!!$    dsetname = 'haloMass'
!!$    Call hdf5_write_mpi_data(file_id, dsetname, haloNB, haloMass, Mpi_Comm_World)
!!$
!!$    ! Close h5 file
!!$    Call hdf5_close_mpi_file(file_id)
!!$
!!$    
!!$  End Subroutine mpih5writehalomass


  
End Module modio
