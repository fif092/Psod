Module modio
  
  Use modparameters
  Use modhdf5
  Use modvariables
  
  Integer :: ioerr

Contains

  !=======================================================================
  !> This subroutine reads the parameters from the file pfof_cone.nml
  Subroutine readparameters()

    Implicit None
    
    Open(Unit=10, file='pfof_cone.nml', iostat=ioerr) 
    If(ioerr>0) Then
       Print *,'** Error opening input file pfof_cone.nml. Please check this file. **'
       Call exit(1)
    End If
    Read(10, nml=shellparameters)
    Read(10, nml=fofparameters)
    Read(10, nml=outputparameters)
    Close(10)
 
    ! nb of shell files 
    shell_nb = shell_lid - shell_fid + 1
    
  End Subroutine readparameters


  !=======================================================================
  !> This subroutine reads the cone size and the number of particles in each
  !> cubic group from each shell file, then perform a reduction to get the 
  !> total number of particles for each cubic group.
  Subroutine h5readshellinfo()
    
    Implicit None

    Character(len=16) :: dataname
    Character(len=5) :: charish
    Character(len=400) :: shellname
    Integer :: ish
    Integer(kind=4) :: isfullsky
    Integer, dimension(:), allocatable :: nptemp
    
    Integer(kind=hid_t) :: file_id
    
    Do ish = shell_fid, shell_lid
       
       Write(charish(1:5),'(I5.5)') ish
       shellname = trim(shell_dir)//trim(shell_filename)//'_'//charish//'.h5'
       Call hdf5_open_file(shellname, file_id)

       fullsky=.false.
       dataname = 'fullsky'
       Call hdf5_read_attr(file_id, dataname, isfullsky)
       If(isfullsky==1) fullsky=.true.
    
       dataname = 'cubesize'
       Call hdf5_read_attr(file_id, dataname, cubesize)
   
       dataname = 'nctab'
       Call hdf5_read_attr(file_id, dataname, 3, nctab)
       
       dataname = 'npartcube'
       ncube = nctab(1)*nctab(2)*nctab(3)
       
       If(.not.Allocated(npartcube)) Then
          Allocate(npartcube(ncube))
          npartcube = 0
       End If
       
       If(.not.Allocated(nptemp)) Allocate(nptemp(ncube))
       Call hdf5_read_data(file_id, dataname, ncube, nptemp)
       
       npartcube = npartcube + nptemp

       Call hdf5_close_file(file_id)
       
    End Do
    
    npart = sum(npartcube)
    Deallocate(nptemp)

  End Subroutine h5readshellinfo


  !=======================================================================
  !> This subroutine writes a map for pfof_cone, this map contains the 
  !> id of the shell cubes and the procID of the neighbours for each 
  !> pfof_cone process.
  Subroutine h5writemap(factor, procNB, dims)

    Implicit None

    Integer(kind=4), intent(in) :: factor
    Integer(kind=4), intent(in) :: procNB
    Integer(kind=4), dimension(3), intent(in) :: dims
    Character(len=400) :: filename
    Character(len=16) :: groupname
    Character(len=16) :: dataname
    Character(len=6) :: charpnb
    Character(len=4) :: charfac
    Character(len=8) :: charpid
    Character(len=50) :: origin
    Integer(kind=hid_t) :: file_id, gr_id
    Integer(kind=4) :: factorcube
    Integer(kind=4) :: ip

    origin='Origintest'
    factorcube = size(pictable,1)

    Write(charfac(1:4),'(I4.4)') factor
    Write(charpnb(1:6),'(I6.6)') procNB
    filename = 'procmap_'//charfac//'_'//charpnb//'.h5'
    Print '(A)','Name of the map file: ', trim(filename)
    Print *,' '

    Call hdf5_create_file(filename, file_id, origin)

    dataname='shellcubes_size'
    Call hdf5_write_attr(file_id, dataname, factorcube)

    dataname='commdims'
    Call hdf5_write_attr(file_id, dataname, 3, dims)
    

    Do ip = 1, procNB
       Write(charpid(1:8),'(I8.8)') ip-1
       groupname = 'process_'//charpid
       Call hdf5_create_group(file_id, groupname, gr_id)

       dataname = 'boundaries'
       Call hdf5_write_attr(gr_id, dataname, 6, boundaries(:,ip))

       dataname='my_npart'
       Call hdf5_write_attr(gr_id, dataname, my_npart_tab(ip))
       
       dataname = 'shellcubes'
       Call hdf5_write_data(gr_id, dataname, factorcube, pictable(:,ip))
       dataname = 'neighbours'
       Call hdf5_write_data(gr_id, dataname, 6, neigh(:,ip))

       Call hdf5_close_group(gr_id)
    End Do

    Call hdf5_close_file(file_id)

  End Subroutine h5writemap

End Module modio
