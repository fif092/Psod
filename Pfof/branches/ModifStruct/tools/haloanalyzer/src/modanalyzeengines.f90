!=======================================================================
! Author: Fabrice Roy (LUTH/CNRS/Observatoire de Paris)
! Fabrice.Roy@obspm.fr
!=======================================================================
Module modanalyzeengines

  Use modvariable
  Use modfunctions
  Use modhdf5
  Implicit None

  Type ListOfVar
     Logical :: pos = .true.
     Logical :: vel = .true.
     Logical :: id = .true.
     Logical :: mp = .false.
     Logical :: tp = .false.
     Logical :: zp = .false.
  End type ListOfVar
  
  Abstract interface 
     Subroutine generic_function()
     End Subroutine generic_function
  End interface

Contains

  !=======================================================================
  !> Analyze a list of haloes defined in the file halolistfilename
  !! Apply the function generic_function to each halo individually
  Subroutine analyzelist(halolistfilename, myfunction, inputvarlist)

    Character(len=400), intent(in) :: halolistfilename
    Procedure(generic_function) :: myfunction
    Type(ListOfVar), intent(in), optional :: inputvarlist
    Type(ListOfVar)      :: varlist
    Character(len=391)    :: basename        ! name of the file without the number of the file and the extension
    Character(len=400)    :: filename        ! name of the file to open
    
    Integer(kind=hid_t)  :: file_id         ! hdf5 id of the opened file
    Integer(kind=hid_t)  :: gr_id           ! hdf5 id of the group
    Integer(kind=hid_t)  :: meta_id           ! hdf5 id of the group metadata
    Integer(kind=hid_t)  :: data_id           ! hdf5 id of the group data

    Character(len=H5STRLEN) :: name               ! name of the dataset/attribute
    Character(len=H5STRLEN) :: groupname          ! name of the group
    Character(len=9) :: filenbchar          ! string used for the id of the file
    Integer(kind=4) :: haloNB               ! number of halos in a file
    Integer(kind=4) :: halopartNB           ! number of particles in a halo
    Integer(kind=PRI), dimension(:), allocatable :: haloID  ! array containing the ID of each halo in the file
    Integer :: ih, ifile                    ! loop indices
    Integer :: nfile                        ! number of halo files written by pFoF 

    Integer :: fh, lh
    Integer :: halolistsize            ! number of halo in the list of halo read from the parameter fil
    Integer :: myhalolistsize

    Integer(kind=PRI), dimension(:), allocatable :: halolist    ! list of the ID of the halo that we must read
    Namelist / Size / halolistsize    ! namelist
    Namelist / File / basename        
    Namelist / List / halolist        ! namelist

    If(Present(inputvarlist)) Then
       varlist = inputvarlist
    End If

#ifdef WITHMPI
    Call Mpi_Comm_Rank(Mpi_Comm_World, procID, mpierr)
    Call Mpi_Comm_Size(Mpi_Comm_World, procNB, mpierr)
    
    If(procID == 0) Then
#endif

       ! We read the list of halo that we want to read and analyze: their number and ID
       Open(Unit=10,file=trim(halolistfilename),status='old')
       Read(10,nml=Size)
       Read(10,nml=File)
       Allocate(halolist(halolistsize+1))
       Read(10,nml=List)
       Close(10)

       ! We add a fake haloID at the end of the halo list
       halolist(halolistsize+1) = 0

       fh = 1
       lh = halolistsize
       myhalolistsize = halolistsize

#ifdef WITHMPI
    End If
    
    Call Mpi_Bcast(basename, len(basename), Mpi_Char, 0, Mpi_Comm_World, mpierr)
    Call Mpi_Bcast(halolistsize, 1, Mpi_Integer, 0, Mpi_Comm_World, mpierr)
    If(procID /= 0) Then
       Allocate(halolist(halolistsize+1))
    End If

    myhalolistsize = halolistsize / procNB
    fh = procID*myhalolistsize + 1
    If( procID < mod(halolistsize, procNB) ) Then
       myhalolistsize = myhalolistsize + 1
       fh = fh + procID
    Else
       fh = fh + mod(halolistsize,procNB)
    End If

    lh = fh + myhalolistsize - 1

    Call Mpi_Bcast(halolist, halolistsize+1, MPI_PRI, 0, Mpi_Comm_World, mpierr)
    
#endif


!    basename = 'test_halo'
    filenbchar = '_00000.h5'
    nfile = 1

    nbhaloanalyzed = myhalolistsize



    ! We assume that there is at least 1 file.
    ! We read the total number of files written by pFoF in this first file.
    ! We read what is to be read in this file, and move to the next and repeat until the last file is read. 
    ifile = 0
    ih = fh 
    Do While(ifile < nfile  .And. ih <= lh )


       Write(filenbchar(2:6),'(I5.5)') ifile
       filename = trim(basename)//trim(filenbchar)
       Call hdf5_open_file(filename, file_id)

       groupname='metadata'
       Call hdf5_open_group(file_id, groupname, meta_id)

       name = 'nfile'
       Call hdf5_read_attr(meta_id, name, nfile)

       name = 'nhalo_file'
       Call hdf5_read_attr(meta_id, name, haloNB)
       Call hdf5_close_group(meta_id)

       groupname='data'
       Call hdf5_open_group(file_id, groupname, data_id)

       name = 'identity_halo'
       Allocate(haloID(haloNB))
       Call hdf5_read_data(data_id, name, haloNB, haloID)

       ! We assume that the halo list is sorted 
       ! If the ID of the 1st halo we are looking for is greater than the ID of the last halo in the file
       ! we move to the next file
       If( halolist(ih) > haloID(haloNB) ) Then
          ifile = ifile + 1
          Deallocate(haloID)
          Cycle
       Else 
          Do While( halolist(ih) <= haloID(haloNB) .And. ih <= lh  )

             currenthaloID = halolist(ih)
             currenthalo = ih - fh + 1
             
             groupname = "halo_0000000000000000000"
             Write(groupname(6:24),'(I19.19)') currenthaloID
             Call hdf5_open_group(data_id, groupname, gr_id)

             name = 'npart_halo'
             Call hdf5_read_attr(gr_id, name, halopartNB)

             
             If(varlist%pos) Then
                Allocate(pos(3,halopartNB))
                name = 'position_part'
                Call hdf5_read_data(gr_id, name, 3, halopartNB, pos)
             End If

             If(varlist%vel) Then
                Allocate(vel(3,halopartNB))
                name = 'velocitiy_part'
                Call hdf5_read_data(gr_id, name, 3, halopartNB, vel)
             End If

             If(varlist%id) Then
                Allocate(id(halopartNB))
                name = 'identity_part'
                Call hdf5_read_data(gr_id, name, halopartNB, id)
             End If

             If(varlist%mp) Then
                Allocate(mp(halopartNB))
                name = 'mass_part'
                Call hdf5_read_data(gr_id, name, halopartNB, mp)
             End If

             Call hdf5_close_group(gr_id)

             Call myfunction()

             If(allocated(id)) Deallocate(id)
             If(allocated(pos)) Deallocate(pos)
             If(allocated(vel)) Deallocate(vel)
             If(allocated(mp)) Deallocate(mp) 
             
             ih = ih + 1

          End Do

       End If


       ! next file 
       ifile = ifile + 1
       Deallocate(haloID)

       Call hdf5_close_group(data_id)

       Call hdf5_close_file(file_id)
    
    End Do



  End Subroutine analyzelist

  !=======================================================================
  !> Analyze each halo contained in the halo file halofilename 
  !! Apply the function generic_function to each halo individually
  Subroutine analyzefile(halofilename, myfunction, inputvarlist)
    
    Character(len=400), intent(in) :: halofilename
    Procedure(generic_function) :: myfunction
    Type(ListOfVar), intent(in), optional :: inputvarlist

    Type(ListOfVar)      :: varlist
    Integer(kind=hid_t)  :: file_id         ! hdf5 id of the opened file
    Integer(kind=hid_t)  :: gr_id           ! hdf5 id of the group

    Character(len=H5STRLEN) :: name               ! name of the dataset/attribute
    Character(len=H5STRLEN) :: groupname          ! name of the group
    Integer(kind=4) :: haloNB               ! number of halos in a file
    Integer(kind=4) :: halopartNB           ! number of particles in a halo
    Integer(kind=PRI), dimension(:), allocatable :: haloID  ! array containing the ID of each halo in the file
    Integer :: ih, ifile                    ! loop index
    Integer :: ff, lf

    Integer :: myfilelistsize
    Integer :: filelistsize            ! number of halo files in the list
    Character(len=400), dimension(:), allocatable :: filelist    ! list of the filenames that we must read
    Namelist / Size / filelistsize    ! namelist
    Namelist / List / filelist        ! namelist
    
    If(Present(inputvarlist)) Then
       varlist = inputvarlist
    End If

#ifdef WITHMPI
    Call Mpi_Comm_Rank(Mpi_Comm_World, procID, mpierr)
    Call Mpi_Comm_Size(Mpi_Comm_World, procNB, mpierr)
    
    If(procID == 0) Then
#endif

       ! We read the list of halo files that we want to read and analyze: their number and name
       Open(Unit=10,file=trim(halofilename),status='old')
       Read(10,nml=Size)
       Allocate(filelist(filelistsize))
       Read(10,nml=List)
       Close(10)
       
       
       ff = 1
       lf = filelistsize
       myfilelistsize = filelistsize

#ifdef WITHMPI
    End If
    
    Call Mpi_Bcast(filelistsize, 1, Mpi_Integer, 0, Mpi_Comm_World, mpierr)
    If(procID /= 0) Then
       Allocate(filelist(filelistsize))
    End If

    myfilelistsize = filelistsize / procNB
    ff = procID*myfilelistsize + 1
    If( procID < mod(filelistsize, procNB) ) Then
       myfilelistsize = myfilelistsize + 1
       ff = ff + procID
    Else
       ff = ff + mod(filelistsize,procNB)
    End If

    lf = ff + myfilelistsize - 1

    Call Mpi_Bcast(filelist, len(filelist(1))*filelistsize, Mpi_Char, 0, Mpi_Comm_World, mpierr)
#endif


    nbhaloanalyzed = 0
    ! We must first read the attribute haloNB from each file
    Do ifile = ff, lf
       Call hdf5_open_file(filelist(ifile), file_id)
       name = 'haloNB'
       Call hdf5_read_attr(file_id, name, haloNB)
       nbhaloanalyzed = nbhaloanalyzed + haloNB
       Call hdf5_close_file(file_id)
    End Do

    currenthalo = 0

    Do ifile = ff, lf
       Call hdf5_open_file(filelist(ifile), file_id)

       name = 'haloNB'
       Call hdf5_read_attr(file_id, name, haloNB)

       name = 'haloID'
       Allocate(haloID(haloNB))
       Call hdf5_read_data(file_id, name, haloNB, haloID)
       
       Do ih = 1, haloNB
          
          currenthaloID = haloID(ih)
          currenthalo = currenthalo + 1
          
          groupname = "halo_00000000000"
          Write(groupname(6:16),'(I11.11)') currenthaloID
          Call hdf5_open_group(file_id, groupname, gr_id)
          
          name = 'halopartNB'
          Call hdf5_read_attr(gr_id, name, halopartNB)
          
          If(varlist%pos) Then
             Allocate(pos(3,halopartNB))
             name = 'pos'
             Call hdf5_read_data(gr_id, name, 3, halopartNB, pos)
          End If
          
          If(varlist%vel) Then
             Allocate(vel(3,halopartNB))
             name = 'vel'
             Call hdf5_read_data(gr_id, name, 3, halopartNB, vel)
          End If
          
          If(varlist%id) Then
             Allocate(id(halopartNB))
             name = 'ID'
             Call hdf5_read_data(gr_id, name, halopartNB, id)
          End If

          If(varlist%mp) Then
             Allocate(mp(halopartNB))
             name = 'mass'
             Call hdf5_read_data(gr_id, name, halopartNB, mp)
          End If
          
          Call hdf5_close_group(gr_id)
          
          Call myfunction()
          
          If(allocated(id)) Deallocate(id)
          If(allocated(pos)) Deallocate(pos)
          If(allocated(vel)) Deallocate(vel)
          If(allocated(mp)) Deallocate(mp)
          
       End Do
       
       Deallocate(haloID)

       Call hdf5_close_file(file_id)

    End Do

  End Subroutine analyzefile



End Module modanalyzeengines
