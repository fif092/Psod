!=======================================================================
! Author: Fabrice Roy (LUTH/CNRS/Observatoire de Paris)
! and Vincent Bouillot (LUTH/Observatoire de Paris)
! Fabrice.Roy@obspm.fr
!=======================================================================

Module modio
  
  Use mpi
  Use modvariable
  Use modparam
  Use modhdf5

  Integer(kind=4), dimension(:), allocatable :: nparttab
  Integer(kind=4) :: ioerr
  Integer(kind=4) :: npartloc

  Character(len=50) :: origin
  
Contains  
  
  Subroutine readparameters
        
    Implicit none

    Character(len=15) :: codeversion
    Character(len=16) :: infoconefile
    Integer :: fileperproc
    Integer(kind=4) :: buffersize, nbbytes, b_pos
    Character, dimension(:),allocatable :: buffer

    If(procID==0) Then
       ! read parameters
       infoconefile = 'infocone.nml'
       Open(Unit=10, File=infoconefile, status='old')
       Read(10, nml=infocone)
       Read(10, nml=infofile)
       Read(10, nml=infooutput)
       Close(10)
    End If

    ! Memory allocation for the input parameters broadcast buffer
    ! strings dirname, filename, shellname + origin 
    ! 4x integer(kind=4) filenb, ffile, lfile, nstride
    ! 1x integer(kind=8) npart
    ! 4x real(kind=8)    rmax, thetay, thetaz, cubesize

    ! buffer size in bits: character
    buffersize = len(dirname) + len(filename) + len(shellname)

    ! size of integer(kind=4): we add the 4x integer(4)
    Call Mpi_Sizeof(filenb, nbbytes, mpierr)
    buffersize = buffersize + 4*nbbytes

    ! size of integer(kind=8): we add the 1x integer(8)
    Call Mpi_Sizeof(npart, nbbytes, mpierr)
    buffersize = buffersize + nbbytes

    ! size of real(kind=8): we add the 4x real(8)
    Call Mpi_Sizeof(rmax, nbbytes, mpierr)
    buffersize = buffersize + 4*nbbytes

    ! size of a logical 
    Call Mpi_Sizeof(fullsky, nbbytes, mpierr)
    buffersize = buffersize + nbbytes

    ! read version of conecreator
    Open(Unit=10, File='conecreator.version', status='old', Iostat=ioerr)
    If(ioerr/=0) Then
       Print *, 'version of conecreator not defined: file not found'
       codeversion = 'undefined'
    Else
       Read(10,*,Iostat=ioerr) codeversion
       If(ioerr/=0) Then
          Print *, 'version of conecreator not defined: file empty'
          codeversion='undefined'
       End If
       Close(10)
    End If
    origin = 'Created with conecreator version '//codeversion

    ! add origin size to buffersize
    buffersize = buffersize + len(origin)

    allocate (buffer(0:buffersize-1))
    b_pos = 0

    If(procID==0) Then
       ! Pack input parameters in the buffer
       Call Mpi_Pack(   dirname,   len(dirname), Mpi_Character, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(  filename,  len(filename), Mpi_Character, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack( shellname, len(shellname), Mpi_Character, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(    origin,    len(origin), Mpi_Character, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr) 
       Call Mpi_Pack(    filenb,              1,   Mpi_Integer, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(     ffile,              1,   Mpi_Integer, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(     lfile,              1,   Mpi_Integer, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(   nstride,              1,   Mpi_Integer, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(     npart,              1,  Mpi_Integer8, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(   fullsky,              1,   Mpi_Logical, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(  cubesize,              1,     Mpi_Real8, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(      rmax,              1,     Mpi_Real8, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(    thetay,              1,     Mpi_Real8, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(    thetaz,              1,     Mpi_Real8, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
    End If

    ! Broadcast of input parameters
    Call Mpi_Bcast(buffer,buffersize,Mpi_Packed,0,Mpi_Comm_World,mpierr)

    ! Processes with ID != 0 unpack the input parameters
    If(procID /= 0) Then
       Call Mpi_Unpack(buffer, buffersize, b_pos,   dirname,   len(dirname), Mpi_Character, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,  filename,  len(filename), Mpi_Character, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos, shellname, len(shellname), Mpi_Character, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,    origin,    len(origin), Mpi_Character, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,    filenb,              1,   Mpi_Integer, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,     ffile,              1,   Mpi_Integer, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,     lfile,              1,   Mpi_Integer, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,   nstride,              1,   Mpi_Integer, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,     npart,              1,  Mpi_Integer8, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,   fullsky,              1,   Mpi_Logical, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,  cubesize,              1,     Mpi_Real8, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,      rmax,              1,     Mpi_Real8, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,    thetay,              1,     Mpi_Real8, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,    thetaz,              1,     Mpi_Real8, Mpi_Comm_World, mpierr)
    End If

    ! broadcast buffer deallocated
    Deallocate(buffer)

    Allocate(filelist(filenb))
    If(procID == 0) Then
       Call retrievefilelist()
    End If

    Call Mpi_Bcast(filelist, len(filelist(1))*filenb, Mpi_Character, 0, Mpi_Comm_World, mpierr)

    fileperproc = filenb / procNB
    ffile = procID * fileperproc + 1

    If( procID < mod(filenb,procNB) ) Then
       fileperproc = fileperproc + 1
       ffile = ffile + procID
    Else If (mod(filenb,procNB) > 0) Then
       ffile = ffile + mod(filenb,procNB) 
    End If

    lfile = ffile + fileperproc - 1


  End Subroutine readparameters

  !=======================================================================

  Subroutine retrievefilelist()

    Character(len=200) :: conefilename
    Character(len=5) :: ifilechar
    Integer :: indfile
    Integer :: ifile
    Logical :: fileexist

    indfile = 0
    Do ifile = ffile, lfile
       Write(ifilechar(1:5),'(I5.5)') ifile
       conefilename = trim(filename)//'_'//ifilechar//'.dat'
       Inquire(File=trim(dirname)//'/'//conefilename, Exist=fileexist)
       If(fileexist) Then
          indfile = indfile + 1
          filelist(indfile) = conefilename
       End If
    End Do

    If(indfile /= filenb) Then 
       Print *, 'Error while retrieving the list of Ramses cone files name.'
    End If

  End Subroutine retrievefilelist

  !=======================================================================

  Subroutine readramsesfiles()

    Implicit None

    Integer :: ifile
    Character(len=401) :: filedat
    Character(len=401) :: filehdr

    Allocate(nparttab(ffile:lfile))
    
    Do ifile = ffile, lfile
       filedat = trim(dirname)//'/'//trim(filelist(ifile))
       filehdr = filedat
       Write(filehdr(len_trim(filedat)-2:len_trim(filedat)),'(A)') 'hdr'
       Call readconehdr(filehdr, npartloc)
       nparttab(ifile) = npartloc
    End Do

    npartloc = sum(nparttab) 

    Allocate(pos(3,npartloc), vel(3,npartloc), id(npartloc))

    Do ifile = ffile, lfile
       filedat = trim(dirname)//'/'//trim(filelist(ifile))
       Call readconedat(filedat, ifile)
    End Do

  End Subroutine readramsesfiles

  !=======================================================================

  Subroutine readconehdr(file, npartloc)

    Character(len=401), intent(in) :: file
    Integer(kind=4), intent(out) :: npartloc
    
    Integer(kind=4) :: nloc

    Open(Unit=10, File=trim(file), Form='unformatted', Status='old', iostat=ioerr)
    Read(10,iostat=ioerr) nloc
    If(ioerr/=0) Print *,'Error while reading nloc in '//trim(file)
    Read(10,iostat=ioerr) nstride
    If(ioerr/=0) Print *,'Error while reading nstride in '//trim(file)
    Read(10,iostat=ioerr) npartloc
    If(ioerr/=0) Print *,'Error while reading npartloc in '//trim(file)
    Close(10)

  End Subroutine readconehdr

  !=======================================================================

  Subroutine readconedat(file,ifile)

    Implicit None

    Character(len=401), intent(in) :: file
    Integer(kind=4), intent(in) :: ifile
    Integer :: pp, fp
    Integer :: npartloc
    Integer :: nblocs, rest
    Integer :: ibloc, idim, i
    Integer :: ndim 
    Real(kind=4), dimension(:), allocatable :: tmpsimple
!!!!!!!!!!!!!!!!!!!
    Integer(kind=8), dimension(:), allocatable :: tmpint
!!!!!!!!!!!!!!!!!!!

    ndim = 3
    npartloc = nparttab(ifile)
    fp = 0
    Do i = ffile, ifile-1
       fp = fp + nparttab(i)
    End Do
    nblocs = npartloc/nstride
    rest = npartloc - nblocs*nstride

    Open(Unit=10, File=trim(file), Status='old', Form='unformatted', iostat=ioerr)
    Read(10)
    Read(10)
    Read(10)
    
    Do ibloc = 0,nblocs-1
       pp = fp + ibloc*nstride + 1
       Allocate(tmpsimple(nstride))
       Do idim = 1, ndim
          Read(10) tmpsimple
          pos(idim,pp:pp+nstride-1) = tmpsimple
          Read(10) tmpsimple
          vel(idim,pp:pp+nstride-1) = tmpsimple
       End Do
       Read(10)


!!!!!!!!!!!!!!!!!!!
       Allocate(tmpint(nstride))
       Read(10) tmpint
       Read(10) tmpsimple
       Do idim = 1,ndim
          Read(10) tmpsimple
       End Do
       Deallocate(tmpint)
!!!!!!!!!!!!!!!!!!!

       Deallocate(tmpsimple)
    End Do
    
    If(rest > 0)  Then
       pp = fp + nblocs*nstride + 1
       Allocate(tmpsimple(rest))
       Do idim = 1,ndim
          Read(10) tmpsimple
          pos(idim,pp:pp+rest-1) = tmpsimple
          Read(10) tmpsimple
          vel(idim,pp:pp+rest-1) = tmpsimple
       End Do
       Read(10) 

!!!!!!!!!!!!!!!!!!!
       Allocate(tmpint(rest))
       Read(10) tmpint
       Read(10) tmpsimple
       Do idim = 1,ndim
          Read(10) tmpsimple
       End Do
       Deallocate(tmpint)
!!!!!!!!!!!!!!!!!!!

       Deallocate(tmpsimple)
    End If

    Close(10)
    
  End Subroutine readconedat

  !=======================================================================

  Subroutine h5writecone()

    Implicit none

    Character(len=8) :: charic
    Character(len=16) :: name
    Integer(kind=hid_t) :: file_id
    Integer(kind=hid_t) :: gr_id
    Integer :: ic
    Integer :: ipointer
    Integer :: np
    Integer(kind=4) :: isfullsky
    Integer(kind=4), dimension(3) :: nctab
    Real(kind=8), dimension(3) :: dimtab

    Call hdf5_create_mpi_file(shellname,Mpi_Comm_World,file_id, origin)

    name = 'npartcube'
    If(procID == 0) Then
       Call hdf5_write_mpi_data(file_id, name, ncube, npartcube, Mpi_Comm_World)
    Else
       Call hdf5_write_mpi_data(file_id, name, 0, npartcube, Mpi_Comm_World)
    End If

    isfullsky = 0
    name = 'fullsky'
    If(fullsky) isfullsky = 1
    Call hdf5_write_attr(file_id, name, isfullsky)

    nctab(1) = ncx
    nctab(2) = ncy
    nctab(3) = ncz
    name = 'nctab'
    Call hdf5_write_attr(file_id, name, 3, nctab)

    name = 'cubesize'
    Call hdf5_write_attr(file_id, name, cubesize)

    dimtab(1) = rmax
    dimtab(2) = hy
    dimtab(3) = hz
    name = 'dimtab'
    Call hdf5_write_attr(file_id, name, 3, dimtab)

    ipointer = 1
    Do ic = 1, ncube
       If(npartcube(ic) /= 0) Then
          Write(charic(1:8),'(I8.8)') ic
          name = 'cube'//charic
          Call hdf5_create_group(file_id, name, gr_id)
          np = npartcubeloc(ic)
          name = 'pos'
          Call hdf5_write_mpi_data(gr_id, name, 3, np, pos(:,ipointer:ipointer+np-1), Mpi_Comm_World)
          name = 'vel'
          Call hdf5_write_mpi_data(gr_id, name, 3, np, vel(:,ipointer:ipointer+np-1), Mpi_Comm_World)
          ipointer = ipointer + np

          Call hdf5_close_group(gr_id)
       End If
    End Do

    Call hdf5_close_mpi_file(file_id)

  End Subroutine h5writecone


End Module modio

