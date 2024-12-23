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

  Integer(kind=4), dimension(:), allocatable :: ncelltab
  Integer(kind=4), dimension(:,:), allocatable :: ncellperleveltab
  Integer(kind=4) :: ioerr
  Integer(kind=4) :: ncellloc

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
       infoconefile = 'infoconegrav.nml'
       Open(Unit=10, File=infoconefile, status='old')
       Read(10, nml=infocone)
       Read(10, nml=infofile)
       Read(10, nml=infooutput)
       Close(10)
    End If

    ! Memory allocation for the input parameters broadcast buffer
    ! strings dirname, filename, shellname + origin 
    ! 6x integer(kind=4) filenb, ffile, lfile, nstride, nlevel, levelmin
    ! 1x integer(kind=8) ncell
    ! 4x real(kind=8)    rmax, thetay, thetaz, cubesize
    ! 1x logical         fullsky

    ! buffer size in bits: character
    buffersize = len(dirname) + len(filename) + len(shellname)

    ! size of integer(kind=4): we add the 6x integer(4) and the 1x logical(4) (same size)
    Call Mpi_Sizeof(filenb, nbbytes, mpierr)
    buffersize = buffersize + 7*nbbytes

    ! size of integer(kind=8): we add the 1x integer(8)
    Call Mpi_Sizeof(ncell, nbbytes, mpierr)
    buffersize = buffersize + nbbytes

    ! size of real(kind=8): we add the 4x real(8)
    Call Mpi_Sizeof(rmax, nbbytes, mpierr)
    buffersize = buffersize + 4*nbbytes

!!$    ! size of a logical : we add the 2x logical
!!$    Call Mpi_Sizeof(fullsky, nbbytes, mpierr)
!!$    buffersize = buffersize + 1*nbbytes

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
       Call Mpi_Pack(    nlevel,              1,   Mpi_Integer, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(  levelmin,              1,   Mpi_Integer, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack(     ncell,              1,  Mpi_Integer8, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
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
       Call Mpi_Unpack(buffer, buffersize, b_pos,    nlevel,              1,   Mpi_Integer, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,  levelmin,              1,   Mpi_Integer, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,     ncell,              1,  Mpi_Integer8, Mpi_Comm_World, mpierr)
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

    Allocate(ncelltab(ffile:lfile))
    Allocate(ncellperleveltab(nlevel,ffile:lfile))
    Allocate(ncellperlevel(nlevel))

    ncellperleveltab = 0
    
    Do ifile = ffile, lfile
       filedat = trim(dirname)//'/'//trim(filelist(ifile))
       filehdr = filedat
       Write(filehdr(len_trim(filedat)-2:len_trim(filedat)),'(A)') 'hdr'
       Call readconegravhdr(filehdr, ncellloc)
       ncelltab(ifile) = ncellloc
       ncellperleveltab(:,ifile) = ncellperlevel(:)
    End Do

    ncellloc = sum(ncelltab)
    ncellperlevel(:) = sum(ncellperleveltab,2)
    If(ncellloc /= sum(ncellperlevel)) Then
       Print *,'Error in ncell for process ',procID
       Print *,ncellloc, ' should be equal to the sum of ',ncellperlevel
       Print *,ncellperleveltab
    End If
       
    If(ncellloc /= 0) Then
       Allocate(pos(3,ncellloc), acc(3,ncellloc))
       Allocate(phi(ncellloc), rho(ncellloc), refined(ncellloc))
    Else
       ! If npartloc == 0 we allocate arrays with size (3,1) to avoid difficulties later with 0-size arrays
       Allocate(pos(3,1), acc(3,1))
       Allocate(phi(1), rho(1), refined(1))
       pos = -1.
       acc = -1.
       phi = -1.
       rho = -1.
       refined = -1
    End If

    Do ifile = ffile, lfile
       filedat = trim(dirname)//'/'//trim(filelist(ifile))
       Call readconegravdat(filedat, ifile)
    End Do

  End Subroutine readramsesfiles

  !=======================================================================

  Subroutine readconegravhdr(file, ncellloc)

    Character(len=401), intent(in) :: file
    Integer(kind=4), intent(out) :: ncellloc
    
    Integer(kind=4) :: nloc

    Open(Unit=10, File=trim(file), Form='unformatted', Status='old', iostat=ioerr)
    Read(10,iostat=ioerr) nloc
    If(ioerr/=0) Print *,'Error while reading nloc in '//trim(file)
    Read(10,iostat=ioerr) nstride
    If(ioerr/=0) Print *,'Error while reading nstride in '//trim(file)
    Read(10,iostat=ioerr) ncellloc
    If(ioerr/=0) Print *,'Error while reading ncell in '//trim(file)
    Read(10,iostat=ioerr) ncellperlevel
    If(ioerr/=0) Print *,'Error while reading ncellperlevel in '//trim(file)
    Close(10)

  End Subroutine readconegravhdr

  !=======================================================================

  Subroutine readconegravdat(file,ifile)

    Implicit None

    Character(len=401), intent(in) :: file
    Integer(kind=4), intent(in) :: ifile
    Integer :: ilvl
    Integer, dimension(:), allocatable :: fp
    Integer(kind=4) :: ncellloc
    Integer :: nblocs, rest
    Integer :: ibloc, idim, i, pp
    Integer :: ndim 
    Real(kind=4), dimension(:), allocatable :: tmpsimple
    Integer(kind=4), dimension(:), allocatable :: tmpint

    ! data is stored per lvl for pos, acc, phi, rho and refined
    ! we use ncellperlevel and ncellperleveltab to point to the right place in the arrays

    Allocate(fp(nlevel))
    ndim = 3
    fp(1) = 0
    Do i = 2, nlevel
       fp(i) = fp(i-1) + ncellperlevel(i-1)
    End Do
    Do i = ffile, ifile-1
       fp(:) = fp(:) + ncellperleveltab(:,i)
    End Do

#ifdef DEBUG
       Print *,'readconegravdat'
       Print *,'process:',procID, ' fp=',fp
#endif


    Open(Unit=10, File=trim(file), Status='old', Form='unformatted', iostat=ioerr)
    Read(10)
    Read(10)
    Read(10)

    level:Do ilvl = 1, nlevel

#ifdef DEBUG
       Print *,'lvl=',ilvl
#endif

       ncellloc = ncellperleveltab(ilvl,ifile)
       nblocs = ncellloc/nstride
       rest = ncellloc - nblocs*nstride

       Allocate(tmpsimple(nstride))
       Allocate(tmpint(nstride))

       Do ibloc = 0,nblocs-1
          pp = fp(ilvl) + ibloc*nstride + 1
          Do idim = 1, ndim
             Read(10) tmpsimple
             pos(idim,pp:pp+nstride-1) = tmpsimple
             Read(10) tmpsimple
             acc(idim,pp:pp+nstride-1) = tmpsimple
          End Do
          Read(10) tmpsimple
          phi(pp:pp+nstride-1) = tmpsimple
          Read(10) tmpsimple
          rho(pp:pp+nstride-1) = tmpsimple
          Read(10) tmpint
          refined(pp:pp+nstride-1) = tmpsimple
       End Do
       Deallocate(tmpsimple)
       Deallocate(tmpint)

       If(rest > 0)  Then
          pp = fp(ilvl) + nblocs*nstride + 1
          Allocate(tmpsimple(rest))
          Allocate(tmpint(rest)) 
         Do idim = 1,ndim
             Read(10) tmpsimple
             pos(idim,pp:pp+rest-1) = tmpsimple
             Read(10) tmpsimple
             acc(idim,pp:pp+rest-1) = tmpsimple
          End Do
          Read(10) tmpsimple
          phi(pp:pp+rest-1) = tmpsimple
          Read(10) tmpsimple
          rho(pp:pp+rest-1) = tmpsimple
          Read(10) tmpint
          refined(pp:pp+rest-1) = tmpsimple
          Deallocate(tmpsimple)
          Deallocate(tmpint)
       End If

    End Do level

    Close(10)
    
  End Subroutine readconegravdat

  !=======================================================================

  Subroutine h5writeconegrav()

    Implicit none

    Character(len=8) :: charic
    Character(len=16) :: name
    Integer(kind=hid_t) :: file_id
    Integer(kind=hid_t) :: gr_id, grlvl_id
    Integer :: ic, ilvl
    Integer :: ipointer
    Integer :: np
    Integer(kind=4) :: isfullsky
    Integer(kind=4), dimension(3) :: nctab
    Real(kind=8), dimension(3) :: dimtab
    Integer(kind=4) :: beglvl, endlvl, nlevelperproc, begc, endc
    Logical :: empty

    Call hdf5_create_mpi_file(shellname,Mpi_Comm_World,file_id, origin)


    name = 'nlevel'
    Call hdf5_write_attr(file_id, name, nlevel)
    
    name = 'levelmin'
    Call hdf5_write_attr(file_id, name, levelmin)

    name = 'ncellcube'
!!$    If(procID == 0) Then
!!$       Call hdf5_write_mpi_data(file_id, name, ncube, npartcube, Mpi_Comm_World)
!!$    Else
!!$       Call hdf5_write_mpi_data(file_id, name, 0, npartcube, Mpi_Comm_World)
!!$    End If

    !! 
    empty = .false.
    If(nlevel >= procNB) Then
       nlevelperproc = nlevel / procNB
       If(procID < modulo(nlevel,procNB)) Then
          nlevelperproc = nlevelperproc + 1
          beglvl = nlevelperproc*procID + 1
          endlvl = beglvl + nlevelperproc - 1
       Else
          beglvl = nlevelperproc*procID + modulo(nlevel,procNB) + 1
          endlvl = beglvl + nlevelperproc - 1
       End If
    Else
       nlevelperproc = 1
       If(procID >= nlevel) Then
          nlevelperproc = 0
          empty = .true.
       End If
    End If
    Call hdf5_write_mpi_data(file_id,name,ncube,nlevelperproc,ncellcube(:,beglvl:endlvl),&
         Mpi_Comm_World,empty)
    !!

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
    
    !    Do ilvl = 1, nlevel
    Do ic = 1, ncube
       If(ncellcubeloc(ic,2)/=0) Print *,'PROC',procID,'  ic',ic, '  ncell',ncellcubeloc(ic,2),ncellcube(ic,2)
    End Do
    !    End Do
    
    
    Print *,'LEVELMIN', levelmin
    ! loop over the levels
    ipointer = 1
    Do ilvl = 1, nlevel
       Write(charic(1:2),'(I2.2)') levelmin+ilvl-1
       name = 'level'//charic(1:2)
       Call hdf5_create_group(file_id, name, grlvl_id)
       Do ic = 1, ncube
          If(ncellcube(ic,ilvl) /= 0) Then
             Write(charic(1:8),'(I8.8)') ic
             name = 'cube'//charic
             Call hdf5_create_group(grlvl_id, name, gr_id)
             np = ncellcubeloc(ic,ilvl)
             If(np == 0) Then
                np = 1
                begc = lbound(pos,2)
                endc = begc
                empty = .true.
             Else
                begc = ipointer
                endc = ipointer+np-1
                empty = .false.
             End If

#ifdef DEBUG
             Print *,'process:',procID,' call hdf5_write_mpi_data with arguments:',gr_id, ilvl,ic, np, begc, endc, empty
!!$             If(empty) Then
!!$                Print *,'pos and vel:',pos(:,begc:endc),acc(:,begc:endc)
!!$             End If
#endif
             name = 'pos'          
             Call hdf5_write_mpi_data(gr_id, name, 3, np, pos(:,begc:endc), Mpi_Comm_World,empty)
             name = 'acc'
             Call hdf5_write_mpi_data(gr_id, name, 3, np, acc(:,begc:endc), Mpi_Comm_World,empty)
             name = 'phi'
             Call hdf5_write_mpi_data(gr_id, name, np, phi(begc:endc), Mpi_Comm_World,empty)
             name = 'rho'
             Call hdf5_write_mpi_data(gr_id, name, np, rho(begc:endc), Mpi_Comm_World,empty)
             name = 'refined'
             Call hdf5_write_mpi_data(gr_id, name, np, refined(begc:endc), Mpi_Comm_World,empty)
             
             If(.not.empty) Then
                ipointer = ipointer + np
             End If
             
             Call hdf5_close_group(gr_id)
          End If
       End Do
       Call hdf5_close_group(grlvl_id)
    End Do

    Call hdf5_close_mpi_file(file_id)

  End Subroutine h5writeconegrav


End Module modio

