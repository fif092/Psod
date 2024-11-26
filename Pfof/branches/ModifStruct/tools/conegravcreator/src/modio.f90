!==============================================================================
! Project: pFoF
! File: tools/conegravcreator/src/modio.f90
! Copyright Fabrice Roy and Vincent Bouillot (2015)
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


!=======================================================================
! Author: Fabrice Roy (LUTH/CNRS/Observatoire de Paris)
! and Vincent Bouillot (LUTH/Observatoire de Paris)
! Fabrice.Roy@obspm.fr
!=======================================================================

Module modio
  
  Use mpi
  Use modvariables
  Use modhdf5

  Integer(kind=4), dimension(:), allocatable :: ncelltab
  Integer(kind=4), dimension(:,:), allocatable :: ncellperleveltab
  Integer(kind=4) :: ioerr
  Integer(kind=4) :: ncellloc
  Integer(kind=4) :: nstride

  Character(len=50) :: origin

  ! these two variables are used in several i/o routines
  Integer(kind=4) :: first_file
  Integer(kind=4) :: last_file

  
Contains  
  
  !=======================================================================
  Subroutine readparameters

    Use modmpicommons
    Implicit none

    Character(len=16) :: infoconefile
    Character(len=200) :: input_path
    Character(len=200) :: simulation_name
    Character(len=200) :: cone_input_file
    Character(len=200) :: info_cone_input_file
    Character(len=200) :: info_ramses_input_file
    Integer(kind=4) :: nfile
    Real(kind=8) :: cone_max_radius
    Real(kind=8) :: cube_size    
    Integer(kind=4) :: nlevel
    Integer(kind=4) :: levelmin
    
    Namelist/input_parameters/input_path, cone_input_file, info_cone_input_file, &
         info_ramses_input_file, nfile, first_file, last_file, cone_max_radius, nlevel, levelmin
    Namelist/output_parameters/simulation_name, cube_size

    Integer :: fileperproc

    If(procID==0) Then
       ! read parameters
       infoconefile = 'infoconegrav.nml'
       Open(Unit=10, File=infoconefile, status='old')
       Read(10, nml=input_parameters)
       Read(10, nml=output_parameters)
       Close(10)

       param%input_path = input_path
       param%cone_input_file = cone_input_file
       param%info_cone_input_file = info_cone_input_file
       param%info_ramses_input_file = info_ramses_input_file
       param%nfile = nfile
       param%first_file = first_file
       param%last_file = last_file
       param%cone_max_radius = cone_max_radius
       param%nlevel = nlevel
       param%levelmin = levelmin
       param%simulation_name = simulation_name
       param%cube_size = cube_size
       
    End If

    Call create_mpi_type_param_cone_grav()

    Call Mpi_Bcast(param, 1, Mpi_Type_parameter_cone_grav, 0, Mpi_Comm_World, mpierr)

    Allocate(filelist(param%nfile))
    If(procID == 0) Then
       Call retrievefilelist()
    End If

    Call Mpi_Bcast(filelist, len(filelist(1))*param%nfile, Mpi_Character, 0, Mpi_Comm_World, mpierr)

    fileperproc = param%nfile / procNB
    first_file = procID * fileperproc + 1

    If( procID < mod(param%nfile,procNB) ) Then
       fileperproc = fileperproc + 1
       first_file = first_file + procID
    Else If (mod(param%nfile,procNB) > 0) Then
       first_file = first_file + mod(param%nfile,procNB) 
    End If

    last_file = first_file + fileperproc - 1


  End Subroutine readparameters

  !=======================================================================

  Subroutine retrievefilelist()

    Character(len=200) :: conefilename
    Character(len=5) :: ifilechar
    Integer :: indfile
    Integer :: ifile
    Logical :: fileexist

    indfile = 0
    Do ifile = param%first_file, param%last_file
       Write(ifilechar(1:5),'(I5.5)') ifile
       conefilename = trim(param%cone_input_file)//'_'//ifilechar//'.dat'
       Inquire(File=trim(param%input_path)//'/'//conefilename, Exist=fileexist)
       If(fileexist) Then
          indfile = indfile + 1
          filelist(indfile) = conefilename
       End If
    End Do

    If(indfile /= param%nfile) Then 
       Print *, 'Error while retrieving the list of Ramses cone files name.'
    End If

  End Subroutine retrievefilelist

  !=======================================================================

  Subroutine readramsesfiles()
    
    Use modmpicommons
    Use modreadinfo
    Use modiocommons
    Implicit None

    Integer(kind=4) :: ifile
    Character(len=401) :: filedat
    Character(len=401) :: filehdr
    Character(len=400) :: filename
    Integer(kind=4) :: ioerr
    Character(len=500) :: errormessage

    If(procID==0) Then
       filename = trim(param%input_path)//'/'//trim(param%info_cone_input_file)
       Call readinfoconegrav(filename, infocone, ioerr, errormessage)

       If(ioerr > 0) Then
          Print *,errormessage
          Call EmergencyStop(errormessage,ioerr)
       End If

       filename = trim(param%input_path)//'/'//trim(param%info_ramses_input_file)
       Call readinforamses(filename, inforamses, ioerr, errormessage)
       
       If(ioerr > 0) Then
          Print *,errormessage
          Call EmergencyStop(errormessage,ioerr)
       End If
    End If

    Call create_mpi_type_info_ramses()
    Call Mpi_Bcast(inforamses, 1, Mpi_Type_info_ramses, 0, Mpi_Comm_World, mpierr)
    Call create_mpi_type_info_cone_part()
    Call Mpi_Bcast(infocone, 1, Mpi_Type_info_cone_grav, 0, Mpi_Comm_World, mpierr)

    Allocate(ncelltab(first_file:last_file))
    Allocate(ncellperleveltab(param%nlevel,first_file:last_file))
    Allocate(ncellperlevel(param%nlevel))

    ncellperleveltab = 0
    
    Do ifile = first_file, last_file
       filedat = trim(param%input_path)//'/'//trim(filelist(ifile))
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

    Do ifile = first_file, last_file
       filedat = trim(param%input_path)//'/'//trim(filelist(ifile))
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

    Allocate(fp(param%nlevel))
    ndim = 3
    fp(1) = 0
    Do i = 2, param%nlevel
       fp(i) = fp(i-1) + ncellperlevel(i-1)
    End Do
    Do i = first_file, ifile-1
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

    level:Do ilvl = 1, param%nlevel

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
          refined(pp:pp+nstride-1) = tmpint
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
          refined(pp:pp+rest-1) = tmpint
          Deallocate(tmpsimple)
          Deallocate(tmpint)
       End If

    End Do level

    Close(10)
    
  End Subroutine readconegravdat

  !=======================================================================

  Subroutine h5writeconegrav()

    Use modiocommons
    Implicit none

    Character(len=8) :: charic
    Character(len=5) :: charncoarse
    Character(len=400) :: shellname
    Character(len=H5STRLEN) :: aname
    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: dsetname
    Character(len=H5STRLEN) :: adata
    Character(len=H5STRLEN) :: codename
    Integer(kind=hid_t) :: file_id
    Integer(kind=hid_t) :: gr_id, grlvl_id
    Integer(kind=4) :: ic, ilvl
    Integer(kind=4) :: ipointer
    Integer(kind=4) :: np
    Integer(kind=4), dimension(3) :: nctab
    Real(kind=8), dimension(3) :: dimtab
    Integer(kind=4) :: beglvl, endlvl, nlevelperproc, begc, endc
    Logical(kind=4) :: empty
    Integer(kind=8),dimension(1) :: ncelllvl
    Integer(kind=4) :: tmpint4
    Integer(kind=8) :: npart8
    Logical(kind=4) :: islast

    islast = .false.

    Write(charncoarse(1:5), '(I5.5)') infocone%nstep_coarse
    shellname = 'shell_cone_grav_data_'//trim(param%simulation_name)//'_'//charncoarse//'.h5'
    
    Call hdf5_create_mpi_file(shellname,Mpi_Comm_World,file_id)

    codename = 'conecreator_grav'
    tmpint4 = (2**inforamses%levelmin)
    npart8 = tmpint4**3
    Call h5write_meta_common(file_id, codename, npart8)
    Call h5write_meta_info_cone(file_id, infocone, islast)
    Call h5write_meta_info_ramses(file_id, inforamses, islast)
    Call h5write_meta_conecreator_parameter(file_id, param)

    ! Write type as attribute    
    aname = 'file_type'
    adata = 'conegrav'
    Call hdf5_write_attr(file_id, aname, adata)
    
    empty = .false.
    If(param%nlevel >= procNB) Then
       nlevelperproc = param%nlevel / procNB
       If(procID < modulo(param%nlevel,procNB)) Then
          nlevelperproc = nlevelperproc + 1
          beglvl = nlevelperproc*procID + 1
          endlvl = beglvl + nlevelperproc - 1
       Else
          beglvl = nlevelperproc*procID + modulo(param%nlevel,procNB) + 1
          endlvl = beglvl + nlevelperproc - 1
       End If
    Else
       nlevelperproc = 1
       If(procID >= param%nlevel) Then
          nlevelperproc = 0
          empty = .true.
       End If
    End If
    dsetname = 'ncell_cube_array'
    Call hdf5_write_mpi_data(file_id,dsetname,ncube,nlevelperproc,ncellcube(:,beglvl:endlvl),&
         Mpi_Comm_World,empty)

    nctab(1) = ncx
    nctab(2) = ncy
    nctab(3) = ncz
    aname = 'ncube_array'
    Call hdf5_write_attr(file_id, aname, 3, nctab)

    dimtab(1) = param%cone_max_radius
    dimtab(2) = hy
    dimtab(3) = hz
    aname = 'dimensions_array'
    Call hdf5_write_attr(file_id, aname, 3, dimtab)
    
        
    ! loop over the levels
    ipointer = 1
    Do ilvl = 1, param%nlevel
       Write(charic(1:2),'(I2.2)') param%levelmin+ilvl-1
       groupname = 'level'//charic(1:2)
       Call hdf5_create_group(file_id, groupname, grlvl_id)

       ncelllvl = 0

       Do ic = 1, ncube
          ncelllvl = ncelllvl + ncellcube(ic, ilvl)
          If(ncellcube(ic,ilvl) /= 0) Then
             Write(charic(1:8),'(I8.8)') ic
             groupname = 'cube'//charic
             Call hdf5_create_group(grlvl_id, groupname, gr_id)
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
#endif
             dsetname = 'position_cell'
             Call hdf5_write_mpi_data(gr_id, dsetname, 3, np, pos(:,begc:endc), Mpi_Comm_World,empty)
             dsetname = 'gravity_field_cell'
             Call hdf5_write_mpi_data(gr_id, dsetname, 3, np, acc(:,begc:endc), Mpi_Comm_World,empty)
             dsetname = 'potential_cell'
             Call hdf5_write_mpi_data(gr_id, dsetname, np, phi(begc:endc), Mpi_Comm_World,empty)
             dsetname = 'density_cell'
             Call hdf5_write_mpi_data(gr_id, dsetname, np, rho(begc:endc), Mpi_Comm_World,empty)
             dsetname = 'refined_bool'
             Call hdf5_write_mpi_data(gr_id, dsetname, np, refined(begc:endc), Mpi_Comm_World,empty)
             
             If(.not.empty) Then
                ipointer = ipointer + np
             End If
             
             Call hdf5_close_group(gr_id)
          End If
       End Do

       dsetname = 'ncell_level'
       empty = .true.
       If(procID==0) Then
          empty = .false.
       End If
       np = 1
       Call hdf5_write_mpi_data(grlvl_id, dsetname, np, ncelllvl(1:1),Mpi_Comm_World,empty)

       Call hdf5_close_group(grlvl_id)
    End Do

    Call hdf5_close_mpi_file(file_id)

  End Subroutine h5writeconegrav


End Module modio

