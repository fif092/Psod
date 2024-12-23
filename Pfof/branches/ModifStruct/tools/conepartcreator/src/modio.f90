!==============================================================================
! Project: pFoF
! File: tools/conepartcreator/src/modio.f90
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

  Integer(kind=4), dimension(:), allocatable :: nparttab
  Integer(kind=4) :: ioerr
  Integer(kind=4) :: npartloc

  Character(len=200), dimension(:), allocatable :: filelist    ! list of the filenames that we must read

Contains  

  !=======================================================================

  Subroutine theend()

    Implicit none
    
    Print *,' '
    Print *,'Run Completed!'
    Print *,' '

  End Subroutine theend

  !=======================================================================
  
  Subroutine readparameters
    
    Use modconstant
    Use modmpicommons
   
    Implicit none

    Character(len=16) :: infoconefile
    Integer :: fileperproc

    Character(len=200) :: input_path
    Character(len=200) :: cone_input_file
    Character(len=200) :: info_cone_input_file
    Character(len=200) :: info_ramses_input_file
    Integer(kind=4) :: nfile
    Integer(kind=4) :: first_file
    Integer(kind=4) :: last_file
    Logical(kind=4) :: do_read_ramses_part_id
    Logical(kind=4) :: do_read_potential
    Logical(kind=4) :: do_read_gravitational_field
    Real(kind=8) :: cone_max_radius
    Character(len=200) :: simulation_name
    Real(kind=8) :: cube_size
    Integer(kind=4) :: mpierr
 
    Namelist/input_parameters/ input_path, cone_input_file, info_cone_input_file, &
         info_ramses_input_file, nfile, first_file, last_file, &
         do_read_ramses_part_id, do_read_potential, do_read_gravitational_field, cone_max_radius
    Namelist/output_parameters/ simulation_name, cube_size
    

    If(procID==0) Then
       ! read parameters
       infoconefile = 'infoconepart.nml'
       Open(Unit=10, File=infoconefile, status='old')
       Read(10, nml=input_parameters)
       Read(10, nml=output_parameters)
       Close(10)
       
       param%input_path=input_path
       param%cone_input_file = cone_input_file
       param%info_cone_input_file = info_cone_input_file
       param%info_ramses_input_file = info_ramses_input_file
       param%nfile = nfile
       param%first_file = first_file
       param%last_file = last_file
       param%do_read_ramses_part_id = do_read_ramses_part_id
       param%do_read_potential = do_read_potential
       param%do_read_gravitational_field = do_read_gravitational_field
       param%cone_max_radius = cone_max_radius
       param%simulation_name = simulation_name
       param%cube_size = cube_size
       
    End If

    Call create_mpi_type_param_cone_part()

    Call Mpi_Bcast(param, 1, Mpi_type_parameter_cone_part, 0, Mpi_Comm_World, mpierr)


    Allocate(filelist(param%nfile))
    If(procID == 0) Then
       Call retrievefilelist()
    End If

    Call Mpi_Bcast(filelist, len(filelist(1))*param%nfile, Mpi_Character, 0, Mpi_Comm_World, mpierr)

    fileperproc = param%nfile / procNB
    param%first_file = procID * fileperproc + 1

    If( procID < mod(param%nfile,procNB) ) Then
       fileperproc = fileperproc + 1
       param%first_file = param%first_file + procID
    Else If (mod(param%nfile,procNB) > 0) Then
       param%first_file = param%first_file + mod(param%nfile,procNB) 
    End If

    param%last_file = param%first_file + fileperproc - 1


  End Subroutine readparameters

  !=======================================================================

  Subroutine retrievefilelist()

    Implicit none

    Character(len=200) :: conefilename
    Character(len=5) :: ifilechar
    Integer :: indfile
    Integer :: ifile
    Logical :: fileexist

    indfile = 0
    Do ifile = param%first_file, param%last_file
       Write(ifilechar(1:5),'(I5.5)') ifile
       conefilename = trim(param%cone_input_file)//'_'//ifilechar//'.dat'
       Inquire(File=trim(param%input_path)//'/'//trim(conefilename), Exist=fileexist)
       If(fileexist) Then
          indfile = indfile + 1
          filelist(indfile) = conefilename
       End If
    End Do

    If(indfile /= param%nfile) Then 
       Print *, 'Error while retrieving the list of Ramses cone files name.'
       Print *,trim(param%input_path)//'/'//trim(conefilename)
    End If

  End Subroutine retrievefilelist

  !=======================================================================

  Subroutine readramsesfiles()

    Use modmpicommons
    Use modreadinfo
    Implicit None

    Integer(kind=4) :: ifile
    Character(len=401) :: filedat
    Character(len=401) :: filehdr
    Character(len=400) :: filename
    Integer(kind=4) :: ioerr
    Character(len=500) :: errormessage
    Integer(kind=4) :: mpierr

!    Call create_mpi_type_info()
    Allocate(nparttab(param%first_file:param%last_file))

    If(procID==0) Then
       filename = trim(param%input_path)//'/'//trim(param%info_cone_input_file)
       Call readinfoconepart(filename, infocone, ioerr, errormessage)

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
    Call Mpi_Bcast(infocone, 1, Mpi_Type_info_cone_part, 0, Mpi_Comm_World, mpierr)
    
    Do ifile = param%first_file, param%last_file
       filedat = trim(param%input_path)//'/'//trim(filelist(ifile))
       filehdr = filedat
       Write(filehdr(len_trim(filedat)-2:len_trim(filedat)),'(A)') 'hdr'
       Call readconehdr(filehdr, npartloc)
       nparttab(ifile) = npartloc
    End Do

    npartloc = sum(nparttab) 
    If(npartloc /= 0) Then
       Allocate(pos(3,npartloc), vel(3,npartloc), id(npartloc))
       If(param%do_read_ramses_part_id) Then
          Allocate(ramsespartid(npartloc))
       End If
       If(param%do_read_potential) Then
          Allocate(pot(npartloc))
       End If
       If(param%do_read_gravitational_field) Then
          Allocate(field(3,npartloc))
       End If
    Else
       ! If npartloc == 0 we allocate arrays with size (3,1) to avoid difficulties later with 0-size arrays
       Allocate(pos(3,1), vel(3,1), id(1))
       If(param%do_read_ramses_part_id) Then
          Allocate(ramsespartid(1))
          ramsespartid = -1
       End If
       If(param%do_read_potential) Then
          Allocate(pot(1))
          pot = -1.
       End If
       If(param%do_read_gravitational_field) Then
          Allocate(field(3,1))
          field = -1.
       End If
       pos = -1.
       vel = -1.
       id = -1
    End If

    Do ifile = param%first_file, param%last_file
       filedat = trim(param%input_path)//'/'//trim(filelist(ifile))
       Call readconedat(filedat, ifile)
    End Do

  End Subroutine readramsesfiles

  !=======================================================================

  Subroutine readconehdr(file, npartloc)

    Implicit none

    Character(len=401), intent(in) :: file
    Integer(kind=4), intent(out) :: npartloc
    Integer(kind=4) :: nstride
    
    Integer(kind=4) :: nloc

    Open(Unit=10, File=trim(file), Form='unformatted', Status='old', iostat=ioerr)
    Read(10,iostat=ioerr) nloc
    If(ioerr/=0) Print *,'Error while reading nloc in '//trim(file)
    Read(10,iostat=ioerr) nstride
    If(ioerr/=0) Print *,'Error while reading nstride in '//trim(file)
    Read(10,iostat=ioerr) npartloc
    If(ioerr/=0) Print *,'Error while reading npartloc in '//trim(file)
    Close(10)

    If(nstride /= infocone%nstride) Then
       Print *,'Error: nstride has different values in Ramses HDR files and Ramses info cone files!'
    End If

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
    Integer(kind=PRI), dimension(:), allocatable :: tmpinteger

    ndim = 3
    npartloc = nparttab(ifile)
    fp = 0
    Do i = param%first_file, ifile-1
       fp = fp + nparttab(i)
    End Do
    nblocs = npartloc/infocone%nstride
    rest = npartloc - nblocs*infocone%nstride

    Open(Unit=10, File=trim(file), Status='old', Form='unformatted', iostat=ioerr)
    Read(10)
    Read(10)
    Read(10)

    Allocate(tmpsimple(infocone%nstride))
    If(param%do_read_ramses_part_id) Allocate(tmpinteger(infocone%nstride))

    Do ibloc = 0,nblocs-1
       pp = fp + ibloc*infocone%nstride + 1
       Do idim = 1, ndim
          Read(10) tmpsimple
          pos(idim,pp:pp+infocone%nstride-1) = tmpsimple
          Read(10) tmpsimple
          vel(idim,pp:pp+infocone%nstride-1) = tmpsimple
       End Do
       Read(10)
       If(param%do_read_ramses_part_id) Then
          Read(10) tmpinteger
          ramsespartid(pp:pp+infocone%nstride-1) = tmpinteger
       End If

       ! Read potential if potential parameter is .true.
       If(param%do_read_potential) Then
          Read(10) tmpsimple
          pot(pp:pp+infocone%nstride-1) = tmpsimple
       End If

       ! Read force if force parameter is .true.
       If(param%do_read_gravitational_field) Then
          Do idim = 1,inforamses%ndim
             Read(10) tmpsimple
             field(idim,pp:pp+infocone%nstride-1) = tmpsimple
          End Do
       End If

    End Do

    Deallocate(tmpsimple)
    If(Allocated(tmpinteger)) Deallocate(tmpinteger)
 
    If(rest > 0)  Then
       pp = fp + nblocs*infocone%nstride + 1
       Allocate(tmpsimple(rest))
       If(param%do_read_ramses_part_id) Allocate(tmpinteger(rest))
       Do idim = 1,ndim
          Read(10) tmpsimple
          pos(idim,pp:pp+rest-1) = tmpsimple
          Read(10) tmpsimple
          vel(idim,pp:pp+rest-1) = tmpsimple
       End Do
       Read(10) 
       If(param%do_read_ramses_part_id) Then
          Read(10) tmpinteger
          ramsespartid(pp:pp+rest-1) = tmpinteger
       End If       
       ! Read potential if potential parameter is .true.
       If(param%do_read_potential) Then
          Read(10) tmpsimple
          pot(pp:pp+rest-1) = tmpsimple
       End If

       ! Read force if force parameter is .true.
       If(param%do_read_gravitational_field) Then
          Do idim = 1,inforamses%ndim
             Read(10) tmpsimple
             field(idim,pp:pp+rest-1) = tmpsimple
          End Do
       End If

       Deallocate(tmpsimple)
       If(Allocated(tmpinteger)) Deallocate(tmpinteger)

    End If

    Close(10)
    
  End Subroutine readconedat

  !=======================================================================

  Subroutine h5writecone()

    Use modiocommons
    Use modmpicommons, only : procNB, procID
    Implicit none

    Character(len=8) :: charic
    Character(len=H5STRLEN) :: name
    Character(len=20) :: adata
    Integer(kind=hid_t) :: file_id
    Integer(kind=hid_t) :: gr_id
    Integer(kind=hid_t) :: gr_data_id
    Integer(kind=4) :: ic
    Integer(kind=4) :: ipointer
    Integer(kind=4) :: np
!    Integer(kind=4) :: isfullsky
    Integer(kind=4), dimension(3) :: nctab
    Real(kind=8), dimension(3) :: dimtab
    Integer(kind=4) :: beg, end, ncubeperproc
    Logical(kind=4) :: empty
    Integer(kind=4) :: tmpint4
    Character(len=400) :: shellname
    Character(len=5) :: charncoarse
    Integer(kind=8) :: npart8
    Character(len=H5STRLEN) :: codename
    Character(len=H5STRLEN) :: groupname
    Logical(kind=4) :: islast    

    islast = .false.

    Write(charncoarse(1:5), '(I5.5)') infocone%nstep_coarse
    shellname = 'pfof_shell_cone_part_data_'//trim(param%simulation_name)//'_'//charncoarse//'.h5'
    
    Call hdf5_create_mpi_file(shellname,Mpi_Comm_World,file_id)

    codename = 'conecreator_part'
    tmpint4 = (2**inforamses%levelmin)
    npart8 = tmpint4**3
    Call h5write_meta_common(file_id, codename, npart8, procID)
    Call h5write_meta_info_cone(file_id, infocone, islast)
    Call h5write_meta_info_ramses(file_id, inforamses, islast)
    Call h5write_meta_conecreator_parameter(file_id, param)


    groupname='metadata'
    Call hdf5_open_group(file_id, groupname, gr_id)

    ! Write type as attribute    
    name = 'file_type'
    adata = 'cone_part'
    Call hdf5_write_attr(gr_id, name, adata)
    
    ncubeperproc = ncube / procNB
    If(procID==0) Then
       ncubeperproc = ncubeperproc + modulo(ncube,procNB)
       beg = 1
       end = ncubeperproc
    Else
       beg = ncubeperproc*procID + modulo(ncube,procNB) + 1
       end = beg + ncubeperproc - 1
    End If

    name = 'npart_cube_array'
    Call hdf5_write_mpi_data(gr_id, name, ncubeperproc, npartcube(beg:end), Mpi_Comm_World)

    nctab(1) = ncx
    nctab(2) = ncy
    nctab(3) = ncz
    name = 'ncube_array'
    Call hdf5_write_attr(gr_id, name, 3, nctab)

    dimtab(1) = param%cone_max_radius
    dimtab(2) = hy
    dimtab(3) = hz
    name = 'dimensions_array'
    Call hdf5_write_attr(gr_id, name, 3, dimtab)

    name = 'npart_file'
    Call hdf5_write_data(gr_id, name, infocone%npart)
    Call hdf5_close_group(gr_id)

    groupname='data'
    Call hdf5_create_group(file_id, groupname, gr_data_id)

    ipointer = 1
    Do ic = 1, ncube
       If(npartcube(ic) /= 0) Then
          Write(charic(1:8),'(I8.8)') ic
          name = 'cube'//charic
          Call hdf5_create_group(gr_data_id, name, gr_id)
          np = npartcubeloc(ic)
          If(np == 0) Then
             np = 1
             beg = lbound(pos,2)
             end = beg
             empty = .true.
          Else
             beg = ipointer
             end = ipointer+np-1
             empty = .false.
          End If

#ifdef DEBUG
          Print *,'process:',procID,' call hdf5_write_mpi_data with arguments:',gr_id, name, 3, np, beg, end, empty
          If(empty) Then
             Print *,'pos and vel:',pos(:,beg:end),vel(:,beg:end)
          End If
          Print *,'process:',procID, ' size of ramsespartid:',size(ramsespartid)
#endif
          name = 'position_part'          
          Call hdf5_write_mpi_data(gr_id, name, 3, np, pos(:,beg:end), Mpi_Comm_World,empty)
          name = 'velocity_part'
          Call hdf5_write_mpi_data(gr_id, name, 3, np, vel(:,beg:end), Mpi_Comm_World,empty)

          If(param%do_read_ramses_part_id) Then
             name = 'identity_part_ramses'
             Call hdf5_write_mpi_data(gr_id, name, np, ramsespartid(beg:end), Mpi_Comm_World,empty)
          End If

          If(param%do_read_potential) Then
             name = 'potential_part'
             Call hdf5_write_mpi_data(gr_id, name, np, pot(beg:end), Mpi_Comm_World,empty)
          End If

          If(param%do_read_gravitational_field) Then
             name = 'gravitational_field_part'
             Call hdf5_write_mpi_data(gr_id, name, 3, np, field(:,beg:end), Mpi_Comm_World,empty)
          End If


          If(.not.empty) Then
             ipointer = ipointer + np
          End If

          Call hdf5_close_group(gr_id)
       End If
    End Do

    Call hdf5_close_group(gr_data_id)
    Call hdf5_close_mpi_file(file_id)

  End Subroutine h5writecone


End Module modio

