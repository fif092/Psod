!==============================================================================
! Project: pFoF
! File: pfof_snap/src/modio.f90
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
!! This file contains subroutines used for I/O.

!> This module contains subroutines used for I/O.
!> Authors: F. Roy, V. Bouillot
Module modio

  Use mpi
  Use modconstant
  Use modvarcommons
  Use modvariables
  Use modreadinfo
 
  Implicit None

  Private 

  Public :: title, &
       theend, &
       readparameters,&
       print_screen_parameters,&
       writetimings, &
       ramses_lecture

Contains

  !=======================================================================
  !> This subroutines prints the title message on screen.
  Subroutine title()

    Use modmpicommons, only : procNB

    Implicit None

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
    Print *,'Snapshot version'
    Print *,'Number of processes: '//trim(adjustl(charpnb))
    Print *,' '

  End Subroutine title

  !=======================================================================
  !> This subroutine print the final completion message on screen.
  Subroutine theend()

    Implicit none
    
    Print *,' '
    Print *,'Run Completed!'
    Print *,' '

  End Subroutine theend

  !=======================================================================
  !> This subroutine reads the input parameters from pfof_snap.nml file (default name) with process 0
  !! and broadcasts them to every processes.
  !! It also writes the parameters in a log file and in a backup .nml file.
  Subroutine readparameters()

    Use modvariables, only : param
    Use modmpicommons
    Use modreadparameters

    Implicit none

    Character(len=400) :: filename
    
    Integer(kind=4) :: ioerr
    Character(len=84) :: filelog 
    Integer(kind=4) :: mpierr
    Character(len=500) :: errormessage
    

    ! The process 0 read the input parameters and pack them for the broadcast.
    If(procID==0) Then
       filename='pfof_snap.nml'
       Call read_pfof_snap_parameters(filename, param, ioerr, errormessage, .true.)
       
       If(ioerr /= 0) Then
          Print *,errormessage
          Call EmergencyStop(errormessage,ioerr)
       End If
       
       Print *,'Parallel FoF snapshot version'
       Print *,procNB,' processes:'

       ! Open log file
       filelog = 'pfof_log_'//trim(param%simulation_name)//'.log'
       Open(Unit=Ulog,file=filelog)
       Write(Ulog,*) 'Parallel FoF'
       Write(Ulog,*) procNB,' processes:'
       
    End If

    Call create_mpi_type_param_pfof_snap()
    
    Call Mpi_Bcast(param, 1, Mpi_Type_parameter_pfof_snap, 0, Mpi_Comm_World, mpierr)
    
  End Subroutine readparameters


  !=======================================================================
  !> This subroutine writes the input parameters in a .nml file,
  !! print them on screen and writes them in the .log file.
  !! It should be called only from process 0.
  Subroutine print_screen_parameters()
    
    Use modvariables, only : param
    Use modmpicommons

    Implicit none
        
    If(procID==0) Then

       ! Print input parameters on screen
       Print *, 'Input parameters:'
       Print *, ' '
       Print *, 'Type of RAMSES input files:                     ',param%code_index
       Print *, 'Path to input files:                            ',trim(param%input_path)
       Print *, 'Particle files base name:                       ',trim(param%part_input_file)
       Print *, 'Info file base name:                            ',trim(param%info_input_file)
       Print *, 'Size of groups of inputs:                       ',param%grpsize
       Print *, 'Were stars written in RAMSES files:             ',param%do_skip_star
       Print *, 'Were metallicities written in RAMSES files:     ',param%do_skip_metal
       Print *, 'Were potentials written in RAMSES files:        ',param%do_read_potential
       Print *, 'Were forces written in RAMSES files:            ',param%do_read_gravitational_field
       Print *, 'Read particles from cube files:                 ',param%do_read_from_cube
       Print *, 'Gather factor for cube input:                   ',param%gatherread_factor
       Print *,' '
       Print *, 'Halo detection parameters:'
       Print *, 'Percolation parameter:                          ',param%percolation_length
       Print *, 'Minimum mass of halo to be analyzed:            ',param%mmin
       Print *, 'Maximum mass of halo to be analyzed:            ',param%mmax
       Print *, 'Perform friends of friends halo detection:      ',param%do_fof
       Print *, 'Perform unbinding:                              ',param%do_unbinding
       Print *, 'Perform subhalo detection:                      ',param%do_subhalo
       Print *,' '
       Print *, 'Output parameters:' 
       Print *, 'Simulation name:                                ',trim(param%simulation_name)
       Print *, 'Snapshot number:                                ',param%snapshot
       Print *, 'Write cubes of particles:                       ',param%do_write_cube
       Print *, 'Gather factor for cube output:                  ',param%gatherwrite_factor
       Print *, 'Sort particles in cube files:                   ',param%do_sort_cube
       Print *, 'Perform timings (imply extra synchronisations): ',param%do_timings       
       Print *, ' '
       
    End If

  End Subroutine print_screen_parameters


  !=======================================================================
  !> This subroutine writes timings in the log file and prints them on screen.
  !! It should be called only by process 0.
  !! It also closes the log file.
  Subroutine writetimings()
    
    Use modtiming

    Implicit none

    tFoF = tFoFinit + tFoFloc + tRaccord
    tOut = tOuthalopart+tOutmass
    
    Print *,''
    Print *,'Timings:'
    Print *,'Input:',tRead
    Print *,'        initialization        :',tInitRead
    Print *,'        read input files      :',tReadFile
    Print *,'        scatter particles     :',tTailPart
    Print *,''
    Print *,'Friend of Friend:',tFoF
    Print *,'        initialization:',tFoFinit
    Print *,'        local FoF     :',tFoFloc
    Print *,'        merge         :',tRaccord
    Print *,''
    Print *,'Observables computation and output:',tObs + tOut + tSort + tGatherhalo + tSelecthalo 
    Print *,'        sort particles following haloID:',tSort
    Print *,'        gather particles following haloID: ',tGatherhalo
    Print *,'        select halo with M > Mmin: ',tSelecthalo
    Print *,'        computation of observables:',tObs
    Print *,'        write files:',tOuthalopart+tOutmass
    
    Write(Ulog,*) ''
    Write(Ulog,*) 'End of pFoF'
    Write(Ulog,*) ''
    Write(Ulog,*) 'Timings:'
    Write(Ulog,*) 'Input:',tRead
    Write(Ulog,*) '        initialization        :',tInitRead
    Write(Ulog,*) '        read input files      :',tReadFile
    Write(Ulog,*) '        scatter particles     :',tTailPart
    Write(Ulog,*) ''
    Write(Ulog,*) 'Friend of Friend:',tFoF
    Write(Ulog,*) '        initialization:',tFoFinit
    Write(Ulog,*) '        local FoF     :',tFoFloc
    Write(Ulog,*) '        merge         :',tRaccord
    Write(Ulog,*) ''
    Write(Ulog,*) 'Observables computation and output:',tObs + tOut+ tSort + tGatherhalo + tSelecthalo
    Write(Ulog,*) '        sort particles following haloID:',tSort
    Write(Ulog,*) '        gather particles following haloID: ',tGatherhalo
    Write(Ulog,*) '        select halo with M > Mmin: ',tSelecthalo
    Write(Ulog,*) '        computation of observables:',tObs
    Write(Ulog,*) '        write files:',tOuthalopart+tOutmass


    ! Close log file
    Close(Ulog)
    
  End Subroutine writetimings


  !=======================================================================
  !> This subroutine reads the particles files created by RAMSES that pFOF has to analyze.
  Subroutine ramses_lecture()

    Use modvariables, only : param 
    Use modmpicom
    Use modmpicommons
    Use modtiming

    Implicit none

    !-----------------------------------------------
    ! Lecture du fichier particules au format Ramses
    !-----------------------------------------------

    ! Local variables
    Character(len=5)               :: ncharcpu
    Character(len=9)               :: tmpstr1, tmpstr2
    Character(len=400)             :: nomfich
    Character(len=11)              :: grpchar

    Integer(kind=4)                :: i, j, icpu, idim   ! loop variables
    Integer(kind=4)                :: destCoord(3)       ! coords of the destination MPI process in MPI process cart
    Integer(kind=4)                :: nrecv              ! number of elements received in a Mpi_Recv
    Integer(kind=4)                :: recvpoint          ! address of the 1st element received in the local vector
    Integer(kind=4)                :: mynbfile           ! number of RAMSES part files read by local process
    Integer(kind=4)                :: nmod               ! 
    Integer(kind=4)                :: firstp, lastp      ! id of 1st and last RAMSES part file read
    Integer(kind=4), allocatable   :: npartvloc(:), npartv(:)  ! temp and global table of particle numbers for each process
    Integer(kind=4)                :: n_i, n_j, n_k, nsd, ind 
    Integer(kind=4)                :: ncpu2       ! process number  read in RAMSES part files
    Integer(kind=4)                :: ndim2       ! dimension       read in RAMSES part files
    Integer(kind=4)                :: npartloc    ! particle number read in RAMSES part files
    Integer(kind=4)                :: prov, dest  ! provenance and destination process number for p2p MPI communications
    Integer(kind=4)                :: mpistat(Mpi_Status_Size)   ! status of MPI communication
    Integer(kind=PRI)              :: tmplongint              ! temp integer8 variable
    Integer(kind=PRI), allocatable :: tmpi(:), tmpsendi(:)    ! TYPE VARIABLE EN FONCTION DU NB DE PART
    Integer(kind=4)                :: grpnb

    Real(kind=4), allocatable     :: tmpsimple(:)            ! temporary variable for Ramses v2 output
    Real(kind=8), allocatable     :: tmpdouble(:)            ! temporary variable for Ramses v3 output
    Real(kind=4), allocatable     :: tmpsendx(:,:),tmpsendv(:,:), tmpsendp(:), tmpsendf(:,:)
    Real(kind=4), allocatable     :: tmpx(:,:), tmpv(:,:), tmpp(:), tmpf(:,:)
    Real(kind=4)                  :: deltasd
    Integer(kind=4)               :: tmpinteger
    Integer(kind=4) :: ioerr
    Character(len=500) :: errormessage

    Integer(kind=4) :: mpierr
    Integer(kind=4) :: mpireqs1,mpireqs2,mpireqs3,mpireqs4,mpireqs5
    Integer(kind=4) :: mpireqr1,mpireqr2,mpireqr3,mpireqr4,mpireqr5
    

    ! Initialisation timer
    time0 = Mpi_Wtime()

    grpchar = 'group_00001'

    ! Lecture parametres et remplissage du tampon pour diffusion des parametres
    If(procID == 0) Then
       If(param%code_index.eq.'RA2') Then 
          Print *,'Reading Ramses v2 output...'
          Write(Ulog,*) 'Reading Ramses v2 output...'
       Else if(param%code_index.eq.'RA3') Then
          Print *,'Reading Ramses v3 output...'
          Write(Ulog,*) 'Reading Ramses v3 output...'
       End If
       
       If( param%grpsize == 0 ) Then
          nomfich = trim(param%input_path)//'/'//trim(param%info_input_file)
       Else
          nomfich = trim(param%input_path)//'/'//trim(grpchar)//'/'//trim(param%info_input_file)
       End If
       Print *,'Reading RAMSES info file:',trim(nomfich)
       Call readinforamses(nomfich, inforamses, ioerr, errormessage)

       If(ioerr > 0) Then
          Print *,errormessage
          Call EmergencyStop(errormessage,ioerr)
       End If
    End If
    Call create_mpi_type_info_ramses()
    
    Call Mpi_Bcast(inforamses, 1, Mpi_Type_info_ramses, 0, Mpi_Comm_World, mpierr)

    nres = 2**inforamses%levelmin
    If(procID==0) Then

       Write(*,*) 'Number of:' 
       Write(*,'(A25,I6)') ' - files for each output:',inforamses%ncpu
       Write(*,'(A25,I6)') ' - dimensions:           ',inforamses%ndim
       Write(*,'(A25,I6)') ' - grid points:          ',nres
       Write(Ulog,*) 'nb_proc = ',inforamses%ncpu,'ndim = ',inforamses%ndim,'nres = ',nres

    End If

    ngrid = int(nres,kind=8)**3
    
    If(procID==0) Print *,'Reading positions...'
    If(procID==0) Write(Ulog,*) 'Reading positions...'

    global_npart = 0

    nmod = mod(inforamses%ncpu,procNB)
    mynbfile = inforamses%ncpu / procNB
    If(procID <= nmod-1) Then
       mynbfile = mynbfile+1
       firstp   = procID * mynbfile + 1
       lastp    = (procID+1) * mynbfile
    Else
       firstp   = procID * mynbfile + 1 + nmod
       lastp    = (procID+1) * mynbfile + nmod
    End If
    
 
    local_npart = 0

    Allocate(npartv(procNB))
    Allocate(npartvloc(procNB))

    npartv = 0
    npartvloc = 0

    nsd = int(procNB**(1./3.))
    deltasd = 1./nsd

    If(procID == 0) Then
       Write(*,*) 'Number of subdomains in each dimension:',nsd
       Write(*,*) 'Size of each subdomain:',deltasd
    End If

    xmin =  info_proc%global_comm%coords(1)      * deltasd
    xmax = (info_proc%global_comm%coords(1) + 1) * deltasd
    ymin =  info_proc%global_comm%coords(2)      * deltasd
    ymax = (info_proc%global_comm%coords(2) + 1) * deltasd
    zmin =  info_proc%global_comm%coords(3)      * deltasd
    zmax = (info_proc%global_comm%coords(3) + 1) * deltasd
    If(info_proc%global_comm%coords(1) == info_proc%global_comm%dims(1) - 1) xmax = 1.e0
    If(info_proc%global_comm%coords(2) == info_proc%global_comm%dims(2) - 1) ymax = 1.e0
    If(info_proc%global_comm%coords(3) == info_proc%global_comm%dims(3) - 1) zmax = 1.e0


    Do icpu = firstp,lastp
       If( param%grpsize == 0 ) Then
          Write(ncharcpu(1:5),'(I5.5)') icpu
          nomfich = trim(param%input_path)//'/'//trim(param%part_input_file)//trim(ncharcpu)
       Else
          Write(ncharcpu(1:5),'(I5.5)') icpu
          grpnb = (icpu-1)/param%grpsize + 1
          Write(grpchar(7:11),'(I5.5)') grpnb
          nomfich = trim(param%input_path)//'/'//trim(grpchar)//'/'//trim(param%part_input_file)//trim(ncharcpu)
       End If

       Open(unit=1,file=nomfich,status='old',form='unformatted')
       Read(1) ncpu2
       Read(1) ndim2
       Read(1) npartloc
       Close(1)

       If((ncpu2/=inforamses%ncpu).Or.(ndim2/=inforamses%ndim)) Then
          Call EmergencyStop('Files'//trim(nomfich)// ' and '//&
               trim(param%info_input_file)//' are not consistent for ncpu and/or ndim',22)
       End If

       local_npart = local_npart + npartloc
    End Do

    tmplongint = local_npart
    Call Mpi_AllReduce(tmplongint,global_npart,1,MPI_PRI,Mpi_Sum,Mpi_Comm_World,mpierr)

    If(procID == 0) Then
       Write(* ,*)'There are ',global_npart,' DM particles'
       Write(Ulog,*)'There are ',global_npart,' DM particles'
    End If

    Allocate(tmpx(3,local_npart))
    Allocate(tmpv(3,local_npart))
    Allocate(tmpi(local_npart))
    If(param%do_read_gravitational_field) Allocate(tmpf(3,local_npart))
    If(param%do_read_potential) Allocate(tmpp(local_npart))

    tmpx=0.
    local_npart = 0

    If(param%do_timings) Then
       Call Mpi_Barrier(info_proc%global_comm%name,mpierr)
       timeInt = Mpi_Wtime()
       tInitRead = timeInt - time0
    End If
    
    Do icpu = firstp,lastp
       If( param%grpsize == 0 ) Then
          Write(ncharcpu(1:5),'(I5.5)') icpu
          nomfich = trim(param%input_path)//'/'//trim(param%part_input_file)//trim(ncharcpu)
       Else
          Write(ncharcpu(1:5),'(I5.5)') icpu
          grpnb = (icpu-1)/param%grpsize + 1
          Write(grpchar(7:11),'(I5.5)') grpnb
          nomfich = trim(param%input_path)//'/'//trim(grpchar)//'/'//trim(param%part_input_file)//trim(ncharcpu)
       End If


       Open(unit=1,file=nomfich,status='old',form='unformatted')
       Read(1) ncpu2
       Read(1) ndim2       
       Read(1) npartloc
       ramsesV2 : If(param%code_index.eq.'RA2') Then

          Allocate(tmpsimple(1:npartloc))
          ! Read positions
          Do idim = 1,inforamses%ndim
             Read(1) tmpsimple
             ! put all positions in tmpx vector
             tmpx(idim,local_npart+1:local_npart+npartloc) = tmpsimple
          End Do
          
          ! Read velocities in a dummy variable
          Do idim = 1,inforamses%ndim
             Read(1) tmpsimple
             ! put all velocities in tmpv vector
             tmpv(idim,local_npart+1:local_npart+npartloc) = tmpsimple
          End Do
          
          ! Read masses in a dummy variable
          Read(1) tmpsimple
          Deallocate(tmpsimple)
          
       End If ramsesV2

       ramsesV3 : If(param%code_index == 'RA3') Then
          Allocate(tmpdouble(1:npartloc))
          
          Read(1)
          Read(1)
          Read(1)
          Read(1)
          Read(1)
          ! Read positions
          Do idim = 1,inforamses%ndim
             Read(1) tmpdouble
             ! put all positions in tmpx vector
             tmpx(idim,local_npart+1:local_npart+npartloc) = real(tmpdouble, kind=4)
          End Do
          
          ! Read velocities in a dummy variable
          Do idim = 1,inforamses%ndim
             Read(1) tmpdouble
             ! put all velocities in tmpv vector
             tmpv(idim,local_npart+1:local_npart+npartloc) = real(tmpdouble,kind=4)
          End Do
          
          ! Read masses in a dummy variable
          Read(1) tmpdouble

          If(param%do_skip_star) Then
             Read(1) tmpdouble
             If(param%do_skip_metal) Read(1) tmpdouble
          End If

       End If ramsesV3

       Deallocate(tmpdouble)

       ! Read particle id
       Read(1) tmpi(local_npart+1:local_npart+npartloc)

       ! Read potential if potential parameter is .true.
       If(param%do_read_potential) Then
          Allocate(tmpsimple(1:npartloc))
          ! First: skip the level (integer 4)
          Read(1) tmpinteger
          ! Then read the potential: same kind as position and velocity
          Read(1) tmpsimple
          tmpp(local_npart+1:local_npart+npartloc) = tmpsimple
       End If

       ! Read force if force parameter is .true.
       If(param%do_read_gravitational_field) Then
          If( .not. param%do_read_potential) Then
             Allocate(tmpsimple(1:npartloc))
             ! If the potential skip the level (integer 4)
             Read(1) tmpinteger
          End If
          ! Then read the force: same kind as position and velocity
          Do idim = 1,inforamses%ndim
             Read(1) tmpsimple
             tmpf(idim,local_npart+1:local_npart+npartloc) = tmpsimple
          End Do
       End If
       If(Allocated(tmpsimple)) Deallocate(tmpsimple)
       
       Do j = local_npart+1,local_npart+npartloc
          Do idim = 1,inforamses%ndim
             If(abs(tmpx(idim,j)-1.0e0) < 1.e-8) tmpx(idim,j) = 0.e0
          End Do
          n_i = int(tmpx(1,j)/deltasd)
          n_j = int(tmpx(2,j)/deltasd)
          n_k = int(tmpx(3,j)/deltasd)
          ind = nsd**2 *n_i + nsd*n_j + n_k + 1
          npartvloc(ind) = npartvloc(ind)+1
       End Do
       
       Close(1)
       local_npart = local_npart+npartloc
    End Do

    Call Mpi_AllReduce(npartvloc,npartv,procNB,Mpi_Integer,Mpi_Sum,Mpi_Comm_World,mpierr)
    tReadfile = Mpi_Wtime() - timeInt
    timeInt = Mpi_Wtime()

    ! ------------------------------------------------
    ! Repartition des particules entre les processeurs
    ! ------------------------------------------------

    Allocate (position(3,npartv(procID+1)))
    Allocate (velocity(3,npartv(procID+1)))
    Allocate (pfof_id(npartv(procID+1)))
    If(param%do_read_gravitational_field) Allocate(field(3, npartv(procID+1)))
    If(param%do_read_potential) Allocate(potential(npartv(procID+1)))

    recvpoint = 1

    processus : Do i = 1,procNB - 1
       dest = mod(procID + i,procNB)
       prov = mod(procID + procNB - i, procNB)

       Call Mpi_Cart_coords(info_proc%global_comm%name,dest,3,destCoord,mpierr)
       Call Mpi_Isend(npartvloc(dest+1),1,Mpi_Integer,dest,procID,&
            info_proc%global_comm%name,mpireqs1,mpierr)
       Call Mpi_Irecv(nrecv,1,Mpi_Integer,prov,prov,&
            info_proc%global_comm%name,mpireqr1,mpierr)
       xmin =  destCoord(1)      * deltasd
       xmax = (destCoord(1) + 1) * deltasd
       ymin =  destCoord(2)      * deltasd
       ymax = (destCoord(2) + 1) * deltasd
       zmin =  destCoord(3)      * deltasd
       zmax = (destCoord(3) + 1) * deltasd
       If(destCoord(1) == info_proc%global_comm%dims(1) - 1) xmax = 1.e0
       If(destCoord(2) == info_proc%global_comm%dims(2) - 1) ymax = 1.e0
       If(destCoord(3) == info_proc%global_comm%dims(3) - 1) zmax = 1.e0
       
       If(npartvloc(dest+1)/=0) Then
          Allocate(tmpsendx(3,npartvloc(dest+1)))
          Allocate(tmpsendv(3,npartvloc(dest+1)))
          Allocate(tmpsendi(  npartvloc(dest+1)))
          If(param%do_read_gravitational_field) Allocate(tmpsendf(3, npartvloc(dest+1)))
          If(param%do_read_potential) Allocate(tmpsendp(npartvloc(dest+1)))
          ind = 1

          Do j=1,local_npart
             If(tmpx(1,j)>=xmin .and. tmpx(1,j)<xmax .and. &
                  tmpx(2,j)>= ymin .and. tmpx(2,j) < ymax .and. &
                  tmpx(3,j)>= zmin .and. tmpx(3,j) < zmax) Then
                tmpsendx(:,ind) = tmpx(:,j)
                tmpsendv(:,ind) = tmpv(:,j)
                If(param%do_read_gravitational_field) tmpsendf(:,ind) = tmpf(:,j)
                If(param%do_read_potential) tmpsendp(  ind) = tmpp(j)
                tmpsendi(  ind) = tmpi(j)
                ind=ind+1
             End If
          End Do

          If(ind/=npartvloc(dest+1)+1) Then
             Call EmergencyStop('Erreur dans la repartition des particules',23)
          End If

       End If

       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

       If(npartvloc(dest+1)/=0) Then
          Call Mpi_Isend(tmpsendx,3*npartvloc(dest+1),Mpi_Real,dest,1,&
               info_proc%global_comm%name,mpireqs1,mpierr)
          Call Mpi_Isend(tmpsendv,3*npartvloc(dest+1),Mpi_Real,dest,2,&
               info_proc%global_comm%name,mpireqs2,mpierr)
          Call Mpi_Isend(tmpsendi,  npartvloc(dest+1), MPI_PRI,dest,3,&
               info_proc%global_comm%name,mpireqs3,mpierr)
          If(param%do_read_potential) Call Mpi_Isend(tmpsendp, npartvloc(dest+1), &
               Mpi_Real,dest,4,info_proc%global_comm%name,mpireqs4,mpierr)
          If(param%do_read_gravitational_field) Call Mpi_Isend(tmpsendf,3*npartvloc(dest+1),&
               Mpi_Real,dest,5, info_proc%global_comm%name,mpireqs5,mpierr)
       End If
       If(nrecv/=0) Then
          Call Mpi_Irecv(position(1,recvpoint),3*nrecv,Mpi_Real,prov,1,&
               info_proc%global_comm%name,mpireqr1,mpierr)
          Call Mpi_Irecv(velocity(1,recvpoint),3*nrecv,Mpi_Real,prov,2,&
               info_proc%global_comm%name,mpireqr2,mpierr)
          Call Mpi_Irecv(pfof_id(recvpoint),  nrecv, MPI_PRI,prov,  3,&
               info_proc%global_comm%name,mpireqr3,mpierr)
          If(param%do_read_potential) Call Mpi_Irecv(potential(recvpoint),nrecv,&
               Mpi_Real,prov,4,info_proc%global_comm%name,mpireqr4,mpierr)
          If(param%do_read_gravitational_field) Call Mpi_Irecv(field(1,recvpoint),3*nrecv,&
               Mpi_Real,prov,5,info_proc%global_comm%name,mpireqr5,mpierr)
       End If
       recvpoint=recvpoint+nrecv

       If(npartvloc(dest+1)/=0) Then
          Call Mpi_Wait(mpireqs1,mpistat,mpierr)
          Deallocate(tmpsendx)
          Call Mpi_Wait(mpireqs2,mpistat,mpierr)
          Deallocate(tmpsendv)
          Call Mpi_Wait(mpireqs3,mpistat,mpierr)
          Deallocate(tmpsendi)
       End If
       If(nrecv/=0) Then
          Call Mpi_Wait(mpireqr1,mpistat,mpierr)
          Call Mpi_Wait(mpireqr2,mpistat,mpierr)
          Call Mpi_Wait(mpireqr3,mpistat,mpierr)
       End If
       
       If(param%do_read_potential) Then
          If(npartvloc(dest+1)/=0) Then
             Call Mpi_Wait(mpireqs4, mpistat, mpierr)
             Deallocate(tmpsendp)
          End If
          If(nrecv/=0) Then
             Call Mpi_Wait(mpireqr4, mpistat, mpierr)
          End If
       End If

       If(param%do_read_gravitational_field) Then
          If(npartvloc(dest+1)/=0) Then
             Call Mpi_Wait(mpireqs5, mpistat, mpierr)
             Deallocate(tmpsendf)
          End If
          If(nrecv/=0) Then
             Call Mpi_Wait(mpireqr5, mpistat, mpierr)
          End If
       End If

    End Do processus

    xmin =  info_proc%global_comm%coords(1)      * deltasd
    xmax = (info_proc%global_comm%coords(1) + 1) * deltasd
    ymin =  info_proc%global_comm%coords(2)      * deltasd
    ymax = (info_proc%global_comm%coords(2) + 1) * deltasd
    zmin =  info_proc%global_comm%coords(3)      * deltasd
    zmax = (info_proc%global_comm%coords(3) + 1) * deltasd
    If(info_proc%global_comm%coords(1) == info_proc%global_comm%dims(1) - 1) xmax = 1.e0
    If(info_proc%global_comm%coords(2) == info_proc%global_comm%dims(2) - 1) ymax = 1.e0
    If(info_proc%global_comm%coords(3) == info_proc%global_comm%dims(3) - 1) zmax = 1.e0

    ind = 0
    Do j=1,local_npart
       If(tmpx(1,j)>=xmin .and. tmpx(1,j)<xmax .and. &
            tmpx(2,j)>= ymin .and. tmpx(2,j) < ymax .and. &
            tmpx(3,j)>= zmin .and. tmpx(3,j) < zmax) Then
          position(:,recvpoint+ind) = tmpx(:,j)
          velocity(:,recvpoint+ind) = tmpv(:,j)
          pfof_id(recvpoint+ind)  = tmpi(j)
          If(param%do_read_potential) potential(recvpoint+ind) = tmpp(j)
          If(param%do_read_gravitational_field) field(:,recvpoint+ind) = tmpf(:,j)
          ind = ind+1
       End If
    End Do

    If(recvpoint+ind /= npartv(procID+1)+1) Then
       Write(tmpstr1,'(I9.9)') recvpoint+ind
       Write(tmpstr2,'(I9.9)') npartv(procID+1)+1
       Call EmergencyStop('Wrong particles number found after send/recv swaps:'//tmpstr1//' ; '//tmpstr2,24)
    End If

    local_npart = npartv(procID+1)

    Deallocate(tmpx,tmpv,tmpi)
    Deallocate(npartv, npartvloc)
    If(param%do_read_potential) Deallocate(tmpp)
    If(param%do_read_gravitational_field) Deallocate(tmpf)
    tTailPart = Mpi_Wtime() - timeInt
    tRead = Mpi_Wtime() - time0

  End Subroutine ramses_lecture


  !=======================================================================


End Module modio
