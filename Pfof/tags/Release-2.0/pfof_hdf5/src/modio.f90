!> @file
!! This file contains subroutines used for I/O.

!> This module contains subroutines used for I/O.
!> Authors: F. Roy, V. Bouillot
Module modio

  Use modconstant
  Use modvariables
  Implicit None

Contains

  !=======================================================================

  Subroutine title()

    Use modmpicom, only : procNB

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
    Print *,'Cone version run'
    Print *,'Number of processes: '//trim(adjustl(charpnb))
    Print *,' '

  End Subroutine title


  !=======================================================================
  !> This subroutine reads the input parameters from pfof.nml file (default name) with process 0
  !! and broadcasts them to every processes.
  Subroutine readparameters()

    Use modparameters
    Use modmpicom

    Implicit none

    Integer :: ioerr

    ! The process 0 read the input parameters and pack them for the broadcast.
    If(procID==0) Then
       
       ! Read input parameters in file 'fof.nml'
       Open(10, file='pfof.nml', iostat=ioerr) 
       
       If(ioerr>0) Then
          Print *,'** Error opening input file pfof.nml. Please check this file. **'
          Call Mpi_Abort(Mpi_Comm_World,ioerr,mpierr)
       End If
       
       Read(10, nml=ramsesinput)
       Read(10, nml=fofparam)
       Read(10, nml=outputparam)
       Close(10)
       
       Print *,'Parallel FoF'
       Print *,procNB,' processes:'
                     
    End If
    
    Call bcastparam()
    
  End Subroutine readparameters


  !=======================================================================
  !> This subroutine writes the input parameters in a .nml file,
  !! print them on screen and writes them in the .log file.
  !! It should be called only from process 0.
  Subroutine writeparameters()
    
    Use modparameters
    Use modmpicom
    Implicit none
    
    Character(len=84) :: fileopa 
    Character(len=84) :: filelog 
    
    If(procID == 0) Then
       ! Open file .nml and write input parameters
       fileopa = trim(output_root)//'.nml'
       Open(Unit=Uopa,file=fileopa)
       Write(Uopa, nml=ramsesinput)
       Write(Uopa, nml=fofparam)
       Write(Uopa, nml=outputparam)
       Close(Uopa)    
       
       ! Print input parameters on screen
       Print *, 'Input parameters:'
       Print *, ' '
       Print *, 'Simulation name:                                ',trim(simulationname)
       Print *, 'Output number:                                  ',outputNB
       Print *, 'Type of RAMSES input files:                     ',code_index
       Print *, 'Path to input files:                            ',trim(pathinput)
       Print *, 'Particle files base name:                       ',trim(namepart)
       Print *, 'Info file base name:                            ',trim(nameinfo)
       Print *, 'Size of groups of inputs:                       ',grpsize
       Print *, 'Were stars written in RAMSES files:             ',star
       Print *, 'Were metallicities written in RAMSES files:     ',metal
       Print *, 'Were potentials written in RAMSES files:        ',potential
       Print *, 'Were forces written in RAMSES files:            ',force
       Print *, 'Read particles from cube files:                 ',readfromcube
       Print *, 'Read particles from sorted cube files:          ',readfromsortedcube
       Print *, 'Gather factor for cube input:                   ',gatherread
       Print *,' '
       Print *, 'Halo detection parameters:'
       Print *, 'Percolation parameter:                          ',perco
       Print *, 'Minimum mass of halo to be analyzed:            ',Mmin
       Print *, 'Maximum mass of halo to be analyzed:            ',Mmax
       Print *, 'Perform friends of friends halo detection:      ',dofof
       Print *, 'Perform unbinding:                              ',doUnbinding
       Print *, 'Perform subhalo detection:                      ',dosubhalo
       Print *,' '
       Print *, 'Output parameters:' 
       Print *, 'Output filename base:                           ',trim(output_root)
       Print *, 'Write cubes of particles:                       ',outcube
       Print *, 'Gather factor for cube output:                  ',gatherwrite
       Print *, 'Sort particles in cube files:                   ',sortcube
       Print *, 'Perform timings (imply extra synchronisations): ',dotimings       
       Print *, ' '
       
       
       ! Open log file
       filelog = trim(output_root)//'.log'
       Open(Unit=Ulog,file=filelog)
       Write(Ulog,*) 'Parallel FoF'
       Write(Ulog,*) procNB,' processes:'
       
       ! Write input parameters in .log file
       Write(Ulog,*) 'Input parameters'
       Write(Ulog, nml=ramsesinput)
       Write(Ulog, nml=fofparam)
       Write(Ulog, nml=outputparam)

    End If

  End Subroutine writeparameters


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
    Use modparameters
    Use modmpicom
    Use modtiming
    Implicit none

    !-----------------------------------------------
    ! Lecture du fichier particules au format Ramses
    !-----------------------------------------------

    ! Local variables
    Character(len=5)               :: ncharcpu
    Character(len=9)               :: tmpstr1, tmpstr2
    Character(len=400)             :: nomfich
    Character(len=13)              :: dumchar
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
    Integer(kind=4)                :: nproc       ! process number  read in RAMSES info file
    Integer(kind=4)                :: ncpu2       ! process number  read in RAMSES part files
    Integer(kind=4)                :: ndim2       ! dimension       read in RAMSES part files
    Integer(kind=4)                :: npartloc    ! particle number read in RAMSES part files
    Integer(kind=4)                :: prov, dest  ! provenance and destination process number for p2p MPI communications
    Integer(kind=4)                :: mpistat(Mpi_Status_Size)   ! status of MPI communication
    Integer(kind=PRI)              :: tmplongint              ! temp integer8 variable
    Integer(kind=PRI), allocatable :: tmpi(:), tmpsendi(:)    ! TYPE VARIABLE EN FONCTION DU NB DE PART
    Integer(kind=4)                :: errcode
    Integer(kind=4)                :: grpnb

    Real(kind=4), allocatable     :: tmpsimple(:)            ! temporary variable for Ramses v2 output
    Real(kind=8), allocatable     :: tmpdouble(:)            ! temporary variable for Ramses v3 output
    Real(kind=4), allocatable     :: tmpsendx(:,:),tmpsendv(:,:), tmpsendp(:), tmpsendf(:,:)
    Real(kind=4), allocatable     :: tmpx(:,:), tmpv(:,:), tmpp(:), tmpf(:,:)
    Real(kind=4)                  :: deltasd
    Integer(kind=4)               :: tmpinteger

    Integer(kind=4) :: buffersize, nbbytes, b_pos
    Character, dimension(:),allocatable :: buffer

    ! tampon pour Mpi_Pack
    Call Mpi_Sizeof(nres, nbbytes, mpierr)
    buffersize = 3* nbbytes
    Allocate (buffer(0:buffersize-1))
    b_pos = 0

    ! Initialisation timer
    time0 = Mpi_Wtime()

    grpchar = 'group_00001'

    ! Lecture parametres et remplissage du tampon pour diffusion des parametres
    If(procID == 0) Then
       If(code_index.eq.'RA2') Then 
          Print *,'Reading Ramses v2 output...'
          Write(Ulog,*) 'Reading Ramses v2 output...'
       Else if(code_index.eq.'RA3') Then
          Print *,'Reading Ramses v3 output...'
          Write(Ulog,*) 'Reading Ramses v3 output...'
       End If
       
       If( grpsize == 0 ) Then
          nomfich = trim(pathinput)//'/'//trim(nameinfo)
       Else
          nomfich = trim(pathinput)//'/'//trim(grpchar)//'/'//trim(nameinfo)
       End If
       Print *,'Reading RAMSES info file:',trim(nomfich)

       Open(11,file=nomfich, form='formatted', Status='Old', Iostat=errcode) 
       If(errcode > 0) Then
          Call EmergencyStop('Error opening '//trim(nameinfo)//' file',errcode)
       End If
       Rewind 11
       
       Read(11,'(A13,I11)') dumchar,nproc  ! number of processes
       Read(11,'(A13,I11)') dumchar,ndim
       Read(11,'(A13,I11)') dumchar,lmin

       Close(11)
       
       nres = 2**lmin

       Write(*,*) 'Number of:' 
       Write(*,'(A25,I6)') ' - files for each output:',nproc
       Write(*,'(A25,I6)') ' - dimensions:           ',ndim
       Write(*,'(A25,I6)') ' - grid points:          ',nres
       Write(Ulog,*) 'nb_proc = ',nproc,'ndim = ',ndim,'nres = ',nres

       Call Mpi_Pack(nproc, 1, Mpi_Integer, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack( ndim, 1, Mpi_Integer, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack( nres, 1, Mpi_Integer, buffer, buffersize, b_pos, Mpi_Comm_World, mpierr)

    End If

    Call Mpi_Bcast(buffer,buffersize,Mpi_Packed,0,Mpi_Comm_World,mpierr)

    If(procID /= 0) Then

       Call Mpi_Unpack(buffer, buffersize, b_pos, nproc, 1, Mpi_Integer, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,  ndim, 1, Mpi_Integer, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(buffer, buffersize, b_pos,  nres, 1, Mpi_Integer, Mpi_Comm_World, mpierr)

    End If

    ngrid = int(nres,kind=8)**3
    
    If(allocated(buffer)) Deallocate(buffer)

    If(procID==0) Print *,'Reading positions...'
    If(procID==0) Write(Ulog,*) 'Reading positions...'

    npart = 0

    nmod = mod(nproc,procNB)
    mynbfile = nproc / procNB
    If(procID <= nmod-1) Then
       mynbfile = mynbfile+1
       firstp   = procID * mynbfile + 1
       lastp    = (procID+1) * mynbfile
    Else
       firstp   = procID * mynbfile + 1 + nmod
       lastp    = (procID+1) * mynbfile + nmod
    End If
    
 
    mynpart = 0

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

    xmin =  CubeCoord(1)      * deltasd
    xmax = (CubeCoord(1) + 1) * deltasd
    ymin =  CubeCoord(2)      * deltasd
    ymax = (CubeCoord(2) + 1) * deltasd
    zmin =  CubeCoord(3)      * deltasd
    zmax = (CubeCoord(3) + 1) * deltasd


    Do icpu = firstp,lastp
       If( grpsize == 0 ) Then
          Write(ncharcpu(1:5),'(I5.5)') icpu
          nomfich = trim(pathinput)//'/'//trim(namepart)//trim(ncharcpu)
       Else
          Write(ncharcpu(1:5),'(I5.5)') icpu
          grpnb = (icpu-1)/grpsize + 1
          Write(grpchar(7:11),'(I5.5)') grpnb
          nomfich = trim(pathinput)//'/'//trim(grpchar)//'/'//trim(namepart)//trim(ncharcpu)
       End If

       Open(unit=1,file=nomfich,status='old',form='unformatted')
       Read(1) ncpu2
       Read(1) ndim2
       Read(1) npartloc
       Close(1)

       If((ncpu2/=nproc).Or.(ndim2/=ndim)) Then
          Call EmergencyStop('Files'//trim(nomfich)// ' and '//trim(nameinfo)//' are not consistent for ncpu and/or ndim',2)
       End If

       mynpart = mynpart + npartloc
    End Do

    tmplongint = mynpart
    Call Mpi_AllReduce(tmplongint,npart,1,MPI_PRI,Mpi_Sum,Mpi_Comm_World,mpierr)

    If(procID == 0) Then
       Write(* ,*)'There are ',npart,' DM particles'
       Write(Ulog,*)'There are ',npart,' DM particles'
    End If

    Allocate(tmpx(3,mynpart))
    Allocate(tmpv(3,mynpart))
    Allocate(tmpi(mynpart))
    If(force) Allocate(tmpf(3,mynpart))
    If(potential) Allocate(tmpp(mynpart))

    tmpx=0.
    mynpart = 0

    If(dotimings) Then
       Call Mpi_Barrier(MPICube,mpierr)
       timeInt = Mpi_Wtime()
       tInitRead = timeInt - time0
    End If

    Do icpu = firstp,lastp
       If( grpsize == 0 ) Then
          Write(ncharcpu(1:5),'(I5.5)') icpu
          nomfich = trim(pathinput)//'/'//trim(namepart)//trim(ncharcpu)
       Else
          Write(ncharcpu(1:5),'(I5.5)') icpu
          grpnb = (icpu-1)/grpsize + 1
          Write(grpchar(7:11),'(I5.5)') grpnb
          nomfich = trim(pathinput)//'/'//trim(grpchar)//'/'//trim(namepart)//trim(ncharcpu)
       End If


       Open(unit=1,file=nomfich,status='old',form='unformatted')
       Read(1) ncpu2
       Read(1) ndim2       
       Read(1) npartloc
       ramsesV2 : If(code_index.eq.'RA2') Then

          Allocate(tmpsimple(1:npartloc))
          ! Read positions
          Do idim = 1,ndim
             Read(1) tmpsimple
             ! put all positions in tmpx vector
             tmpx(idim,mynpart+1:mynpart+npartloc) = tmpsimple
          End Do
          
          ! Read velocities in a dummy variable
          Do idim = 1,ndim
             Read(1) tmpsimple
             ! put all velocities in tmpv vector
             tmpv(idim,mynpart+1:mynpart+npartloc) = tmpsimple
          End Do
          
          ! Read masses in a dummy variable
          Read(1) tmpsimple
          Deallocate(tmpsimple)
          
       End If ramsesV2

       ramsesV3 : If(code_index.eq.'RA3') Then 
          Allocate(tmpdouble(1:npartloc))
          
          Read(1)
          Read(1)
          Read(1)
          Read(1)
          Read(1)
          
          ! Read positions
          Do idim = 1,ndim
             Read(1) tmpdouble
             ! put all positions in tmpx vector
             tmpx(idim,mynpart+1:mynpart+npartloc) = real(tmpdouble, kind=4)
          End Do
          
          ! Read velocities in a dummy variable
          Do idim = 1,ndim
             Read(1) tmpdouble
             ! put all velocities in tmpv vector
             tmpv(idim,mynpart+1:mynpart+npartloc) = real(tmpdouble,kind=4)
          End Do
          
          ! Read masses in a dummy variable
          Read(1) tmpdouble

          If(star) Then
             Read(1) tmpdouble
             If(metal) Read(1) tmpdouble
          End If

       End If ramsesV3

       Deallocate(tmpdouble)

       ! Read particle id
       Read(1) tmpi(mynpart+1:mynpart+npartloc)

       ! Read potential if potential parameter is .true.
       If(potential) Then
          ! First: skip the level (integer 4)
          Read(1) tmpinteger
          ! Then read the potential: same kind as pos and vel
          Read(1) tmpp(mynpart+1:mynpart+npartloc)
       End If

       ! Read force if force parameter is .true.
       If(force) Then
          If( .not. potential) Then
             ! If the potential skip the level (integer 4)
             Read(1) tmpinteger
          End If
          ! Then read the force: same kind as pos and vel
          Do idim = 1,ndim
             Read(1) tmpf(idim,mynpart+1:mynpart+npartloc)
          End Do
       End If

       
       Do j = mynpart+1,mynpart+npartloc
          Do idim = 1,ndim
             If(tmpx(idim,j)==1.0) tmpx(idim,j) = 0.0
          End Do
          n_i = int(tmpx(1,j)/deltasd)
          n_j = int(tmpx(2,j)/deltasd)
          n_k = int(tmpx(3,j)/deltasd)
          ind = nsd**2 *n_i + nsd*n_j + n_k + 1
          npartvloc(ind) = npartvloc(ind)+1
       End Do
       
       Close(1)
       mynpart = mynpart+npartloc
    End Do


    Call Mpi_AllReduce(npartvloc,npartv,procNB,Mpi_Integer,Mpi_Sum,Mpi_Comm_World,mpierr)
    tReadfile = Mpi_Wtime() - timeInt
    timeInt = Mpi_Wtime()

    ! ------------------------------------------------
    ! Repartition des particules entre les processeurs
    ! ------------------------------------------------

    Allocate (pos(3,npartv(procID+1)))
    Allocate (vel(3,npartv(procID+1)))
    Allocate (id(npartv(procID+1)))
    If(force) Allocate(for(3, npartv(procID+1)))
    If(potential) Allocate(pot(npartv(procID+1)))

    recvpoint = 1

    processus : Do i = 1,procNB - 1
       dest = mod(procID + i,procNB)
       prov = mod(procID + procNB - i, procNB)

       Call Mpi_Cart_coords(MPICube,dest,3,destCoord,mpierr)
       Call Mpi_Isend(npartvloc(dest+1),1,Mpi_Integer,dest,procID,MpiCube,mpireqs1,mpierr)
       Call Mpi_Irecv(nrecv,1,Mpi_Integer,prov,prov,MpiCube,mpireqr1,mpierr)
       xmin =  destCoord(1)      * deltasd
       xmax = (destCoord(1) + 1) * deltasd
       ymin =  destCoord(2)      * deltasd
       ymax = (destCoord(2) + 1) * deltasd
       zmin =  destCoord(3)      * deltasd
       zmax = (destCoord(3) + 1) * deltasd

       If(npartvloc(dest+1)/=0) Then
          Allocate(tmpsendx(3,npartvloc(dest+1)))
          Allocate(tmpsendv(3,npartvloc(dest+1)))
          Allocate(tmpsendi(  npartvloc(dest+1)))
          If(force) Allocate(tmpsendf(3, npartvloc(dest+1)))
          If(potential) Allocate(tmpsendp(npartvloc(dest+1)))
          ind = 1

          Do j=1,mynpart
             If(tmpx(1,j)>=xmin .and. tmpx(1,j)<xmax .and. &
                  tmpx(2,j)>= ymin .and. tmpx(2,j) < ymax .and. &
                  tmpx(3,j)>= zmin .and. tmpx(3,j) < zmax) Then
                tmpsendx(:,ind) = tmpx(:,j)
                tmpsendv(:,ind) = tmpv(:,j)
                If(force) tmpsendf(:,ind) = tmpf(:,j)
                If(potential) tmpsendp(  ind) = tmpp(j)
                tmpsendi(  ind) = tmpi(j)
                ind=ind+1
             End If
          End Do

          If(ind/=npartvloc(dest+1)+1) Then
             Call EmergencyStop('Erreur dans la repartition des particules',2)
          End If

       End If

       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

       If(npartvloc(dest+1)/=0) Then
          Call Mpi_Isend(tmpsendx,3*npartvloc(dest+1),Mpi_Real,dest,1,MpiCube,mpireqs1,mpierr)
          Call Mpi_Isend(tmpsendv,3*npartvloc(dest+1),Mpi_Real,dest,2,MpiCube,mpireqs2,mpierr)
          Call Mpi_Isend(tmpsendi,  npartvloc(dest+1), MPI_PRI,dest,3,MpiCube,mpireqs3,mpierr)
          If(potential) Call Mpi_Isend(tmpp, npartvloc(dest+1), Mpi_Real,dest,4, MpiCube,mpireqs4,mpierr)
          If(force) Call Mpi_Isend(tmpsendv,3*npartvloc(dest+1),Mpi_Real,dest,5, MpiCube,mpireqs5,mpierr)
       End If
       If(nrecv/=0) Then
          Call Mpi_Irecv(pos(1,recvpoint),3*nrecv,Mpi_Real,prov,1,MpiCube,mpireqr1,mpierr)
          Call Mpi_Irecv(vel(1,recvpoint),3*nrecv,Mpi_Real,prov,2,MpiCube,mpireqr2,mpierr)
          Call Mpi_Irecv( id(recvpoint),  nrecv, MPI_PRI,prov,  3,MpiCube,mpireqr3,mpierr)
          If(potential) Call Mpi_Irecv(pot(recvpoint),nrecv,Mpi_Real,prov,4,MpiCube,mpireqr4,mpierr)
          If(force) Call Mpi_Irecv(for(1,recvpoint),3*nrecv,Mpi_Real,prov,5,MpiCube,mpireqr5,mpierr)
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
       
       If(potential) Then
          If(npartvloc(dest+1)/=0) Then
             Call Mpi_Wait(mpireqs4, mpistat, mpierr)
             Deallocate(tmpsendp)
          End If
          If(nrecv/=0) Then
             Call Mpi_Wait(mpireqr4, mpistat, mpierr)
          End If
       End If

       If(force) Then
          If(npartvloc(dest+1)/=0) Then
             Call Mpi_Wait(mpireqs5, mpistat, mpierr)
             Deallocate(tmpsendf)
          End If
          If(nrecv/=0) Then
             Call Mpi_Wait(mpireqr5, mpistat, mpierr)
          End If
       End If

    End Do processus

    xmin =  CubeCoord(1)      * deltasd
    xmax = (CubeCoord(1) + 1) * deltasd
    ymin =  CubeCoord(2)      * deltasd
    ymax = (CubeCoord(2) + 1) * deltasd
    zmin =  CubeCoord(3)      * deltasd
    zmax = (CubeCoord(3) + 1) * deltasd

    ind = 0
    Do j=1,mynpart
       If(tmpx(1,j)>=xmin .and. tmpx(1,j)<xmax .and. &
            tmpx(2,j)>= ymin .and. tmpx(2,j) < ymax .and. &
            tmpx(3,j)>= zmin .and. tmpx(3,j) < zmax) Then
          pos(:,recvpoint+ind) = tmpx(:,j)
          vel(:,recvpoint+ind) = tmpv(:,j)
          id(recvpoint+ind)  = tmpi(j)
          If(potential) pot(recvpoint+ind) = tmpp(j)
          If(force) for(:,recvpoint+ind) = tmpf(:,j)
          ind = ind+1
       End If
    End Do

    If(recvpoint+ind /= npartv(procID+1)+1) Then
       Write(tmpstr1,'(I9.9)') recvpoint+ind
       Write(tmpstr2,'(I9.9)') npartv(procID+1)+1
       Call EmergencyStop('Wrong particles number found after send/recv swaps:'//tmpstr1//' ; '//tmpstr2,2)
    End If

    mynpart = npartv(procID+1)


    Deallocate(tmpx,tmpv,tmpi)
    Deallocate(npartv, npartvloc)
    If(potential) Deallocate(tmpp)
    If(force) Deallocate(tmpf)
    tTailPart = Mpi_Wtime() - timeInt
    tRead = Mpi_Wtime() - time0

  End Subroutine ramses_lecture


  !=======================================================================
  !> This subroutine writes the position, velocity and id of each particle on the process in a hdf5 file.
  !! One file is written per MPI process
  Subroutine h5writecube()

    Use modhdf5
    Use modparameters
    Use modmpicom
    Use modxdmf
    Implicit none

    Character(len=400) :: filecube
    Character(len=391) :: filebase
    Character(len=5)  :: pid_char

    Character(len=16) :: dsetname                           ! Dataset name
    Character(len=16) :: aname                              ! Attribute name
    Character(len=20) :: adata
    Character(len=16) :: groupname

    Integer(hid_t) :: file_id                               ! File identifier
    Integer(hid_t) :: gr_root_id                            ! Group identifier

    Integer(kind=4) :: potattr
    Integer(kind=4) :: forattr
    Real(kind=4), dimension(6) :: boundaries
    Integer(kind=4), dimension(procNB) :: nparttab

#ifdef DEBUG
    Print *,"Enter h5writecube on process ",procID
#endif    
    
    Write(pid_char(1:5),'(I5.5)') procID
    filebase = trim(output_root)//'_cube'
    filecube = trim(output_root)//'_cube_'//pid_char//'.h5'
    
    Call Mpi_Gather(mynpart, 1, Mpi_Integer, nparttab, 1, Mpi_Integer, 0, MpiCube, mpierr)
!    If(procID==0) Then
!       Call writesinglecubexdmf(procNB,filebase,nparttab)
!    End If

    ! create the hdf5 file
    Call hdf5_create_file(filecube, file_id, origin)

    ! open the root group
    groupname = '/'
    Call hdf5_open_group(file_id,groupname, gr_root_id)

    ! Write type as attribute    
    aname = 'Type'
    adata = 'cube format'
    Call hdf5_write_attr(gr_root_id, aname, adata)
    
    ! Write the number of particles as an attribute
    aname = 'partNB'
    Call hdf5_write_attr(gr_root_id, aname, mynpart)

    ! Write the process ID as an attribute
    aname = 'procID'
    Call hdf5_write_attr(gr_root_id, aname, procID)

    aname = 'nfile'
    Call hdf5_write_attr(gr_root_id, aname, procNB)
    
    ! Write nres as an attribute
    aname = 'nres'
    Call hdf5_write_attr(gr_root_id, aname, nres)

    ! Write potential logical as an integer attribute (1=true, 0=false)
    aname = 'potential'
    If(potential) Then
       potattr=1
    Else
       potattr=0
    End If
    Call hdf5_write_attr(gr_root_id, aname, potattr)

    ! Write force logical as an integer attribute (1=true, 0=false)
    aname = 'force'
    If(force) Then
       forattr=1
    Else
       forattr=0
    End If
    Call hdf5_write_attr(gr_root_id, aname, forattr)

    ! Write the boundaries of the cube as an attribute
    boundaries(1) = xmin
    boundaries(2) = xmax
    boundaries(3) = ymin
    boundaries(4) = ymax
    boundaries(5) = zmin
    boundaries(6) = zmax
    aname = 'boundaries'

    Call hdf5_write_attr(gr_root_id, aname, 6, boundaries)

    ! Write the position of the particles
    dsetname='pos'
    Call hdf5_write_data(gr_root_id, dsetname, 3, mynpart, pos)

    ! Write the velocity of the particles
    dsetname='vel'
    Call hdf5_write_data(gr_root_id, dsetname, 3, mynpart, vel) 

    ! If we use potential, write potential of the particles
    If(force) Then
       dsetname = 'for'
       Call hdf5_write_data(gr_root_id, dsetname, 3, mynpart, for)
    End If   

    ! If we use potential, write potential of the particles
    If(potential) Then
       dsetname = 'pot'
       Call hdf5_write_data(gr_root_id, dsetname, mynpart, pot)
    End If

    ! Write the ID of the particles
    dsetname='ID'
    Call hdf5_write_data(gr_root_id, dsetname, mynpart, id) 

    ! Close the root group.
    Call hdf5_close_group(gr_root_id)

    Call hdf5_close_file(file_id)

#ifdef DEBUG
    Print *,'Exit h5writecube on process ',procID
#endif    

  End Subroutine h5writecube



  !=======================================================================
  !> This subroutine reads hdf5 cube files created by pFOF.
  Subroutine h5readcube()

    Use modhdf5
    Use modparameters
    Use modmpicom


    Character(len=400) :: filename                          ! File name
    Integer(hid_t) :: file_id                               ! File identifier
    Integer(hid_t) :: gr_root_id                            ! Group identifier
    Character(len=5)  :: pid_char
    Character(len=16) :: aname                              ! Attribute name
    Character(len=16) :: dsetname                           ! Dataset name
    Character(len=16) :: groupname
    Real(kind=4), dimension(6) :: boundaries
    Integer(kind=4) :: potattr
    Integer(kind=4) :: forattr

#ifdef DEBUG
    Print *,"Enter h5readcube on process ",procID
#endif

    Write(pid_char(1:5),'(I5.5)') procID
    filename = trim(output_root)//'_cube_'//pid_char//'.h5'

    Print *,'filename : ', filename
    ! open the file
    Call hdf5_open_file(filename, file_id)


    ! open the root group
    groupname = '/'
    Call hdf5_open_group(file_id,groupname, gr_root_id)

    ! read attribute partNB
    aname = 'partNB'
    Call hdf5_read_attr(gr_root_id, aname, mynpart)

    ! read attribute nres
    aname = 'nres'
    Call hdf5_read_attr(gr_root_id, aname, nres)

    ngrid = int(nres,kind=8)**3
    npart = ngrid

    ! read attribute boundaries
    aname = 'boundaries'
    Call hdf5_read_attr(gr_root_id, aname, 6, boundaries)

    xmin=boundaries(1)
    xmax=boundaries(2)
    ymin=boundaries(3)
    ymax=boundaries(4)
    zmin=boundaries(5)
    zmax=boundaries(6)

    potattr=0
    ! if requested, try to read the potential attribute
    If(potential) Then
       aname = 'potential'
       Call hdf5_read_attr(gr_root_id, aname, potattr)
       If(potattr==0) Then ! potential is asked by the user but not found in the cube file: abort
          If(procID==0) Then
             Print *,'The potential was not written in the cube files used for this analysis.'
             Print *,'You should change the potential parameter to false.'
          End If
          Call EmergencyStop('Potential missing in the cube files: abort.',100)
       End If
       If(.not.Allocated(pot)) Allocate(pot(mynpart))
       dsetname = 'pot'
       Call hdf5_read_data(gr_root_id, dsetname, mynpart, pot)
    End If

    forattr=0
    ! if requested, try to read the force attribute
    If(force) Then
       aname = 'force'
       Call hdf5_read_attr(gr_root_id, aname, forattr)
       If(forattr==0) Then ! force is asked by the user but not found in the cube file: abort
          If(procID==0) Then
             Print *,'The force was not written in the cube files used for this analysis.'
             Print *,'You should change the force parameter to false.'
          End If
          Call EmergencyStop('Force missing in the cube files: abort.',100)
       End If
       If(.not.Allocated(for)) Allocate(for(3,mynpart))
       dsetname = 'for'
       Call hdf5_read_data(gr_root_id, dsetname, 3, mynpart, for)
    End If

    If(.not.Allocated(id)) Allocate(id(mynpart))
    If(.not.Allocated(pos))  Allocate(pos(3,mynpart))
    If(.not.Allocated(vel))  Allocate(vel(3,mynpart))

    ! read position of the particles
    dsetname = 'pos'
    Call hdf5_read_data(gr_root_id, dsetname, 3, mynpart, pos)

    ! read velocity of the particles
    dsetname = 'vel'
    Call hdf5_read_data(gr_root_id, dsetname, 3, mynpart, vel)

    ! read id of the particles
    dsetname = 'ID'
    Call hdf5_read_data(gr_root_id, dsetname, mynpart, id)


    ! Close the root group.
    Call hdf5_close_group(gr_root_id)

    Call hdf5_close_file(file_id)


#ifdef DEBUG
    Print *,"Exit h5readcube on process ",procID
#endif

  End Subroutine h5readcube



  !=======================================================================
  !> This subroutine reads hdf5 sorted cube files created by pFOF.
  Subroutine h5readsortedcube()

    Use modhdf5
    Use modparameters
    Use modmpicom


    Character(len=400) :: filename                          ! File name
    Integer(hid_t) :: file_id                               ! File identifier
    Integer(hid_t) :: gr_root_id                            ! Group identifier
    Integer(hid_t) :: gr_id                                 ! Group identifier
    Character(len=5)  :: pid_char
    Character(len=8)  :: gid_char
    Character(len=16) :: aname                              ! Attribute name
    Character(len=16) :: dsetname                           ! Dataset name
    Character(len=16) :: groupname
    Real(kind=4), dimension(6) :: boundaries
    Integer(kind=4) :: potattr
    Integer(kind=4) :: forattr
    Integer(kind=4) :: ngroup
    Integer(kind=4), dimension(:), allocatable :: npartpergroup
    Integer(kind=4) :: igroup
    Integer(kind=4) :: indbeg, indend

#ifdef DEBUG
    Print *,"Enter h5readsortedcube on process ",procID
#endif

    Write(pid_char(1:5),'(I5.5)') procID
    filename = trim(output_root)//'_sortedcube_'//pid_char//'.h5'

    ! open the file
    Call hdf5_open_file(filename, file_id)

    ! open the root group
    groupname = '/'
    Call hdf5_open_group(file_id,groupname, gr_root_id)

    ! read attribute partNB
    aname = 'partNB'
    Call hdf5_read_attr(gr_root_id, aname, mynpart)

    ! read attribute nres
    aname = 'nres'
    Call hdf5_read_attr(gr_root_id, aname, nres)

    ngrid = int(nres,kind=8)**3
    npart = ngrid

    ! read attribute boundaries
    aname = 'boundaries'
    Call hdf5_read_attr(gr_root_id, aname, 6, boundaries)

    xmin=boundaries(1)
    xmax=boundaries(2)
    ymin=boundaries(3)
    ymax=boundaries(4)
    zmin=boundaries(5)
    zmax=boundaries(6)

    ! read attribute ngroup
    aname = 'ngroup'
    Call hdf5_read_attr(gr_root_id, aname, ngroup)

    Allocate(npartpergroup(ngroup))
    aname = 'npartpergroup'
    Call hdf5_read_data(gr_root_id, aname, ngroup, npartpergroup)

    potattr=0
    ! if requested, try to read the potential attribute
    If(potential) Then
       aname = 'potential'
       Call hdf5_read_attr(gr_root_id, aname, potattr)
       If(potattr==0) Then ! potential is asked by the user but not found in the cube file: abort
          If(procID==0) Then
             Print *,'The potential was not written in the cube files used for this analysis.'
             Print *,'You should change the potential parameter to false.'
          End If
          Call EmergencyStop('Potential missing in the cube files: abort.',100)
       End If
       If(.not.Allocated(pot)) Allocate(pot(mynpart))
    End If

    forattr=0
    ! if requested, try to read the potential attribute
    If(force) Then
       aname = 'force'
       Call hdf5_read_attr(gr_root_id, aname, forattr)
       If(forattr==0) Then ! force is asked by the user but not found in the cube file: abort
          If(procID==0) Then
             Print *,'The force was not written in the cube files used for this analysis.'
             Print *,'You should change the force parameter to false.'
          End If
          Call EmergencyStop('force missing in the cube files: abort.',100)
       End If
       If(.not.Allocated(for)) Allocate(for(3,mynpart))
    End If

    If(.not.Allocated(id)) Allocate(id(mynpart))
    If(.not.Allocated(pos))  Allocate(pos(3,mynpart))
    If(.not.Allocated(vel))  Allocate(vel(3,mynpart))

    indbeg = 1
    indend = 0
    ! loop over the groups
    Do igroup=1, ngroup
       If(npartpergroup(igroup) /= 0) Then ! there is at least one part. to read in this group
          Write(gid_char(1:8),'(I8.8)') igroup
          groupname='group'//gid_char
          Call hdf5_open_group(gr_root_id,groupname,gr_id)
          
          indend = indbeg + npartpergroup(igroup) - 1
          ! read position of the particles
          dsetname = 'pos'
          Call hdf5_read_data(gr_id, dsetname, 3, mynpart, pos(:,indbeg:indend))
          
          ! read velocity of the particles
          dsetname = 'vel'
          Call hdf5_read_data(gr_id, dsetname, 3, mynpart, vel(:,indbeg:indend))
          
          ! read id of the particles
          dsetname = 'ID'
          Call hdf5_read_data(gr_id, dsetname, mynpart, id(indbeg:indend))
          
          ! read potential if requested
          If(potential) Then
             dsetname = 'pot'
             Call hdf5_read_data(gr_id, dsetname, mynpart, pot(indbeg:indend))
          End If

          ! read force if requested
          If(force) Then
             dsetname = 'for'
             Call hdf5_read_data(gr_id, dsetname, 3, mynpart, for(:,indbeg:indend))
          End If
          
          indbeg = indend + 1

          Call hdf5_close_group(gr_id)
       End If
    End Do

    If(indend /= mynpart) Then
       Print *,'Error while reading particles from file ',filename
       Call EmergencyStop('Error in h5readsortedcube',100)
    End If

    ! Close the root group.
    Call hdf5_close_group(gr_root_id)

    Call hdf5_close_file(file_id)

    Deallocate(npartpergroup)

#ifdef DEBUG
    Print *,"Exit h5readsortedcube on process ",procID
#endif

  End Subroutine h5readsortedcube



  !=======================================================================
  !> This subroutine writes the position, the velocity and the id of each particle 
  !! on the process in a hdf5 file. 
  !! The particles are gathered 
  Subroutine mpih5writecube()

    Use modhdf5
    Use modparameters
    Use modmpicom
    Implicit none

    Character(len=400) :: filecube
    Character(len=5)  :: pid_char

    Character(len=16) :: dsetname                           ! Dataset name
    Character(len=16) :: aname                              ! Attribute name
    Character(len=16) :: groupname
    Character(len=20) :: adata

    Integer(hid_t) :: file_id                               ! File identifier
    Integer(hid_t) :: gr_root_id                            ! Group identifier

    Integer(kind=4) :: procperfile
    Integer(kind=4) :: npart
    Integer(kind=4) :: nfile
    Integer(kind=4) :: potattr
    Integer(kind=4) :: forattr
    Integer(kind=4), dimension(:), allocatable :: partnb_tab

    Real(kind=4), dimension(6) :: boundaries
    Real(kind=4), dimension(:,:), allocatable :: bound_tab

#ifdef DEBUG
    Print *,"Enter mpih5writecube on process ",procID
#endif    

    ! number of processes writing in the same file
    procperfile = gatherwrite**3
        
    Allocate(partnb_tab(procperfile))
    Allocate(bound_tab(6,procperfile))

    Write(pid_char(1:5),'(I5.5)') commcolorWrite
    filecube = trim(output_root)//'_mpicube_'//pid_char//'.h5'

    ! create the hdf5 file
    Call hdf5_create_mpi_file(filecube, MpiSubCubeWrite, file_id, origin)

    ! open the root group
    groupname = '/'
    Call hdf5_open_group(file_id,groupname, gr_root_id)

    ! Offset and global dimensions have to be computed
    Call Mpi_Allgather(mynpart, 1, Mpi_Integer, partnb_tab, 1, Mpi_Integer, MpiSubCubeWrite, mpierr)

    aname='partNB_tab'
    Call hdf5_write_attr(gr_root_id, aname, procperfile, partnb_tab)

    ! Write type as attribute    
    aname = 'Type'
    adata = 'mpicube format'
    Call hdf5_write_attr(gr_root_id, aname, adata)

    npart = sum(partnb_tab)
    aname = 'partNB'
    Call hdf5_write_attr(gr_root_id, aname, npart)
    
    aname='nres'
    Call hdf5_write_attr(gr_root_id, aname, nres)

    aname = 'nfile'
    nfile = procNB / procperfile
    Call hdf5_write_attr(gr_root_id, aname, nfile)

    aname = 'procID'
    Call hdf5_write_attr(gr_root_id, aname, commcolorWrite) 

    ! Write the boundaries of the cube as an attribute
    boundaries(1) = xmin
    boundaries(2) = xmax
    boundaries(3) = ymin
    boundaries(4) = ymax
    boundaries(5) = zmin
    boundaries(6) = zmax

    Call Mpi_Allgather(boundaries,6,Mpi_Real, bound_tab, 6, Mpi_Real, MpiSubCubeWrite, mpierr)

    aname = 'boundaries_tab'
    Call hdf5_write_attr(gr_root_id, aname, 6, procperfile, bound_tab)

    boundaries(1) = minval(bound_tab(1,:))
    boundaries(2) = maxval(bound_tab(2,:))
    boundaries(3) = minval(bound_tab(3,:))
    boundaries(4) = maxval(bound_tab(4,:))
    boundaries(5) = minval(bound_tab(5,:))
    boundaries(6) = maxval(bound_tab(6,:))
    aname = 'boundaries'
    Call hdf5_write_attr(gr_root_id, aname, 6, boundaries)

    aname = 'potential'
    If(potential) Then
       potattr = 1
    Else
       potattr = 0
    End If
    Call hdf5_write_attr(gr_root_id, aname, potattr) 

    aname = 'force'
    If(force) Then
       forattr = 1
    Else
       forattr = 0
    End If
    Call hdf5_write_attr(gr_root_id, aname, forattr) 

    dsetname = 'pos'
    Call hdf5_write_mpi_data(gr_root_id, dsetname, 3, mynpart, pos, MpiSubCubeWrite)

    dsetname = 'vel'
    Call hdf5_write_mpi_data(gr_root_id, dsetname, 3, mynpart, vel, MpiSubCubeWrite)

    dsetname = 'ID'
    Call hdf5_write_mpi_data(gr_root_id, dsetname, mynpart, id, MpiSubCubeWrite)

    If(potential) Then
       dsetname = 'pot'
       Call hdf5_write_mpi_data(gr_root_id, dsetname, mynpart, pot, MpiSubCubeWrite)
    End If

    If(force) Then
       dsetname = 'for'
       Call hdf5_write_mpi_data(gr_root_id, dsetname, 3, mynpart, for, MpiSubCubeWrite)
    End If

    ! Close the root group.
    Call hdf5_close_group(gr_root_id)

    Deallocate(partnb_tab, bound_tab)

    ! Close h5 file
    Call hdf5_close_mpi_file(file_id)

#ifdef DEBUG
    Print *,"Exit mpih5writecube on process ",procID
#endif    

  End Subroutine mpih5writecube


  !=======================================================================
 
  Subroutine mpih5readcube()
    
    Use modhdf5
    Use modparameters
    Use modmpicom
    Implicit none

    Character(len=400) :: filecube
    Character(len=5)  :: pid_char

    Character(len=16) :: dsetname                           ! Dataset name
    Character(len=16) :: aname                              ! Attribute name
    Character(len=16) :: groupname

    Integer(hid_t) :: file_id                               ! File identifier
    Integer(hid_t) :: gr_root_id                            ! Group identifier
    
    Integer(kind=4) :: procperfile
    Integer(kind=4), dimension(:), allocatable :: partnb_tab
    Integer(kind=4) :: potattr
    Integer(kind=4) :: forattr

    Real(kind=4), dimension(:,:), allocatable :: bound_tab

    
#ifdef DEBUG
    Print *,"Enter mpih5readcube on process ",procID
#endif    
    
    ! number of processes writing in the same file
    procperfile = gatherread**3
    
    Allocate(partnb_tab(procperfile))
    Allocate(bound_tab(6,procperfile))
        
    Write(pid_char(1:5),'(I5.5)') commcolorRead
    filecube = trim(output_root)//'_mpicube_'//pid_char//'.h5'

    Call hdf5_open_mpi_file(filecube, MpiSubCubeRead, file_id)

     ! open the root group
    groupname = '/'
    Call hdf5_open_group(file_id,groupname, gr_root_id)
    ! read attribute partNB
    aname = 'nres'
    Call hdf5_read_attr(gr_root_id, aname, nres)

#ifdef DEBUG
    Print *,"Value of nres read = ",nres, " on process ",procID
    Call Mpi_Barrier(Mpi_Comm_World, mpierr)
#endif

    ngrid = int(nres,kind=8)**3
    npart = ngrid
    
    ! read attribute partNB
    aname = 'partNB_tab'
    Call hdf5_read_attr(gr_root_id, aname, 6, partnb_tab)

#ifdef DEBUG
    Print *,"partNB_tab read on process ",procID
    Call Mpi_Barrier(Mpi_Comm_World, mpierr)
#endif    

    mynpart = partnb_tab(scprocIDRead+1)

    ! read attribute boundaries
    aname = 'boundaries_tab'
    Call hdf5_read_attr(gr_root_id, aname, 6, procperfile, bound_tab)

#ifdef DEBUG
    Print *,"boundaries_tab read on process ",procID
    Call Mpi_Barrier(Mpi_Comm_World, mpierr)
#endif    


    xmin=bound_tab(1,scprocIDRead+1)
    xmax=bound_tab(2,scprocIDRead+1)
    ymin=bound_tab(3,scprocIDRead+1)
    ymax=bound_tab(4,scprocIDRead+1)
    zmin=bound_tab(5,scprocIDRead+1)
    zmax=bound_tab(6,scprocIDRead+1)

    potattr=0
    ! if requested, try to read the potential attribute
    If(potential) Then
       aname = 'potential'
       Call hdf5_read_attr(gr_root_id, aname, potattr)
       If(potattr==0) Then ! potential is asked by the user but not found in the cube file: abort
          If(procID==0) Then
             Print *,'The potential was not written in the cube files used for this analysis.'
             Print *,'You should change the potential parameter to false.'
          End If
          Call EmergencyStop('Potential missing in the cube files: abort.',100)
       End If
       If(.not.Allocated(pot)) Allocate(pot(mynpart))
       dsetname = 'pot'
       Call hdf5_read_mpi_data(gr_root_id, dsetname, mynpart, pot, MpiSubCubeRead)
    End If

    forattr=0
    ! if requested, try to read the force attribute
    If(force) Then
       aname = 'force'
       Call hdf5_read_attr(gr_root_id, aname, forattr)
       If(forattr==0) Then ! force is asked by the user but not found in the cube file: abort
          If(procID==0) Then
             Print *,'The force was not written in the cube files used for this analysis.'
             Print *,'You should change the force parameter to false.'
          End If
          Call EmergencyStop('Force missing in the cube files: abort.',100)
       End If
       If(.not.Allocated(for)) Allocate(for(3,mynpart))
       dsetname = 'for'
       Call hdf5_read_mpi_data(gr_root_id, dsetname, 3, mynpart, for, MpiSubCubeRead)
    End If

    If(.not.Allocated(id)) Allocate(id(mynpart))
    If(.not.Allocated(pos))  Allocate(pos(3,mynpart))
    If(.not.Allocated(vel))  Allocate(vel(3,mynpart))
    
    dsetname = 'pos'
    Call hdf5_read_mpi_data(gr_root_id, dsetname, 3, mynpart, pos, MpiSubCubeRead)

#ifdef DEBUG
    Print *,"pos read on process ",procID
    Call Mpi_Barrier(Mpi_Comm_World, mpierr)
#endif    

    dsetname = 'vel'
    Call hdf5_read_mpi_data(gr_root_id, dsetname, 3, mynpart, vel, MpiSubCubeRead)

#ifdef DEBUG
    Print *,"vel read on process ",procID
    Call Mpi_Barrier(Mpi_Comm_World, mpierr)
#endif    

    dsetname = 'ID'
    Call hdf5_read_mpi_data(gr_root_id, dsetname, mynpart, id, MpiSubCubeRead)

#ifdef DEBUG
    Print *,"ID read on process ",procID
    Call Mpi_Barrier(Mpi_Comm_World, mpierr)
#endif    
  
    ! Close the root group.
    Call hdf5_close_group(gr_root_id)

    Call hdf5_close_file(file_id)
    
    Deallocate(partnb_tab)
    Deallocate(bound_tab)


#ifdef DEBUG
    Print *,"Exit mpih5readcube on process ",procID
#endif    
    
    
  End Subroutine mpih5readcube


  !=======================================================================

  Subroutine h5writesortedcube()

    Use modhdf5
    Use modparameters
    Use modmpicom
    Use modxdmf
    Use modsortinterf
    Implicit none

    Character(len=400) :: filecube
    Character(len=391) :: filebase
    Character(len=5)  :: pid_char
    Character(len=8) :: charic

    Character(len=16) :: dsetname                           ! Dataset name
    Character(len=16) :: aname                              ! Attribute name
    Character(len=20) :: adata
    Character(len=16) :: groupname

    Integer(hid_t) :: file_id                               ! File identifier
    Integer(hid_t) :: gr_root_id                            ! Group identifier
    Integer(hid_t) :: gr_id

    Real(kind=4), dimension(6) :: boundaries

    Integer(kind=4) :: potattr
    Integer(kind=4) :: forattr
    Integer(kind=4) :: ic, ix, iy, iz, ip, deb, fin, nc, ncdim, deltam1
    Integer(kind=4), dimension(:), allocatable :: npartcube
    Integer(kind=8), dimension(:), allocatable :: ictable
    Logical(kind=4) :: fileexist, fileopened

#ifdef DEBUG
    Print *,"Enter h5writesortedcube on process ",procID
#endif    

    ! each FoF cube is divided into nc=512=8x8x8 groups, with ncdim=8
    ncdim = 8
    nc = ncdim*ncdim*ncdim
    Allocate(ictable(mynpart))
    Allocate(npartcube(nc))
    npartcube = 0
    ! Ramses coarse grid is composed of nres^3 cells 
    ! on a process: nres^3 / procNB cells 
    ! => there is (nres/(ncdim*dims(1))^3 coarse cells in each group
    ! the size of a group is nres/(ncdim*dims(1))*(1/nres) where 1/nres is the size of 1 Ramses coarse cell
    ! => delta = 1/(ncdim*dims(1)) => 1/delta = ncdim*dims(1)
    deltam1 = ncdim*dims(1)

    ! We compute the "group" index of each particle and the number of particle in each group
    Do ip=1, mynpart
       ix = int((pos(1,ip) - xmin)*deltam1 + 1)
       iy = int((pos(2,ip) - ymin)*deltam1 + 1)
       iz = int((pos(3,ip) - zmin)*deltam1 + 1)
       ! rounding issue
       If(ix>ncdim) ix=ncdim
       If(iy>ncdim) iy=ncdim
       If(iz>ncdim) iz=ncdim
       ic = ix + (iy-1)*ncdim + (iz-1)*ncdim*ncdim
       ictable(ip) = ic
       npartcube(ic) = npartcube(ic) + 1
    End Do

    ! We sort the particles along their group id
    If(potential .and. force) Then
!       Call quicksort(1,mynpart,ictable,pos,vel,pot,id)
       Call heapsort(mynpart,ictable,pos,vel,for,pot,id)
    Else If (potential .and. .not. force) Then
       !       Call quicksort(1,mynpart,ictable,pos,vel,id)
       Call heapsort(mynpart,ictable,pos,vel,pot,id)
    Else If (force .and. .not. potential) Then
       Call heapsort(mynpart,ictable,pos,vel,for,id)
    Else
       Call heapsort(mynpart,ictable,pos,vel,id)
    End If

    ! We open a sortedcube file and write the groups into it
    Write(pid_char(1:5),'(I5.5)') procID
    filebase = trim(output_root)//'_sortedcube'
    filecube = trim(output_root)//'_sortedcube_'//pid_char//'.h5'

    Inquire(File=filecube,exist=fileexist, opened=fileopened)
    If(fileopened) Then
       Print *,'Error on process ',procID,' : File ',filecube,' already opened'
    End If

    ! create the hdf5 file
    Call hdf5_create_file(filecube, file_id, origin)

    ! open the root group
    groupname = '/'
    Call hdf5_open_group(file_id,groupname, gr_root_id)

    ! Write type as attribute    
    aname = 'Type'
    adata = 'sortedcube format'
    Call hdf5_write_attr(gr_root_id, aname, adata)
    
    ! Write the number of particles as an attribute
    aname = 'partNB'
    Call hdf5_write_attr(gr_root_id, aname, mynpart)

    ! Write the process ID as an attribute
    aname = 'procID'
    Call hdf5_write_attr(gr_root_id, aname, procID)

    aname = 'nfile'
    Call hdf5_write_attr(gr_root_id, aname, procNB)
    
    ! Write nres as an attribute
    aname = 'nres'
    Call hdf5_write_attr(gr_root_id, aname, nres)

    ! Write the boundaries of the cube as an attribute
    boundaries(1) = xmin
    boundaries(2) = xmax
    boundaries(3) = ymin
    boundaries(4) = ymax
    boundaries(5) = zmin
    boundaries(6) = zmax
    aname = 'boundaries'
    Call hdf5_write_attr(gr_root_id, aname, 6, boundaries)

    aname = 'potential'
    If(potential) Then
       potattr = 1
    Else
       potattr = 0
    End If
    Call hdf5_write_attr(gr_root_id, aname, potattr) 

    aname = 'force'
    If(force) Then
       forattr = 1
    Else
       forattr = 0
    End If
    Call hdf5_write_attr(gr_root_id, aname, forattr) 

    aname='ngroup'
    Call hdf5_write_attr(gr_root_id, aname, nc)

    aname='1/groupsize'
    Call hdf5_write_attr(gr_root_id, aname, deltam1)

    dsetname = 'npartpergroup'
    Call hdf5_write_data(gr_root_id, dsetname, nc, npartcube)

    deb = 1
    ! For each non empty group we create an HDF5 group and write dataset into it
    Do ic = 1, nc
       If(npartcube(ic) /= 0) Then
          Write(charic(1:8),'(I8.8)') ic
          groupname = 'group'//charic
          
          Call hdf5_create_group(gr_root_id,groupname,gr_id)
       
          fin = deb + npartcube(ic) - 1

          ! Write the position of the particles
          dsetname='pos'
          Call hdf5_write_data(gr_id, dsetname, 3, npartcube(ic), pos(:,deb:fin))
          
          ! Write the velocity of the particles
          dsetname='vel'
          Call hdf5_write_data(gr_id, dsetname, 3, npartcube(ic), vel(:,deb:fin)) 
          
          ! Write the ID of the particles
          dsetname='ID'
          Call hdf5_write_data(gr_id, dsetname, npartcube(ic), id(deb:fin)) 
          
          ! Write potential if it is used
          If(potential) Then
             dsetname = 'pot'
             Call hdf5_write_data(gr_id, dsetname, npartcube(ic), pot(deb:fin))
          End If

          ! Write force if it is used
          If(force) Then
             dsetname = 'for'
             Call hdf5_write_data(gr_id, dsetname, 3, npartcube(ic), for(:,deb:fin))
          End If
          
          Call hdf5_close_group(gr_id)
          deb = fin + 1

       End If
    End Do

    ! Close the root group.
    Call hdf5_close_group(gr_root_id)
          
    Call hdf5_close_file(file_id)

#ifdef DEBUG
    Print *,'Exit h5writesortedcube on process ',procID
#endif    

  End Subroutine h5writesortedcube




  !=======================================================================
  !> This subroutine writes the position, the velocity and the id of each particle 
  !! on the process in a hdf5 file. 
  !! The particles are gathered 
  Subroutine mpih5writesortedcube()

    Use modhdf5
    Use modparameters
    Use modmpicom
    Use modsortinterf
    Implicit none

    Character(len=400) :: filecube
    Character(len=5)  :: pid_char
    Character(len=8) :: charic

    Character(len=16) :: dsetname                           ! Dataset name
    Character(len=16) :: aname                              ! Attribute name
    Character(len=16) :: groupname
    Character(len=20) :: adata

    Integer(hid_t) :: file_id                               ! File identifier
    Integer(hid_t) :: gr_root_id                            ! Group identifier
    Integer(hid_t) :: gr_id

    Integer(kind=4) :: procperfile
    Integer(kind=4) :: npart
    Integer(kind=4) :: nfile
    Integer(kind=4) :: ic, ix, iy, iz, ip, deb, fin, nc, ncdim, deltam1, fic
    Integer(kind=4) :: potattr
    Integer(kind=4) :: forattr
    Integer(kind=4), dimension(:), allocatable :: npartcube
    Integer(kind=8), dimension(:), allocatable :: ictable

    Integer(kind=4), dimension(:), allocatable :: partnb_tab

    Real(kind=4), dimension(6) :: boundaries
    Real(kind=4), dimension(:,:), allocatable :: bound_tab

    Logical :: empty

#ifdef DEBUG
    Print *,"Enter mpih5writesortedcube on process ",procID
#endif    


    ! each FoF cube is divided into nc=512=8x8x8 groups, with ncdim=8
    ncdim = 8
    nc = ncdim*ncdim*ncdim
    Allocate(ictable(mynpart))
    Allocate(npartcube(nc))
    npartcube = 0
    ! Ramses coarse grid is composed of nres^3 cells 
    ! on a process: nres^3 / procNB cells 
    ! => there is (nres/(ncdim*dims(1))^3 coarse cells in each group
    ! the size of a group is nres/(ncdim*dims(1))*(1/nres) where 1/nres is the size of 1 Ramses coarse cell
    ! => delta = 1/(ncdim*dims(1)) => 1/delta = ncdim*dims(1)
    deltam1 = ncdim*dims(1)

    ! We compute the "group" index of each particle and the number of particle in each group
    Do ip=1, mynpart
       ix = int((pos(1,ip) - xmin)*deltam1 + 1)
       iy = int((pos(2,ip) - ymin)*deltam1 + 1)
       iz = int((pos(3,ip) - zmin)*deltam1 + 1)
       ! rounding issue
       If(ix>ncdim) ix=ncdim
       If(iy>ncdim) iy=ncdim
       If(iz>ncdim) iz=ncdim
       ic = ix + (iy-1)*ncdim + (iz-1)*ncdim*ncdim
       ictable(ip) = ic
       npartcube(ic) = npartcube(ic) + 1
    End Do

    ! We sort the particles along their group id
    If(potential .and. force) Then
!       Call quicksort(1,mynpart,ictable,pos,vel,pot,id)
       Call heapsort(mynpart,ictable,pos,vel,for,pot,id)
    Else If (potential .and. .not. force) Then
       !       Call quicksort(1,mynpart,ictable,pos,vel,id)
       Call heapsort(mynpart,ictable,pos,vel,pot,id)
    Else If (force .and. .not. potential) Then
       Call heapsort(mynpart,ictable,pos,vel,for,id)
    Else
       Call heapsort(mynpart,ictable,pos,vel,id)
    End If


#ifdef DEBUG
    Print *,'Particles sorted'
#endif


    ! number of processes writing in the same file
    procperfile = gatherwrite**3
        
    Allocate(partnb_tab(procperfile))
    Allocate(bound_tab(6,procperfile))

    Write(pid_char(1:5),'(I5.5)') commcolorWrite
    filecube = trim(output_root)//'_mpisortedcube_'//pid_char//'.h5'

    ! create the hdf5 file
    Call hdf5_create_mpi_file(filecube, MpiSubCubeWrite, file_id, origin)

    ! open the root group
    groupname = '/'
    Call hdf5_open_group(file_id,groupname, gr_root_id)

    ! Offset and global dimensions have to be computed
    Call Mpi_Allgather(mynpart, 1, Mpi_Integer, partnb_tab, 1, Mpi_Integer, MpiSubCubeWrite, mpierr)

    aname='partNB_tab'
    Call hdf5_write_attr(gr_root_id, aname, procperfile, partnb_tab)

    aname='ngroup'
    Call hdf5_write_attr(gr_root_id, aname, nc*procperfile)

    aname='1/groupsize'
    Call hdf5_write_attr(gr_root_id, aname, deltam1)

    dsetname = 'npartpergroup'
    Call hdf5_write_mpi_data(gr_root_id, dsetname, nc, npartcube, MpiSubCubeWrite)


    ! Write type as attribute    
    aname = 'Type'
    adata = 'mpisortedcube format'
    Call hdf5_write_attr(gr_root_id, aname, adata)

    npart = sum(partnb_tab)
    aname = 'partNB'
    Call hdf5_write_attr(gr_root_id, aname, npart)
    
    aname='nres'
    Call hdf5_write_attr(gr_root_id, aname, nres)

    aname = 'nfile'
    nfile = procNB / procperfile
    Call hdf5_write_attr(gr_root_id, aname, nfile)

    aname = 'procID'
    Call hdf5_write_attr(gr_root_id, aname, commcolorWrite) 

    ! Write the boundaries of the cube as an attribute
    boundaries(1) = xmin
    boundaries(2) = xmax
    boundaries(3) = ymin
    boundaries(4) = ymax
    boundaries(5) = zmin
    boundaries(6) = zmax

    Call Mpi_Allgather(boundaries,6,Mpi_Real, bound_tab, 6, Mpi_Real, MpiSubCubeWrite, mpierr)

    aname = 'boundaries_tab'
    Call hdf5_write_attr(gr_root_id, aname, 6, procperfile, bound_tab)

    boundaries(1) = minval(bound_tab(1,:))
    boundaries(2) = maxval(bound_tab(2,:))
    boundaries(3) = minval(bound_tab(3,:))
    boundaries(4) = maxval(bound_tab(4,:))
    boundaries(5) = minval(bound_tab(5,:))
    boundaries(6) = maxval(bound_tab(6,:))
    aname = 'boundaries'
    Call hdf5_write_attr(gr_root_id, aname, 6, boundaries)

    aname = 'potential'
    If(potential) Then
       potattr = 1
    Else
       potattr = 0
    End If
    Call hdf5_write_attr(gr_root_id, aname, potattr) 

    aname = 'force'
    If(force) Then
       forattr = 1
    Else
       forattr = 0
    End If
    Call hdf5_write_attr(gr_root_id, aname, forattr) 

    deb = 1
    fin = 1
    fic = nc*scprocIDWrite

    ! For each non empty group we create an HDF5 group and write dataset into it
    Do ic = 1, nc*procperfile

       Write(charic(1:8),'(I8.8)') ic
       groupname = 'group'//charic
       Call hdf5_create_group(gr_root_id, groupname, gr_id)

       fic = ic - nc*scprocIDWrite
       If(fic>=1 .and. fic<=nc) Then
          If(npartcube(fic) == 0) Then
             empty = .true.
             If(deb>mynpart) Then
                deb=1
                fin=1
             End If
          Else
             empty = .false.
             fin = deb + npartcube(fic) - 1
#ifdef DEBUG
             Print *,'group ',groupname,' opened:', procID, fic, deb, fin
#endif
          End If
       Else
          empty = .true.
       End If

       ! Write the position of the particles
       dsetname="pos"
       Call hdf5_write_mpi_data(gr_id, dsetname, 3, fin-deb+1, pos(:,deb:fin), &
            MpiSubCubeWrite, empty)
       
       ! Write the velocity of the particles
       dsetname="vel"
       Call hdf5_write_mpi_data(gr_id, dsetname, 3, fin-deb+1, vel(:,deb:fin),&
            MpiSubCubeWrite, empty) 
       
       ! Write the ID of the particles
       dsetname="ID"
       Call hdf5_write_mpi_data(gr_id, dsetname, fin-deb+1, id(deb:fin),&
            MpiSubCubeWrite,empty) 

       ! Write potential if it is used
       If(potential) Then
          dsetname = 'pot'
          Call hdf5_write_mpi_data(gr_id, dsetname, fin-deb+1, pot(deb:fin),&
               MpiSubCubeWrite,empty)
       End If

       ! Write force if it is used
       If(force) Then
          dsetname = 'for'
          Call hdf5_write_mpi_data(gr_id, dsetname, 3, fin-deb+1, for(:,deb:fin),&
               MpiSubCubeWrite,empty)
       End If
       
       Call hdf5_close_group(gr_id)

       If(fic>=1 .and. fic<=nc) Then          
          If(npartcube(fic) /= 0) Then
             deb = fin + 1
#ifdef DEBUG
             Print *,'group ',groupname,' closed'
#endif
          End If
       End If
          
    End Do

    ! Close the root group.
    Call hdf5_close_group(gr_root_id)

    Deallocate(partnb_tab, bound_tab)

    ! Close h5 file
    Call hdf5_close_mpi_file(file_id)

#ifdef DEBUG
    Print *,"Exit mpih5writesortedcube on process ",procID
#endif    

  End Subroutine mpih5writesortedcube


  !=======================================================================
  !!! A terminer
  Subroutine mpih5readsortedcube()
    
    Use modhdf5
    Use modparameters
    Use modmpicom
    Implicit none

    Character(len=400) :: filecube
    Character(len=5)  :: pid_char
    Character(len=8)  :: gid_char

    Character(len=16) :: dsetname                           ! Dataset name
    Character(len=16) :: aname                              ! Attribute name
    Character(len=16) :: groupname

    Integer(hid_t) :: file_id                               ! File identifier
    Integer(hid_t) :: gr_root_id                            ! Group identifier
    Integer(hid_t) :: gr_id                                 ! Group identifier
    
    Integer(kind=4) :: procperfile
    Integer(kind=4), dimension(:), allocatable :: partnb_tab
    Integer(kind=4) :: potattr
    Integer(kind=4) :: forattr
    Integer(kind=4) :: igroup, firstgroup, lastgroup
    Integer(kind=4) :: indbeg, indend
    Integer(kind=4) :: nc, ngroup
    
    Integer(kind=4), dimension(:), allocatable :: npartpergroup

    Real(kind=4), dimension(:,:), allocatable :: bound_tab

    
#ifdef DEBUG
    Print *,"Enter mpih5readsortedcube on process ",procID
#endif    
    
    ! number of processes reading in the same file
    procperfile = gatherread**3
    
    Allocate(partnb_tab(procperfile))
    Allocate(bound_tab(6,procperfile))
        
    Write(pid_char(1:5),'(I5.5)') commcolorRead
    filecube = trim(output_root)//'_mpisortedcube_'//pid_char//'.h5'

    ! we open the file for a serial read
    Call hdf5_open_file(filecube, file_id)

     ! open the root group
    groupname = '/'
    Call hdf5_open_group(file_id,groupname, gr_root_id)
    ! read attribute partNB
    aname = 'nres'
    Call hdf5_read_attr(gr_root_id, aname, nres)

#ifdef DEBUG
    Print *,"Value of nres read = ",nres, " on process ",procID
    Call Mpi_Barrier(Mpi_Comm_World, mpierr)
#endif

    ngrid = int(nres,kind=8)**3
    npart = ngrid
    
    ! read attribute partNB
    aname = 'partNB_tab'
    Call hdf5_read_attr(gr_root_id, aname, 6, partnb_tab)

#ifdef DEBUG
    Print *,"partNB_tab read on process ",procID
    Call Mpi_Barrier(Mpi_Comm_World, mpierr)
#endif    

    mynpart = partnb_tab(scprocIDRead+1)

    ! read attribute boundaries
    aname = 'boundaries_tab'
    Call hdf5_read_attr(gr_root_id, aname, 6, procperfile, bound_tab)

#ifdef DEBUG
    Print *,"boundaries_tab read on process ",procID
    Call Mpi_Barrier(Mpi_Comm_World, mpierr)
#endif    

    xmin=bound_tab(1,scprocIDRead+1)
    xmax=bound_tab(2,scprocIDRead+1)
    ymin=bound_tab(3,scprocIDRead+1)
    ymax=bound_tab(4,scprocIDRead+1)
    zmin=bound_tab(5,scprocIDRead+1)
    zmax=bound_tab(6,scprocIDRead+1)

    ! read attribute ngroup
    aname = 'ngroup'
    Call hdf5_read_attr(gr_root_id, aname, ngroup)

    Allocate(npartpergroup(ngroup))
    aname = 'npartpergroup'
    Call hdf5_read_data(gr_root_id, aname, ngroup, npartpergroup)

    nc = ngroup / procperfile
    firstgroup = nc * scprocIDRead + 1
    lastgroup  = nc *(scprocIDRead+1)

    potattr=0
    ! if requested, try to read the potential attribute
    If(potential) Then
       aname = 'potential'
       Call hdf5_read_attr(gr_root_id, aname, potattr)
       If(potattr==0) Then ! potential is asked by the user but not found in the cube file: abort
          If(procID==0) Then
             Print *,'The potential was not written in the cube files used for this analysis.'
             Print *,'You should change the potential parameter to false.'
          End If
          Call EmergencyStop('Potential missing in the cube files: abort.',100)
       End If
       If(.not.Allocated(pot)) Allocate(pot(mynpart))
    End If

    forattr=0
    ! if requested, try to read the force attribute
    If(force) Then
       aname = 'force'
       Call hdf5_read_attr(gr_root_id, aname, forattr)
       If(forattr==0) Then ! force is asked by the user but not found in the cube file: abort
          If(procID==0) Then
             Print *,'The force was not written in the cube files used for this analysis.'
             Print *,'You should change the force parameter to false.'
          End If
          Call EmergencyStop('Force missing in the cube files: abort.',100)
       End If
       If(.not.Allocated(for)) Allocate(for(3,mynpart))
    End If

    If(.not.Allocated(id)) Allocate(id(mynpart))
    If(.not.Allocated(pos))  Allocate(pos(3,mynpart))
    If(.not.Allocated(vel))  Allocate(vel(3,mynpart))



    indbeg = 1
    indend = 0
    ! loop over the groups
    Do igroup=firstgroup, lastgroup
       If(npartpergroup(igroup) /= 0) Then ! there is at least one part. to read in this group
          Write(gid_char(1:8),'(I8.8)') igroup
          groupname='group'//gid_char
          Call hdf5_open_group(gr_root_id,groupname,gr_id)
          
          indend = indbeg + npartpergroup(igroup) - 1
          ! read position of the particles
          dsetname = 'pos'
          Call hdf5_read_data(gr_id, dsetname, 3, mynpart, pos(:,indbeg:indend))
          
          ! read velocity of the particles
          dsetname = 'vel'
          Call hdf5_read_data(gr_id, dsetname, 3, mynpart, vel(:,indbeg:indend))
          
          ! read id of the particles
          dsetname = 'ID'
          Call hdf5_read_data(gr_id, dsetname, mynpart, id(indbeg:indend))
          
          ! read potential if requested
          If(potential) Then
             dsetname = 'pot'
             Call hdf5_read_data(gr_id, dsetname, mynpart, pot(indbeg:indend))
          End If

          ! read force if requested
          If(force) Then
             dsetname = 'for'
             Call hdf5_read_data(gr_id, dsetname, 3, mynpart, for(:,indbeg:indend))
          End If
          
          indbeg = indend + 1

          Call hdf5_close_group(gr_id)
       End If
    End Do


    ! Close the root group.
    Call hdf5_close_group(gr_root_id)

    Call hdf5_close_file(file_id)
    
    Deallocate(partnb_tab)
    Deallocate(bound_tab)


#ifdef DEBUG
    Print *,"Exit mpih5readsortedcube on process ",procID
#endif    
    
    
  End Subroutine mpih5readsortedcube


  !=======================================================================


End Module modio
