!> @file
!! This file contains subroutines used for I/O.

!> This module contains subroutines used for I/O.
Module modio

  Use modconst

  Integer(kind=4) :: mynpart       ! local particle number
  Integer(kind=4) :: ndim          ! number of dimensions
  Integer(kind=4) :: lmin          ! minimum mesh refinement level
  Integer(kind=4) :: nres          ! 1-D resolution: number of grid points in each dimension ; nres = 2 ** lmin
  Integer(kind=PRI) :: nptot       ! total particle number: nptot = nres ** 3
  Integer(kind=PRI) :: ngrid       ! total number of grid points: ngrid = nres ** 3
  Real   (kind=SP)  :: xmin, xmax, ymin, ymax, zmin, zmax  ! min and max (x,y,z) for each process

  Character(len=50) :: origin

Contains

  !=======================================================================

  Subroutine title()

    Use modmpicom, only : procNB

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

    Use modparam
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
              
#ifndef WITHHDF5
       If(usehdf5) Then
          Print *,'** pFoF was compiled without HDF5 support. ** '
          Print *,'** hdf5 parameter should be .false. or pFoF should be compiled with HDF5 support. **'
          Call Mpi_Abort(Mpi_Comm_World,ioerr,mpierr)
       End If
#endif
       
    End If
    
    Call bcastparam()
    
  End Subroutine readparameters


  !=======================================================================
  !> This subroutine writes the input parameters in a .nml file,
  !! print them on screen and writes them in the .log file.
  !! It should be called only from process 0.
  Subroutine writeparameters()
    
    Use modparam
    Use modmpicom
    Implicit none
    
    Character(len=84) :: fileopa 
    Character(len=84) :: filelog 
    
    If(procID == 0) Then
       ! Open file .nml and write input parameters
       fileopa = trim(root)//'.nml'
       Open(Unit=Uopa,file=fileopa)
       Write(Uopa, nml=ramsesinput)
       Write(Uopa, nml=fofparam)
       Write(Uopa, nml=outputparam)
       Close(Uopa)    
       
       ! Print input parameters on screen
       Print *,'Input parameters:'
       Print *, " "
       Print *, "Root:                                           ",trim(root)
       Print *, "Type of input files:                            ",code_index
       Print *, "Path to input files:                            ",trim(pathinput)
       Print *, "Particle files base name:                       ",trim(namepart)
       Print *, "Info file base name:                            ",trim(nameinfo)
       Print *, "Size of groups of inputs:                       ",grpsize
       Print *, "Percolation parameter:                          ",perco
       Print *, "Minimum mass of halo to be analyzed:            ",Mmin
       Print *, "Maximum mass of halo to be analyzed:            ",Mmax
       Print *, "Were stars written in RAMSES files:             ",star
       Print *, "Were metallicities written in RAMSES files:     ",metal
       Print *, "Write cubes of particles:                       ",outcube
       Print *, "Perform friends of friends halo detection:      ",dofof
       Print *, "Read particles from cube files:                 ",readfromcube
       Print *, "Perform timings (imply extra synchronisations): ",dotimings
       Print *, 'Use hdf5 i/o file format:                       ', usehdf5 
       Print *, " "
       
       
       ! Open log file
       filelog = trim(root)//'.log'
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

    Print *,'Friend of Friend termine'
    Print *,''
    Print *,'Temps d''execution:'
    Print *,'Lecture:',tReadRA,' dont'
    Print *,'        initialisation        :',tInitRead
    Print *,'        lecture des fichiers  :',tReadFile
    Print *,'        partage des particules:',tTailPart
    Print *,''
    Print *,'Friend of Friend:',tFoF,'dont'
    Print *,'        initialisatin:',tFoFinit
    Print *,'        FoF local    :',tFoFloc
    Print *,'        raccordement :',tRaccord
    Print *,'        calcul d''observables:',tObs
    Print *,'        sorties:',tOut
    
    Write(Ulog,*) 'Friend of Friend termine'
    Write(Ulog,*) ''
    Write(Ulog,*) 'Temps d''execution:'
    Write(Ulog,*) 'Lecture:',tReadRA,' dont'
    Write(Ulog,*) '        initialisation        :',tInitRead
    Write(Ulog,*) '        lecture des fichiers  :',tReadFile
    Write(Ulog,*) '        partage des particules:',tTailPart
    Write(Ulog,*) ''
    Write(Ulog,*) 'Friend of Friend:',tFoF,'dont'
    Write(Ulog,*) '        initialisatin:',tFoFinit
    Write(Ulog,*) '        FoF local    :',tFoFloc
    Write(Ulog,*) '        raccordement :',tRaccord
    Write(Ulog,*) '        calcul d''observables:',tObs
    Write(Ulog,*) '        sorties:',tOut

    ! Close log file
    Close(Ulog)
    
  End Subroutine writetimings


  !=======================================================================
  !> This subroutine reads the particles files created by RAMSES that pFOF has to analyze.
  Subroutine ramses_lecture()
    Use modvariable
    Use modparam
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

    Real(kind=SP), allocatable     :: tmpsimple(:)            ! temporary variable for Ramses v2 output
    Real(kind=DP), allocatable     :: tmpdouble(:)            ! temporary variable for Ramses v3 output
    Real(kind=SP), allocatable     :: tmpsendx(:,:),tmpsendv(:,:)
    Real(kind=SP), allocatable     :: tmpx(:,:), tmpv(:,:)
    Real(kind=SP)                  :: deltasd

    Integer(kind=4) :: buffersize, nbbytes, b_pos
    Character, dimension(:),allocatable :: buffer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(kind=4) :: tmpinteger


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

    nptot = 0

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
    Call Mpi_AllReduce(tmplongint,nptot,1,MPI_PRI,Mpi_Sum,Mpi_Comm_World,mpierr)

    If(procID == 0) Then
       Write(* ,*)'There are ',nptot,' DM particles'
       Write(Ulog,*)'There are ',nptot,' DM particles'
    End If

    Allocate(tmpx(3,mynpart))
    Allocate(tmpv(3,mynpart))
    Allocate(tmpi(mynpart))

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
             tmpx(idim,mynpart+1:mynpart+npartloc) = real(tmpdouble, kind=SP)
          End Do
          
          ! Read velocities in a dummy variable
          Do idim = 1,ndim
             Read(1) tmpdouble
             ! put all velocities in tmpv vector
             tmpv(idim,mynpart+1:mynpart+npartloc) = real(tmpdouble,kind=SP)
          End Do
          
          ! Read masses in a dummy variable
          Read(1) tmpdouble

          If(star) Then
             Read(1) tmpdouble
             If(metal) Read(1) tmpdouble
          End If

          Deallocate(tmpdouble)

       End If ramsesV3

       ! Read particle id
       Read(1) tmpi(mynpart+1:mynpart+npartloc)

!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!     AJOUT TEMPORAIRE POUR TEST
!!$       ! Read potential if potential parameter is .true.
!!$       ! First: skip the level (integer 4)
!!$       Read(1) tmpinteger
!!$       ! Then read the potential: same kind as pos and vel
!!$       Allocate(tmpsimple(1:npartloc))
!!$       Read(1) tmpsimple
!!$       tmpp(mynpart+1:mynpart+npartloc) = tmpsimple
!!$
!!$       ! Then read the force: same kind as pos and vel
!!$       Do idim = 1,ndim
!!$          Read(1) tmpsimple
!!$          tmpf(idim,mynpart+1:mynpart+npartloc) = tmpsimple
!!$       End Do
!!$       Deallocate(tmpsimple)
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       
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

    Allocate (x(3,npartv(procID+1)))
    Allocate (v(3,npartv(procID+1)))
    Allocate (id(npartv(procID+1)))

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
          ind = 1

          Do j=1,mynpart
             If(tmpx(1,j)>=xmin .and. tmpx(1,j)<xmax .and. &
                  tmpx(2,j)>= ymin .and. tmpx(2,j) < ymax .and. &
                  tmpx(3,j)>= zmin .and. tmpx(3,j) < zmax) Then
                tmpsendx(:,ind) = tmpx(:,j)
                tmpsendv(:,ind) = tmpv(:,j)
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
          Call Mpi_Isend(tmpsendx,3*npartvloc(dest+1),Mpi_Real,dest,procID,MpiCube,mpireqs1,mpierr)
          Call Mpi_Isend(tmpsendv,3*npartvloc(dest+1),Mpi_Real,dest,  dest,MpiCube,mpireqs2,mpierr)
          Call Mpi_Isend(tmpsendi,  npartvloc(dest+1), MPI_PRI,dest,procID,MpiCube,mpireqs3,mpierr)
       End If
       If(nrecv/=0) Then
          Call Mpi_Irecv(x(1,recvpoint),3*nrecv,Mpi_Real,prov,  prov,MpiCube,mpireqr1,mpierr)
          Call Mpi_Irecv(v(1,recvpoint),3*nrecv,Mpi_Real,prov,procID,MpiCube,mpireqr2,mpierr)
          Call Mpi_Irecv( id(recvpoint),  nrecv, MPI_PRI,prov,  prov,MpiCube,mpireqr3,mpierr)
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
          x(:,recvpoint+ind) = tmpx(:,j)
          v(:,recvpoint+ind) = tmpv(:,j)
          id(recvpoint+ind)  = tmpi(j)
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
    tTailPart = Mpi_Wtime() - timeInt
    tReadRA = Mpi_Wtime() - time0

  End Subroutine ramses_lecture


  !=======================================================================
  !> This subroutine writes binary cube files.
  Subroutine outputcube()

    Use modparam
    Use modvariable
    Use modmpicom
    Implicit none

    Integer(kind=PRI) :: i, j
    Character(len=400) :: filecube
    Character(len=5)  :: pid_char
    
    Write(pid_char(1:5),'(I5.5)') procID
    filecube = trim(root)//"_cube_"//pid_char
    Open(Unit=Ucub,file=filecube,Form='Unformatted')

    Write(Ucub) int(mynpart,kind=4)
    Write(Ucub) procID
    Write(Ucub) xmin,xmax,ymin,ymax,zmin,zmax
    Write(Ucub) ((x(j,i),j=1,3),i=1,mynpart)
    Write(Ucub) ((v(j,i),j=1,3),i=1,mynpart)
    Write(Ucub) (id(i),i=1,mynpart)
    
    Close(Ucub)

  End Subroutine outputcube

  !=======================================================================
  !> This subroutine writes binary files containing the mass, the id and the position of the center
  !! of mass for each halo found by pFOF. One file is written per process.
  Subroutine writehalomass(haloNB, haloMass, halocomPos, haloID)

    Use modparam
    Use modvariable
    Use modmpicom
    Implicit none

    
    Integer(kind=4),                      intent(in) :: haloNB
    Integer(kind=4),  dimension(haloNB),  intent(in) :: haloMass
    Real(kind=DP),    dimension(3,haloNB),intent(in) :: halocomPos
    Integer(kind=PRI), dimension(haloNB), intent(in) :: haloID

    Character(len=400) :: fileamas
    Character(len=5)  :: pid_char
    Integer(kind=4) :: i

    
    Write(pid_char(1:5),'(I5.5)') procID
    fileamas = trim(root)//"_masst_"//pid_char

    Open(Umas, File=trim(fileamas), Status='Unknown', Form='Unformatted')
    Write(Umas) haloNB
    Do i=1,haloNB
       Write(Umas) haloID(i), haloMass(i), halocomPos(:,i)
    End Do
    Close(Umas)
    
  End Subroutine writehalomass

  !=======================================================================
  !> This subroutine writes binary files containing the position, the velocity and the id of each particle 
  !! in each halo found by pFOF. One file is written per process.
  Subroutine writehalopart(haloNB, halopartNB, haloMass, halopartPos, halopartVel, halopartID)

    Use modparam
    Use modvariable
    Use modmpicom
    Implicit none

    Integer(kind=4),                          intent(in) :: haloNB, halopartNB
    Integer(kind=PRI),dimension(halopartNB),  intent(in) :: halopartID
    Integer(kind=4),  dimension(haloNB),      intent(in) :: haloMass
    Real   (kind=SP), dimension(3,halopartNB),intent(in) :: halopartPos, halopartVel

    Integer(kind=4) :: h, j, k, b
    Character(len=400) :: filestrct
    Character(len=5)  :: pid_char
    
    Write(pid_char(1:5),'(I5.5)') procID
    filestrct = trim(root)//'_strct_'//pid_char
    
    Open(Ustr, File=trim(filestrct), Status='Unknown',Form='Unformatted')

    Write(Ustr) haloNB

    b = 1
    Do h = 1, haloNB
       Write(Ustr) haloMass(h)
       Write(Ustr) ((halopartPos(k,j),k=1,3),j=b,b+haloMass(h)-1)
       Write(Ustr) ((halopartVel(k,j),k=1,3),j=b,b+haloMass(h)-1)
       Write(Ustr) (halopartID(j),j=b,b+haloMass(h)-1)
       b = b + haloMass(h)
    End Do
       
    Close(Ustr)

  End Subroutine writehalopart

  !=======================================================================
  
#ifdef WITHHDF5
  !> This subroutine writes the position, velocity and id of each particle on the process in a hdf5 file.
  !! One file is written per MPI process
  Subroutine h5writecube()

    Use modhdf5
    Use modparam
    Use modvariable
    Use modmpicom
    Use modxdmf
    Implicit none

    Character(len=400) :: filecube
    Character(len=391) :: filebase
    Character(len=5)  :: pid_char

    Integer(kind=4), dimension(2) :: nbelem
    Character(len=16) :: dsetname                           ! Dataset name
    Character(len=16) :: aname                              ! Attribute name
    Character(len=20) :: adata
    Character(len=16) :: groupname

    Integer(hid_t) :: file_id                               ! File identifier
    Integer(hid_t) :: gr_root_id                            ! Group identifier

    Real(kind=SP), dimension(6) :: boundaries
    Integer(kind=4), dimension(procNB) :: nparttab

#ifdef DEBUG
    Print *,"Enter h5writecube on process ",procID
#endif    
    
    Write(pid_char(1:5),'(I5.5)') procID
    filebase = trim(root)//'_cube'
    filecube = trim(root)//'_cube_'//pid_char//'.h5'
    
    Call Mpi_Gather(mynpart, 1, Mpi_Integer, nparttab, 1, Mpi_Integer, 0, MpiCube, mpierr)
    If(procID==0) Then
       Call writesinglecubexdmf(procNB,filebase,nparttab)
    End If

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

    ! Write the boundaries of the cube as an attribute
    boundaries(1) = xmin
    boundaries(2) = xmax
    boundaries(3) = ymin
    boundaries(4) = ymax
    boundaries(5) = zmin
    boundaries(6) = zmax
    aname = 'boundaries'

    Call hdf5_write_attr(gr_root_id, aname, 6, boundaries)

    nbelem(1) = 3
    nbelem(2) = mynpart
    ! Write the position of the particles
    dsetname="pos"
    Call hdf5_write_data(gr_root_id, dsetname, 3, mynpart, x)

    ! Write the velocity of the particles
    dsetname="vel"
    Call hdf5_write_data(gr_root_id, dsetname, 3, mynpart, v) 

    ! Write the ID of the particles
    dsetname="ID"
    Call hdf5_write_data(gr_root_id, dsetname, mynpart, id) 

    ! Close the root group.
    Call hdf5_close_group(gr_root_id)

    Call hdf5_close_file(file_id)

#ifdef DEBUG
    Print *,"Exit h5writecube on process ",procID
#endif    

  End Subroutine h5writecube


  !=======================================================================
  !> This subroutine writes for each halo the position, the velocity and the id of each particle in this halo
  !! in a hdf5 file.
  !! One file is written per MPI process.
  !! There is one group per halo, the name of the group is the ID of the halo.
  Subroutine h5writehalopart(haloNB, halopartNB, haloMass, haloID, halopartPos, halopartVel, halopartID)

    Use modhdf5
    Use modparam
    Use modmpicom
    Use modvariable
    Use modxdmf
    Implicit none
    
    Integer(kind=4),                            intent(in) :: haloNB, halopartNB
    Integer(kind=PRI), dimension(haloNB),       intent(in) :: haloID
    Integer(kind=PRI), dimension(halopartNB),   intent(in), target :: halopartID
    Integer(kind=4),   dimension(haloNB),       intent(in) :: haloMass
    Real   (kind=SP),  dimension(3,halopartNB), intent(in), target :: halopartPos, halopartVel

    Integer(kind=4) :: ih, hptr
    Character(len=400) :: filestrct
    Character(len=391) :: filebase
    Character(len=5)  :: pid_char
    Character(len=16) :: groupname
    Character(len=16) :: dsetname                           ! Dataset name
    Character(len=16) :: aname                              ! Attribute name

    Integer(hid_t) :: file_id                               ! File identifier
    Integer(hid_t) :: gr_root_id                            ! Group identifier
    Integer(hid_t) :: gr_halo_id                            ! Group identifier

    Integer(kind=4), dimension(procNB) :: haloNBtab, dspl
    Integer(kind=4), dimension(:), allocatable :: haloMasstab
    Integer(kind=PRI), dimension(:), allocatable :: haloIDtab
    Integer(kind=PRI), dimension(2) :: IDminmax
    Integer(kind=4) :: haloNBall, p

    
#ifdef DEBUG
    Print *,"Enter h5writehalopart on process ",procID
#endif    

    Write(pid_char(1:5),'(I5.5)') procID
    filestrct = trim(root)//'_halo_'//pid_char//'.h5'
    filebase = trim(root)//'_halo'

    Call Mpi_Gather(haloNB, 1, Mpi_Integer, haloNBtab, 1, Mpi_Integer, 0, MpiCube, mpierr)
    If(procID==0) Then
       haloNBall = haloNBtab(1)
       dspl(1) = 0
       Do p=2, procNB 
          dspl(p) = haloNBall
          haloNBall= haloNBall + haloNBtab(p)
       End Do
       Allocate(haloMasstab(haloNBall), haloIDtab(haloNBall))
    End If
    Call Mpi_Gatherv(haloMass, haloNB, Mpi_Integer, haloMasstab, haloNBtab, dspl, Mpi_Integer, 0, MpiCube, mpierr)
    Call Mpi_Gatherv(haloID, haloNB, Mpi_PRI, haloIDtab, haloNBtab, dspl, Mpi_PRI, 0, MpiCube, mpierr)

    If(procID==0) Then
       Call writesinglehaloxdmf(procNB,filebase,haloNBtab,haloMasstab,haloIDtab)
    End If

    ! create the hdf5 file
    Call hdf5_create_file(filestrct, file_id, origin)

    ! open the root group
    groupname = '/'
    Call hdf5_open_group(file_id,groupname, gr_root_id)

    aname = 'nfile'
    Call hdf5_write_attr(gr_root_id, aname, procNB)
    
    ! write the number of haloes as an attribute
    aname = "haloNB"
    Call hdf5_write_attr(gr_root_id, aname, haloNB)

    !! write the halo ID as data and not attribute: it seems that we cannot write integer(kind=8) attribute 
    aname = 'haloID'
    Call hdf5_write_data(gr_root_id, aname, haloNB, haloID)

    aname = 'haloIDminmax'
    IDminmax(1) = haloID(1)
    IDminmax(2) = haloID(haloNB)
    Call hdf5_write_data(gr_root_id, aname, 2, IDminmax)

    ! pointer to the current halo
    hptr = 1
    
    Do ih = 1, haloNB
       ! create a group for each halo
       groupname = "halo_00000000000"
       Write(groupname(6:16),'(I11.11)') haloID(ih)
       Call hdf5_create_group(gr_root_id, groupname, gr_halo_id)
       ! create an attribute containing the number of particles in the halo
       aname = "halopartNB"
       Call hdf5_write_attr(gr_halo_id, aname, haloMass(ih))
       
       dsetname="pos"
       Call hdf5_write_data(gr_halo_id, dsetname, 3, haloMass(ih), halopartPos(:,hptr:hptr+haloMass(ih)-1)) 
       
       dsetname="vel"
       Call hdf5_write_data(gr_halo_id, dsetname, 3, haloMass(ih), halopartVel(:,hptr:hptr+haloMass(ih)-1)) 
       
       dsetname = "ID"
       Call hdf5_write_data(gr_halo_id, dsetname, haloMass(ih), halopartID(hptr:hptr+haloMass(ih)-1)) 
       
       ! Close the halo group.
       Call hdf5_close_group(gr_halo_id)
       
       ! move the pointer to the next halo
       hptr = hptr + haloMass(ih)
    End Do
    
    ! Close the root group.
    Call hdf5_close_group(gr_root_id)
       
    Call hdf5_close_file(file_id)

#ifdef DEBUG
    Print *,"Exit h5writehalopart on process ",procID
#endif

  End Subroutine h5writehalopart


  !=======================================================================
  !> This subroutine writes, for each halo, its mass (number of particles), 
  !! the position of its center of mass and its ID
  !! in only one hdf5 file using parallel HDF5.
  Subroutine mpih5writehalomass(haloNB, haloMass, halocomPos, halocomVel, haloID)

    Use modparam
    Use modvariable
    Use modmpicom
    Use modhdf5

    Implicit none
    
    Integer(kind=4),                        intent(in) :: haloNB
    Integer(kind=4),   dimension(haloNB),   intent(in), target :: haloMass
    Real(kind=DP),     dimension(3,haloNB), intent(in), target :: halocomPos
    Real(kind=DP),     dimension(3,haloNB), intent(in), target :: halocomVel
    Integer(kind=PRI), dimension(haloNB),   intent(in), target :: haloID

    Integer(kind=4) :: haloNB_all
    Character(len=400) :: filename
    Character(len=16) :: dsetname                        ! Dataset name
    Integer(hid_t) :: file_id                            ! File identifier
    

    filename = trim(root)//'_halomass.h5'

    ! Create h5 parallel file
    Call hdf5_create_mpi_file(filename, MpiCube, file_id, origin)

    Call Mpi_Allreduce(haloNB, haloNB_all, 1, Mpi_Integer, Mpi_Sum, Mpi_Comm_World, mpierr)

    dsetname = 'haloNB'
    Call hdf5_write_attr(file_id, dsetname, haloNB_all)
    
    dsetname = 'halocomPos'
    Call hdf5_write_mpi_data(file_id, dsetname, 3, haloNB, halocomPos, MpiCube)

    dsetname = 'halocomVel'
    Call hdf5_write_mpi_data(file_id, dsetname, 3, haloNB, halocomVel, MpiCube)

    dsetname = 'haloID'
    Call hdf5_write_mpi_data(file_id, dsetname, haloNB, haloID, MpiCube)

    dsetname = 'haloMass'
    Call hdf5_write_mpi_data(file_id, dsetname, haloNB, haloMass, MpiCube)

    ! Close h5 file
    Call hdf5_close_mpi_file(file_id)

    
  End Subroutine mpih5writehalomass


  !=======================================================================
  !> This subroutine reads hdf5 cube files created by pFOF.
  Subroutine h5readcube()

    Use modhdf5
    Use modparam
    Use modmpicom
    Use modvariable


    Character(len=400) :: filename      ! File name
    Integer(hid_t) :: file_id                               ! File identifier
    Integer(hid_t) :: gr_root_id                            ! Group identifier
    Character(len=5)  :: pid_char
    Character(len=16) :: aname                              ! Attribute name
    Character(len=16) :: dsetname                           ! Dataset name
    Character(len=16) :: groupname
    Real(kind=SP), dimension(6) :: boundaries


#ifdef DEBUG
    Print *,"Enter h5readcube on process ",procID
#endif

    Write(pid_char(1:5),'(I5.5)') procID
    filename = trim(root)//'_cube_'//pid_char//'.h5'

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
    nptot = ngrid

    ! read attribute boundaries
    aname = 'boundaries'
    Call hdf5_read_attr(gr_root_id, aname, 6, boundaries)

    xmin=boundaries(1)
    xmax=boundaries(2)
    ymin=boundaries(3)
    ymax=boundaries(4)
    zmin=boundaries(5)
    zmax=boundaries(6)

    If(.not.Allocated(id)) Allocate(id(mynpart))
    If(.not.Allocated(x))  Allocate(x(3,mynpart))
    If(.not.Allocated(v))  Allocate(v(3,mynpart))

    ! read position of the particles
    dsetname = 'pos'
    Call hdf5_read_data(gr_root_id, dsetname, 3, mynpart, x)

    ! read velocity of the particles
    dsetname = 'vel'
    Call hdf5_read_data(gr_root_id, dsetname, 3, mynpart, v)

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
  !> This subroutine writes the position, the velocity and the id of each particle 
  !! on the process in a hdf5 file. 
  !! The particles are gathered 
  Subroutine mpih5writecube()

    Use modhdf5
    Use modparam
    Use modvariable
    Use modmpicom
    Implicit none

    Character(len=400) :: filecube
    Character(len=5)  :: pid_char

    Character(len=16) :: dsetname                           ! Dataset name
    Character(len=16) :: aname                              ! Attribute name
    Character(len=16) :: groupname

    Integer(hid_t) :: file_id                               ! File identifier
    Integer(hid_t) :: gr_root_id                            ! Group identifier

    Integer :: procperfile
    Integer :: npart
    Integer :: nfile
    Integer, dimension(:), allocatable :: partnb_tab

    Real(kind=SP), dimension(6) :: boundaries
    Real(kind=SP), dimension(:,:), allocatable :: bound_tab

#ifdef DEBUG
    Print *,"Enter mpih5writecube on process ",procID
#endif    

    ! number of processes writing in the same file
    procperfile = gatherwrite**3
        
    Allocate(partnb_tab(procperfile))
    Allocate(bound_tab(6,procperfile))

    Write(pid_char(1:5),'(I5.5)') commcolorWrite
    filecube = trim(root)//'_mpicube_'//pid_char//'.h5'

    ! create the hdf5 file
    Call hdf5_create_mpi_file(filecube, MpiSubCubeWrite, file_id, origin)

    ! open the root group
    groupname = '/'
    Call hdf5_open_group(file_id,groupname, gr_root_id)

    ! Offset and global dimensions have to be computed
    Call Mpi_Allgather(mynpart, 1, Mpi_Integer, partnb_tab, 1, Mpi_Integer, MpiSubCubeWrite, mpierr)

    aname='partNB_tab'
    Call hdf5_write_attr(gr_root_id, aname, procperfile, partnb_tab)

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

    dsetname = 'pos'
    Call hdf5_write_mpi_data(gr_root_id, dsetname, 3, mynpart, x, MpiSubCubeWrite)

    dsetname = 'vel'
    Call hdf5_write_mpi_data(gr_root_id, dsetname, 3, mynpart, v, MpiSubCubeWrite)

    dsetname = 'ID'
    Call hdf5_write_mpi_data(gr_root_id, dsetname, mynpart, id, MpiSubCubeWrite)

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

  ! Write for each halo the position, velocity and id of each particle in this halo
  ! using hdf5 format.
  ! One file is written per communicator, several processes write in the same file
  ! There is one group per halo, the name of the group is the ID of the halo
!!$  Subroutine mpih5writehalopart(haloNB, halopartNB, haloMass, haloID, halopartPos, halopartVel, halopartID)
!!$
!!$    Use modhdf5
!!$    Use modparam
!!$    Use modmpicom
!!$    Use modvariable
!!$    Use modxdmf
!!$    Implicit none
!!$    
!!$    Integer(kind=4),                            intent(in) :: haloNB, halopartNB
!!$    Integer(kind=PRI), dimension(haloNB),       intent(in) :: haloID
!!$    Integer(kind=PRI), dimension(halopartNB),   intent(in), target :: halopartID
!!$    Integer(kind=4),   dimension(haloNB),       intent(in) :: haloMass
!!$    Real   (kind=SP),  dimension(3,halopartNB), intent(in), target :: halopartPos, halopartVel
!!$
!!$    Integer(kind=4), dimension(2) :: nbelem
!!$    Integer(kind=4) :: ih, hptr
!!$    Character(len=85) :: filestrct
!!$    Character(len=76) :: filebase
!!$    Character(len=5)  :: pid_char
!!$    Character(len=16) :: groupname
!!$    Character(len=16) :: dsetname                           ! Dataset name
!!$    Character(len=16) :: aname                              ! Attribute name
!!$
!!$    Integer(hid_t) :: file_id                               ! File identifier
!!$    Integer(hid_t) :: gr_root_id                            ! Group identifier
!!$    Integer(hid_t) :: gr_halo_id                            ! Group identifier
!!$
!!$    Integer :: h5err                                        ! Error flag
!!$    Integer(kind=4), dimension(procNB) :: haloNBtab, dspl
!!$    Integer(kind=4), dimension(:), allocatable :: haloMasstab
!!$    Integer(kind=PRI), dimension(:), allocatable :: haloIDtab
!!$    Integer(kind=4) :: haloNBall, p
!!$
!!$    
!!$#ifdef DEBUG
!!$    Print *,"Enter h5writehalopart on process ",procID
!!$#endif    
!!$
!!$    Write(pid_char(1:5),'(I5.5)') procID
!!$    filestrct = trim(root)//'_halo_'//pid_char//'.h5'
!!$    filebase = trim(root)//'_halo'
!!$
!!$    Call Mpi_Gather(haloNB, 1, Mpi_Integer, haloNBtab, 1, Mpi_Integer, 0, MpiCube, mpierr)
!!$    If(procID==0) Then
!!$       haloNBall = haloNBtab(1)
!!$       dspl(1) = 0
!!$       Do p=2, procNB 
!!$          dspl(p) = haloNBall
!!$          haloNBall= haloNBall + haloNBtab(p)
!!$       End Do
!!$       Allocate(haloMasstab(haloNBall), haloIDtab(haloNBall))
!!$    End If
!!$    Call Mpi_Gatherv(haloMass, haloNB, Mpi_Integer, haloMasstab, haloNBtab, dspl, Mpi_Integer, 0, MpiCube, mpierr)
!!$    Call Mpi_Gatherv(haloID, haloNB, Mpi_PRI, haloIDtab, haloNBtab, dspl, Mpi_PRI, 0, MpiCube, mpierr)
!!$
!!$    If(procID==0) Then
!!$       Call writesinglehaloxdmf(procNB,filebase,haloNBtab,haloMasstab,haloIDtab)
!!$    End If
!!$
!!$    ! create the hdf5 file
!!$    Call hdf5_create_file(filestrct, file_id)
!!$
!!$    ! open the root group
!!$    Call h5gopen_f(file_id, "/", gr_root_id, h5err)
!!$        
!!$    aname = "haloNB"
!!$    Call hdf5_write_attr(gr_root_id, aname, haloNB)
!!$
!!$    
!!$    ! pointer to the current halo
!!$    hptr = 1
!!$    
!!$    Do ih = 1, haloNB
!!$       ! create a group for each halo
!!$       groupname = "halo_00000000000"
!!$       Write(groupname(6:16),'(I11.11)') haloID(ih)
!!$       Call h5gcreate_f(gr_root_id, groupname, gr_halo_id, h5err)
!!$       
!!$       ! create an attribute containing the number of particles in the halo
!!$       aname = "halopartNB"
!!$       Call hdf5_write_attr(gr_halo_id, aname, haloMass(ih))
!!$       
!!$       nbelem(1) = 3
!!$       nbelem(2) = haloMass(ih)
!!$       
!!$       dsetname="pos"
!!$       Call hdf5_write_real2D(gr_halo_id, C_LOC(halopartPos(1,hptr)), nbelem, 4, dsetname) 
!!$       
!!$       dsetname="vel"
!!$       Call hdf5_write_real2D(gr_halo_id, C_LOC(halopartVel(1,hptr)), nbelem, 4, dsetname) 
!!$       
!!$       dsetname = "ID"
!!$       Call hdf5_write_int1D(gr_halo_id, C_LOC(halopartID(hptr)), haloMass(ih), 8, dsetname)
!!$       
!!$       ! Close the halo group.
!!$       Call h5gclose_f(gr_halo_id, h5err)
!!$       
!!$       ! move the pointer to the next halo
!!$       hptr = hptr + haloMass(ih)
!!$    End Do
!!$    
!!$    ! Close the root group.
!!$    Call h5gclose_f(gr_root_id, h5err)
!!$       
!!$    Call hdf5_close_file(file_id)
!!$
!!$#ifdef DEBUG
!!$    Print *,"Exit h5writehalopart on process ",procID
!!$#endif
!!$
!!$  End Subroutine mpih5writehalopart




  !=======================================================================
 
  Subroutine mpih5readcube()
    
    Use modhdf5
    Use modparam
    Use modvariable
    Use modmpicom
    Implicit none

    Character(len=400) :: filecube
    Character(len=5)  :: pid_char

    Character(len=16) :: dsetname                           ! Dataset name
    Character(len=16) :: aname                              ! Attribute name
    Character(len=16) :: groupname

    Integer(hid_t) :: file_id                               ! File identifier
    Integer(hid_t) :: gr_root_id                            ! Group identifier
    
    Integer :: procperfile
    Integer, dimension(:), allocatable :: partnb_tab

    Real(kind=SP), dimension(:,:), allocatable :: bound_tab

    
#ifdef DEBUG
    Print *,"Enter mpih5readcube on process ",procID
#endif    
    
    ! number of processes writing in the same file
    procperfile = gatherread**3
    
    Allocate(partnb_tab(procperfile))
    Allocate(bound_tab(6,procperfile))
        
    Write(pid_char(1:5),'(I5.5)') commcolorRead
    filecube = trim(root)//'_mpicube_'//pid_char//'.h5'

    Call hdf5_open_mpi_file(filecube, MpiSubCubeRead, file_id)

     ! open the root group
    groupname = '/'
    Call hdf5_open_group(file_id,groupname, gr_root_id)
    ! read attribute partNB
    aname = 'nres'
    Call hdf5_read_attr(gr_root_id, aname, nres)

    ngrid = int(nres,kind=8)**3
    nptot = ngrid
    
    ! read attribute partNB
    aname = 'partNB_tab'
    Call hdf5_read_attr(gr_root_id, aname, 6, partnb_tab)

    mynpart = partnb_tab(scprocIDRead+1)

    ! read attribute boundaries
    aname = 'boundaries_tab'
    Call hdf5_read_attr(gr_root_id, aname, 6, procperfile, bound_tab)

    xmin=bound_tab(1,scprocIDRead+1)
    xmax=bound_tab(2,scprocIDRead+1)
    ymin=bound_tab(3,scprocIDRead+1)
    ymax=bound_tab(4,scprocIDRead+1)
    zmin=bound_tab(5,scprocIDRead+1)
    zmax=bound_tab(6,scprocIDRead+1)

    If(.not.Allocated(id)) Allocate(id(mynpart))
    If(.not.Allocated(x))  Allocate(x(3,mynpart))
    If(.not.Allocated(v))  Allocate(v(3,mynpart))
    
    dsetname = 'pos'
    Call hdf5_read_mpi_data(gr_root_id, dsetname, 3, mynpart, x, MpiSubCubeRead)

    dsetname = 'vel'
    Call hdf5_read_mpi_data(gr_root_id, dsetname, 3, mynpart, v, MpiSubCubeRead)

    dsetname = 'ID'
    Call hdf5_read_mpi_data(gr_root_id, dsetname, mynpart, id, MpiSubCubeRead)
  
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



#endif


End Module modio
