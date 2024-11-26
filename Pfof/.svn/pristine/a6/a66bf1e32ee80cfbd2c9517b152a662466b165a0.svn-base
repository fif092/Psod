Module modio

  Use modconst

  Integer(kind=4) :: mynpart       ! local particle number
  Integer(kind=4) :: ndim          ! number of dimensions
  Integer(kind=4) :: lmin          ! minimum mesh refinement level
  Integer(kind=4) :: nres          ! 1-D resolution: number of grid points in each dimension ; nres = 2 ** lmin
  Integer(kind=PRI) :: nptot       ! total particle number: nptot = nres ** 3
  Integer(kind=PRI) :: ngrid       ! total number of grid points: ngrid = nres ** 3
  Real   (kind=SP)  :: xmin, xmax, ymin, ymax, zmin, zmax  ! min and max (x,y,z) for each process

Contains


  ! Write input parameters in file .inp
  Subroutine outputparameters()
    Use modparam
    Use modmpicom
    Implicit none

    Character(len=84) :: fileopa 

    ! Open file .inp and write input parameters
    fileopa = trim(root)//'.inp'
    Open(Unit=Uopa,file=fileopa)
    Write(Uopa,*) 'Nb of processes:'
    Write(Uopa,*) procNB
    Write(Uopa,*) ''
    Write(Uopa,*) 'Input parameters:'
    Write(Uopa,*) root
    Write(Uopa,*) code_index
    Write(Uopa,*) pathinput
    Write(Uopa,*) namepart
    Write(Uopa,*) nameinfo
    Write(Uopa,*) grpsize
    Write(Uopa,*) perco
    Write(Uopa,*) Mmin
    Write(Uopa,*) Mmax
    Write(Uopa,*) star
    Write(Uopa,*) metal
    Write(Uopa,*) outcube
    Write(Uopa,*) dofof
    Write(Uopa,*) readfromcube
    Write(Uopa,*) dotimings
     
    Close(Uopa)    

  End Subroutine outputparameters


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
    Character(len=100)             :: nomfich
    Character(len=13)              :: dumchar
    Character(len=11)              :: grpchar

    Integer(kind=4)                :: i, j, icpu, idim   ! loop variables
    Integer(kind=4)                :: destCoord(3)       ! coords of the destination MPI process in MPI process cart
    Integer(kind=4)                :: nrecv              ! number of elements received in a Mpi_Recv
    Integer(kind=4)                :: recvpoint          ! address of the 1st element received in the local vector
    Integer(kind=4)                :: mynproc            ! number of RAMSES part files read by local process
    Integer(kind=4)                :: firstp, lastp      ! id of 1st and last RAMSES part file read
    Integer(kind=4), allocatable :: npartvloc(:), npartv(:)  ! temp and global table of particle numbers for each process
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

    ! tampon pour Mpi_Pack
    h_length = 3* bit_size(nres)/8
    Allocate (header(0:h_length-1))
    h_pos = 0

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

       Call Mpi_Pack(nproc, 1, Mpi_Integer, header, h_length, h_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack( ndim, 1, Mpi_Integer, header, h_length, h_pos, Mpi_Comm_World, mpierr)
       Call Mpi_Pack( nres, 1, Mpi_Integer, header, h_length, h_pos, Mpi_Comm_World, mpierr)

    End If

    Call Mpi_Bcast(header,h_length,Mpi_Packed,0,Mpi_Comm_World,mpierr)

    If(procID /= 0) Then

       Call Mpi_Unpack(header, h_length, h_pos, nproc, 1, Mpi_Integer, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(header, h_length, h_pos,  ndim, 1, Mpi_Integer, Mpi_Comm_World, mpierr)
       Call Mpi_Unpack(header, h_length, h_pos,  nres, 1, Mpi_Integer, Mpi_Comm_World, mpierr)

    End If

    ngrid = int(nres,kind=8)**3
    
    If(allocated(header)) Deallocate(header)

    If(procID==0) Print *,'Reading positions...'
    If(procID==0) Write(Ulog,*) 'Reading positions...'

    nptot = 0
    mynproc = nproc / procNB
    firstp  = procID * mynproc + 1
    lastp   = (procID+1) * mynproc
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

    Call Mpi_Barrier(MPICube,mpierr)
    timeInt = Mpi_Wtime()
    tInitRead = timeInt - time0


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
             tmpx(idim,mynpart+1:mynpart+npartloc) = tmpdouble
          End Do
          
          ! Read velocities in a dummy variable
          Do idim = 1,ndim
             Read(1) tmpdouble
             ! put all velocities in tmpv vector
             tmpv(idim,mynpart+1:mynpart+npartloc) = tmpdouble
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
             Print *,'Erreur dans la repartition des particules'
             Call Mpi_Finalize(mpierr)
             Stop
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

  Subroutine outputcube()

    Use modparam
    Use modvariable
    Use modmpicom
    Implicit none

    Integer(kind=PRI) :: i, j
    Character(len=80) :: filecube
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

  Subroutine outputmass(ns,smin,nbs,massamas,cdmamas)

    Use modparam
    Use modvariable
    Use modmpicom
    Implicit none

    
    Integer(kind=4),                   intent(in) :: ns
    Integer(kind=PRI),                 intent(in) :: smin
    Integer(kind=4)                               :: nbs
    Integer(kind=4),  dimension(nbs),  intent(in) :: massamas
    Real(kind=DP),    dimension(3,nbs),intent(in) :: cdmamas

    Character(len=85) :: fileamas
    Character(len=5)  :: pid_char
    Integer(kind=4) :: i

    Write(pid_char(1:5),'(I5.5)') procID
    fileamas = trim(root)//"_masst_"//pid_char

    Open(Umas, File=trim(fileamas), Status='Unknown', Form='Unformatted')
    Write(Umas) int(ns,kind=4)
    Do i=1,nbs
       If(massamas(i) >= Mmin ) Write(Umas) int(i+smin-1,kind=8), massamas(i), cdmamas(:,i)
    End Do
    Close(Umas)



  End Subroutine outputmass

  !=======================================================================

  Subroutine outputstruct(np,ns,smin,nbs,idf,xf,vf,massamas)

    Use modparam
    Use modvariable
    Use modmpicom
    Implicit none

    Integer(kind=4),                  intent(in) :: np, ns, nbs
    Integer(kind=PRI),                intent(in) :: smin
    Integer(kind=PRI),dimension(np),  intent(in) :: idf
    Integer(kind=4),  dimension(nbs), intent(in) :: massamas
    Real   (kind=SP), dimension(3,np),intent(in) :: xf, vf


    Integer(kind=4) :: i
    Integer(kind=4) :: j, k, b
    Character(len=85) :: filestrct
    Character(len=5)  :: pid_char
    
    Write(pid_char(1:5),'(I5.5)') procID
    filestrct = trim(root)//"_strct_"//pid_char

    Open(Ustr, File=trim(filestrct), Status='Unknown',Form='Unformatted')

    Write(Ustr) int(ns,kind=4)
        
    b = 1
    Do i=1, nbs
       If(massamas(i) >= Mmin) Then
          Write(Ustr) massamas(i)
          Write(Ustr) ((xf(k,j),k=1,3),j=b,b+massamas(i)-1)
          Write(Ustr) ((vf(k,j),k=1,3),j=b,b+massamas(i)-1)
          Write(Ustr) (idf(j),j=b,b+massamas(i)-1)
          b = b+massamas(i)
       Else
          b = b+massamas(i)
       End If
    End Do

    Close(Ustr)

  End Subroutine outputstruct

End Module modio
