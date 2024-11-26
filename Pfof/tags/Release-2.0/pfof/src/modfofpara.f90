Module modfofpara

  Use modconst
  Use modmpicom

Contains

  Subroutine parafof()

    ! Halo detection is performed locally in each subdomains (by each process).
    ! Structures which are cut by a border between two subdomains are gathered.
    ! The criterium to gather two parts of a halo is simple: if one particle from in a halo is seperated from a particle in 
    ! another halo (on the other side of a border) by a distance smaller than the percolation parameter, then these two halos 
    ! are two parts of a single halo.


    Use modparam
    Use modio
    Use modtiming
    Use modsort
    Use modvariable

    Implicit none

    ! Local variables
    !    Integer(kind=1), dimension(mynpart) :: border       ! if border(i)!=0, then part. i is close to a border
    Integer(kind=1), dimension(:), allocatable :: border
    Integer(kind=4), dimension(:),allocatable :: adfirst,npcase
    !    Integer(kind=4), dimension(mynpart) :: lamas,pile,indl
    Integer(kind=4), dimension(:), allocatable :: lamas, pile, indl
    Integer(kind=4) :: i, j, k

    Integer(kind=PRI) :: strNBbeg, strNBend
    Integer(kind=PRI) :: smin,smax
    Integer(kind=PRI) :: tmpdi

    Integer(kind=4) :: nbridgearr, nbridgeava, nbridgegau
    Integer(kind=4) :: nbridgedro, nbridgebas, nbridgehau
    Integer(kind=PRI), dimension(:,:), allocatable :: bridgeArr, bridgeAva, bridgeGau
    Integer(kind=PRI), dimension(:,:), allocatable :: bridgeDro, bridgeBas, bridgeHau
    Integer(kind=PRI) :: li, ri
    Integer(kind=4) :: ib, f, fi, fp
    Integer(kind=4) :: fpID, fiID
    Integer(kind=4), dimension(6) :: mpireqsf, mpireqrf
    Integer(kind=4) :: mpistatf(MPI_STATUS_SIZE,6)

    Integer(kind=4) :: nbs
    Integer(kind=4) :: index, ipart
    Integer(kind=4) :: ipermut, sttmp, ind, ipile, ipilepr
    Integer(kind=4) :: ix, iy, iz, ic
    Integer(kind=4) :: i1, i2, i3, j1, j2, j3
    Integer(kind=4) :: nbborderloc, nbborder
    Integer(kind=4) :: deba,debb,debd,debg,debh,debr!!,nba,nbb,nbd,nbg,nbh,nbr
    Integer(kind=4), dimension(:), allocatable :: massamas
    Integer(kind=4) :: NPparProc
    Integer(kind=4) :: iproc
    Integer(kind=4) :: ngpp,nsd,nrespp
    Integer(kind=4) :: strPID
    Integer(kind=4) :: myNpartW
    Integer(kind=4), dimension(:),allocatable :: strPIDvec, strPIDvecloc  ! array of particle nb by process after
    ! distribution by structure ID
    Integer(kind=4) :: sendID, recvID   ! ID for process to send to and to receive from
    Integer(kind=4) :: nbrec, nbsend    ! nb of elements to recv and to send
    Integer(kind=4) :: mpistat(MPI_STATUS_SIZE)   ! MPI comm. status
    Integer(kind=4) :: nbpassage        ! nb of iteration in structure connection phase
    Integer(kind=4), dimension(6) :: nflagloc, nflagrecv
    Integer(kind=4) :: allStat
    Integer(kind=4) :: recvpoint
    Integer(kind=4) :: NstrW,NstrW_all
    Integer(kind=4) :: sid      ! tmp structure ID with 1 < sid < Nb of structure in the process (=> kind=4)

    Integer(kind=PRI) :: namas      ! tmp structure ID
    Integer(kind=PRI), dimension(:), allocatable :: strAva,strArr,strBas,strDro,strGau,strHau,strRec
    !  structure ID for particles close to a border and checked during the connection phase
    Integer(kind=PRI), dimension(:), allocatable :: stf, strSend ! structure ID array, and tmp array used for comm
    Integer(kind=PRI), dimension(:), allocatable :: idf, idSend  ! particle ID array and tmp array used for comm
    !    Integer(kind=PRI), dimension(mynpart) :: structure  ! structure ID for each particle in the process
    Integer(kind=PRI), dimension(:), allocatable :: structure

    Logical :: nopermut, finished

    Real(kind=SP), dimension(:,:), allocatable :: xf, posSend, vf, velSend
    Real(kind=DP), dimension(:,:), allocatable :: cdmamas,speedamas
    Real(kind=SP), dimension(:,:), allocatable :: posAva,posArr,posBas,posDro,posGau,posHau,posRec
    Real(kind=DP) :: delta
    Real(kind=SP) :: r, d, taille, xx, yy, zz
    Real(kind=SP) :: vxx, vyy, vzz                      ! temporary position
    Real(kind=SP) :: dx, dy, dz
    !Real(kind=SP) :: dvx, dvy, dvz, modv
    Real(kind=SP) :: r2
    Real(kind=SP) :: rT

    Real(kind=DP) :: ttmp0, ttmp1

!!$    ! Quantitites to read ramses_input file
!!$    Integer(kind=4) :: cnt, errcode
!!$    Character(len=200)             :: nomfich
!!$    Real(kind=SP)   :: dumbreal 

    ! INITIALIZATION


    ! Initialize timings if necessary
    If(dotimings) Then
       time0 = Mpi_Wtime()
    End If

    
    nflagloc = 0
    
    nbborderloc = 0

    nsd = int(procNB**(1./3.))
    tmpdi = ngrid / procNB
    ngpp = int(tmpdi,kind=4)

    nrespp = nres/nsd

    taille = float(nres)

    Allocate(adfirst(ngpp),npcase(ngpp))

    Allocate(border(mynpart), lamas(mynpart), indl(mynpart), pile(mynpart), structure(mynpart))

    Call init(nrespp,mynpart,ngpp,x,pile,adfirst,npcase,xmin,ymin,zmin,taille)

    r = perco
    rT = r/taille
    r2 = rT*rT

    Do i = 1, mynpart
       lamas(i) = i
       indl(i)  = i
    End Do

!!!namas = id(1)
    if(mynpart> 0) then
       namas = id(1) !!!!!RY
    else
       namas=0
    endif
    ipermut  = 2

    border = 0

    Call Mpi_Barrier(Mpi_Comm_World,mpierr)
    timeInt = Mpi_Wtime()
    tFoFinit = timeInt - time0

    If(procID==0) Print *,'Debut de FOF Local'

    !-----------!
    ! Local Friends of Friends
    !-----------!

    particules : Do i=1, mynpart

       If(procID==0 .and. mod(i,100000)==0) Print *,'Particule num.',i
       ipart = lamas(i)

       If (x(1,ipart) == 1.0) x(1,ipart) = 0.00000
       If (x(2,ipart) == 1.0) x(2,ipart) = 0.00000
       If (x(3,ipart) == 1.0) x(3,ipart) = 0.00000


       If(abs(x(1,ipart)-xmin) <= rT ) Then
          border(ipart) = border(ipart) + 4
          nflagloc(1) = nflagloc(1) + 1          
       End If
       If(abs(x(1,ipart)-xmax) <= rT ) Then
          border(ipart) = border(ipart) + 8
          nflagloc(2) = nflagloc(2) + 1
       End If
       If(abs(x(2,ipart)-ymin) <= rT ) Then
          border(ipart) = border(ipart) + 1
          nflagloc(3) = nflagloc(3) + 1
       End If
       If(abs(x(2,ipart)-ymax) <= rT ) Then
          border(ipart) = border(ipart) + 2
          nflagloc(4) = nflagloc(4) + 1
       End If
       If(abs(x(3,ipart)-zmin) <= rT ) Then
          border(ipart) = border(ipart) + 16
          nflagloc(5) = nflagloc(5) + 1
       End If
       If(abs(x(3,ipart)-zmax) <= rT ) Then
          border(ipart) = border(ipart) + 32
          nflagloc(6) = nflagloc(6) + 1
       End If
       If(border(ipart) /= 0 ) nbborderloc = nbborderloc + 1

       xx = x(1,ipart)
       yy = x(2,ipart)
       zz = x(3,ipart)

       vxx= v(1,ipart)
       vyy= v(2,ipart)
       vzz= v(3,ipart)

       ix = int((x(1,ipart)-xmin)*taille)
       iy = int((x(2,ipart)-ymin)*taille)
       iz = int((x(3,ipart)-zmin)*taille)
       structure(ipart) = namas

       dimx : Do i1 = -1,1
          j1 = ix + i1
          If(j1 >= nrespp) j1 = 0
          If(j1 <= -1) j1 = nrespp-1
          dimy : Do i2 = -1,1
             j2 = iy + i2
             If(j2 >= nrespp) j2 = 0
             If(j2 <= -1) j2 = nrespp-1
             dimz : Do i3 = -1,1
                j3 = iz + i3
                If(j3 >= nrespp) j3 = 0
                If(j3 <= -1) j3 = nrespp-1

                index = 1 + j1 + j2*nrespp + j3*nrespp*nrespp
                ic = 0

                nonvide : If (npcase(index) >= 1) Then  
                   ipile   = adfirst(index)
                   ipilepr = adfirst(index)
                   feuille : Do ind = 1,npcase(index)
                      nonself : If (ipart /= ipile) Then
                         dx = abs(xx-x(1,ipile))
                         dy = abs(yy-x(2,ipile))
                         dz = abs(zz-x(3,ipile))

                         d  = dx*dx+dy*dy+dz*dz

                         !  dvx = abs(vxx-v(1,ipile))
                         !  dvy = abs(vyy-v(2,ipile))
                         !  dvz = abs(vzz-v(3,ipile))

                         !                         modv= dvx*dvx+dvy*dvy+dvz*dvz
                         !write(*,*) procID,modv
                         !inhalo : If(d <= r2.and. modv<= 0.0005) Then
                         inhalo : If(d <= r2) Then
                            sttmp              = lamas(ipermut)
                            lamas(ipermut)     = ipile
                            lamas(indl(ipile)) = sttmp
                            indl(sttmp)        = indl(ipile)
                            indl(ipile)        = ipermut
                            ipermut            = ipermut+1
                            isfirst : If(ipile == adfirst(index))Then
                               adfirst(index) = pile(ipile)
                               ipilepr        = pile(ipile)
                            Else
                               pile(ipilepr)=pile(ipile)
                            End If isfirst
                            ic = ic + 1
                         Else
                            ipilepr = ipile
                         End If inhalo
                      Else
                         isfirst2 : If(ipile == adfirst(index))Then
                            adfirst(index) = pile(ipile)
                            ipilepr        = pile(ipile)
                         Else
                            pile(ipilepr) = pile(ipile)
                         End If isfirst2
                         ic = ic + 1
                      End If nonself
                      ipile = pile(ipile)
                   End Do feuille
                End If nonvide
                npcase(index) = npcase(index) - ic
             End Do dimz
          End Do dimy
       End Do dimx

       !-----------------------------------------------!
       ! Si l'amas est fini, on passe a l'amas suivant !
       !-----------------------------------------------!

       If (ipermut == (i+1)) Then

          ipermut  = ipermut + 1
          If(i<mynpart) namas = id(lamas(i+1))
       End If

    End Do particules

    If(procID==0) Print *,'FoF local termine'

    Print *,'Process ',procID,' termine apres ',Mpi_Wtime()-timeInt,'s dans fof local.'

    Deallocate(lamas,pile,indl)
    Deallocate(adfirst,npcase)

    ! call writememfilestep(1, 0, 2, 0, ' After fof local ', 'memall.txt',974, 0.0, 0.0)

    !-------------------!
    ! FOF local termine !
    !-------------------!

    Call Mpi_Reduce(nbborderloc, nbborder, 1, Mpi_Integer,Mpi_Sum,0,Mpi_Comm_World,mpierr)


    tFoFloc = Mpi_Wtime() - timeInt
    timeInt = Mpi_Wtime()

    If(procID==0) Print *,'Initialisation du raccordement'

    !----------------------------------------------------------------------!
    !Boucle sur les particules : si au bord, remplissage des tableaux face !
    !----------------------------------------------------------------------!

    debh = 1
    debb = 1
    deba = 1
    debr = 1
    debg = 1
    debd = 1  

    If(nflagloc(1) /= 0) Allocate (posArr(3,nflagloc(1)), strArr(nflagloc(1)))
    If(nflagloc(2) /= 0) Allocate (posAva(3,nflagloc(2)), strAva(nflagloc(2)))
    If(nflagloc(3) /= 0) Allocate (posGau(3,nflagloc(3)), strGau(nflagloc(3)))
    If(nflagloc(4) /= 0) Allocate (posDro(3,nflagloc(4)), strDro(nflagloc(4)))
    If(nflagloc(5) /= 0) Allocate (posBas(3,nflagloc(5)), strBas(nflagloc(5)))
    If(nflagloc(6) /= 0) Allocate (posHau(3,nflagloc(6)), strHau(nflagloc(6)))

    bufferface : Do i = 1, mynpart
       bord : If(border(i) /= 0) Then
          haut : If(border(i) >= 32) Then
             posHau (:,debh) = x(:,i)
             strHau (debh) = structure(i)
             debh = debh+1
             border(i) = border(i) - 32
          End If haut
          bas : If(border(i) >= 16) Then
             posBas (:,debb) = x(:,i)
             strBas (debb) = structure(i)
             debb = debb+1
             border(i) = border(i) - 16
          End If bas
          avant : If(border(i) >= 8) Then
             posAva (:,deba) = x(:,i)
             strAva (deba) = structure(i)
             deba = deba+1
             border(i) = border(i) - 8
          End If avant
          arriere : If(border(i) >=4 ) Then
             posArr (:,debr) = x(:,i)
             strArr (debr) = structure(i)
             debr = debr+1
             border(i) = border(i) - 4
          End If arriere
          droite : If(border(i) >=2 ) Then
             posDro (:,debd) = x(:,i)
             strDro (debd) = structure(i)
             debd = debd+1
             border(i) = border(i) - 2
          End If droite
          gauche : If(border(i) >= 1) Then
             posGau (:,debg) = x(:,i)
             strGau (debg) = structure(i)
             debg = debg+1
          End If gauche
       End If bord

    End Do bufferface


!!!!!! MODIF VINCENT BOUILLOT !!!!!!!!!
    Deallocate(border)
!!!!!! MODIF VINCENT BOUILLOT !!!!!!!!!

    ! Echange des positions, calcul des distances, determination des ponts de raccord 
    ! ---------------------
    
    Do f=1,3 ! boucle sur les faces en les prenant 2 par 2
       
       fi = 2*f-1
       fp = 2*f

       fiID = voisin(fi)
       fpID = voisin(fp) 

       Call Mpi_Irecv(nflagrecv(fi),1,Mpi_Integer,fiID,fiID,MPICube,mpireqrf(fi),mpierr)
       Call Mpi_Irecv(nflagrecv(fp),1,Mpi_Integer,fpID,fpID,MPICube,mpireqrf(fp),mpierr)
       Call Mpi_Isend(nflagloc(fp) ,1,Mpi_Integer,fpID,procID,MPICube,mpireqsf(fp),mpierr)
       Call Mpi_Isend(nflagloc(fi) ,1,Mpi_Integer,fiID,procID,MPICube,mpireqsf(fi),mpierr)

    End Do

    Call Mpi_Waitall(6,mpireqrf,mpistatf,mpierr)
    Call Mpi_Waitall(6,mpireqsf,mpistatf,mpierr)


    ! Envoi vers l'avant, reception venant de l'arriere
    Allocate (posRec(3,nflagrecv(1)))
    Call Mpi_Irecv(posRec,3*nflagrecv(1),Mpi_Real,voisin(1),voisin(1),MPICube,mpireqr1,mpierr)
    Call Mpi_ISend(posAva,3*nflagloc(2), Mpi_Real,voisin(2),procID,MPICube,mpireqs1,mpierr)

    Call Mpi_Wait(mpireqs1,mpistat,mpierr)
    Call Mpi_Wait(mpireqr1,mpistat,mpierr)
    
    nbridgearr = 0
    Do li = 1, nflagloc(1)
       Do ri = 1, nflagrecv(1)
          dx = abs(posArr(1,li)-posRec(1,ri))
          dx = min(dx,1.0-dx)
          dy = abs(posArr(2,li)-posRec(2,ri))
          dy = min(dy,1.0-dy)
          dz = abs(posArr(3,li)-posRec(3,ri))
          dz = min(dz,1.0-dz)
          
          d  = dx*dx+dy*dy+dz*dz
          If(d <= r2 ) Then
             nbridgearr = nbridgearr + 1
          End If
       End Do
    End Do

    Allocate(bridgeArr(2,nbridgearr))

    ind = 1
    Do li = 1, nflagloc(1)
       Do ri = 1, nflagrecv(1)
          dx = abs(posArr(1,li)-posRec(1,ri))
          dx = min(dx,1.0-dx)
          dy = abs(posArr(2,li)-posRec(2,ri))
          dy = min(dy,1.0-dy)
          dz = abs(posArr(3,li)-posRec(3,ri))
          dz = min(dz,1.0-dz)
          
          d  = dx*dx+dy*dy+dz*dz
          If(d <= r2 ) Then
             bridgeArr(1,ind) = li
             bridgeArr(2,ind) = ri
             ind = ind + 1
          End If
       End Do
    End Do

    Deallocate(posRec)


    ! Envoi vers l'arriere, reception venant de l'avant
    Allocate (posRec(3,nflagrecv(2)))
    Call Mpi_Irecv(posRec,3*nflagrecv(2),Mpi_Real,voisin(2),voisin(2),MPICube,mpireqr1,mpierr)
    Call Mpi_ISend(posArr,3*nflagloc(1), Mpi_Real,voisin(1),procID,MPICube,mpireqs1,mpierr)

    Call Mpi_Wait(mpireqs1,mpistat,mpierr)
    Call Mpi_Wait(mpireqr1,mpistat,mpierr)
    
    nbridgeava = 0
    Do li = 1, nflagloc(2)
       Do ri = 1, nflagrecv(2)
          dx = abs(posAva(1,li)-posRec(1,ri))
          dx = min(dx,1.0-dx)
          dy = abs(posAva(2,li)-posRec(2,ri))
          dy = min(dy,1.0-dy)
          dz = abs(posAva(3,li)-posRec(3,ri))
          dz = min(dz,1.0-dz)
          
          d  = dx*dx+dy*dy+dz*dz
          If(d <= r2 ) Then
             nbridgeava = nbridgeava+1
          End If
       End Do
    End Do

    Allocate(bridgeAva(2,nbridgeava))

    ind = 1
    Do li = 1, nflagloc(2)
       Do ri = 1, nflagrecv(2)
          dx = abs(posAva(1,li)-posRec(1,ri))
          dx = min(dx,1.0-dx)
          dy = abs(posAva(2,li)-posRec(2,ri))
          dy = min(dy,1.0-dy)
          dz = abs(posAva(3,li)-posRec(3,ri))
          dz = min(dz,1.0-dz)
          
          d  = dx*dx+dy*dy+dz*dz
          If(d <= r2 ) Then
             bridgeAva(1,ind) = li
             bridgeAva(2,ind) = ri
             ind = ind + 1
          End If
       End Do
    End Do

    Deallocate(posRec)
    If(Allocated(posAva)) Deallocate(posAva)
    If(Allocated(posArr)) Deallocate(posArr)


    ! Envoi vers la droite, reception venant de la gauche
    Allocate (posRec(3,nflagrecv(3)))
    Call Mpi_Irecv(posRec,3*nflagrecv(3),Mpi_Real,voisin(3),voisin(3),MPICube,mpireqr1,mpierr)
    Call Mpi_ISend(posDro,3*nflagloc(4), Mpi_Real,voisin(4),procID,MPICube,mpireqs1,mpierr)

    Call Mpi_Wait(mpireqs1,mpistat,mpierr)
    Call Mpi_Wait(mpireqr1,mpistat,mpierr)
    
    nbridgegau = 0
    Do li = 1, nflagloc(3)
       Do ri = 1, nflagrecv(3)
          dx = abs(posGau(1,li)-posRec(1,ri))
          dx = min(dx,1.0-dx)
          dy = abs(posGau(2,li)-posRec(2,ri))
          dy = min(dy,1.0-dy)
          dz = abs(posGau(3,li)-posRec(3,ri))
          dz = min(dz,1.0-dz)
          
          d  = dx*dx+dy*dy+dz*dz
          If(d <= r2 ) Then
             nbridgegau = nbridgegau+1
          End If
       End Do
    End Do

    Allocate(bridgeGau(2,nbridgegau))

    ind = 1
    Do li = 1, nflagloc(3)
       Do ri = 1, nflagrecv(3)
          dx = abs(posGau(1,li)-posRec(1,ri))
          dx = min(dx,1.0-dx)
          dy = abs(posGau(2,li)-posRec(2,ri))
          dy = min(dy,1.0-dy)
          dz = abs(posGau(3,li)-posRec(3,ri))
          dz = min(dz,1.0-dz)
          
          d  = dx*dx+dy*dy+dz*dz
          If(d <= r2 ) Then
             bridgeGau(1,ind) = li
             bridgeGau(2,ind) = ri
             ind = ind + 1
          End If
       End Do
    End Do

    Deallocate(posRec)


    ! Envoi vers la gauche, reception venant de la droite
    Allocate (posRec(3,nflagrecv(4)))
    Call Mpi_Irecv(posRec,3*nflagrecv(4),Mpi_Real,voisin(4),voisin(4),MPICube,mpireqr1,mpierr)
    Call Mpi_ISend(posGau,3*nflagloc(3), Mpi_Real,voisin(3),procID   ,MPICube,mpireqs1,mpierr)

    Call Mpi_Wait(mpireqs1,mpistat,mpierr)
    Call Mpi_Wait(mpireqr1,mpistat,mpierr)
    
    nbridgedro = 0
    Do li = 1, nflagloc(4)
       Do ri = 1, nflagrecv(4)
          dx = abs(posDro(1,li)-posRec(1,ri))
          dx = min(dx,1.0-dx)
          dy = abs(posDro(2,li)-posRec(2,ri))
          dy = min(dy,1.0-dy)
          dz = abs(posDro(3,li)-posRec(3,ri))
          dz = min(dz,1.0-dz)
          
          d  = dx*dx+dy*dy+dz*dz
          If(d <= r2 ) Then
             nbridgedro = nbridgedro+1
          End If
       End Do
    End Do

    Allocate(bridgeDro(2,nbridgedro))

    ind = 1
    Do li = 1, nflagloc(4)
       Do ri = 1, nflagrecv(4)
          dx = abs(posDro(1,li)-posRec(1,ri))
          dx = min(dx,1.0-dx)
          dy = abs(posDro(2,li)-posRec(2,ri))
          dy = min(dy,1.0-dy)
          dz = abs(posDro(3,li)-posRec(3,ri))
          dz = min(dz,1.0-dz)
          
          d  = dx*dx+dy*dy+dz*dz
          If(d <= r2 ) Then
             bridgeDro(1,ind) = li
             bridgeDro(2,ind) = ri
             ind = ind + 1
          End If
       End Do
    End Do

    Deallocate(posRec)
    If(Allocated(posGau)) Deallocate(posGau)
    If(Allocated(posDro)) Deallocate(posDro)



    ! Envoi vers le haut, reception venant du bas
    Allocate (posRec(3,nflagrecv(5)))
    Call Mpi_Irecv(posRec,3*nflagrecv(5),Mpi_Real,voisin(5),voisin(5),MPICube,mpireqr1,mpierr)
    Call Mpi_ISend(posHau,3*nflagloc(6), Mpi_Real,voisin(6),procID,   MPICube,mpireqs1,mpierr)

    Call Mpi_Wait(mpireqs1,mpistat,mpierr)
    Call Mpi_Wait(mpireqr1,mpistat,mpierr)
    
    nbridgebas = 0
    Do li = 1, nflagloc(5)
       Do ri = 1, nflagrecv(5)
          dx = abs(posBas(1,li)-posRec(1,ri))
          dx = min(dx,1.0-dx)
          dy = abs(posBas(2,li)-posRec(2,ri))
          dy = min(dy,1.0-dy)
          dz = abs(posBas(3,li)-posRec(3,ri))
          dz = min(dz,1.0-dz)
          
          d  = dx*dx+dy*dy+dz*dz
          If(d <= r2 ) Then
             nbridgebas = nbridgebas+1
          End If
       End Do
    End Do

    Allocate(bridgeBas(2,nbridgebas))

    ind = 1
    Do li = 1, nflagloc(5)
       Do ri = 1, nflagrecv(5)
          dx = abs(posBas(1,li)-posRec(1,ri))
          dx = min(dx,1.0-dx)
          dy = abs(posBas(2,li)-posRec(2,ri))
          dy = min(dy,1.0-dy)
          dz = abs(posBas(3,li)-posRec(3,ri))
          dz = min(dz,1.0-dz)
          
          d  = dx*dx+dy*dy+dz*dz
          If(d <= r2 ) Then
             bridgeBas(1,ind) = li
             bridgeBas(2,ind) = ri
             ind = ind + 1
          End If
       End Do
    End Do

    Deallocate(posRec)

    ! Envoi vers le bas, reception venant du haut
    Allocate (posRec(3,nflagrecv(6)))
    Call Mpi_Irecv(posRec,3*nflagrecv(6),Mpi_Real,voisin(6),voisin(6),MPICube,mpireqr1,mpierr)
    Call Mpi_ISend(posBas,3*nflagloc(5), Mpi_Real,voisin(5),procID   ,MPICube,mpireqs1,mpierr)

    Call Mpi_Wait(mpireqs1,mpistat,mpierr)
    Call Mpi_Wait(mpireqr1,mpistat,mpierr)
    
    nbridgehau = 0
    Do li = 1, nflagloc(6)
       Do ri = 1, nflagrecv(6)
          dx = abs(posHau(1,li)-posRec(1,ri))
          dx = min(dx,1.0-dx)
          dy = abs(posHau(2,li)-posRec(2,ri))
          dy = min(dy,1.0-dy)
          dz = abs(posHau(3,li)-posRec(3,ri))
          dz = min(dz,1.0-dz)
          
          d  = dx*dx+dy*dy+dz*dz
          If(d <= r2 ) Then
             nbridgehau = nbridgehau+1
          End If
       End Do
    End Do

    Allocate(bridgeHau(2,nbridgehau))

    ind = 1
    Do li = 1, nflagloc(6)
       Do ri = 1, nflagrecv(6)
          dx = abs(posHau(1,li)-posRec(1,ri))
          dx = min(dx,1.0-dx)
          dy = abs(posHau(2,li)-posRec(2,ri))
          dy = min(dy,1.0-dy)
          dz = abs(posHau(3,li)-posRec(3,ri))
          dz = min(dz,1.0-dz)
          
          d  = dx*dx+dy*dy+dz*dz
          If(d <= r2 ) Then
             bridgeHau(1,ind) = li
             bridgeHau(2,ind) = ri
             ind = ind + 1
          End If
       End Do
    End Do

    Deallocate(posRec)
    If(Allocated(posBas)) Deallocate(posBas)
    If(Allocated(posHau)) Deallocate(posHau)



    nbpassage = 0

    If(procID==0) Print *,'Debut du raccordement'

    raccordement : Do
       nopermut = .true.
       nbpassage = nbpassage + 1
       If(ProcID == 0) Print *,'Passage numero',nbpassage

       !------------------------------------------------------!
       ! Envoi du bas vers le bas et reception venant du haut !
       !------------------------------------------------------!

       sendID = voisin(5)
       recvID = voisin(6)
      
       Allocate (strRec(nflagrecv(6)))

       Call Mpi_ISend(strBas,nflagloc(5), MPI_PRI,sendID,sendID,MPICube,mpireqs1,mpierr)
       Call Mpi_IRecv(strRec,nflagrecv(6),MPI_PRI,recvID,procID,MPICube,mpireqr1,mpierr)
       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

       Do ib = 1, nbridgehau
          li = bridgeHau(1,ib)
          ri = bridgeHau(2,ib)
          raccordh : If(strRec(ri)< strHau(li) ) Then

             nopermut = .false.
             renumallh : Do k=1,mynpart
                If(structure(k) == strHau(li)) structure(k) = strRec(ri)
             End Do renumallh
             Do k=1, nflagloc(4)
                If(strDro(k) == strHau(li)) strDro(k) = strRec(ri)
             End Do
             Do k=1, nflagloc(3)
                If(strGau(k) == strHau(li)) strGau(k) = strRec(ri)
             End Do
             Do k=1, nflagloc(2)
                If(strAva(k) == strHau(li)) strAva(k) = strRec(ri)
             End Do
             Do k=1, nflagloc(1)
                If(strArr(k) == strHau(li)) strArr(k) = strRec(ri)
             End Do
             Do k=1, nflagloc(6)
                If(strHau(k) == strHau(li) .and. k/= li) strHau(k) = strRec(ri)
             End Do
             strHau(li) = strRec(ri)
          End If raccordh
       End Do

       Deallocate(strRec)


       !-------------------------------------------------------------------
       ! Envoi de ma frontiere haut vers le haut et reception venant du bas
       !-------------------------------------------------------------------

       sendID = voisin(6)
       recvID = voisin(5)
       Allocate (strRec(nflagrecv(5)))

       Call Mpi_ISend(strHau,nflagloc(6), MPI_PRI,sendID,sendID,MPICube,mpireqs1,mpierr)
       Call Mpi_IRecv(strRec,nflagrecv(5),MPI_PRI,recvID,procID,MPICube,mpireqr1,mpierr)
       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

       Do ib = 1, nbridgebas
          li = bridgeBas(1,ib)
          ri = bridgeBas(2,ib)

          raccordb : If(strRec(ri)< strBas(li) ) Then
             
             nopermut = .false.
             renumallb : Do k = 1,mynpart
                If(structure(k) == strBas(li)) structure(k) = strRec(ri)
             End Do renumallb
             Do k=1, nflagloc(4)
                If(strDro(k) == strBas(li)) strDro(k) = strRec(ri)
             End Do
             Do k=1, nflagloc(3)
                If(strGau(k) == strBas(li)) strGau(k) = strRec(ri)
             End Do
             Do k=1, nflagloc(2)
                If(strAva(k) == strBas(li)) strAva(k) = strRec(ri)
             End Do
             Do k=1, nflagloc(1)
                If(strArr(k) == strBas(li)) strArr(k) = strRec(ri)
             End Do
             Do k=1, nflagloc(5)
                If(strBas(k) == strBas(li) .and. k/= li) strBas(k) = strRec(ri)
             End Do
             strBas(li) = strRec(ri)
          End If raccordb
       End Do

       Deallocate(strRec)


       !--------------------------------------------------------------------------
       ! Envoi de ma frontiere droite vers la droite et reception venant de gauche
       !--------------------------------------------------------------------------

       sendID = voisin(4)
       recvID = voisin(3)

       Allocate (strRec(nflagrecv(3)))

       Call Mpi_ISend(strDro,nflagloc(4), MPI_PRI,sendID,sendID,MPICube,mpireqs1,mpierr)
       Call Mpi_IRecv(strRec,nflagrecv(3),MPI_PRI,recvID,procID,MPICube,mpireqr1,mpierr)
       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

       ! comparaison des part. proche de ma face 'gauche' avec les part. venant de mon voisin de gauche
       Do ib = 1, nbridgegau
          li = bridgeGau(1,ib)
          ri = bridgeGau(2,ib)
      
          raccordg : If(strRec(ri)< strGau(li) ) Then

             nopermut = .false.
             renumallg : Do k = 1,mynpart
                If(structure(k) == strGau(li)) structure(k) = strRec(ri)
             End Do renumallg
             Do k=1, nflagloc(5)
                If(strBas(k) == strGau(li)) strBas(k) = strRec(ri)
             End Do
             Do k=1, nflagloc(3)
                If(strGau(k) == strGau(li) .and. k/= li) strGau(k) = strRec(ri)
             End Do
             Do k=1, nflagloc(2)
                If(strAva(k) == strGau(li)) strAva(k) = strRec(ri)
             End Do
             Do k=1, nflagloc(1)
                If(strArr(k) == strGau(li)) strArr(k) = strRec(ri)
             End Do
             Do k=1, nflagloc(6)
                If(strHau(k) == strGau(li)) strHau(k) = strRec(ri)
             End Do
             strGau(li) = strRec(ri)
          End If raccordg
       End Do

       Deallocate(strRec)


       !---------------------------------------------------------------------
       ! Envoi de ma face gauche vers la gauche et reception venant de droite
       !---------------------------------------------------------------------

       sendID = voisin(3)
       recvID = voisin(4)

       Allocate (strRec(nflagrecv(4)))

       Call Mpi_ISend(strGau,nflagloc(3), MPI_PRI,sendID,sendID,MPICube,mpireqs1,mpierr)
       Call Mpi_IRecv(strRec,nflagrecv(4),MPI_PRI,recvID,procID,MPICube,mpireqr1,mpierr)
       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

       ! comparaison des part. proche de ma face 'droite' avec les part. venant de mon voisin de droite
       Do ib = 1, nbridgedro
          li = bridgeDro(1,ib)
          ri = bridgeDro(2,ib)
      
          raccordd : If(strRec(ri)< strDro(li) ) Then

             nopermut = .false.
             renumalld : Do k = 1,mynpart
                If(structure(k) == strDro(li)) structure(k) = strRec(ri)
             End Do renumalld
             Do k=1, nflagloc(4)
                If(strDro(k) == strDro(li) .and. k/= li) strDro(k) = strRec(ri)
             End Do
             Do k=1, nflagloc(5)
                If(strBas(k) == strDro(li)) strBas(k) = strRec(ri)
             End Do
             Do k=1, nflagloc(2)
                If(strAva(k) == strDro(li)) strAva(k) = strRec(ri)
             End Do
             Do k=1, nflagloc(1)
                If(strArr(k) == strDro(li)) strArr(k) = strRec(ri)
             End Do
             Do k=1, nflagloc(6)
                If(strHau(k) == strDro(li)) strHau(k) = strRec(ri)
             End Do
             strDro(li) = strRec(ri)
          End If raccordd
       End Do
    
       Deallocate(strRec)


       !--------------------------------------------------------------------------
       ! Envoi de ma frontiere avant vers l'avant et reception venant de l'arriere
       !--------------------------------------------------------------------------

       sendID = voisin(2)
       recvID = voisin(1)
       
       Allocate (strRec(nflagrecv(1)))


       Call Mpi_ISend(strAva,nflagloc(2), MPI_PRI,sendID,sendID,MPICube,mpireqs1,mpierr)
       Call Mpi_IRecv(strRec,nflagrecv(1),MPI_PRI,recvID,procID,MPICube,mpireqr1,mpierr)
       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

       ! comparaison des part. proche de ma face 'arriere' avec les part. venant de mon voisin de l'arriere
       Do ib = 1, nbridgearr
          li = bridgeArr(1,ib)
          ri = bridgeArr(2,ib)

          raccordr : If(strRec(ri)< strArr(li) ) Then

             nopermut = .false.
             renumallr : Do k = 1,mynpart
                If(structure(k) == strArr(li)) structure(k) = strRec(ri)
             End Do renumallr
             Do k=1, nflagloc(5)
                If(strBas(k) == strArr(li)) strBas(k) = strRec(ri)
             End Do
             Do k=1, nflagloc(3)
                If(strGau(k) == strArr(li)) strGau(k) = strRec(ri)
             End Do
             Do k=1, nflagloc(4)
                If(strDro(k) == strArr(li)) strDro(k) = strRec(ri)
             End Do
             Do k=1, nflagloc(1)
                If(strArr(k) == strArr(li) .and. k/= li) strArr(k) = strRec(ri)
             End Do
             Do k=1, nflagloc(6)
                If(strHau(k) == strArr(li)) strHau(k) = strRec(ri)
             End Do
             strArr(li) = strRec(ri)
          End If raccordr
       End Do

       Deallocate(strRec)


       !---------------------------------------------------------------------
       ! Envoi de ma face arriere vers l'arriere et reception venant de l'avant
       !---------------------------------------------------------------------

       sendID = voisin(1)
       recvID = voisin(2)

       Allocate (strRec(nflagrecv(2)))


       Call Mpi_ISend(strArr,nflagloc(1), MPI_PRI,sendID,sendID,MPICube,mpireqs1,mpierr)
       Call Mpi_IRecv(strRec,nflagrecv(2),MPI_PRI,recvID,procID,MPICube,mpireqr1,mpierr)
       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

       ! comparaison des part. proche de ma face 'avant' avec les part. venant de mon voisin de l'avant
       Do ib = 1, nbridgeava
          li = bridgeAva(1,ib)
          ri = bridgeAva(2,ib)

          raccorda : If(strRec(ri)< strAva(li) ) Then

             nopermut = .false.
             renumalla : Do k = 1,mynpart
                If(structure(k) == strAva(li)) structure(k) = strRec(ri)
             End Do renumalla
             Do k=1, nflagloc(4)
                If(strDro(k) == strAva(li)) strDro(k) = strRec(ri)
             End Do
             Do k=1, nflagloc(5)
                If(strBas(k) == strAva(li)) strBas(k) = strRec(ri)
             End Do
             Do k=1, nflagloc(2)
                If(strAva(k) == strAva(li) .and. k/= li) strAva(k) = strRec(ri)
             End Do
             Do k=1, nflagloc(3)
                If(strGau(k) == strAva(li)) strGau(k) = strRec(ri)
             End Do
             Do k=1, nflagloc(6)
                If(strHau(k) == strAva(li)) strHau(k) = strRec(ri)
             End Do
             strAva(li) = strRec(ri)
          End If raccorda
       End Do

       Deallocate(strRec)

       Call Mpi_AllReduce(nopermut,finished,1,Mpi_Logical,Mpi_Land,MPICube,mpierr)
       If(procID==0) Print *,'permutations terminees'
       If(finished) Exit
       
    End Do raccordement

    If(Allocated(strHau)) Deallocate(strHau)
    If(Allocated(strBas)) Deallocate(strBas)
    If(Allocated(strDro)) Deallocate(strDro)
    If(Allocated(strGau)) Deallocate(strGau)
    If(Allocated(strAva)) Deallocate(strAva)
    If(Allocated(strArr)) Deallocate(strArr) 

    Call Mpi_Barrier(Mpi_Comm_World,mpierr)
    tRaccord = Mpi_Wtime() - timeInt
    timeInt = Mpi_Wtime()

    If(procID==0) Print *,'Fin du raccordement'

    ! call writememfilestep(1, 0, 2, 0, ' End of the raccordement ', 'memall.txt',974, 0.0, 0.0)
    !---------------------!
    ! FIN du raccordement !
    !---------------------!

    If(procID==0) Print *,'Debut du calcul de la masse et de la position du cdm des stuctures'

    NPparProc = nptot/procNB

    strNBbeg = int(NPparProc,kind=8) * int(procID,kind=8) + 1
    strNBend = int(NPparProc,kind=8) * int((procID + 1),kind=8)
!    if(procID==procNB-1)strNBend=nptot!!!!!Y

    Allocate(strPIDvecloc(procNB), strPIDvec(procNB))
    strPIDvecloc = 0
    strPIDvec = 0

    Do i = 1,mynpart
       strPID = (structure(i)-1) / NPparProc + 1
       strPID = min(strPID,procNB)!!!!Y
       strPIDvecloc(strPID) = strPIDvecloc(strPID) + 1
    End Do

    Call Mpi_Allreduce(strPIDvecloc,strPIDvec,procNB,Mpi_Integer,Mpi_Sum,Mpi_Comm_World,mpierr)

    myNpartW = strPIDvec(procID+1)
    !    Print *,'P',procID,'MNP',myNpartW
    Allocate(xf(3,myNpartW),STAT=allStat)
    If(allStat > 0) Call EmergencyStop('Allocate failed for xf in amidami',2)
    Allocate(vf(3,myNpartW),STAT=allStat)
    If(allStat > 0) Call EmergencyStop('Allocate failed for vf in amidami',2)
    Allocate(stf(myNpartW),STAT=allStat)
    If(allStat > 0) Call EmergencyStop('Allocate failed for stf in amidami',2)
    Allocate(idf(myNpartW),STAT=allStat)
    If(allStat > 0) Call EmergencyStop('Allocate failed for idf in amidami',2)

    recvpoint = 1
    procMasse : Do iproc = 1,procNB-1
       If(procID==0) Print *,'Permutation ',iproc
       sendID = mod(procID + iproc,procNB)
       recvID = mod(procID + procNB - iproc,procNB)
       nbsend = strPIDvecloc(sendID+1)

       Call Mpi_ISend(nbsend,1,Mpi_Integer,sendID,sendID,MPICube,mpireqs1,mpierr)
       Call Mpi_IRecv(nbrec, 1,Mpi_Integer,recvID,procID,MPICube,mpireqr1,mpierr)

       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

       If(nbsend /= 0) Then
          Allocate(strSend(nbsend),posSend(3,nbsend),velSend(3,nbsend),idSend(nbsend))
          smin = int(NPparProc,kind=8) * int(sendID,kind=8) + 1
          smax = int(NPparProc,kind=8) * int((sendID+1),kind=8)
!          if(sendID==procNB-1)smax=nptot!!!!Y

          ind=1
          Do i=1, mynpart
             If(structure(i)>= smin .and. structure(i)<= smax) Then
                strSend(ind) = structure(i)
                posSend(:,ind) = x(:,i)
                velSend(:,ind) = v(:,i)
                idSend(ind) = id(i)
                ind = ind+1
             End If
          End Do

          If(ind /= nbsend +1 ) Then
             write(*,*)ind,nbsend,sendID,smin,smax,strPIDvec
             Call EmergencyStop('Error  1 while sharing structures for output.',2)
          End If

          Call Mpi_Isend(strSend,nbsend,  MPI_PRI, sendID,1,MPICube,mpireqs1,mpierr)
          Call Mpi_Isend(posSend,3*nbsend,Mpi_Real,sendID,2,MPICube,mpireqs2,mpierr)
          Call Mpi_Isend(velSend,3*nbsend,Mpi_Real,sendID,3,MPICube,mpireqs3,mpierr)
          Call Mpi_Isend(idSend, nbsend,  MPI_PRI, sendID,4,MPICube,mpireqs4,mpierr)

       End If

       If(nbrec /= 0) Then
          Call Mpi_IRecv(stf(recvpoint), nbrec,  MPI_PRI, recvID,1,MPICube,mpireqr1,mpierr)
          Call Mpi_IRecv(xf(1,recvpoint),3*nbrec,Mpi_Real,recvID,2,MPICube,mpireqr2,mpierr)
          Call Mpi_IRecv(vf(1,recvpoint),3*nbrec,Mpi_Real,recvID,3,MPICube,mpireqr3,mpierr)
          Call Mpi_IRecv(idf(recvpoint), nbrec,  MPI_PRI, recvID,4,MPICube,mpireqr4,mpierr)
          recvpoint=recvpoint+nbrec
       End If


       If(nbsend/=0) Then
          Call Mpi_Wait(mpireqs1,mpistat,mpierr)
          Call Mpi_Wait(mpireqs2,mpistat,mpierr)
          Call Mpi_Wait(mpireqs3,mpistat,mpierr)
          Call Mpi_Wait(mpireqs4,mpistat,mpierr)
          Deallocate(strSend,posSend,idSend,velSend)
       End If
       If(nbrec/=0) Then
          Call Mpi_Wait(mpireqr1,mpistat,mpierr)
          Call Mpi_Wait(mpireqr2,mpistat,mpierr)
          Call Mpi_Wait(mpireqr3,mpistat,mpierr)
          Call Mpi_Wait(mpireqr4,mpistat,mpierr)
       End If

    End Do procMasse


    ind= recvpoint
    Do i=1, mynpart
       If(structure(i)>= strNBbeg .and. structure(i)<= strNBend) Then
          stf(recvpoint)  = structure(i)
          xf(:,recvpoint) = x(:,i)
          vf(:,recvpoint) = v(:,i)
          idf(recvpoint)  = id(i)
          recvpoint = recvpoint+1
       End If
    End Do

    Deallocate(structure)
    Deallocate(strPIDvecloc,strPIDvec)
    Deallocate(x, v, id)

    If( recvpoint/= myNpartW +1) Then
       write(*,*)procID,' recvpoint=',recvpoint
       write(*,*)procID,' mynpartw=',mynpartw
       Call EmergencyStop('Error 2 while sharing structures for output.',2)
    End If

    ttmp0 = Mpi_Wtime()
    Print *,'Process ',procID,' repartition:',ttmp0-timeInt

    !!Call trirapide(1,myNpartW,stf,idf,xf,vf)
    Call tritas(stf,idf,xf,vf,myNpartW)

    ttmp1 = Mpi_Wtime()
    Print *,'Process ',procID,' sort :', ttmp1-ttmp0

    if(myNpartW.ne.0) then
       smin = stf(1)
       smax = stf(myNpartW)
    else
       smin=0
       smax=-1
    endif
    nbs = int(smax - smin + 1, kind=4)

    Allocate(massamas(nbs))

    massamas = 0

    Do i=1, myNpartW
       sid = int(stf(i) - smin + 1,kind=4)
       massamas(sid) = massamas(sid) + 1
    End Do

    NstrW = 0

    Do i = 1, nbs
       If(massamas(i) >= Mmin) NstrW = NstrW + 1
    End Do

    Call outputstruct(myNpartW,NstrW,smin,nbs,idf,xf,vf,massamas)

    Deallocate(idf)

    Allocate(cdmamas(3,nbs))

    massamas = 0
    cdmamas = 0.0d0

    Do i=1, myNpartW
       sid = int(stf(i) - smin + 1,kind=4)

       Do j=1,3
          If(massamas(sid) /= 0) Then
             delta = cdmamas(j,sid)/massamas(sid) - xf(j,i)
             If(abs(delta) > 0.5d0) Then
                If(delta > 0d0) cdmamas(j,sid) = cdmamas(j,sid) + 1.0d0
                If(delta < 0d0) cdmamas(j,sid) = cdmamas(j,sid) - 1.0d0
             End If
          End If
          cdmamas(j,sid) = cdmamas(j,sid) + xf(j,i)

       End Do
       massamas(sid) = massamas(sid) + 1

    End Do

    Deallocate(xf)

    Allocate(speedamas(3,nbs))

    massamas = 0

    Do i=1, myNpartW
       sid = int(stf(i) - smin + 1,kind=4)

       Do j=1,3
          If(massamas(sid) /= 0) speedamas(j,sid)=speedamas(j,sid)+vf(j,i)
       End Do
       massamas(sid) = massamas(sid) + 1

    End Do

    Deallocate(vf,stf)

    !    NstrW = 0

    Do i = 1, nbs
       !      If(massamas(i) >= Mmin) NstrW = NstrW + 1

       Do j=1,3
          If(massamas(i)  > 0  ) cdmamas(j,i) = cdmamas(j,i) / float(massamas(i))
          If(massamas(i)  > 0  ) speedamas(j,i) = speedamas(j,i) / float(massamas(i))
          If(cdmamas(j,i) < 0.0d0) cdmamas(j,i) = cdmamas(j,i) + 1.0d0
          If(cdmamas(j,i) > 1.0d0) cdmamas(j,i) = cdmamas(j,i) - 1.0d0
       End Do

    End Do

    Print *,'Process ',procID,' termine apres ',Mpi_Wtime()-timeInt,' s dans calcul des observables'

    Call Mpi_Reduce(NstrW,NstrW_all,1,Mpi_Integer,Mpi_Sum,0,Mpi_Comm_World,mpierr)

    If(procID==0) Print *,'Fin du calcul de la masse et des centres de masse des structures'
    If(procID==0) Print *,'Nombre de halos de masse superieure a', Mmin,':',NstrW_all
    If(procID==0) Write(Ulog,*) 'Nombre de halos de masse superieure a', Mmin,':',NstrW_all

    Call Mpi_Barrier(Mpi_Comm_World,mpierr)
    tObs = Mpi_Wtime() - timeInt
    timeInt = Mpi_Wtime()

    Call outputmass(NstrW,smin,nbs,massamas,cdmamas)

    ! call writememfilestep(1, 0, 2, 0, ' Before obs output ', 'memall.txt',974, 0.0, 0.0)

    Call Mpi_Barrier(Mpi_Comm_World,mpierr)
    tOut = Mpi_Wtime() - timeInt
    tFoF = Mpi_Wtime() - time0


    Deallocate(massamas, cdmamas,speedamas)


  End Subroutine parafof


  ! ======================================================================
  !       Initialize the stacks of pointers for local friends of friends

  ! The whole domain is divided into t^3 cubes, where t is the resolution of the simulation analized.
  ! Each process will implicitly considere only its subdomains because every particles treated by the process is located
  ! in the same subdomain.
  ! This routine initialize one stack of "pointers" to the index of the particles located in each cube.
  Subroutine init(n,np,ngpp,x,pile,adfirst,npcase,xmin,ymin,zmin,t)

    Implicit none

    ! Input parameters
    Integer(kind=4), Intent(in) :: n                    ! number of "cube" in each dimension in a subdomain, i.e. there are
    ! n^3 cubes in a subdomain
    Integer(kind=4), Intent(in) :: np                   ! number of particles
    Integer(kind=4), Intent(in) :: ngpp                 ! number of cubes in the subdomain
    Real(kind=SP), dimension(3,np), Intent(in)   :: x   ! positions of the particles
    Real(kind=SP), Intent(in) :: t                      ! resolution of the cosmological simulation

    ! Output parameters
    Integer(kind=4), dimension(ngpp), Intent(out) :: npcase    ! number of particles in a cube
    Integer(kind=4), dimension(ngpp), Intent(out) :: adfirst   !index of the first particle to examine in the cube
    Integer(kind=4), dimension(np), Intent(out) :: pile        ! stack of indices

    ! Local variables
    Integer(kind=4), dimension(:), allocatable :: adlast          ! index of last particle "added" to a cube
    Integer(kind=4) :: ipart, icase, ix, iy, iz         ! temporary indices
    Real(kind=SP) :: xx, yy, zz                         ! temporary position
    Real(kind=SP) :: xmin,ymin,zmin                     ! minimum x,y and z of the subdomain

    ! Initialization

    Allocate(adlast(ngpp))

    pile = 0
    adfirst = 0
    adlast = 0
    npcase = 0

    ! loop over the particles
    Do ipart = 1,np
       ! positions scaled to the resolution, i.e. between 0 and 512 for a 512^3 simulation for example
       xx = (x(1,ipart)-xmin)*t
       yy = (x(2,ipart)-ymin)*t
       zz = (x(3,ipart)-zmin)*t
       ! coordinates of the cube containing the particle: this cube is always in the subdomain of the process
       ix = int(xx)
       iy = int(yy)
       iz = int(zz)

       ! there may be some rounding problem because we switch from double precision to simple precision positions
       ! so in case we have som positions = 1.0 we place the particle in the "last" cube
       If(ix == n) ix = ix-1
       If(iy == n) iy = iy-1
       If(iz == n) iz = iz-1

       ! index of the cube
       icase = 1 + ix + iy*n + iz*n*n
       ! is there a bug somewhere with the index?
       If(icase>ngpp .or. icase<=0) Then
          Call EmergencyStop('Probleme in initialization of friends of friends ',2)
       End If
       ! there is one more particle in this cube
       npcase(icase) = npcase(icase) + 1
       ! if this particle is the first added to the cube then
       If(adfirst(icase) == 0) Then
          ! we set the adfirst array to the particle index
          adfirst(icase) = ipart
          ! and the adlast as well because it's the last particle added a this time
          adlast(icase)  = ipart
       Else
          ! else the stack points from the previous particle to this one
          pile(adlast(icase)) = ipart
          ! and we set the adlast array
          adlast(icase) = ipart
       End If

    End Do

    Deallocate(adlast)

  End Subroutine init


End Module modfofpara
