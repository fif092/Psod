!> @file
!!This file contains the serial FOF algorithm and the communications that merge the halo parts of haloes that extend across several processes.

!> This module contains the serial FOF algorithm and the communications that merge the halo parts of haloes that extend across several processes.
Module modfofpara

  Use modparameters
  Use modmpicom  

Contains

  !> Parallel FOF subroutine
  !! Finds haloes on each process then apply the merging process for haloes that extend across several processes.
  Subroutine fofparacone()
    ! Halo detection is performed locally in each subdomains (by each process).
    ! Structures which are cut by a border between two subdomains are gathered.
    ! The criterium to gather two parts of a halo is simple: if one particle from in a halo is seperated from a particle in 
    ! another halo (on the other side of a border) by a distance smaller than the percolation parameter, then these two halos 
    ! are two parts of a single halo.


    Use modio
    Use modsort
    Use modvariables

    Implicit none

    ! Local variables
    Integer(kind=1), dimension(:), allocatable :: border
    Integer(kind=4), dimension(:),allocatable :: adfirst,npcase
    Integer(kind=4), dimension(:), allocatable :: lamas, pile, indl
    Integer(kind=4) :: i, k, h

    Integer(kind=PRI) :: strNBbeg, strNBend
    Integer(kind=PRI) :: smin,smax
    Integer(kind=PRI) :: tmpdi

    Integer(kind=4) :: nbridgearr, nbridgeava, nbridgegau
    Integer(kind=4) :: nbridgedro, nbridgebas, nbridgehau
    Integer(kind=PRI), dimension(:,:), allocatable :: bridgeArr, bridgeAva, bridgeGau
    Integer(kind=PRI), dimension(:,:), allocatable :: bridgeDro, bridgeBas, bridgeHau
    Integer(kind=PRI) :: il, ir
    Integer(kind=4) :: ib, f, fi, li, fp, lp
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
    Integer(kind=4), dimension(:), allocatable :: haloMass
    Integer(kind=4) :: NPparProc
    Integer(kind=4) :: iproc
    Integer(kind=4) :: ngrid_fof,res_fof
    Integer(kind=4) :: refine_fof
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
    Integer(kind=4) :: haloNB,haloNB_all, halopartNB
    Integer(kind=4) :: sid      ! tmp structure ID with 1 < sid < Nb of structure in the process (=> kind=4)
    Integer(kind=4) :: halom
    Integer(kind=4) :: signx(3), signy(3), signz(3)
    Integer(kind=4) :: nsign

    Integer(kind=PRI) :: namas      ! tmp structure ID
    Integer(kind=PRI), dimension(:), allocatable :: strAva,strArr,strBas,strDro,strGau,strHau,strRec
    !  structure ID for particles close to a border and checked during the connection phase
    Integer(kind=PRI), dimension(:), allocatable :: stf, strSend ! structure ID array, and tmp array used for comm
    Integer(kind=PRI), dimension(:), allocatable :: idf, idSend  ! particle ID array and tmp array used for comm
    !    Integer(kind=PRI), dimension(mynpart) :: structure  ! structure ID for each particle in the process
    Integer(kind=PRI), dimension(:), allocatable :: structure
    Integer(kind=PRI), dimension(:), allocatable :: haloID, halopartID

    Logical :: nopermut, finished

    Real(kind=4), dimension(:,:), allocatable :: xf, posSend, vf, velSend
    Real(kind=4), dimension(:,:), allocatable :: halopartPos, halopartVel
    Real(kind=8), dimension(:,:), allocatable :: halocomPos, halocomVel
    Real(kind=4), dimension(:,:), allocatable :: posAva,posArr,posBas,posDro,posGau,posHau
    Real(kind=4) :: r, d, size, size_fof, xx, yy, zz
    Real(kind=4) :: dx, dy, dz
    Real(kind=8) :: r2
    Real(kind=4) :: rT
    Real(kind=4) :: xmin, xmax, ymin, ymax, zmin, zmax
    Real(kind=4) :: maxdim, percolim
    Integer(kind=4) :: log2maxdim

    Real(kind=4) :: cubeedge
    Real(kind=4) :: edgem1

    ! INITIALIZATION

    ! Initialize timings if necessary
!!$    If(dotimings) Then
!!$       time0 = Mpi_Wtime()
!!$    End If

    If(procID==0) Then
       Print *,'Beginning friend of friend halo detection'
       Print *,'...'
    End If

    ! defining boundaries of the cubic domain
    xmin = Real(boundaries(1),kind=4)
    xmax = Real(boundaries(2),kind=4)
    ymin = Real(boundaries(3),kind=4)
    ymax = Real(boundaries(4),kind=4)
    zmin = Real(boundaries(5),kind=4)
    zmax = Real(boundaries(6),kind=4)


!    cubeedge = maxval(partmax-partmin)

    cubeedge=xmax-xmin
    edgem1 = 1.0/cubeedge

    nflagloc = 0
    nbborderloc = 0

    ! we used a grid fof FoF
    ! this grid is more refined than the coarse grid used for the dynamic evolution
    ! the refinement factor is set to 2
    refine_fof = 2

    ! Changement : a verifier !!!
    ! Demander a Yann comment on fait ici !!!
    ! on va essayer un truc complique...
    ! on prend le nombre de process et on fait la racine cubique
    maxdim = procNB**(1./3.)
    ! on prend la valeur entiere du log2 de ce max
    log2maxdim = int(log(maxdim)/log(2.))
    ! ca nous donne le facteur par lequel on divise la res. Ramses pour avoir la res. FoF
    ! on fait ca pour maintenir le nombre de cellules fof a un niveau raisonnable 
    ! pour un souci d'occupation memoire
    res_fof = refine_fof * nres / (2**log2maxdim)
    ngrid_fof = res_fof**3 
    size = float(nres)
    size_fof = real(res_fof)

    Allocate(adfirst(ngrid_fof),npcase(ngrid_fof))

    Allocate(border(mynpart), lamas(mynpart), indl(mynpart), pile(mynpart), structure(mynpart))

    Call init(res_fof,mynpart,ngrid_fof,pos,pile,adfirst,npcase,xmin,ymin,zmin,size_fof,edgem1)

    percolim = 2**log2maxdim*cubeedge / real(refine_fof)

    nsign = 2
    !    If( perco*refine_fof >= 0.5 ) Then
    If( perco >= percolim/2.d0 ) Then
       nsign = 3
    End If

    If( perco >= percolim ) Then 
       If(procID == 0) Then
          Print *, ' ' 
          Print *, '*********************************************************************************'
          Print *, '*** perco * refine_fof should be < 2^log2maxdim*cubeedge                      ***'
          Print *, '*** refine_fof is set and used in modfofpara.f90, perco is an input parameter ***'
          Print *, '*** please contact Fabrice Roy (fabrice.roy@obspm.fr) for more information    ***'
          Print *, '*** Pfof is exiting                                                           ***'
          Print *, '*********************************************************************************'
          Print *, ' ' 
       End If
       Call Mpi_Barrier(Mpi_Comm_World, mpierr)
       Call Mpi_Abort(Mpi_Comm_World, 2, mpierr)
    End If

    r = perco
    rT = r / size
    r2 = rT*rT

    If(procID==0) Then
       Print *, '***********************'
       Print *, 'resolution fof',res_fof
       Print *, 'r/taille', rT
    End If

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

!!$    If(dotimings) Then
!!$       Call Mpi_Barrier(Mpi_Comm_World,mpierr)
!!$       timeInt = Mpi_Wtime()
!!$       tFoFinit = timeInt - time0
!!$    End If

    If(procID==0) Then
       Print *,'Local FoF halo detection'
    End If

    !-----------!
    ! Local Friends of Friends
    !-----------!

    signx(1) = 0
    signy(1) = 0
    signz(1) = 0
    signx(3) = 1
    signy(3) = 1
    signz(3) = 1

    particules : Do i=1, mynpart

       If(procID==0 .and. mod(i,100000)==0) Print *,'Particle n.',i
       ipart = lamas(i)

!!$       If (x(1,ipart) == 1.0) x(1,ipart) = 0.00000
!!$       If (x(2,ipart) == 1.0) x(2,ipart) = 0.00000
!!$       If (x(3,ipart) == 1.0) x(3,ipart) = 0.00000

       If(abs(pos(1,ipart)-xmin) <= rT ) Then
          border(ipart) = border(ipart) + 4_1
          nflagloc(1) = nflagloc(1) + 1          
       End If
       If(abs(pos(1,ipart)-xmax) <= rT ) Then
          border(ipart) = border(ipart) + 8_1
          nflagloc(2) = nflagloc(2) + 1
       End If
       If(abs(pos(2,ipart)-ymin) <= rT ) Then
          border(ipart) = border(ipart) + 1_1
          nflagloc(3) = nflagloc(3) + 1
       End If
       If(abs(pos(2,ipart)-ymax) <= rT ) Then
          border(ipart) = border(ipart) + 2_1
          nflagloc(4) = nflagloc(4) + 1
       End If
       If(abs(pos(3,ipart)-zmin) <= rT ) Then
          border(ipart) = border(ipart) + 16_1
          nflagloc(5) = nflagloc(5) + 1
       End If
       If(abs(pos(3,ipart)-zmax) <= rT ) Then
          border(ipart) = border(ipart) + 32_1
          nflagloc(6) = nflagloc(6) + 1
       End If
       If(border(ipart) /= 0 ) nbborderloc = nbborderloc + 1


       xx = pos(1,ipart)
       yy = pos(2,ipart)
       zz = pos(3,ipart)

       ix = int((pos(1,ipart)-xmin)*size_fof*edgem1)
       iy = int((pos(2,ipart)-ymin)*size_fof*edgem1)
       iz = int((pos(3,ipart)-zmin)*size_fof*edgem1)

       signx(2) = -1
       signy(2) = -1
       signz(2) = -1
       If(nsign == 2) Then
          !          If( (pos(1,ipart)-xmin)*size_fof - ix > 0.5 ) Then
          If( (pos(1,ipart)-xmin)*size_fof*edgem1 - ix > 0.5 ) Then
             signx(2) = 1
          End If
          !          If( (pos(2,ipart)-ymin)*size_fof - iy > 0.5 ) Then
          If( (pos(2,ipart)-ymin)*size_fof*edgem1 - iy > 0.5 ) Then
             signy(2) = 1
          End If
          !          If( (pos(3,ipart)-zmin)*size_fof - iz > 0.5 ) Then
          If( (pos(3,ipart)-zmin)*size_fof*edgem1 - iz > 0.5 ) Then
             signz(2) = 1
          End If
       End If
       structure(ipart) = namas

       ! nsign depends on perco * refine_fof: see above
       dimx : Do i1 = 1, nsign
          j1 = ix + signx(i1)
          If( (j1>= res_fof) .Or. (j1 < 0) ) Cycle

          dimy : Do i2 = 1, nsign
             j2 = iy + signy(i2)
             If( (j2>= res_fof) .Or. (j2 < 0) ) Cycle

             dimz : Do i3 = 1, nsign
                j3 = iz + signz(i3)
                If( (j3>= res_fof) .Or. (j3 < 0) ) Cycle

                index = 1 + j1 + j2*res_fof + j3*res_fof*res_fof
                ic = 0

                nonvide : If (npcase(index) >= 1) Then  
                   ipile   = adfirst(index)
                   ipilepr = adfirst(index)
                   feuille : Do ind = 1,npcase(index)
                      nonself : If (ipart /= ipile) Then
                         dx = abs(xx-pos(1,ipile))
                         dy = abs(yy-pos(2,ipile))
                         dz = abs(zz-pos(3,ipile))

                         d  = dx*dx+dy*dy+dz*dz

                         inhalo : If(d <= r2) Then
                            !                            Print *, 'D', procID, d, r2
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

    If(procID==0) Then
       Print *,'Local FoF halo detection completed'
    End If
!!$    If(dotimings) Then
!!$       Print *,'Process ',procID,' termine apres ',Mpi_Wtime()-timeInt,'s dans fof local.'
!!$    End If

    Deallocate(lamas,pile,indl)
    Deallocate(adfirst,npcase)



!!$    Allocate(massamas(mynpart))
!!$    massamas = 0
!!$    Do i=1,mynpart
!!$       ic = structure(i) - mynpart*procID
!!$       massamas(ic) = massamas(ic) + 1
!!$    End Do
!!$    
!!$    Do i=1,mynpart
!!$       If(massamas(i) > 10) Then
!!$          Print *,'MASS',procID,i,massamas(i)
!!$       End If
!!$    End Do
!!$    Deallocate(massamas)



    ! call writememfilestep(1, 0, 2, 0, ' After fof local ', 'memall.txt',974, 0.0, 0.0)

    !-------------------!
    ! FOF local termine !
    !-------------------!

    Call Mpi_Reduce(nbborderloc, nbborder, 1, Mpi_Integer,Mpi_Sum,0,Mpi_Comm_World,mpierr)


!!$    tFoFloc = Mpi_Wtime() - timeInt
!!$    timeInt = Mpi_Wtime()

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

    !If(nflagloc(1) /= 0)
    Allocate (posArr(3,nflagloc(1)), strArr(nflagloc(1)))
    !If(nflagloc(2) /= 0) 
    Allocate (posAva(3,nflagloc(2)), strAva(nflagloc(2)))
    !If(nflagloc(3) /= 0) 
    Allocate (posGau(3,nflagloc(3)), strGau(nflagloc(3)))
    !If(nflagloc(4) /= 0) 
    Allocate (posDro(3,nflagloc(4)), strDro(nflagloc(4)))
    !If(nflagloc(5) /= 0) 
    Allocate (posBas(3,nflagloc(5)), strBas(nflagloc(5)))
    !If(nflagloc(6) /= 0) 
    Allocate (posHau(3,nflagloc(6)), strHau(nflagloc(6)))

    bufferface : Do i = 1, mynpart
       bord : If(border(i) /= 0) Then
          haut : If(border(i) >= 32) Then
             posHau (:,debh) = pos(:,i)
             strHau (debh) = structure(i)
             debh = debh+1
             border(i) = border(i) - 32_1
          End If haut
          bas : If(border(i) >= 16) Then
             posBas (:,debb) = pos(:,i)
             strBas (debb) = structure(i)
             debb = debb+1
             border(i) = border(i) - 16_1
          End If bas
          avant : If(border(i) >= 8) Then
             posAva (:,deba) = pos(:,i)
             strAva (deba) = structure(i)
             deba = deba+1
             border(i) = border(i) - 8_1
          End If avant
          arriere : If(border(i) >=4 ) Then
             posArr (:,debr) = pos(:,i)
             strArr (debr) = structure(i)
             debr = debr+1
             border(i) = border(i) - 4_1
          End If arriere
          droite : If(border(i) >=2 ) Then
             posDro (:,debd) = pos(:,i)
             strDro (debd) = structure(i)
             debd = debd+1
             border(i) = border(i) - 2_1
          End If droite
          gauche : If(border(i) >= 1) Then
             posGau (:,debg) = pos(:,i)
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

    nflagrecv=0
    Do f=1,3 ! boucle sur les faces en les prenant 2 par 2

       fi = 2*f-1
       fp = 2*f


       fiID = neighbours(fi)
       fpID = neighbours(fp) 

       Call Mpi_Irecv(nflagrecv(fi),1,Mpi_Integer,fiID,1,Mpi_Comm_World,mpireqrf(fi),mpierr)
       Call Mpi_Irecv(nflagrecv(fp),1,Mpi_Integer,fpID,2,Mpi_Comm_World,mpireqrf(fp),mpierr)
       Call Mpi_Isend(nflagloc(fp) ,1,Mpi_Integer,fpID,1,Mpi_Comm_World,mpireqsf(fp),mpierr)
       Call Mpi_Isend(nflagloc(fi) ,1,Mpi_Integer,fiID,2,Mpi_Comm_World,mpireqsf(fi),mpierr)

    End Do

    Call Mpi_Waitall(6,mpireqrf,mpistatf,mpierr)
    Call Mpi_Waitall(6,mpireqsf,mpistatf,mpierr)

    If(.not.Allocated(posArr)) Allocate(posArr(3,0))
    If(.not.Allocated(posAva)) Allocate(posAva(3,0))
    Call findbridge(nflagloc, nflagrecv, posArr, posAva, nbridgearr, bridgeArr, nbridgeava, bridgeAva, r2, 1, 2)
    If(Allocated(posAva)) Deallocate(posAva)
    If(Allocated(posArr)) Deallocate(posArr)

    If(.not.Allocated(posGau)) Allocate(posGau(3,0))
    If(.not.Allocated(posDro)) Allocate(posDro(3,0))
    Call findbridge(nflagloc, nflagrecv, posGau, posDro, nbridgegau, bridgeGau, nbridgedro, bridgeDro, r2, 3, 4)
    If(Allocated(posGau)) Deallocate(posGau)
    If(Allocated(posDro)) Deallocate(posDro)

    If(.not.Allocated(posHau)) Allocate(posHau(3,0))
    If(.not.Allocated(posBas)) Allocate(posBas(3,0))
    Call findbridge(nflagloc, nflagrecv, posBas, posHau, nbridgebas, bridgeBas, nbridgehau, bridgeHau, r2, 5, 6)
    If(Allocated(posBas)) Deallocate(posBas)
    If(Allocated(posHau)) Deallocate(posHau)



    nbpassage = 0

    If(procID==0) Print *,'Debut du raccordement'

    raccordement : Do
       nopermut = .true.
       nbpassage = nbpassage + 1
       If(ProcID == 0) Print *,'Passage numero',nbpassage

#ifdef DEBUG
       Print *,'Process ',procID, ' : merge number ',nbpassage
       Print *,'Process ',procID, ' : nb of particles : ',mynpart
#endif

       !------------------------------------------------------!
       ! Envoi du bas vers le bas et reception venant du haut !
       !------------------------------------------------------!

       sendID = neighbours(5)
       recvID = neighbours(6)

       Allocate (strRec(nflagrecv(6)))

       Call Mpi_ISend(strBas,nflagloc(5), MPI_PRI,sendID,5,Mpi_Comm_World,mpireqs1,mpierr)
       Call Mpi_IRecv(strRec,nflagrecv(6),MPI_PRI,recvID,5,Mpi_Comm_World,mpireqr1,mpierr)
       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

#ifdef DEBUG
       Print *,'Process ',procID, ' : up begins'
#endif

       Do ib = 1, nbridgehau
          il = bridgeHau(1,ib)
          ir = bridgeHau(2,ib)
          raccordh : If(strRec(ir)< strHau(il) ) Then

             nopermut = .false.
             renumallh : Do k=1,mynpart
                If(structure(k) == strHau(il)) structure(k) = strRec(ir)
             End Do renumallh
             Do k=1, nflagloc(4)
                If(strDro(k) == strHau(il)) strDro(k) = strRec(ir)
             End Do
             Do k=1, nflagloc(3)
                If(strGau(k) == strHau(il)) strGau(k) = strRec(ir)
             End Do
             Do k=1, nflagloc(2)
                If(strAva(k) == strHau(il)) strAva(k) = strRec(ir)
             End Do
             Do k=1, nflagloc(1)
                If(strArr(k) == strHau(il)) strArr(k) = strRec(ir)
             End Do
             Do k=1, nflagloc(6)
                If(strHau(k) == strHau(il) .and. k/= il) strHau(k) = strRec(ir)
             End Do
             strHau(il) = strRec(ir)
          End If raccordh
       End Do

       Deallocate(strRec)

#ifdef DEBUG
       Print *,'Process ',procID, ' : up ends'
#endif

       !-------------------------------------------------------------------
       ! Envoi de ma frontiere haut vers le haut et reception venant du bas
       !-------------------------------------------------------------------

       sendID = neighbours(6)
       recvID = neighbours(5)
       Allocate (strRec(nflagrecv(5)))

       Call Mpi_ISend(strHau,nflagloc(6), MPI_PRI,sendID,6,Mpi_Comm_World,mpireqs1,mpierr)
       Call Mpi_IRecv(strRec,nflagrecv(5),MPI_PRI,recvID,6,Mpi_Comm_World,mpireqr1,mpierr)
       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

#ifdef DEBUG
       Print *,'Process ',procID, ' : down begins'
#endif

       Do ib = 1, nbridgebas
          il = bridgeBas(1,ib)
          ir = bridgeBas(2,ib)

          raccordb : If(strRec(ir)< strBas(il) ) Then

             nopermut = .false.
             renumallb : Do k = 1,mynpart
                If(structure(k) == strBas(il)) structure(k) = strRec(ir)
             End Do renumallb
             Do k=1, nflagloc(4)
                If(strDro(k) == strBas(il)) strDro(k) = strRec(ir)
             End Do
             Do k=1, nflagloc(3)
                If(strGau(k) == strBas(il)) strGau(k) = strRec(ir)
             End Do
             Do k=1, nflagloc(2)
                If(strAva(k) == strBas(il)) strAva(k) = strRec(ir)
             End Do
             Do k=1, nflagloc(1)
                If(strArr(k) == strBas(il)) strArr(k) = strRec(ir)
             End Do
             Do k=1, nflagloc(5)
                If(strBas(k) == strBas(il) .and. k/= il) strBas(k) = strRec(ir)
             End Do
             strBas(il) = strRec(ir)
          End If raccordb
       End Do

       Deallocate(strRec)

#ifdef DEBUG
       Print *,'Process ',procID, ' : down ends'
#endif

       !--------------------------------------------------------------------------
       ! Envoi de ma frontiere droite vers la droite et reception venant de gauche
       !--------------------------------------------------------------------------

       sendID = neighbours(4)
       recvID = neighbours(3)

       Allocate (strRec(nflagrecv(3)))

       Call Mpi_ISend(strDro,nflagloc(4), MPI_PRI,sendID,4,Mpi_Comm_World,mpireqs1,mpierr)
       Call Mpi_IRecv(strRec,nflagrecv(3),MPI_PRI,recvID,4,Mpi_Comm_World,mpireqr1,mpierr)
       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

#ifdef DEBUG
       Print *,'Process ',procID, ' : left begins'
#endif

       ! comparaison des part. proche de ma face 'gauche' avec les part. venant de mon neighbours de gauche
       Do ib = 1, nbridgegau
          il = bridgeGau(1,ib)
          ir = bridgeGau(2,ib)

          raccordg : If(strRec(ir)< strGau(il) ) Then

             nopermut = .false.
             renumallg : Do k = 1,mynpart
                If(structure(k) == strGau(il)) structure(k) = strRec(ir)
             End Do renumallg
             Do k=1, nflagloc(5)
                If(strBas(k) == strGau(il)) strBas(k) = strRec(ir)
             End Do
             Do k=1, nflagloc(3)
                If(strGau(k) == strGau(il) .and. k/= il) strGau(k) = strRec(ir)
             End Do
             Do k=1, nflagloc(2)
                If(strAva(k) == strGau(il)) strAva(k) = strRec(ir)
             End Do
             Do k=1, nflagloc(1)
                If(strArr(k) == strGau(il)) strArr(k) = strRec(ir)
             End Do
             Do k=1, nflagloc(6)
                If(strHau(k) == strGau(il)) strHau(k) = strRec(ir)
             End Do
             strGau(il) = strRec(ir)
          End If raccordg
       End Do

       Deallocate(strRec)

#ifdef DEBUG
       Print *,'Process ',procID, ' : left ends'
#endif


       !---------------------------------------------------------------------
       ! Envoi de ma face gauche vers la gauche et reception venant de droite
       !---------------------------------------------------------------------

       sendID = neighbours(3)
       recvID = neighbours(4)

       Allocate (strRec(nflagrecv(4)))

       Call Mpi_ISend(strGau,nflagloc(3), MPI_PRI,sendID,3,Mpi_Comm_World,mpireqs1,mpierr)
       Call Mpi_IRecv(strRec,nflagrecv(4),MPI_PRI,recvID,3,Mpi_Comm_World,mpireqr1,mpierr)
       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

#ifdef DEBUG
       Print *,'Process ',procID, ' : right begins'
#endif

       ! comparaison des part. proche de ma face 'droite' avec les part. venant de mon voisin de droite
       Do ib = 1, nbridgedro
          il = bridgeDro(1,ib)
          ir = bridgeDro(2,ib)

          raccordd : If(strRec(ir)< strDro(il) ) Then

             nopermut = .false.
             renumalld : Do k = 1,mynpart
                If(structure(k) == strDro(il)) structure(k) = strRec(ir)
             End Do renumalld
             Do k=1, nflagloc(4)
                If(strDro(k) == strDro(il) .and. k/= il) strDro(k) = strRec(ir)
             End Do
             Do k=1, nflagloc(5)
                If(strBas(k) == strDro(il)) strBas(k) = strRec(ir)
             End Do
             Do k=1, nflagloc(2)
                If(strAva(k) == strDro(il)) strAva(k) = strRec(ir)
             End Do
             Do k=1, nflagloc(1)
                If(strArr(k) == strDro(il)) strArr(k) = strRec(ir)
             End Do
             Do k=1, nflagloc(6)
                If(strHau(k) == strDro(il)) strHau(k) = strRec(ir)
             End Do
             strDro(il) = strRec(ir)
          End If raccordd
       End Do

       Deallocate(strRec)

#ifdef DEBUG
       Print *,'Process ',procID, ' : right ends'
#endif

       !--------------------------------------------------------------------------
       ! Envoi de ma frontiere avant vers l'avant et reception venant de l'arriere
       !--------------------------------------------------------------------------

       sendID = neighbours(2)
       recvID = neighbours(1)

       Allocate (strRec(nflagrecv(1)))


       Call Mpi_ISend(strAva,nflagloc(2), MPI_PRI,sendID,2,Mpi_Comm_World,mpireqs1,mpierr)
       Call Mpi_IRecv(strRec,nflagrecv(1),MPI_PRI,recvID,2,Mpi_Comm_World,mpireqr1,mpierr)
       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

#ifdef DEBUG
       Print *,'Process ',procID, ' : back begins'
#endif

       ! comparaison des part. proche de ma face 'arriere' avec les part. venant de mon voisin de l'arriere
       Do ib = 1, nbridgearr
          il = bridgeArr(1,ib)
          ir = bridgeArr(2,ib)

          raccordr : If(strRec(ir)< strArr(il) ) Then

             nopermut = .false.
             renumallr : Do k = 1,mynpart
                If(structure(k) == strArr(il)) structure(k) = strRec(ir)
             End Do renumallr
             Do k=1, nflagloc(5)
                If(strBas(k) == strArr(il)) strBas(k) = strRec(ir)
             End Do
             Do k=1, nflagloc(3)
                If(strGau(k) == strArr(il)) strGau(k) = strRec(ir)
             End Do
             Do k=1, nflagloc(4)
                If(strDro(k) == strArr(il)) strDro(k) = strRec(ir)
             End Do
             Do k=1, nflagloc(1)
                If(strArr(k) == strArr(il) .and. k/= il) strArr(k) = strRec(ir)
             End Do
             Do k=1, nflagloc(6)
                If(strHau(k) == strArr(il)) strHau(k) = strRec(ir)
             End Do
             strArr(il) = strRec(ir)
          End If raccordr
       End Do

       Deallocate(strRec)

#ifdef DEBUG
       Print *,'Process ',procID, ' : back ends'
#endif

       !---------------------------------------------------------------------
       ! Envoi de ma face arriere vers l'arriere et reception venant de l'avant
       !---------------------------------------------------------------------

       sendID = neighbours(1)
       recvID = neighbours(2)

       Allocate (strRec(nflagrecv(2)))

       Call Mpi_ISend(strArr,nflagloc(1), MPI_PRI,sendID,1,Mpi_Comm_World,mpireqs1,mpierr)
       Call Mpi_IRecv(strRec,nflagrecv(2),MPI_PRI,recvID,1,Mpi_Comm_World,mpireqr1,mpierr)
       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

#ifdef DEBUG
       Print *,'Process ',procID, ' : front begins'
#endif

       ! comparaison des part. proche de ma face 'avant' avec les part. venant de mon voisin de l'avant
       Do ib = 1, nbridgeava
          il = bridgeAva(1,ib)
          ir = bridgeAva(2,ib)

          raccorda : If(strRec(ir)< strAva(il) ) Then

             nopermut = .false.
             renumalla : Do k = 1,mynpart
                If(structure(k) == strAva(il)) structure(k) = strRec(ir)
             End Do renumalla
             Do k=1, nflagloc(4)
                If(strDro(k) == strAva(il)) strDro(k) = strRec(ir)
             End Do
             Do k=1, nflagloc(5)
                If(strBas(k) == strAva(il)) strBas(k) = strRec(ir)
             End Do
             Do k=1, nflagloc(2)
                If(strAva(k) == strAva(il) .and. k/= il) strAva(k) = strRec(ir)
             End Do
             Do k=1, nflagloc(3)
                If(strGau(k) == strAva(il)) strGau(k) = strRec(ir)
             End Do
             Do k=1, nflagloc(6)
                If(strHau(k) == strAva(il)) strHau(k) = strRec(ir)
             End Do
             strAva(il) = strRec(ir)
          End If raccorda
       End Do

       Deallocate(strRec)

#ifdef DEBUG
       Print *,'Process ',procID, ' : front ends'
#endif

       Call Mpi_AllReduce(nopermut,finished,1,Mpi_Logical,Mpi_Land,Mpi_Comm_World,mpierr)
       If(procID==0) Print *,'permutations terminees'
       If(finished) Exit

    End Do raccordement

    If(Allocated(strHau)) Deallocate(strHau)
    If(Allocated(strBas)) Deallocate(strBas)
    If(Allocated(strDro)) Deallocate(strDro)
    If(Allocated(strGau)) Deallocate(strGau)
    If(Allocated(strAva)) Deallocate(strAva)
    If(Allocated(strArr)) Deallocate(strArr) 

!!$    If(dotimings) Then
!!$       Call Mpi_Barrier(Mpi_Comm_World,mpierr)
!!$       tRaccord = Mpi_Wtime() - timeInt
!!$       timeInt = Mpi_Wtime()
!!$    End If

    If(procID==0) Print *,'Fin du raccordement'

    ! call writememfilestep(1, 0, 2, 0, ' End of the raccordement ', 'memall.txt',974, 0.0, 0.0)
    !---------------------!
    ! FIN du raccordement !
    !---------------------!

    If(procID==0) Print *,'Debut du calcul de la masse et de la position du cdm des stuctures'

    tmpdi = npart/procNB
    NPparProc = int(tmpdi, kind=4)

    strNBbeg = int(NPparProc,kind=8) * int(procID,kind=8) + 1
    strNBend = int(NPparProc,kind=8) * int((procID + 1),kind=8)
    if(procID==procNB-1)strNBend=npart!!!!!Y   FR : Pas necessaire !

    Allocate(strPIDvecloc(procNB), strPIDvec(procNB))
    strPIDvecloc = 0
    strPIDvec = 0

    Do i = 1,mynpart
       strPID = int((structure(i)-1) / NPparProc, kind=4) + 1
       strPID = min(strPID,procNB)!!!!Y   FR : Pas necessaire !
       strPIDvecloc(strPID) = strPIDvecloc(strPID) + 1
    End Do

    Call Mpi_Allreduce(strPIDvecloc,strPIDvec,procNB,Mpi_Integer,Mpi_Sum,Mpi_Comm_World,mpierr)

    myNpartW = strPIDvec(procID+1)
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

       Call Mpi_ISend(nbsend,1,Mpi_Integer,sendID,1,Mpi_Comm_World,mpireqs1,mpierr)
       Call Mpi_IRecv(nbrec, 1,Mpi_Integer,recvID,1,Mpi_Comm_World,mpireqr1,mpierr)

       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

       If(nbsend /= 0) Then
          Allocate(strSend(nbsend),posSend(3,nbsend),velSend(3,nbsend),idSend(nbsend))
          smin = int(NPparProc,kind=8) * int(sendID,kind=8) + 1
          smax = int(NPparProc,kind=8) * int((sendID+1),kind=8)
          if(sendID==procNB-1)smax=npart!!!!Y

          ind=1
          Do i=1, mynpart
             If(structure(i)>= smin .and. structure(i)<= smax) Then
                strSend(ind) = structure(i)
                posSend(:,ind) = pos(:,i)
                velSend(:,ind) = vel(:,i)
                idSend(ind) = id(i)
                ind = ind+1
             End If
          End Do

          If(ind /= nbsend +1 ) Then
             write(*,*)ind,nbsend,sendID,smin,smax,strPIDvec
             Call EmergencyStop('Error  1 while sharing structures for output.',2)
          End If

          Call Mpi_Isend(strSend,nbsend,  MPI_PRI, sendID,1,Mpi_Comm_World,mpireqs1,mpierr)
          Call Mpi_Isend(posSend,3*nbsend,Mpi_Real,sendID,2,Mpi_Comm_World,mpireqs2,mpierr)
          Call Mpi_Isend(velSend,3*nbsend,Mpi_Real,sendID,3,Mpi_Comm_World,mpireqs3,mpierr)
          Call Mpi_Isend(idSend, nbsend,  MPI_PRI, sendID,4,Mpi_Comm_World,mpireqs4,mpierr)

       End If

       If(nbrec /= 0) Then
          Call Mpi_IRecv(stf(recvpoint), nbrec,  MPI_PRI, recvID,1,Mpi_Comm_World,mpireqr1,mpierr)
          Call Mpi_IRecv(xf(1,recvpoint),3*nbrec,Mpi_Real,recvID,2,Mpi_Comm_World,mpireqr2,mpierr)
          Call Mpi_IRecv(vf(1,recvpoint),3*nbrec,Mpi_Real,recvID,3,Mpi_Comm_World,mpireqr3,mpierr)
          Call Mpi_IRecv(idf(recvpoint), nbrec,  MPI_PRI, recvID,4,Mpi_Comm_World,mpireqr4,mpierr)
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
          xf(:,recvpoint) = pos(:,i)
          vf(:,recvpoint) = vel(:,i)
          idf(recvpoint)  = id(i)
          recvpoint = recvpoint+1
       End If
    End Do

    Deallocate(structure)
    Deallocate(strPIDvecloc,strPIDvec)
    Deallocate(pos, vel, id)

    If( recvpoint/= myNpartW +1) Then
       write(*,*)procID,' recvpoint=',recvpoint
       write(*,*)procID,' mynpartw=',mynpartw
       Call EmergencyStop('Error 2 while sharing structures for output.',2)
    End If

!!$    If(dotimings) Then
!!$       ttmp0 = Mpi_Wtime()
!!$       Print *,'Process ',procID,' repartition:',ttmp0-timeInt
!!$    End If

    !!Call trirapide(1,myNpartW,stf,idf,xf,vf)
    Call tritas(stf,idf,xf,vf,myNpartW)

!!$    If(dotimings) Then
!!$       ttmp1 = Mpi_Wtime()
!!$       Print *,'Process ',procID,' sort :', ttmp1-ttmp0
!!$    End If

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


    haloNB = 0
    halopartNB = 0

    ! Compute total nb of particles in halo with M >= Mmin and nb of halos with M >= Mmin
    Do i = 1, nbs
       If(massamas(i) >= Mmin) Then
          haloNB = haloNB + 1
          halopartNB = halopartNB + massamas(i)
       End If
    End Do

    ! Keep positions, velocities and id for particles in halo with M >= Mmin
    Allocate(halopartPos(3,halopartNB))
    Allocate(halopartVel(3,halopartNB))
    Allocate(halopartID(halopartNB))
    ! Keep mass and id for halos with M >= Mmin
    Allocate(haloMass(haloNB))
    Allocate(haloID(haloNB))
    fp = 1
    lp = 0
    fi = 1
    li = 0
    h = 1
    Do i = 1, nbs
       If(massamas(i) >= Mmin) Then
          lp = fp + massamas(i) - 1
          li = fi + massamas(i) - 1
          halopartPos(:,fp:lp) = xf(:,fi:li)
          halopartVel(:,fp:lp) = vf(:,fi:li)
          halopartID(fp:lp) = idf(fi:li)
          haloMass(h) = massamas(i)
          haloID(h) = stf(fi)
          fp = lp + 1
          fi = li + 1
          h = h + 1
       Else
          li = fi + massamas(i) - 1
          fi = li + 1
       End If
    End Do
    If(li /= myNpartW) Then
       Print *, 'Error when keeping particles positions in halo with M >= Mmin on process ', procID
    End If
    If(lp /= halopartNB) Then
       Print *, 'Error when keeping particles positions in halo with M >= Mmin on process ', procID
    End If

    Call h5writehalopart(haloNB, halopartNB, haloMass, haloID, halopartPos, halopartVel, halopartID)


    Deallocate(xf, vf, idf)
    Deallocate(massamas)
    Deallocate(halopartID)

    Allocate(halocomPos(3,haloNB))
    Allocate(halocomVel(3,haloNB))

    halocomPos = 0.d0
    halocomVel = 0.d0

    fi = 1
    Do h = 1, haloNB
       halom = 1
       li = fi + haloMass(h) - 1
       halocomPos(:,h) = halopartPos(:,fi)
       halocomVel(:,h) = halopartVel(:,fi)
       Do i = fi+1, li
!!$          delta(:) = halocomPos(:,h) / halom - halopartPos(:,i)
!!$          Do j = 1, 3
!!$             If(abs(delta(j)) > 0.5d0) Then
!!$                If(delta(j) > 0d0) halocomPos(j,h) = halocomPos(j,h) + 1.0d0
!!$                If(delta(j) < 0d0) halocomPos(j,h) = halocomPos(j,h) - 1.0d0
!!$             End If
!!$          End Do
          halom = halom + 1
          halocomPos(:,h) = halocomPos(:,h) + halopartPos(:,i)
          halocomVel(:,h) = halocomVel(:,h) + halopartVel(:,i)
       End Do
       fi = li + 1

    End Do

    Deallocate(halopartPos, halopartVel)

    Do h = 1, haloNB
       halocomPos(:,h) = halocomPos(:,h) / haloMass(h)
       halocomVel(:,h) = halocomVel(:,h) / haloMass(h)
!!$       Do j = 1, 3
!!$          If(halocomPos(j,h) > 1.0d0) halocomPos(j,h) = halocomPos(j,h) - 1.0d0
!!$          If(halocomPos(j,h) < 0.0d0) halocomPos(j,h) = halocomPos(j,h) + 1.0d0
!!$       End Do
    End Do


!!$    If(dotimings) Then
!!$       Print *,'Process ',procID,' termine apres ',Mpi_Wtime()-timeInt,' s dans calcul des observables'
!!$    End If

    Call Mpi_Reduce(haloNB,haloNB_all,1,Mpi_Integer,Mpi_Sum,0,Mpi_Comm_World,mpierr)

    If(procID==0) Print *,'Fin du calcul de la masse et des centres de masse des structures'
    If(procID==0) Print *,'Nombre de halos de masse superieure a', Mmin,':',haloNB_all
    !    If(procID==0) Write(Ulog,*) 'Nombre de halos de masse superieure a', Mmin,':',haloNB_all

!!$    If(dotimings) Then
!!$       Call Mpi_Barrier(Mpi_Comm_World,mpierr)
!!$       tObs = Mpi_Wtime() - timeInt
!!$       timeInt = Mpi_Wtime()
!!$    End If

    Call mpih5writehalomass(haloNB, haloMass, halocomPos, halocomVel, haloID)

    Deallocate(haloMass, halocomPos, halocomVel, haloID)

!!$    If(dotimings) Then
!!$       Call Mpi_Barrier(Mpi_Comm_World,mpierr)
!!$       tOut = Mpi_Wtime() - timeInt
!!$       tFoF = Mpi_Wtime() - time0
!!$    End If


  End Subroutine fofparacone


  ! ======================================================================
  !> Initialize the stacks of pointers for local friends of friends.<br>
  !! The whole domain is divided into t^3 cubes, where t is the resolution of the simulation analized times the refine_fof factor.<br>
  !! Each process will implicitly considere only its subdomains because every particles treated by the process is located
  !! in the same subdomain.<br>
  !! This routine initializes one stack of "pointers" to the index of the particles located in each cube.
  Subroutine init(n,np,ngpp,x,pile,adfirst,npcase,xmin,ymin,zmin,t,edgem1)


    Implicit none

    ! Input parameters
    Integer(kind=4), Intent(in) :: n                    !< number of "cube" in each dimension in a subdomain, i.e. there are
    !! n^3 cubes in a subdomain
    Integer(kind=4), Intent(in) :: np                   !< number of particles in the subdomain
    Integer(kind=4), Intent(in) :: ngpp                 !< number of "FOF" cubes in the subdomain
    Real(kind=4), dimension(3,np), Intent(in)   :: x   !< positions of the particles
    Real(kind=4), Intent(in) :: t                      !< resolution of the cosmological simulation times refine_fof factor
    Real(kind=4), Intent(in) :: edgem1

    ! Output parameters
    Integer(kind=4), dimension(ngpp), Intent(out) :: npcase    !< number of particles in each cube
    Integer(kind=4), dimension(ngpp), Intent(out) :: adfirst   !< index of the first particle to examine in the cube
    Integer(kind=4), dimension(np), Intent(out) :: pile        !< stack of indices

    ! Local variables
    Integer(kind=4), dimension(:), allocatable :: adlast          ! index of last particle "added" to a cube
    Integer(kind=4) :: ipart, icase, ix, iy, iz         ! temporary indices
    Real(kind=4) :: xx, yy, zz                         ! temporary position
    Real(kind=4) :: xmin,ymin,zmin                     ! minimum x,y and z of the subdomain

    ! Initialization

    Allocate(adlast(ngpp))

    pile = 0
    adfirst = 0
    adlast = 0
    npcase = 0

!    edgem1 = 1.0/real(boundaries(2)-boundaries(1),kind=4)

    ! loop over the particles
    Do ipart = 1,np
       ! positions scaled to the resolution, i.e. between 0 and 512 for a 512^3 simulation for example
       xx = (x(1,ipart)-xmin)*t*edgem1
       yy = (x(2,ipart)-ymin)*t*edgem1
       zz = (x(3,ipart)-zmin)*t*edgem1
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
          Print *,'Prob. pour une particule:',icase, ngpp,t, ix, iy, iz, n, xx, yy, zz,x(:,ipart), xmin, ymin, zmin
          Call EmergencyStop('Problem in initialization of friends of friends',2)
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

  !------------------------------------------------------------------------------------------------------------------------------------
  !> Finds the pairs of particles that should be linked through the merging procedure. <br>
  !! two directions are considered simultaneously by the local process: <br>
  !! 1= back  2= front <br>
  !! 1= left  2= right <br>
  !! 1= bottom  2= top <br>
  !! meaning that the local process exchanges information with the processes located in these two directions.
  Subroutine findbridge(nloc, nrecv, pos1, pos2, nbridge1, bridge1, nbridge2, bridge2, r2, v1, v2)

    Use modvariables, only : PRI, neighbours

    Implicit None
    ! input
    Integer(kind=4), intent(in), dimension(6) :: nloc   !< local number of particles flagged for link detection in direction 1
    Integer(kind=4), intent(in), dimension(6) :: nrecv  !< local number of particles flagged for link detection in direction 2
    Integer(kind=4), intent(in) :: v1                   !< neighbour index in direction 1
    Integer(kind=4), intent(in) :: v2                   !< neighbour index in direction 2
    Real(kind=4), intent(in), dimension(3,*) :: pos1   !< positions of local particles flagged for link detection in direction 1
    Real(kind=4), intent(in), dimension(3,*) :: pos2   !< positions of distant particles flagged for link detection in direction 2
    Real(kind=8), intent(in) :: r2                     !< percolation length

    ! output
    Integer(kind=4), intent(inout) :: nbridge1                                !< number of "bridges" found, e.g. pairs of particles that should be linked, in direction 1
    Integer(kind=4), intent(inout) :: nbridge2                                !< number of "bridges" found, e.g. pairs of particles that should be linked, in direction 2
    Integer(kind=PRI), dimension(:,:), allocatable, intent(inout) :: bridge1  !< array of indices of the particles forming a bridge in direction 1
    Integer(kind=PRI), dimension(:,:), allocatable, intent(inout) :: bridge2  !< array of indices of the particles forming a bridge in direction 2

    ! local
    Integer(kind=4) :: mpistat(MPI_STATUS_SIZE)   ! MPI comm. status
    Integer(kind=4) :: li, ri, ind
    Real(kind=4) :: d, dx, dy, dz
    Real(kind=4), dimension(:,:), allocatable :: posRec


    ! Envoi vers l'avant, reception venant de l'arriere
    Allocate (posRec(3,nrecv(v1)))
    Call Mpi_Irecv(posRec,3*nrecv(v1),Mpi_Real,neighbours(v1),1,Mpi_Comm_World,mpireqr1,mpierr)
    Call Mpi_ISend(pos2,3*nloc(v2), Mpi_Real,neighbours(v2),1,Mpi_Comm_World,mpireqs1,mpierr)

    Call Mpi_Wait(mpireqs1,mpistat,mpierr)
    Call Mpi_Wait(mpireqr1,mpistat,mpierr)


    nbridge1 = 0
    Do li = 1, nloc(v1)
       Do ri = 1, nrecv(v1)
          dx = abs(pos1(1,li)-posRec(1,ri))
          !          dx = min(dx,1.0-dx)
          dy = abs(pos1(2,li)-posRec(2,ri))
          !          dy = min(dy,1.0-dy)
          dz = abs(pos1(3,li)-posRec(3,ri))
          !          dz = min(dz,1.0-dz)

          d  = dx*dx+dy*dy+dz*dz
          If(d <= r2 ) Then
             nbridge1 = nbridge1 + 1
          End If
       End Do
    End Do

    Allocate(bridge1(2,nbridge1))

    ind = 1
    Do li = 1, nloc(v1)
       Do ri = 1, nrecv(v1)
          dx = abs(pos1(1,li)-posRec(1,ri))
          !          dx = min(dx,1.0-dx)
          dy = abs(pos1(2,li)-posRec(2,ri))
          !          dy = min(dy,1.0-dy)
          dz = abs(pos1(3,li)-posRec(3,ri))
          !          dz = min(dz,1.0-dz)

          d  = dx*dx+dy*dy+dz*dz
          If(d <= r2 ) Then
             bridge1(1,ind) = li
             bridge1(2,ind) = ri
             ind = ind + 1
          End If
       End Do
    End Do

    Deallocate(posRec)


    ! Envoi vers l'arriere, reception venant de l'avant
    Allocate (posRec(3,nrecv(v2)))
    Call Mpi_Irecv(posRec,3*nrecv(v2),Mpi_Real,neighbours(v2),2,Mpi_Comm_World,mpireqr1,mpierr)
    Call Mpi_ISend(pos1,3*nloc(v1), Mpi_Real,neighbours(v1),2,Mpi_Comm_World,mpireqs1,mpierr)

    Call Mpi_Wait(mpireqs1,mpistat,mpierr)
    Call Mpi_Wait(mpireqr1,mpistat,mpierr)

    nbridge2 = 0
    Do li = 1, nloc(v2)
       Do ri = 1, nrecv(v2)
          dx = abs(pos2(1,li)-posRec(1,ri))
          !          dx = min(dx,1.0-dx)
          dy = abs(pos2(2,li)-posRec(2,ri))
          !          dy = min(dy,1.0-dy)
          dz = abs(pos2(3,li)-posRec(3,ri))
          !          dz = min(dz,1.0-dz)

          d  = dx*dx+dy*dy+dz*dz
          If(d <= r2 ) Then
             nbridge2 = nbridge2+1
          End If
       End Do
    End Do

    Allocate(bridge2(2,nbridge2))

    ind = 1
    Do li = 1, nloc(v2)
       Do ri = 1, nrecv(v2)
          dx = abs(pos2(1,li)-posRec(1,ri))
          !          dx = min(dx,1.0-dx)
          dy = abs(pos2(2,li)-posRec(2,ri))
          !          dy = min(dy,1.0-dy)
          dz = abs(pos2(3,li)-posRec(3,ri))
          !          dz = min(dz,1.0-dz)

          d  = dx*dx+dy*dy+dz*dz
          If(d <= r2 ) Then
             bridge2(1,ind) = li
             bridge2(2,ind) = ri
             ind = ind + 1
          End If
       End Do
    End Do

    Deallocate(posRec)


  End Subroutine findbridge


  Subroutine computeminmax()
    
    Use modvariables, only : PR, pos, partmin, partmax
    Implicit none

    Integer(kind=4) :: i
    Real(kind=PR), dimension(3) :: locmin, locmax


    locmin(:) = 1.e10
    locmax(:) = -1.e10

!!$    Do i=1, mynpart
!!$       Do j=1, 3
!!$          If(pos(j,i) < locmin(j)) locmin(j) = pos(j,i)
!!$          If(pos(j,i) > locmax(j)) locmax(j) = pos(j,i)
!!$       End Do
!!$    End Do

    Do i=1, 3
       locmin(i) = minval(pos(i,:))
       locmax(i) = maxval(pos(i,:))
    End Do
      
    Call Mpi_Allreduce(locmin, partmin, 3, Mpi_Real, Mpi_Min, Mpi_Comm_World, mpierr)
    Call Mpi_Allreduce(locmax, partmax, 3, Mpi_Real, Mpi_Max, Mpi_Comm_World, mpierr)

  End Subroutine computeminmax

End Module modfofpara
