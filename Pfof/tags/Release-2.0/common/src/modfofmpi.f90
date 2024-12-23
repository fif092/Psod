!> @file
!!This file contains the initialisation and the communication functions that merge the halo parts of haloes that extend across several processes.

!> This module contains the initialisation and the communication functions that merge the halo parts of haloes that extend across several processes.
!>
!> Authors: F. Roy, V. Bouillot

! ======================================================================
Module modfofmpi

  Use modconstant, only : PRI

  Implicit none

  Real   (kind=4),   dimension(:,:), allocatable :: posAva, posArr, posBas, posDro, posGau, posHau
  Integer(kind=PRI), dimension(:),   allocatable :: strAva, strArr, strBas, strDro, strGau, strHau, strRec
  Integer(kind=4)                                :: nbridgearr, nbridgeava, nbridgegau
  Integer(kind=4)                                :: nbridgedro, nbridgebas, nbridgehau
  Integer(kind=PRI), dimension(:,:), allocatable :: bridgeArr, bridgeAva, bridgeGau
  Integer(kind=PRI), dimension(:,:), allocatable :: bridgeDro, bridgeBas, bridgeHau
  Integer(kind=4),   dimension(6)                :: nflagrecv
  Integer(kind=4),   dimension(6)                :: nflagloc
  Integer(kind=1),   dimension(:),   allocatable :: border

  Private

  Public :: border, nflagloc

  Public :: mergehaloes
  

Contains

  !> This routine merges haloes that extends across several process by setting their halo ID to the same value.
  Subroutine mergehaloes(perco2, mpicomm, periodic)
    
    Use modmpicom

    Implicit none

    Real(kind=4), intent(in) :: perco2  !< Value of the percolation length to the square in normalized units (where the box length in equal to 1.0) 
    Integer, intent(in) :: mpicomm      !< MPI communicator used in this subroutine
    Logical, intent(in) :: periodic     !< .true. for standard pFoF and .false. for cone version

    Integer(kind=4) :: f, fi, fp
    Integer(kind=4) :: fiID, fpID
    Integer(kind=4) :: mpistatf(MPI_STATUS_SIZE,6)
    Integer(kind=4), dimension(6) :: mpireqsf, mpireqrf

    Call initmerging()

    ! Processes exchange positions of the particles close to the faces
    Do f=1,3 ! loop over the faces, we deal with 2 faces at a time
       fi = 2*f-1
       fp = 2*f

       fiID = neighbours(fi)
       fpID = neighbours(fp) 

       Call Mpi_Irecv(nflagrecv(fi), 1, Mpi_Integer, fiID, 1,   mpicomm, mpireqrf(fi), mpierr)
       Call Mpi_Irecv(nflagrecv(fp), 1, Mpi_Integer, fpID, 2,   mpicomm, mpireqrf(fp), mpierr)
       Call Mpi_Isend(nflagloc(fp) , 1, Mpi_Integer, fpID, 1, mpicomm, mpireqsf(fp), mpierr)
       Call Mpi_Isend(nflagloc(fi) , 1, Mpi_Integer, fiID, 2, mpicomm, mpireqsf(fi), mpierr)

    End Do

    Call Mpi_Waitall(6,mpireqrf,mpistatf,mpierr)
    Call Mpi_Waitall(6,mpireqsf,mpistatf,mpierr)

    ! find the pairs of particles separated by a face that are distant by less than the percolation length from each other
    Call findbridge(nflagloc, nflagrecv, posArr, posAva, nbridgearr, bridgeArr, nbridgeava, bridgeAva, &
         perco2, 1, 2, mpicomm, periodic)
    If(Allocated(posAva)) Deallocate(posAva)
    If(Allocated(posArr)) Deallocate(posArr)

    Call findbridge(nflagloc, nflagrecv, posGau, posDro, nbridgegau, bridgeGau, nbridgedro, bridgeDro, &
         perco2, 3, 4, mpicomm, periodic)
    If(Allocated(posGau)) Deallocate(posGau)
    If(Allocated(posDro)) Deallocate(posDro)

    Call findbridge(nflagloc, nflagrecv, posBas, posHau, nbridgebas, bridgeBas, nbridgehau, bridgeHau, &
         perco2, 5, 6, mpicomm, periodic)
    If(Allocated(posBas)) Deallocate(posBas)
    If(Allocated(posHau)) Deallocate(posHau)


    Call setcommonhaloid(mpicomm)

  End Subroutine mergehaloes


  ! ======================================================================
  !> This routine initializes the merging of halo parts. If the flag 
  !! value is not 0, the position and halo ID of the particle are saved
  !! in specific arrays (one array per face).
  Subroutine initmerging()

    Use modio, only : mynpart
    Use modvariables, only : pos, structure

    Implicit None

    Integer(kind=4) :: debh, debb, deba, debr, debg, debd
    Integer(kind=4) :: i

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

    ! loop over the particles
    ! if the particle is located near a face, its position and its halo ID are saved in specific arrays
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

    Deallocate(border)

  End Subroutine initmerging


  ! ======================================================================
  !> Finds the pairs of particles that should be linked through the merging procedure. <br>
  !! two directions are considered simultaneously by the local process: <br>
  !! 1= back  2= front <br>
  !! 1= left  2= right <br>
  !! 1= bottom  2= top <br>
  !! meaning that the local process exchanges information with the processes located in these two directions.
  Subroutine findbridge(nloc, nrecv, pos1, pos2, nbridge1, bridge1, nbridge2, bridge2, r2, v1, v2, mpicomm, periodic)

    Use modmpicom
    Use modconstant
    Implicit None
    ! input
    Integer(kind=4), intent(in), dimension(6) :: nloc   !< local number of particles flagged for link detection in direction 1
    Integer(kind=4), intent(in), dimension(6) :: nrecv  !< local number of particles flagged for link detection in direction 2
    Integer(kind=4), intent(in) :: v1                   !< neighbour index in direction 1
    Integer(kind=4), intent(in) :: v2                   !< neighbour index in direction 2
    Real(kind=4), intent(in), dimension(3,*) :: pos1    !< positions of local particles flagged for link detection in direction 1
    Real(kind=4), intent(in), dimension(3,*) :: pos2    !< positions of distant particles flagged for link detection in direction 2
    Real(kind=4), intent(in) :: r2                      !< percolation length
    Integer, intent(in) :: mpicomm                      !< MPI communicator used in this subroutine
    Logical, intent(in) :: periodic                     !< .true. for standard pFoF and .false. for cone version

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


    ! send flagged positions to direction 2 and receive from direction 1
    Allocate (posRec(3,nrecv(v1)))
    Call Mpi_Irecv(posRec,3*nrecv(v1),Mpi_Real,neighbours(v1),1,mpicomm,mpireqr1,mpierr)
    Call Mpi_ISend(pos2,3*nloc(v2), Mpi_Real,neighbours(v2),1,mpicomm,mpireqs1,mpierr)

    Call Mpi_Wait(mpireqs1,mpistat,mpierr)
    Call Mpi_Wait(mpireqr1,mpistat,mpierr)
    
    ! loop over local flagged positions for direction 1 and received positions from direction 1
    nbridge1 = 0
    Do li = 1, nloc(v1)
       Do ri = 1, nrecv(v1)
          dx = abs(pos1(1,li)-posRec(1,ri))
          dy = abs(pos1(2,li)-posRec(2,ri))
          dz = abs(pos1(3,li)-posRec(3,ri))
          If(periodic) Then
             dx = min(dx,1.0-dx)
             dy = min(dy,1.0-dy)
             dz = min(dz,1.0-dz)
          End If
          d  = dx*dx+dy*dy+dz*dz
          ! count the number of bridges
          If(d <= r2 ) Then
             nbridge1 = nbridge1 + 1
          End If
       End Do
    End Do

    Allocate(bridge1(2,nbridge1))

    ! we have to loop again to save the local ID of the particles that constitute the bridges
    ind = 1
    Do li = 1, nloc(v1)
       Do ri = 1, nrecv(v1)
          dx = abs(pos1(1,li)-posRec(1,ri))
          dy = abs(pos1(2,li)-posRec(2,ri))
          dz = abs(pos1(3,li)-posRec(3,ri))
          If(periodic) Then
             dx = min(dx,1.0-dx)
             dy = min(dy,1.0-dy)
             dz = min(dz,1.0-dz)
          End If
          d  = dx*dx+dy*dy+dz*dz
          If(d <= r2 ) Then
             bridge1(1,ind) = li
             bridge1(2,ind) = ri
             ind = ind + 1
          End If
       End Do
    End Do

    Deallocate(posRec)


    ! send position to direction 1 and receive from direction 2
    Allocate (posRec(3,nrecv(v2)))
    Call Mpi_Irecv(posRec,3*nrecv(v2),Mpi_Real,neighbours(v2),2,mpicomm,mpireqr1,mpierr)
    Call Mpi_ISend(pos1,3*nloc(v1), Mpi_Real,neighbours(v1),2,mpicomm,mpireqs1,mpierr)

    Call Mpi_Wait(mpireqs1,mpistat,mpierr)
    Call Mpi_Wait(mpireqr1,mpistat,mpierr)

    ! same loops with positions corresponding with the second direction
    nbridge2 = 0
    Do li = 1, nloc(v2)
       Do ri = 1, nrecv(v2)
          dx = abs(pos2(1,li)-posRec(1,ri))
          dy = abs(pos2(2,li)-posRec(2,ri))
          dz = abs(pos2(3,li)-posRec(3,ri))
          If(periodic) Then
             dx = min(dx,1.0-dx)
             dy = min(dy,1.0-dy)
             dz = min(dz,1.0-dz)
          End If
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
          dy = abs(pos2(2,li)-posRec(2,ri))
          dz = abs(pos2(3,li)-posRec(3,ri))
          If(periodic) Then
             dx = min(dx,1.0-dx)
             dy = min(dy,1.0-dy)
             dz = min(dz,1.0-dz)
          End If
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


  ! ======================================================================
  !> Exchanges the haloID of the particles located near the faces and edges and set the same haloID
  !! to each particles that are in the same halo.
  Subroutine setcommonhaloid(mpicomm)

    Use modconstant
    Use modmpicom
    Use modvariables, only : structure
    Use modio, only : mynpart

    Implicit none

    Integer, intent(in) :: mpicomm   !< MPI communicator used in this subroutine

    Integer(kind=4) :: nbpassage
    Logical :: nopermut, finished
    Integer(kind=4) :: sendID, recvID
    Integer(kind=4) :: ib, k
    Integer(kind=PRI) :: il, ir
    Integer(kind=4) :: mpistat(MPI_STATUS_SIZE)   ! MPI comm. status


    nbpassage = 0

    If(procID==0) Print *,'Beginning of the merging phase...'

    raccordement : Do
       nopermut = .true.
       nbpassage = nbpassage + 1
       If(ProcID == 0) Print *,'Loop number ',nbpassage

       !------------------------------------------------------!
       ! Envoi du bas vers le bas et reception venant du haut !
       !------------------------------------------------------!

       sendID = neighbours(5)
       recvID = neighbours(6)
      
       Allocate (strRec(nflagrecv(6)))

       Call Mpi_ISend(strBas,nflagloc(5), MPI_PRI,sendID,1,mpicomm,mpireqs1,mpierr)
       Call Mpi_IRecv(strRec,nflagrecv(6),MPI_PRI,recvID,1,mpicomm,mpireqr1,mpierr)
       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

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


       !-------------------------------------------------------------------
       ! Envoi de ma frontiere haut vers le haut et reception venant du bas
       !-------------------------------------------------------------------

       sendID = neighbours(6)
       recvID = neighbours(5)
       Allocate (strRec(nflagrecv(5)))

       Call Mpi_ISend(strHau,nflagloc(6), MPI_PRI,sendID,2,mpicomm,mpireqs1,mpierr)
       Call Mpi_IRecv(strRec,nflagrecv(5),MPI_PRI,recvID,2,mpicomm,mpireqr1,mpierr)
       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

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


       !--------------------------------------------------------------------------
       ! Envoi de ma frontiere droite vers la droite et reception venant de gauche
       !--------------------------------------------------------------------------

       sendID = neighbours(4)
       recvID = neighbours(3)

       Allocate (strRec(nflagrecv(3)))

       Call Mpi_ISend(strDro,nflagloc(4), MPI_PRI,sendID,3,mpicomm,mpireqs1,mpierr)
       Call Mpi_IRecv(strRec,nflagrecv(3),MPI_PRI,recvID,3,mpicomm,mpireqr1,mpierr)
       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

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


       !---------------------------------------------------------------------
       ! Envoi de ma face gauche vers la gauche et reception venant de droite
       !---------------------------------------------------------------------

       sendID = neighbours(3)
       recvID = neighbours(4)

       Allocate (strRec(nflagrecv(4)))

       Call Mpi_ISend(strGau,nflagloc(3), MPI_PRI,sendID,4,mpicomm,mpireqs1,mpierr)
       Call Mpi_IRecv(strRec,nflagrecv(4),MPI_PRI,recvID,4,mpicomm,mpireqr1,mpierr)
       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

       ! comparaison des part. proche de ma face 'droite' avec les part. venant de mon neighbours de droite
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


       !--------------------------------------------------------------------------
       ! Envoi de ma frontiere avant vers l'avant et reception venant de l'arriere
       !--------------------------------------------------------------------------

       sendID = neighbours(2)
       recvID = neighbours(1)
       
       Allocate (strRec(nflagrecv(1)))


       Call Mpi_ISend(strAva,nflagloc(2), MPI_PRI,sendID,5,mpicomm,mpireqs1,mpierr)
       Call Mpi_IRecv(strRec,nflagrecv(1),MPI_PRI,recvID,5,mpicomm,mpireqr1,mpierr)
       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

       ! comparaison des part. proche de ma face 'arriere' avec les part. venant de mon neighbours de l'arriere
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


       !---------------------------------------------------------------------
       ! Envoi de ma face arriere vers l'arriere et reception venant de l'avant
       !---------------------------------------------------------------------

       sendID = neighbours(1)
       recvID = neighbours(2)

       Allocate (strRec(nflagrecv(2)))


       Call Mpi_ISend(strArr,nflagloc(1), MPI_PRI,sendID,6,mpicomm,mpireqs1,mpierr)
       Call Mpi_IRecv(strRec,nflagrecv(2),MPI_PRI,recvID,6,mpicomm,mpireqr1,mpierr)
       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

       ! comparaison des part. proche de ma face 'avant' avec les part. venant de mon neighbours de l'avant
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

       Call Mpi_AllReduce(nopermut,finished,1,Mpi_Logical,Mpi_Land,mpicomm,mpierr)
       If(finished) Exit
       
    End Do raccordement

    If(procID==0) Then
       Print *,''
       Print *,'End of the permutations'
       Print *,''
    End If

    If(Allocated(strHau)) Deallocate(strHau)
    If(Allocated(strBas)) Deallocate(strBas)
    If(Allocated(strDro)) Deallocate(strDro)
    If(Allocated(strGau)) Deallocate(strGau)
    If(Allocated(strAva)) Deallocate(strAva)
    If(Allocated(strArr)) Deallocate(strArr) 


  End Subroutine setcommonhaloid


End Module modfofmpi
