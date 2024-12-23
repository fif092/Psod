!> @file
!! This file contains subroutines used to gather haloes on the differents processes and compute some observables.

!> This module contains subroutines used to gather haloes on the differents processes and compute some observables.
!>
!> Authors: F. Roy, V. Bouillot

Module modhalo

  Use modconstant, only : PRI

  Implicit none

  Integer(kind=4) :: haloNB        !< local number of haloes
  Integer(kind=4) :: haloNB_all    !< global number of haloes
  Integer(kind=4) :: halopartNB    !< number of particles belonging to one of the local haloes
  Integer(kind=4) :: myfinalnpart  !< number of particles after the particles have been exchanged to gather them by haloes
  Integer(kind=4),   dimension(:), allocatable :: haloMass       !< mass of the haloes
  Integer(kind=PRI), dimension(:), allocatable :: haloID         !< ID of the haloes
  Integer(kind=PRI), dimension(:), allocatable :: halopartID     !< ID of the particles belonging to a halo
  Integer(kind=PRI), dimension(:), allocatable :: halopartramsesID     !< RAMSES ID of the particles belonging to a halo, used only for lightcone halo
  Integer(kind=4),   dimension(:), allocatable :: haloSubHaloNB  !< number of subhalo for each halo

  Real(kind=4), dimension(:,:), allocatable :: halopartPos    !< position of the particles belonging to one of the local haloes
  Real(kind=4), dimension(:,:), allocatable :: halopartVel    !< velocity of the particles belonging to one of the local haloes
  Real(kind=4), dimension(:),   allocatable :: halopartPot    !< potential of the particles belonging to one of the local haloes
  Real(kind=4), dimension(:,:), allocatable :: halopartFor    !< force on the particles belonging to one of the local haloes
  Real(kind=8), dimension(:,:), allocatable :: halocomPos     !< position of the center of mass of each halo
  Real(kind=8), dimension(:,:), allocatable :: halocomVel     !< velocity of the center of mass of each halo
  Real(kind=8), dimension(:),   allocatable :: haloRadius     !< radius of each halo

Contains

  !=======================================================================
  !> Exchange the particles so that particles belonging to one halo are gathered on the same process.
  Subroutine gatherhaloes(mpicomm, param)

    Use modvariables
    Use modmpicom
    Use modtiming

    Implicit None

    Integer, intent(in) :: mpicomm   !< MPI communicator used for the communications
    Class(Type_parameter_pfof), intent(in) :: param

    Integer(kind=PRI), dimension(:), allocatable :: strSend ! structure ID array, and tmp array used for comm
    Integer(kind=PRI), dimension(:), allocatable :: idSend  ! particle ID array and tmp array used for comm
    Integer(kind=PRI), dimension(:), allocatable :: ramsesidSend
    Real(kind=4), dimension(:,:), allocatable :: posSend, velSend, forSend
    Real(kind=4), dimension(:), allocatable :: potSend
    Integer(kind=4) :: strPID
    Integer(kind=4), dimension(:),allocatable :: strPIDvec, strPIDvecloc  ! array of particle nb by process after
    Integer(kind=4) :: sendID, recvID   ! ID for process to send to and to receive from
    Integer(kind=4) :: i, ind, iproc
    Integer(kind=4) :: allocStat
    Integer(kind=4) :: mpistat(MPI_STATUS_SIZE)   ! MPI comm. status
    Integer(kind=4) :: nbrec, nbsend    ! nb of elements to recv and to send
    Integer(kind=PRI) :: strNBbeg, strNBend
    Integer(kind=PRI) :: smin,smax
    Integer(kind=PRI) :: tmpdi
    Integer(kind=4) :: NPparProc
    Integer(kind=4) :: recvpoint

    Logical(kind=4) :: do_read_ramses_part_id

    Select Type (param)
    Type is (Type_parameter_pfof_snap)
       do_read_ramses_part_id = .false.
    Type is (Type_parameter_pfof_cone)
       do_read_ramses_part_id = param%do_read_ramses_part_id
    End Select

    tmpdi = npart/procNB
    NPparProc = int(tmpdi, kind=PRI)
    If(mod(npart,int(procNB,kind=PRI)) /= 0) Then
       NPparProc = NPparProc + 1
    End If
       
    strNBbeg = int(NPparProc,kind=8) * int(procID,kind=8) + 1    ! id min de halo sur le process courant
    strNBend = int(NPparProc,kind=8) * int((procID + 1),kind=8)  ! id max de halo sur le process courant
    
    Allocate(strPIDvecloc(procNB), strPIDvec(procNB))
    strPIDvecloc = 0  ! nb de particules par process
    strPIDvec = 0

    Do i = 1,mynpart
       strPID = int((structure(i)-1) / NPparProc, kind=4) + 1  ! id du process ou se trouvera la particule
       strPIDvecloc(strPID) = strPIDvecloc(strPID) + 1         ! on incremente le nb de particules sur ce process
    End Do
    
    Call Mpi_Allreduce(strPIDvecloc,strPIDvec,procNB,Mpi_Integer,Mpi_Sum,Mpi_Comm_World,mpierr)

    myfinalnpart = strPIDvec(procID+1)             ! nb de particules sur le process courant apres redistribution selon les halos
    Allocate(posf(3,myfinalnpart),STAT=allocStat)      
    If(allocStat > 0) Call EmergencyStop('Allocate failed for posf in gatherhaloes',2)
    Allocate(velf(3,myfinalnpart),STAT=allocStat)
    If(allocStat > 0) Call EmergencyStop('Allocate failed for velf in gatherhaloes',2)
    Allocate(stf(myfinalnpart),STAT=allocStat)
    If(allocStat > 0) Call EmergencyStop('Allocate failed for stf in gatherhaloes',2)
    Allocate(idf(myfinalnpart),STAT=allocStat)
    If(allocStat > 0) Call EmergencyStop('Allocate failed for idf in gatherhaloes',2)
    If(param%do_read_potential) Then
       Allocate(potf(myfinalnpart),STAT=allocStat)
       If(allocStat > 0) Call EmergencyStop('Allocate failed for potf in gatherhaloes',2)
    End If
    If(param%do_read_gravitational_field) Then
       Allocate(forf(3,myfinalnpart),STAT=allocStat)
       If(allocStat > 0) Call EmergencyStop('Allocate failed for forf in gatherhaloes',2)
    End If
    If(do_read_ramses_part_id) Then
       Allocate(ramsesidf(myfinalnpart),STAT=allocStat)
       If(allocStat > 0) Call EmergencyStop('Allocate failed for ramsesidf in gatherhaloes',2)
    End If

    recvpoint = 1
    procMasse : Do iproc = 1,procNB-1
       If(procID==0) Write(*,'(A,I5,A)') 'Gatherhalo permutation ',iproc,'...'

       sendID = mod(procID + iproc,procNB)
       recvID = mod(procID + procNB - iproc,procNB)
       nbsend = strPIDvecloc(sendID+1)

       Call Mpi_ISend(nbsend,1,Mpi_Integer,sendID,sendID,mpicomm,mpireqs1,mpierr)
       Call Mpi_IRecv(nbrec, 1,Mpi_Integer,recvID,procID,mpicomm,mpireqr1,mpierr)

       Call Mpi_Wait(mpireqs1,mpistat,mpierr)
       Call Mpi_Wait(mpireqr1,mpistat,mpierr)

       If(nbsend /= 0) Then
          Allocate(strSend(nbsend),posSend(3,nbsend),velSend(3,nbsend),idSend(nbsend))
          If(param%do_read_potential) Allocate(potSend(nbsend))
          If(param%do_read_gravitational_field) Allocate(forSend(3,nbsend))
          If(do_read_ramses_part_id) Allocate(ramsesidSend(nbsend))

          smin = int(NPparProc,kind=8) * int(sendID,kind=8) + 1
          smax = int(NPparProc,kind=8) * int((sendID+1),kind=8)

          ind=1
          Do i=1, mynpart
             If(structure(i)>= smin .and. structure(i)<= smax) Then
                strSend(ind) = structure(i)
                posSend(:,ind) = pos(:,i)
                velSend(:,ind) = vel(:,i)
                If(param%do_read_potential) potSend(ind) = pot(i)
                If(param%do_read_gravitational_field) forSend(:,ind) = for(:,i)
                If(do_read_ramses_part_id) ramsesidSend(ind) = ramsesid(i)
                idSend(ind) = id(i)
                ind = ind+1
             End If
          End Do

          If(ind /= nbsend +1 ) Then
             write(*,*)ind,nbsend,sendID,smin,smax,strPIDvec
             Call EmergencyStop('Error  1 while sharing structures for output.',2)
          End If

          Call Mpi_Isend(strSend,nbsend,  MPI_PRI, sendID,1,mpicomm,mpireqs1,mpierr)
          Call Mpi_Isend(posSend,3*nbsend,Mpi_Real,sendID,2,mpicomm,mpireqs2,mpierr)
          Call Mpi_Isend(velSend,3*nbsend,Mpi_Real,sendID,3,mpicomm,mpireqs3,mpierr)
          Call Mpi_Isend(idSend, nbsend,  MPI_PRI, sendID,4,mpicomm,mpireqs4,mpierr)
          If(param%do_read_potential) Call Mpi_Isend(potSend, nbsend, Mpi_Real, sendID, 5,mpicomm,mpireqs5,mpierr)
          If(param%do_read_gravitational_field) Call Mpi_Isend(forSend,3*nbsend,Mpi_Real,sendID,6,mpicomm,mpireqs6,mpierr)
          If(do_read_ramses_part_id) Call Mpi_Isend(ramsesidSend,nbsend,MPI_PRI,&
               sendID,7,mpicomm,mpireqs7,mpierr)

       End If

       If(nbrec /= 0) Then
          Call Mpi_IRecv(stf(recvpoint), nbrec,  MPI_PRI, recvID,1,mpicomm,mpireqr1,mpierr)
          Call Mpi_IRecv(posf(1,recvpoint),3*nbrec,Mpi_Real,recvID,2,mpicomm,mpireqr2,mpierr)
          Call Mpi_IRecv(velf(1,recvpoint),3*nbrec,Mpi_Real,recvID,3,mpicomm,mpireqr3,mpierr)
          Call Mpi_IRecv(idf(recvpoint), nbrec,  MPI_PRI, recvID,4,mpicomm,mpireqr4,mpierr)
          If(param%do_read_potential) Call Mpi_IRecv(potf(recvpoint), nbrec,  &
               Mpi_Real, recvID,5,mpicomm,mpireqr5,mpierr)
          If(param%do_read_gravitational_field) Call Mpi_IRecv(forf(1,recvpoint),3*nbrec,&
               Mpi_Real,recvID,6,mpicomm,mpireqr6,mpierr)
          If(do_read_ramses_part_id) Call Mpi_Irecv(ramsesidf(recvpoint), nbrec, MPI_PRI, &
               recvID, 7, mpicomm, mpireqr7, mpierr)
          recvpoint=recvpoint+nbrec
       End If


       If(nbsend/=0) Then
          Call Mpi_Wait(mpireqs1,mpistat,mpierr)
          Call Mpi_Wait(mpireqs2,mpistat,mpierr)
          Call Mpi_Wait(mpireqs3,mpistat,mpierr)
          Call Mpi_Wait(mpireqs4,mpistat,mpierr)
          Deallocate(strSend,posSend,idSend,velSend)

          If(param%do_read_potential) Then
             Call Mpi_Wait(mpireqs5,mpistat,mpierr)
             Deallocate(potSend)
          End If
          If(param%do_read_gravitational_field) Then
             Call Mpi_Wait(mpireqs6,mpistat,mpierr)
             Deallocate(forSend)
          End If
          If(do_read_ramses_part_id) Then
             Call Mpi_Wait(mpireqs7, mpistat, mpierr)
             Deallocate(ramsesidSend)
          End If
          
       End If
       If(nbrec/=0) Then
          Call Mpi_Wait(mpireqr1,mpistat,mpierr)
          Call Mpi_Wait(mpireqr2,mpistat,mpierr)
          Call Mpi_Wait(mpireqr3,mpistat,mpierr)
          Call Mpi_Wait(mpireqr4,mpistat,mpierr)
          If(param%do_read_potential) Call Mpi_Wait(mpireqr5,mpistat,mpierr)
          If(param%do_read_gravitational_field) Call Mpi_Wait(mpireqr6,mpistat,mpierr)
          If(do_read_ramses_part_id) Call Mpi_Wait(mpireqr7, mpistat, mpierr)
       End If

    End Do procMasse


    ind= recvpoint
    Do i=1, mynpart
       If(structure(i)>= strNBbeg .and. structure(i)<= strNBend) Then
          stf(recvpoint)  = structure(i)
          posf(:,recvpoint) = pos(:,i)
          velf(:,recvpoint) = vel(:,i)
          idf(recvpoint)  = id(i)
          If(param%do_read_potential) potf(recvpoint) = pot(i)
          If(param%do_read_gravitational_field) forf(:,recvpoint) = for(:,i)
          If(do_read_ramses_part_id) ramsesidf(recvpoint) = ramsesid(i)
          recvpoint = recvpoint+1
       End If
    End Do

    Deallocate(structure)
    Deallocate(strPIDvecloc,strPIDvec)
    Deallocate(pos, vel, id)
    If(param%do_read_potential) Deallocate(pot)
    If(param%do_read_gravitational_field) Deallocate(for)
    If(do_read_ramses_part_id) Deallocate(ramsesid)

    If( recvpoint/= myfinalnpart +1 ) Then
       Print *, procID,' recvpoint=',recvpoint
       Print *, procID,' mynpartw=',myfinalnpart
       Call EmergencyStop('Error while gathering haloes in gatherhaloes.',2)
    End If


  End Subroutine gatherhaloes


  ! ======================================================================
  ! Select haloes whose mass is >= Mmin
  Subroutine selecthaloes(param)

    Use modvariables
    Use modmpicom, only : procID

    Implicit none

    Class(Type_parameter_pfof), intent(in) :: param

    Integer(kind=8) :: hidmin, hidmax
    Integer(kind=4) :: nbhid
    Integer(kind=4) :: hindex
    Integer(kind=4) :: li, lp, fi, fp, h, i
    Integer(kind=4), dimension(:), allocatable :: halomasstmp

    Logical(kind=4) :: do_read_ramses_part_id

    Select Type (param)
    Class is (Type_parameter_pfof_snap)
       do_read_ramses_part_id = .false.
    Class is (Type_parameter_pfof_cone)
       do_read_ramses_part_id = param%do_read_ramses_part_id
    End Select

#ifdef DEBUG
    Print *, 'Process ',procID, ' enters selecthaloes'
#endif
    
    If(myfinalnpart /= 0) Then
       hidmin = stf(1)
       hidmax = stf(myfinalnpart)
    Else
       hidmin=0
       hidmax=-1
    End If
    nbhid = int(hidmax - hidmin + 1, kind=4)

#ifdef DEBUG
    Print *, 'hidmax=',hidmax,' ; hidmin=',hidmin,' ; nbhid=',nbhid, ' ; myfinalnpart=',myfinalnpart 
#endif

    Allocate(halomasstmp(nbhid))

    halomasstmp = 0

    Do i=1, myfinalnpart
       hindex = int(stf(i) - hidmin + 1,kind=4)
       halomasstmp(hindex) = halomasstmp(hindex) + 1
    End Do

    haloNB = 0
    halopartNB = 0

    ! Compute total nb of particles in halo with M >= Mmin and nb of halos with M >= Mmin
    Do i = 1, nbhid
       If(halomasstmp(i) >= param%mmin) Then
          haloNB = haloNB + 1
          halopartNB = halopartNB + halomasstmp(i)
       End If
    End Do

    ! Keep positions, velocities and id for particles in halo with M >= Mmin, and potential if requested
    Allocate(halopartPos(3,halopartNB))
    Allocate(halopartVel(3,halopartNB))
    Allocate(halopartID(halopartNB))
    If(param%do_read_potential) Allocate(halopartPot(halopartNB))
    If(param%do_read_gravitational_field) Allocate(halopartFor(3,halopartNB))
    If(do_read_ramses_part_id) Allocate(halopartramsesid(halopartNB))
    ! Keep mass and id for halos with M >= Mmin
    If(haloNB==0) Then
       Allocate(haloMass(1))
       Allocate(haloID(1))
    Else
       Allocate(haloMass(haloNB))
       Allocate(haloID(haloNB))
    End If
    ! Sub-halo detection is not implemented yet
    ! We allocate halosubhaloNB with a size=1
    Allocate(halosubhaloNB(1))


    fp = 1
    lp = 0
    fi = 1
    li = 0
    h = 1
    Do i = 1, nbhid
       If(halomasstmp(i) >= param%mmin) Then
          lp = fp + halomasstmp(i) - 1
          li = fi + halomasstmp(i) - 1
          halopartPos(:,fp:lp) = posf(:,fi:li)
          halopartVel(:,fp:lp) = velf(:,fi:li)
          If(param%do_read_potential) halopartPot(fp:lp) = potf(fi:li)
          If(param%do_read_gravitational_field) halopartFor(:,fp:lp) = forf(:,fi:li)
          If(do_read_ramses_part_id) halopartramsesid(fp:lp) = ramsesidf(fi:li)
          halopartID(fp:lp) = idf(fi:li)
          haloMass(h) = halomasstmp(i)
          haloID(h) = stf(fi)
          fp = lp + 1
          fi = li + 1
          h = h + 1
       Else
          li = fi + halomasstmp(i) - 1
          fi = li + 1
       End If
    End Do
    If(li /= myfinalnpart) Then
       Print *, 'Error when keeping particles positions in halo with M >= Mmin on process ', procID
    End If
    If(lp /= halopartNB) Then
       Print *, 'Error when keeping particles positions in halo with M >= Mmin on process ', procID
    End If


    If(Allocated(posf)) Deallocate(posf)
    If(Allocated(velf)) Deallocate(velf)
    If(Allocated(idf)) Deallocate(idf)
    If(Allocated(halomasstmp)) Deallocate(halomasstmp)
    If(Allocated(potf)) Deallocate(potf)
    If(Allocated(forf)) Deallocate(forf)
    If(Allocated(ramsesidf)) Deallocate(ramsesidf)
  End Subroutine selecthaloes


  ! ======================================================================
  ! Computes the position and the velocity of the center of mass for each halo
  Subroutine computecom(periodic)

    Use modmpicom

    Implicit none

    Logical, intent(in) :: periodic

    Integer(kind=4) :: fi, li, h, i, j
    Integer(kind=4) :: halom
    Real(kind=8), dimension(3) :: delta

    If(haloNB==0) Then
       Allocate(halocomPos(3,1))
       Allocate(halocomVel(3,1))
    Else
       Allocate(halocomPos(3,haloNB))
       Allocate(halocomVel(3,haloNB))
    End If

    halocomPos = 0.d0
    halocomVel = 0.d0

    fi = 1
    Do h = 1, haloNB
       halom = 1
       li = fi + haloMass(h) - 1
       halocomPos(:,h) = halopartPos(:,fi)
       halocomVel(:,h) = halopartVel(:,fi)
       Do i = fi+1, li
          delta(:) = halocomPos(:,h) / halom - halopartPos(:,i)
          If(periodic) Then
             Do j = 1, 3
                If(abs(delta(j)) > 0.5d0) Then
                   If(delta(j) > 0d0) halocomPos(j,h) = halocomPos(j,h) + 1.0d0
                   If(delta(j) < 0d0) halocomPos(j,h) = halocomPos(j,h) - 1.0d0
                End If
             End Do
          End If
          halom = halom + 1
          halocomPos(:,h) = halocomPos(:,h) + halopartPos(:,i)
          halocomVel(:,h) = halocomVel(:,h) + halopartVel(:,i)
       End Do
       fi = li + 1
       
    End Do

    Do h = 1, haloNB
       halocomPos(:,h) = halocomPos(:,h) / haloMass(h)
       halocomVel(:,h) = halocomVel(:,h) / haloMass(h)
       If(periodic) Then
          Do j = 1, 3
             If(halocomPos(j,h) > 1.0d0) halocomPos(j,h) = halocomPos(j,h) - 1.0d0
             If(halocomPos(j,h) < 0.0d0) halocomPos(j,h) = halocomPos(j,h) + 1.0d0
          End Do
       End If
    End Do
    
    Call Mpi_AllReduce(haloNB,haloNB_all,1,Mpi_Integer,Mpi_Sum,Mpi_Comm_World,mpierr)

  End Subroutine computecom

  
  ! ======================================================================
  Subroutine computeradius(periodic)

    Implicit None
    
    Logical, intent(in) :: periodic

    Integer(kind=4) :: ib
    Integer(kind=4) :: id
    Integer(kind=4) :: ih
    Integer(kind=4) :: ip

    Real(kind=8) :: rmax2
    Real(kind=8) :: d2
    Real(kind=8), dimension(3) :: delta
    
    If(haloNB==0) Then
       Allocate(haloRadius(1))
    Else
       Allocate(haloRadius(haloNB))
    End If

    ib = 0
    Do ih = 1, haloNB
       rmax2 = 0.d0
       Do ip = ib+1, ib+haloMass(ih)
          delta(:) = halopartPos(:,ip) - haloComPos(:,ih)
          If(periodic) Then
             Do id = 1, 3
                If(abs(delta(id)) > 0.5d0) Then
                   delta(id) = 1.d0 - abs(delta(id))
                End If
             End Do
          End If
          d2 = delta(1)*delta(1) + delta(2)*delta(2) + delta(3)*delta(3) 
          If(d2 > rmax2) rmax2 = d2
       End Do
       ib = ib + haloMass(ih)
       haloRadius(ih) = sqrt(rmax2)
    End Do

    Deallocate(halopartPos, halopartVel)


  End Subroutine computeradius

End Module modhalo
