! Copyright (c) 2007-2016 CNRS, Fabrice Roy, Vincent Bouillot, Edouard Audit
! Author: Fabrice Roy (LUTH/CNRS/PSL), fabrice.roy@obspm.fr
! Vincent Bouillot (LUTH/CNRS/PSL)
! Edouard Audit (Maison de la Simulation) 
! Distributed under the CeCILL 2.0 license. For full terms see the file LICENSE.


!> @file
!! Contains the serial FOF algorithm. 
!!
!! @author Edouard Audit
!! @author Fabrice Roy
!! @author Vincent Bouillot
  
!> Contains the serial FOF algorithm. 
Module modfofpara

  Use modmpicom  

Contains

  !=======================================================================
  !> Serial FOF subroutine
  !> @brief
  !!
  !> Halo detection is performed locally in each subdomains (by each process).<br>
  !! Structures which are cut by a border between two subdomains are gathered.<br>
  !! The criterium to gather two parts of a halo is simple: if one particle in a halo is seperated
  !! from a particle in another halo (on the other side of a border) by a distance smaller than the
  !! percolation parameter, then these two halos are two parts of a single halo.
  Subroutine fofpara()

    Use modfofmpi, only : mergehaloes, border, nflagloc
    Use modparameters
    Use modhalo
    Use modio
    Use modwritehalo
    Use modtiming
    Use modsortinterf
    Use modvariables

    Implicit none

    ! Local variables
    Integer(kind=4), dimension(:),allocatable :: adfirst,npcase
    Integer(kind=4), dimension(:), allocatable :: lamas, pile, indl
    Integer(kind=4) :: i

    Integer(kind=4) :: index, ipart
    Integer(kind=4) :: ipermut, sttmp, ind, ipile, ipilepr
    Integer(kind=4) :: ix, iy, iz, ic
    Integer(kind=4) :: i1, i2, i3, j1, j2, j3
    Integer(kind=4) :: nbborderloc, nbborder
    Integer(kind=4) :: ngrid_fof,res_fof
    Integer(kind=4) :: refine_fof
    Integer(kind=4) :: signx(3), signy(3), signz(3)
    Integer(kind=4) :: nsign
    
    Integer(kind=PRI) :: namas      ! tmp structure ID

    Logical :: periodic

    Real(kind=4) :: r, d, size, size_fof, xx, yy, zz
    Real(kind=4) :: dx, dy, dz
    Real(kind=4) :: r2
    Real(kind=4) :: rT
    Integer(kind=4) :: nh

    ! INITIALIZATION

    ! Initialize timings if necessary
    If(param%do_timings) Then
       Call Mpi_Barrier(Mpi_Comm_World,mpierr)
       time0 = Mpi_Wtime()
    End If

    nflagloc = 0
    nbborderloc = 0

    ! we used a grid fof FoF
    ! this grid is more refined than the coarse grid used for the dynamic evolution
    ! the refinement factor is set to 2
    refine_fof = 2

    res_fof = refine_fof * nres / dims(1)
    ngrid_fof = res_fof**3 
    size = float(nres)
    size_fof = refine_fof * size

    Allocate(adfirst(ngrid_fof),npcase(ngrid_fof))

    Allocate(border(mynpart), lamas(mynpart), indl(mynpart), pile(mynpart), structure(mynpart))

    Call init(res_fof,mynpart,ngrid_fof,pos,pile,adfirst,npcase,xmin,ymin,zmin,size_fof)

    nsign = 2
    If( param%percolation_length*refine_fof >= 0.5 ) Then
       nsign = 3
    End If

    If( param%percolation_length*refine_fof >= 1 ) Then 
       If(procID == 0) Then
          Print *, ' ' 
          Print *, '*********************************************************************************'
          Print *, '*** percolation_length * refine_fof should be < 1                                          ***'
          Print *, '*** refine_fof is set and used in modfofpara.f90, percolation_length is an input parameter ***'
          Print *, '*** Pfof is exiting                                                           ***'
          Print *, '*********************************************************************************'
          Print *, ' ' 
       End If
       Call Mpi_Barrier(Mpi_Comm_World, mpierr)
       Call Mpi_Abort(Mpi_Comm_World, 11, mpierr)
    End If

    r = param%percolation_length
    rT = r / size
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

    If(param%do_timings) Then
       Call Mpi_Barrier(Mpi_Comm_World,mpierr)
       timeInt = Mpi_Wtime()
       tFoFinit = timeInt - time0
    End If

    If(procID==0) Then
       Print *,' '
       Print *,'Beginning of local Friend of Friend...'
       Print *,' '
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

       If(procID==0 .and. mod(i,100000)==0) Print *,'Particule num.',i
       ipart = lamas(i)

       If (abs(pos(1,ipart)-1.0)<1.e-8) pos(1,ipart) = 0.e0
       If (abs(pos(2,ipart)-1.0)<1.e-8) pos(2,ipart) = 0.e0
       If (abs(pos(3,ipart)-1.0)<1.e-8) pos(3,ipart) = 0.e0

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

       ix = int((pos(1,ipart)-xmin)*size_fof)
       iy = int((pos(2,ipart)-ymin)*size_fof)
       iz = int((pos(3,ipart)-zmin)*size_fof)

       signx(2) = -1
       signy(2) = -1
       signz(2) = -1
       If(nsign == 2) Then
          If( (pos(1,ipart)-xmin)*size_fof - ix > 0.5 ) Then
             signx(2) = 1
          End If
          If( (pos(2,ipart)-ymin)*size_fof - iy > 0.5 ) Then
             signy(2) = 1
          End If
          If( (pos(3,ipart)-zmin)*size_fof - iz > 0.5 ) Then
             signz(2) = 1
          End If
       End If
       structure(ipart) = namas

       ! nsign depends on percolation_length * refine_fof: see above
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

    If(procID==0) Print *,'Local FoF ends'

    If(param%do_timings) Then
       Call Mpi_Barrier(Mpi_Comm_World,mpierr)
       tFoFloc = Mpi_Wtime() - timeInt
       timeInt = Mpi_Wtime()
       If(procID==0) Print *,'Process 0 spent ',tFoFloc,'s performing the local FoF halo detection.'
    End If

    Deallocate(lamas,pile,indl)
    Deallocate(adfirst,npcase)

    !-------------------!
    ! FOF local termine !
    !-------------------!

    Call Mpi_Reduce(nbborderloc, nbborder, 1, Mpi_Integer,Mpi_Sum,0,Mpi_Comm_World,mpierr)


    If(procID==0) Then
       Print *, ' '
       Print *, 'Merging procedure...'
       Print *,' '
    End If

    periodic = .true.
    !! Call the merging procedure with the r**2 as an argument
    Call mergehaloes(r2, MPICube, periodic)

    If(param%do_timings) Then
       Call Mpi_Barrier(Mpi_Comm_World,mpierr)
       tRaccord = Mpi_Wtime() - timeInt
       timeInt = Mpi_Wtime()
    End If

    If(procID==0) Print *,'End of the merging procedure.'
    

    !---------------------!
    ! FIN du raccordement !
    !---------------------!

    If(procID==0) Print *,'pFoF gathers particles across the processes according their halo ID'

    Call gatherhaloes(MPICube, param)

    If(param%do_timings) Then
       Call Mpi_Barrier(Mpi_Comm_World,mpierr)
       tGatherhalo = Mpi_Wtime() - timeInt
       If(procID==0) Then
          Print *,' ' 
          Print *,'Process 0 spent ',tGatherhalo,' s gathering particles across the processes'
          Print *,'according their halo ID.'
          Print *,' '
       End If
       timeInt = Mpi_Wtime()
    End If

    If(param%do_read_potential .and. param%do_read_gravitational_field) Then
       Call heapsort(myfinalnpart,stf,posf,velf,forf,potf,idf)
    Else If(param%do_read_potential .and. .not. param%do_read_gravitational_field) Then
       Call heapsort(myfinalnpart,stf,posf,velf,potf,idf)
    Else If(param%do_read_gravitational_field .and. .not. param%do_read_potential) Then
       Call heapsort(myfinalnpart,stf,posf,velf,forf,idf)
    Else
       Call heapsort(myfinalnpart,stf,posf,velf,idf)
    End If

    If(param%do_timings) Then
       Call Mpi_Barrier(Mpi_Comm_World,mpierr)
       tSort = Mpi_Wtime() - timeInt
       If(procID==0) Print *,'Process 0 spent ', tSort,' s sorting particles according their halo ID'
       timeInt = Mpi_Wtime()
    End If

    !! select haloes with mass > Mmin
    Call selecthaloes(param)

    If(param%do_timings) Then
       Call Mpi_Barrier(Mpi_Comm_World,mpierr)
       tSelectHalo = Mpi_Wtime() - timeInt
       If(procID==0) Print *,'Process 0 spent ', tSelectHalo, &
            ' s selecting haloes with Mass > ',param%Mmin
       timeInt = Mpi_Wtime()
    End If

    If(procID==0) Then
       Print *, ' '
       Print *,'pFoF writes halo particles files.'
       Print *, ' '
    End If

    If(param%gatherwrite_factor <= 1) Then
       If(param%do_read_potential .and. param%do_read_gravitational_field) Then
          Call h5writehalopart(param, haloNB, halopartNB, haloMass, haloID, halopartPos, &
               halopartVel, halopartID, halopartFor=halopartFor, halopartPot=halopartPot, &
               inforamses=inforamses)
       Else If(param%do_read_potential .and. .not. param%do_read_gravitational_field ) Then
          Call h5writehalopart(param, haloNB, halopartNB, haloMass, haloID, halopartPos, &
               halopartVel, halopartID, halopartPot=halopartPot,inforamses=inforamses)
       Else If(param%do_read_gravitational_field .and. .not. param%do_read_potential) Then
          Call h5writehalopart(param, haloNB, halopartNB, haloMass, haloID, halopartPos, &
               halopartVel, halopartID, halopartFor=halopartFor, inforamses=inforamses)
       Else
          Call h5writehalopart(param, haloNB, halopartNB, haloMass, haloID, halopartPos, &
               halopartVel, halopartID, inforamses=inforamses)
       End If

    Else
       If(param%do_read_potential .and. param%do_read_gravitational_field) Then
          Call mpih5writehalopart(MpiSubCubeWrite, param, haloNB, halopartNB, haloMass, haloID, halopartPos, &
               halopartVel, halopartID, halopartFor=halopartFor, halopartPot=halopartPot, &
               inforamses=inforamses)
       Else If(param%do_read_potential .and. .not. param%do_read_gravitational_field ) Then
          Call mpih5writehalopart(MpiSubCubeWrite, param, haloNB, halopartNB, haloMass, haloID, halopartPos, &
               halopartVel, halopartID, halopartPot=halopartPot,inforamses=inforamses)
       Else If(param%do_read_gravitational_field .and. .not. param%do_read_potential) Then
          Call mpih5writehalopart(MpiSubCubeWrite, param, haloNB, halopartNB, haloMass, haloID, halopartPos, &
               halopartVel, halopartID, halopartFor=halopartFor, inforamses=inforamses)
       Else
          Call mpih5writehalopart(MpiSubCubeWrite, param, haloNB, halopartNB, haloMass, haloID, halopartPos, &
               halopartVel, halopartID, inforamses=inforamses)
       End If
    End If

    
    If(param%do_timings) Then
       Call Mpi_Barrier(Mpi_Comm_World,mpierr)
       tOuthalopart = Mpi_Wtime() - timeInt
       If(procID==0) Print *,'Process 0 spent ',tOuthalopart,' s for the output of halo particles.'
       timeInt = Mpi_Wtime()
    End If

    Deallocate(halopartID)
    If(Allocated(halopartPot)) Deallocate(halopartPot)
    If(Allocated(halopartFor)) Deallocate(halopartFor)
    !! compute position and velocity of the center of mass for each halo

    If(procID==0) Then
       Print *,' '
       Print *,'pFoF computes observables and writes them.'
       Print *,' '
    End If

    Call computecom(periodic)
    Call computeradius(periodic)

    If(param%do_timings) Then
       Call Mpi_Barrier(Mpi_Comm_World,mpierr)
       tObs = Mpi_Wtime() - timeInt
       If(procID==0) Print *,'Process 0 spent',tObs,' s for the computation of observables.'
       timeInt = Mpi_Wtime()
    End If

    nh = ubound(haloMass,1)
    Call mpih5writehalomass(MpiCube, param, haloNB_all, haloNB, nh, haloMass, halocomPos, &
         halocomVel, haloID, haloRadius, haloSubHaloNB, inforamses=inforamses)

    If(param%do_timings) Then
       Call Mpi_Barrier(Mpi_Comm_World,mpierr)
       tOutmass = Mpi_Wtime() - timeInt
       If(procID==0) Print *,'Process 0 spent',tOutmass,' s to write the observables.'
    End If

    If(procID==0) Then
       Print *,' '
       Print *,'******************************************************'
       Print *,'End of pFoF!'
       Print *,'Number of haloes with mass > ', param%Mmin,':',haloNB_all
       Print *,'******************************************************'
       Print *,' '

       Write(Ulog,*) 'Number of haloes with mass > ', param%Mmin,':',haloNB_all
       
    End If
    
    Deallocate(haloMass, halocomPos, halocomVel, haloID)

  End Subroutine fofpara


  !--------------------------------------------------------------------------------------------------
  !> Initializes the stacks of pointers for local friends of friends.
  !> @brief
  !!
  !> The whole domain is divided into t^3 cubes, where t is the resolution of the simulation analized times the refine_fof factor.<br>
  !! Each process will implicitly considere only its subdomains because every particles treated by the process is located
  !! in the same subdomain.<br>
  !! This routine initializes one stack of "pointers" to the index of the particles located in each cube.
  Subroutine init(n,np,ngpp,x,pile,adfirst,npcase,xmin,ymin,zmin,t)

    Implicit none

    ! Input parameters
    Integer(kind=4), Intent(in) :: n                   !< number of "cubes" in each dimension in a subdomain, i.e. there are
                                                       !! n^3 cubes in a subdomain
    Integer(kind=4), Intent(in) :: np                  !< number of particles in the subdomain
    Integer(kind=4), Intent(in) :: ngpp                !< number of "FOF" cubes in the subdomain
    Real(kind=4), dimension(3,np), Intent(in)   :: x   !< positions of the particles
    Real(kind=4), Intent(in) :: xmin,ymin,zmin         !< minimum x,y and z of the subdomain
    Real(kind=4), Intent(in) :: t                      !< resolution of the cosmological simulation times refine_fof factor

    ! Output parameters
    Integer(kind=4), dimension(ngpp), Intent(out) :: npcase    !< number of particles in each cube
    Integer(kind=4), dimension(ngpp), Intent(out) :: adfirst   !< index of the first particle to examine in the cube
    Integer(kind=4), dimension(np), Intent(out) :: pile        !< stack of indices

    ! Local variables
    Integer(kind=4), dimension(:), allocatable :: adlast          ! index of last particle "added" to a cube
    Integer(kind=4) :: ipart, icase, ix, iy, iz         ! temporary indices
    Real(kind=4) :: xx, yy, zz                         ! temporary position

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
          print*,icase,ix,iy,iz,n,xmin,ymin,zmin,xx,yy,zz,ngpp,t
          print*,x(1,ipart),x(2,ipart),x(3,ipart)
          Call EmergencyStop('Probleme in initialization of friends of friends ',12)
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
