!> @file
!!This file contains the serial FOF algorithm. 

!> This module contains the serial FOF algorithm. 
!>
!> Authors: E. Audit, F. Roy, V. Bouillot

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
    Use modwritehalo
    Use modsortinterf
    Use modvariables
    Use modhalo
    Use modtiming
    Use modfofmpi

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

    Real(kind=4) :: r, d, size, size_fof, xx, yy, zz
    Real(kind=4) :: dx, dy, dz
    Real(kind=4) :: r2
    Real(kind=4) :: rT
    Real(kind=4) :: xmin, xmax, ymin, ymax, zmin, zmax
    Real(kind=4) :: maxdim, percolim
    Integer(kind=4) :: log2maxdim

    Real(kind=4) :: cubeedge
    Real(kind=4) :: edgem1
    Integer(kind=4) :: nh

    Logical(kind=4) :: periodic

    Real(kind=8) :: ttmp0, ttmp1

    Logical(kind=4) :: gather

    
    ! INITIALIZATION

    ! Initialize timings if necessary
    If(param%do_timings) Then
       time0 = Mpi_Wtime()
    End If

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
    If( param%percolation_length >= percolim/2.d0 ) Then
       nsign = 3
    End If

    If( param%percolation_length >= percolim ) Then 
       If(procID == 0) Then
          Print *, ' ' 
          Print *, '*********************************************************************************'
          Print *, '*** percolation_length * refine_fof should be < 2^log2maxdim*cubeedge         ***'
          Print *, '*** refine_fof is set and used in modfofpara.f90, percolation_length is an    ***'
          Print *, '*** input parameter                                                           ***'
          Print *, '*** please contact Fabrice Roy (fabrice.roy@obspm.fr) for more information    ***'
          Print *, '*** Pfof is exiting                                                           ***'
          Print *, '*********************************************************************************'
          Print *, ' ' 
       End If
       Call Mpi_Barrier(Mpi_Comm_World, mpierr)
       Call Mpi_Abort(Mpi_Comm_World, 2, mpierr)
    End If

    r = param%percolation_length
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

    If(mynpart> 0) Then
       namas = id(1) !!!!!RY
    Else
       namas = 0
    Endif
    ipermut  = 2

    border = 0

    If(param%do_timings) Then
       Call Mpi_Barrier(Mpi_Comm_World,mpierr)
       timeInt = Mpi_Wtime()
       tFoFinit = timeInt - time0
    End If

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

       ! nsign depends on param%percolation_length * refine_fof: see above
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
    If(param%do_timings) Then
       Print *,'Process ',procID,' termine apres ',Mpi_Wtime()-timeInt,'s dans fof local.'
    End If

    Deallocate(lamas,pile,indl)
    Deallocate(adfirst,npcase)


    !-------------------!
    ! FOF local termine !
    !-------------------!

    Call Mpi_Reduce(nbborderloc, nbborder, 1, Mpi_Integer,Mpi_Sum,0,Mpi_Comm_World,mpierr)


    If(procID==0) Print *,'Initialisation du raccordement'

#ifdef DEBUG
    Call Mpi_Barrier(Mpi_Comm_World, mpierr)
#endif

    periodic = .false.
    !! Call the merging procedure with the r**2 as an argument
    Call mergehaloes(r2, Mpi_Comm_World, periodic)


    If(param%do_timings) Then
       Call Mpi_Barrier(Mpi_Comm_World,mpierr)
       tRaccord = Mpi_Wtime() - timeInt
       timeInt = Mpi_Wtime()
    End If

    If(procID==0) Print *,'Fin du raccordement'

    ! call writememfilestep(1, 0, 2, 0, ' End of the raccordement ', 'memall.txt',974, 0.0, 0.0)
    !---------------------!
    ! FIN du raccordement !
    !---------------------!

    If(procID==0) Print *,'Debut du calcul de la masse et de la position du cdm des stuctures'

    Call gatherhaloes(Mpi_Comm_World, param)

    If(param%do_timings) Then
       ttmp0 = Mpi_Wtime()
       Print *,'Process ',procID,' repartition:',ttmp0-timeInt
    End If

    If(param%do_read_ramses_part_id) Then
       If(param%do_read_potential .and. param%do_read_gravitational_field) Then
          Call heapsort(myfinalnpart,stf,posf,velf,forf,potf,idf,ramsesidf)
       Else If(param%do_read_potential .and. .not. param%do_read_gravitational_field) Then
          Call heapsort(myfinalnpart,stf,posf,velf,potf,idf,ramsesidf)
       Else If(param%do_read_gravitational_field .and. .not. param%do_read_potential) Then
          Call heapsort(myfinalnpart,stf,posf,velf,forf,idf,ramsesidf)
       Else
          Call heapsort(myfinalnpart,stf,posf,velf,idf,ramsesidf)
       End If
    Else
       If(param%do_read_potential .and. param%do_read_gravitational_field) Then
          Call heapsort(myfinalnpart,stf,posf,velf,forf,potf,idf)
       Else If(param%do_read_potential .and. .not. param%do_read_gravitational_field) Then
          Call heapsort(myfinalnpart,stf,posf,velf,potf,idf)
       Else If(param%do_read_gravitational_field .and. .not. param%do_read_potential) Then
          Call heapsort(myfinalnpart,stf,posf,velf,forf,idf)
       Else
          Call heapsort(myfinalnpart,stf,posf,velf,idf)
       End If
    End If

    If(param%do_timings) Then
       ttmp1 = Mpi_Wtime()
       Print *,'Process ',procID,' sort :', ttmp1-ttmp0
    End If


    !! select haloes with mass > Mmin
    Call selecthaloes(param)

    dogather: If(param%do_gather_halo) Then 
       If(param%do_read_ramses_part_id) Then
          If(param%do_read_potential .and. param%do_read_gravitational_field) Then
             Call mpih5writehalopart(CommGatherWrite, param, haloNB, halopartNB, haloMass, haloID, &
                  halopartPos, halopartVel, halopartID, halopartFor=halopartFor, &
                  halopartPot=halopartPot, halopartRamsesID=halopartRamsesID, &
                  inforamses=inforamseslast, infocone=infoconelast)
             
          Else If(param%do_read_potential .and. .not. param%do_read_gravitational_field ) Then
             Call mpih5writehalopart(CommGatherWrite, param, haloNB, halopartNB, haloMass, haloID, &
                  halopartPos, halopartVel, halopartID, halopartPot=halopartPot, &
                  halopartRamsesID=halopartRamsesID, inforamses=inforamseslast, infocone=infoconelast)
             
          Else If(param%do_read_gravitational_field .and. .not. param%do_read_potential) Then
             Call mpih5writehalopart(CommGatherWrite, param, haloNB, halopartNB, haloMass, haloID, &
                  halopartPos, halopartVel, halopartID, halopartFor=halopartFor, &
                  halopartRamsesID=halopartRamsesID, inforamses=inforamseslast, infocone=infoconelast)
             
          Else
             Call mpih5writehalopart(CommGatherWrite, param, haloNB, halopartNB, haloMass, haloID, &
                  halopartPos, halopartVel, halopartID, halopartRamsesID=halopartRamsesID,&
                  inforamses=inforamseslast, infocone=infoconelast)
          End If
          
       Else
          If(param%do_read_potential .and. param%do_read_gravitational_field) Then
             Call mpih5writehalopart(CommGatherWrite, param, haloNB, halopartNB, haloMass, haloID, &
                  halopartPos, halopartVel, halopartID, halopartFor=halopartFor, &
                  halopartPot=halopartPot, inforamses=inforamseslast, infocone=infoconelast)
             
          Else If(param%do_read_potential .and. .not. param%do_read_gravitational_field ) Then
             Call mpih5writehalopart(CommGatherWrite, param, haloNB, halopartNB, haloMass, haloID, &
                  halopartPos, halopartVel, halopartID, halopartPot=halopartPot, &
                  inforamses=inforamseslast, infocone=infoconelast)
             
          Else If(param%do_read_gravitational_field .and. .not. param%do_read_potential) Then
             Call mpih5writehalopart(CommGatherWrite, param, haloNB, halopartNB, haloMass, haloID, &
                  halopartPos, halopartVel, halopartID, halopartFor=halopartFor, &
                  inforamses=inforamseslast, infocone=infoconelast)
          Else
             Call mpih5writehalopart(CommGatherWrite, param, haloNB, halopartNB, haloMass, haloID, &
                  halopartPos, halopartVel, halopartID, inforamses=inforamseslast, infocone=infoconelast)
          End If
       End If

    Else
       
       If(param%do_read_ramses_part_id) Then
          If(param%do_read_potential .and. param%do_read_gravitational_field) Then
             Call h5writehalopart(param, haloNB, halopartNB, haloMass, haloID, &
                  halopartPos, halopartVel, halopartID, halopartFor=halopartFor, &
                  halopartPot=halopartPot, halopartRamsesID=halopartRamsesID, &
                  inforamses=inforamseslast, infocone=infoconelast)
             
          Else If(param%do_read_potential .and. .not. param%do_read_gravitational_field ) Then
             Call h5writehalopart(param, haloNB, halopartNB, haloMass, haloID, &
                  halopartPos, halopartVel, halopartID, halopartPot=halopartPot, &
                  halopartRamsesID=halopartRamsesID, inforamses=inforamseslast, infocone=infoconelast)
             
          Else If(param%do_read_gravitational_field .and. .not. param%do_read_potential) Then
             Call h5writehalopart(param, haloNB, halopartNB, haloMass, haloID, &
                  halopartPos, halopartVel, halopartID, halopartFor=halopartFor, &
                  halopartRamsesID=halopartRamsesID, inforamses=inforamseslast, infocone=infoconelast)
             
          Else
             Call h5writehalopart(param, haloNB, halopartNB, haloMass, haloID, &
                  halopartPos, halopartVel, halopartID, halopartRamsesID=halopartRamsesID,&
                  inforamses=inforamseslast, infocone=infoconelast)
          End If
          
       Else
          If(param%do_read_potential .and. param%do_read_gravitational_field) Then
             Call h5writehalopart(param, haloNB, halopartNB, haloMass, haloID, &
                  halopartPos, halopartVel, halopartID, halopartFor=halopartFor, &
                  halopartPot=halopartPot, inforamses=inforamseslast, infocone=infoconelast)
             
          Else If(param%do_read_potential .and. .not. param%do_read_gravitational_field ) Then
             Call h5writehalopart(param, haloNB, halopartNB, haloMass, haloID, &
                  halopartPos, halopartVel, halopartID, halopartPot=halopartPot, &
                  inforamses=inforamseslast, infocone=infoconelast)
             
          Else If(param%do_read_gravitational_field .and. .not. param%do_read_potential) Then
             Call h5writehalopart(param, haloNB, halopartNB, haloMass, haloID, &
                  halopartPos, halopartVel, halopartID, halopartFor=halopartFor, &
                  inforamses=inforamseslast, infocone=infoconelast)
          Else
             Call h5writehalopart(param, haloNB, halopartNB, haloMass, haloID, &
                  halopartPos, halopartVel, halopartID, inforamses=inforamseslast, infocone=infoconelast)
          End If
       End If
    End If dogather
    

    
    Deallocate(halopartID)

    !! compute position and velocity of the center of mass for each halo
    Call computecom(periodic)

    Call computeradius(periodic)

    If(param%do_timings) Then
       Print *,'Process ',procID,' termine apres ',Mpi_Wtime()-timeInt,' s dans calcul des observables'
    End If

    If(param%do_timings) Then
       Print *,'Process ',procID,' termine apres ',Mpi_Wtime()-timeInt,' s dans calcul des observables'
    End If


    If(procID==0) Print *,'Fin du calcul de la masse et des centres de masse des structures'
    If(procID==0) Print *,'Nombre de halos de masse superieure a', param%mmin,':',haloNB_all

    If(param%do_timings) Then
       Call Mpi_Barrier(Mpi_Comm_World,mpierr)
       tObs = Mpi_Wtime() - timeInt
       timeInt = Mpi_Wtime()
    End If

    nh = ubound(haloMass,1)
    Call mpih5writehalomass(Mpi_Comm_World,param, haloNB_all, haloNB, nh, haloMass, halocomPos, & 
         halocomVel, haloID, haloRadius, haloSubHaloNB, inforamses=inforamseslast, infocone=infoconelast)

    Deallocate(haloMass, halocomPos, halocomVel, haloID)

    If(param%do_timings) Then
       Call Mpi_Barrier(Mpi_Comm_World,mpierr)
       tOut = Mpi_Wtime() - timeInt
       tFoF = Mpi_Wtime() - time0
    End If


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


  ! ======================================================================

  Subroutine computeminmax()
    
    Use modvariables, only : PR, pos, partmin, partmax
    Implicit none

    Integer(kind=4) :: i
    Real(kind=PR), dimension(3) :: locmin, locmax


    locmin(:) = 1.e10
    locmax(:) = -1.e10

    Do i=1, 3
       locmin(i) = minval(pos(i,:))
       locmax(i) = maxval(pos(i,:))
    End Do
      
    Call Mpi_Allreduce(locmin, partmin, 3, Mpi_Real, Mpi_Min, Mpi_Comm_World, mpierr)
    Call Mpi_Allreduce(locmax, partmax, 3, Mpi_Real, Mpi_Max, Mpi_Comm_World, mpierr)

  End Subroutine computeminmax

End Module modfofpara
