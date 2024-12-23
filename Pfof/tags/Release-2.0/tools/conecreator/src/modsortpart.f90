!=======================================================================
! Author: Fabrice Roy (LUTH/CNRS/Observatoire de Paris)
! Fabrice.Roy@obspm.fr
!=======================================================================

Module modsortpart

  Use mpi
  Use modvariable
  Use modparam
  Use modio

  Implicit none


Contains

  Subroutine dividespace()
    
    Implicit none 

    Integer :: ip, ix, iy, iz, icube
    Real(kind=8) :: ihx, ihy, ihz
    Real(kind=8) :: pi
    Real(kind=8) :: cubesizeinv

    If(fullsky) Then
       ! space is divided into cube ; the size of the cube ridge is in Ramses units.
       cubesizeinv = 1.d0/cubesize

       ihx = (int(rmax * cubesizeinv) + 1)*cubesize
       ncx = (int(rmax * cubesizeinv) + 1)*2 
       ihy = ihx
       ncy = ncx
       ihz = ihx
       ncz = ncx
       ncube = ncx*ncy*ncz

       Allocate( npartcubeloc(ncube), npartcube(ncube) )
       Allocate( idcube(npartloc) )
       npartcubeloc = 0

       Do ip = 1, npartloc
          ix = int( (pos(1,ip)+ihx) * cubesizeinv) + 1
          iy = int( (pos(2,ip)+ihy) * cubesizeinv) + 1
          iz = int( (pos(3,ip)+ihz) * cubesizeinv) + 1
          icube = ix + (iy-1)*ncx + (iz-1)*ncx*ncy !!+ 1
          npartcubeloc(icube) = npartcubeloc(icube) + 1
          idcube(ip) = icube
       End Do

    Else
       ! compute size of the cuboid encompassing the cone : 2hy * 2hz * rmax
       pi = 4.0d0* atan(1.d0)
       hy = rmax * tan(thetay*pi/180.d0)
       hz = rmax * tan(thetaz*pi/180.d0)
       
       ! space is divided into cube ; the size of the cube ridge is in Ramses units.
       cubesizeinv = 1.d0/cubesize
       ! max number of cubes in each direction:
       ncx = int(rmax * cubesizeinv) + 1
       ihy = (int(hy * cubesizeinv) + 1)*cubesize
       ncy = (int(hy * cubesizeinv) + 1)*2
       ihz = (int(hz * cubesizeinv) + 1)*cubesize
       ncz = (int(hz * cubesizeinv) + 1)*2
       ncube = ncx*ncy*ncz
    
       Allocate( npartcubeloc(ncube), npartcube(ncube) )
       Allocate( idcube(npartloc) )
       npartcubeloc = 0
       
       Do ip = 1, npartloc
          ix = int(pos(1,ip) * cubesizeinv) + 1
          iy = int( (pos(2,ip)+ihy) * cubesizeinv) + 1
          iz = int( (pos(3,ip)+ihz) * cubesizeinv) + 1
          icube = ix + (iy-1)*ncx + (iz-1)*ncx*ncy !!+ 1
          npartcubeloc(icube) = npartcubeloc(icube) + 1
          idcube(ip) = icube
       End Do

    End If

#ifdef WITHMPI3
    Call Mpi_Allreduce(npartcubeloc, npartcube, ncube, Mpi_Integer, Mpi_Sum, Mpi_Comm_World, req_sumnpc, mpierr)
#else
    Call Mpi_Allreduce(npartcubeloc, npartcube, ncube, Mpi_Integer, Mpi_Sum, Mpi_Comm_World, mpierr)
#endif

  End Subroutine dividespace



  Subroutine tritas(n,tref,tx,tv,tid)

    Implicit None

    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=4), Intent(inout),dimension(*)   :: tref 
    Integer(kind=PRI), Intent(inout), dimension(*), optional :: tid

    ! Input variable
    Integer(kind=4),Intent(in) :: n

    ! Local variables
    Integer(kind=4) :: i,size

    Call construction(size,n,tref,tx,tv,tid)

    Do i = n, 2, -1
       Call echanger(1,i,tref,tx,tv,tid)
       size = size - 1
       Call entasser(size,1,tref,tx,tv,tid)
    End Do

  End Subroutine tritas

  
  Subroutine construction(size,n,tref,tx,tv,tid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=4), Intent(inout),dimension(*)   :: tref 
    Integer(kind=PRI), Intent(inout), dimension(*), optional :: tid

    !Output variable
    Integer(kind=4),Intent(out)                :: size

    ! Input variable
    Integer(kind=4),Intent(in) :: n

    ! Local variables
    Integer(kind=4) :: i

    size = n
    Do i = n/2, 1, -1
       Call entasser(size,i,tref,tx,tv,tid)
    End Do

  End Subroutine construction


  Recursive Subroutine entasser(size,i,tref,tx,tv,tid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=4), Intent(inout),dimension(*)   :: tref
    Integer(kind=PRI), Intent(inout), dimension(*), optional :: tid

    ! Input variable
    Integer(kind=4),Intent(in) :: size,i

    ! Local variables
    Integer(kind=4) :: l, r, max


    l = 2*i
    r = 2*i+1
    
    If( (l <= size) .and. tref(l) > tref(i) ) Then
       max = l
    Else
       max = i
    End If
    
    If( (r <= size) .and. tref(r) > tref(max) ) Then
       max = r
    End If
    
    If( max /= i ) Then
       Call echanger(i,max,tref,tx,tv,tid)
       Call entasser(size,max,tref,tx,tv,tid)
    End If

  End Subroutine entasser


  Subroutine echanger(i,j,tref,tx,tv,tid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=4), Intent(inout),dimension(*)   :: tref
    Integer(kind=PRI), Intent(inout), dimension(*), optional :: tid

    ! Input variable
    Integer(kind=4),Intent(in) :: i,j

    ! Local variables
    Integer(kind=4) :: tmpi
    Real(kind=4 ), dimension(3) :: tmpr
    Integer(kind=PRI) :: tmpipri

    tmpi = tref(i)
    tref(i) = tref(j)
    tref(j) = tmpi

    tmpr(:) = tx(:,i)
    tx(:,i) = tx(:,j)
    tx(:,j) = tmpr(:)

    tmpr(:) = tv(:,i)
    tv(:,i) = tv(:,j)
    tv(:,j) = tmpr(:)

    If(present(tid)) Then
       tmpipri = tid(i)
       tid(i) = tid(j)
       tid(j) = tmpipri
    End If

  End Subroutine echanger




End Module modsortpart
