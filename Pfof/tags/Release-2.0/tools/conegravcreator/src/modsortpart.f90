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

    Integer :: ic, ix, iy, iz, icube, ilvl
    Integer :: icell, fclevel
    Real(kind=8) :: ihx, ihy, ihz
    Real(kind=8) :: pi
    Real(kind=8) :: cubesizeinv
    Integer :: ncellloc

    ncellloc = sum(ncellperlevel)

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

       Allocate( ncellcubeloc(ncube,nlevel), ncellcube(ncube,nlevel) )
       Allocate( idcube(ncellloc) )
       ncellcubeloc = 0

       fclevel = 0
       Do ilvl=1,nlevel
          Do icell = 1, ncellperlevel(ilvl)
             ic = fclevel + icell
             ix = int( (pos(1,ic)+ihx) * cubesizeinv) + 1
             iy = int( (pos(2,ic)+ihy) * cubesizeinv) + 1
             iz = int( (pos(3,ic)+ihz) * cubesizeinv) + 1
             icube = ix + (iy-1)*ncx + (iz-1)*ncx*ncy !!+ 1
             ncellcubeloc(icube,ilvl) = ncellcubeloc(icube,ilvl) + 1
             idcube(ic) = icube
          End Do
          fclevel = fclevel + ncellperlevel(ilvl)
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
    
       Allocate( ncellcubeloc(ncube,nlevel), ncellcube(ncube,nlevel) )
       Allocate( idcube(ncellloc) )
       ncellcubeloc = 0
       
       fclevel = 0
       Do ilvl = 1, nlevel
          Do icell = 1, ncellperlevel(ilvl)
             ic = fclevel + icell
             ix = int(pos(1,ic) * cubesizeinv) + 1
             iy = int( (pos(2,ic)+ihy) * cubesizeinv) + 1
             iz = int( (pos(3,ic)+ihz) * cubesizeinv) + 1
             icube = ix + (iy-1)*ncx + (iz-1)*ncx*ncy !!+ 1
             ncellcubeloc(icube,ilvl) = ncellcubeloc(icube,ilvl) + 1
             idcube(ic) = icube
          End Do
          fclevel = fclevel + ncellperlevel(ilvl)
       End Do

    End If

#ifdef WITHMPI3
    Call Mpi_Allreduce(ncellcubeloc, ncellcube, ncube*nlevel, Mpi_Integer, Mpi_Sum, Mpi_Comm_World, req_sumnpc, mpierr)
#else
    Call Mpi_Allreduce(ncellcubeloc, ncellcube, ncube*nlevel, Mpi_Integer, Mpi_Sum, Mpi_Comm_World, mpierr)
#endif

  End Subroutine dividespace



  Subroutine tritas(n,tref,tx,ta,tp,tr,ti)

    Implicit None

    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,ta     ! positions and accelerations
    Real   (kind=4), Intent(inout),dimension(*)   :: tr, tp
    Integer(kind=4), Intent(inout),dimension(*)   :: tref 
    Integer(kind=4), Intent(inout), dimension(*) :: ti

    ! Input variable
    Integer(kind=4),Intent(in) :: n

    ! Local variables
    Integer(kind=4) :: i,size

    Call construction(size,n,tref,tx,ta,tp,tr,ti)

    Do i = n, 2, -1
       Call echanger(1,i,tref,tx,ta,tp,tr,ti)
       size = size - 1
       Call entasser(size,1,tref,tx,ta,tp,tr,ti)
    End Do

  End Subroutine tritas

  
  Subroutine construction(size,n,tref,tx,ta,tp,tr,ti)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,ta     ! positions and accelerations
    Real   (kind=4), Intent(inout),dimension(*)   :: tr, tp
    Integer(kind=4), Intent(inout),dimension(*)   :: tref 
    Integer(kind=4), Intent(inout), dimension(*) :: ti

    !Output variable
    Integer(kind=4),Intent(out)                :: size

    ! Input variable
    Integer(kind=4),Intent(in) :: n

    ! Local variables
    Integer(kind=4) :: i

    size = n
    Do i = n/2, 1, -1
       Call entasser(size,i,tref,tx,ta,tp,tr,ti)
    End Do

  End Subroutine construction


  Recursive Subroutine entasser(size,i,tref,tx,ta,tp,tr,ti)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,ta     ! positions and accelerations
    Real   (kind=4), Intent(inout),dimension(*)   :: tr, tp
    Integer(kind=4), Intent(inout),dimension(*)   :: tref 
    Integer(kind=4), Intent(inout), dimension(*) :: ti

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
       Call echanger(i,max,tref,tx,ta,tp,tr,ti)
       Call entasser(size,max,tref,tx,ta,tp,tr,ti)
    End If

  End Subroutine entasser


  Subroutine echanger(i,j,tref,tx,ta,tp,tr,ti)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,ta     ! positions and accelerations
    Real   (kind=4), Intent(inout),dimension(*)   :: tr, tp
    Integer(kind=4), Intent(inout),dimension(*)   :: tref 
    Integer(kind=4), Intent(inout), dimension(*) :: ti

    ! Input variable
    Integer(kind=4),Intent(in) :: i,j

    ! Local variables
    Integer(kind=4) :: tmpi
    Real(kind=4 ), dimension(3) :: tmpr

    tmpi = tref(i)
    tref(i) = tref(j)
    tref(j) = tmpi

    tmpr(:) = tx(:,i)
    tx(:,i) = tx(:,j)
    tx(:,j) = tmpr(:)

    tmpr(:) = ta(:,i)
    ta(:,i) = ta(:,j)
    ta(:,j) = tmpr(:)

    tmpr(1) = tp(i)
    tp(i) = tp(j)
    tp(j) = tmpr(1)

    tmpr(1) = tr(i)
    tr(i) = tr(j)
    tr(j) = tmpr(1)

    tmpi = ti(i)
    ti(i) = ti(j)
    ti(j) = tmpi

  End Subroutine echanger




End Module modsortpart
