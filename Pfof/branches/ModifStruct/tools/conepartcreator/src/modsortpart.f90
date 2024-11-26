!==============================================================================
! Project: pFoF
! File: tools/conepartcreator/src/modsortpart.f90
! Copyright Fabrice Roy (2015)
! Fabrice.Roy@obspm.fr
! 
! This software is a computer program whose purpose is to detect dark matter
! haloes in cosmological simulation with a parallel Friend of Friend algorithm. 
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can  use, 
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info". 
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability. 
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or 
! data to be ensured and, more generally, to use and operate it in the 
! same conditions as regards security. 
!
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL license and that you accept its terms.
!==============================================================================

!=======================================================================
! Author: Fabrice Roy (LUTH/CNRS/Observatoire de Paris)
! Fabrice.Roy@obspm.fr
!=======================================================================

Module modsortpart

  Use mpi
  Use modvariables
  Use modio
  Use modmpicommons, only : procID, procNB

  Implicit none


Contains

  Subroutine dividespace()
    
    Implicit none 

    Integer :: ip, ix, iy, iz, icube
    Real(kind=8) :: ihx, ihy, ihz
    Real(kind=8) :: pi
    Real(kind=8) :: cubesizeinv
    Integer(kind=4) :: mpierr

    If(infocone%isfullsky==1) Then
       ! space is divided into cube ; the size of the cube ridge is in Ramses units.
       cubesizeinv = 1.d0/param%cube_size

       ihx = (int(param%cone_max_radius * cubesizeinv) + 1)*param%cube_size
       ncx = (int(param%cone_max_radius * cubesizeinv) + 1)*2 
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
       ! compute size of the cuboid encompassing the cone : 2hy * 2hz * param%cone_max_radius
       pi = 4.0d0* atan(1.d0)
       hy = param%cone_max_radius * tan(infocone%thetay*pi/180.d0)
       hz = param%cone_max_radius * tan(infocone%thetaz*pi/180.d0)
       
       ! space is divided into cube ; the size of the cube ridge is in Ramses units.
       cubesizeinv = 1.d0/param%cube_size
       ! max number of cubes in each direction:
       ncx = int(param%cone_max_radius * cubesizeinv) + 1
       ihy = (int(hy * cubesizeinv) + 1)*param%cube_size
       ncy = (int(hy * cubesizeinv) + 1)*2
       ihz = (int(hz * cubesizeinv) + 1)*param%cube_size
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


  !=======================================================================

  Subroutine heapsort(n,tref,tx,tv,tf,tp,tid)

    Implicit None

    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=4), Intent(inout),dimension(*)   :: tref 
    Integer(kind=PRI), Intent(inout), dimension(*), optional :: tid
    Real(kind=4), Intent(inout), dimension(*), optional :: tp
    Real(kind=4), Intent(inout), dimension(3,*), optional :: tf
    
    ! Input variable
    Integer(kind=4),Intent(in) :: n

    ! Local variables
    Integer(kind=4) :: i,size

    If(present(tid)) Then
       If(present(tp)) Then
          If(present(tf)) Then
             Call construction(size,n,tref,tx,tv,tf=tf,tp=tp,tid=tid)
             Do i = n, 2, -1
                Call echanger(1,i,tref,tx,tv,tf=tf,tp=tp,tid=tid)
                size = size - 1
                Call entasser(size,1,tref,tx,tv,tf=tf,tp=tp,tid=tid)
             End Do

          Else
             Call construction(size,n,tref,tx,tv,tp=tp,tid=tid)
             Do i = n, 2, -1
                Call echanger(1,i,tref,tx,tv,tp=tp,tid=tid)
                size = size - 1
                Call entasser(size,1,tref,tx,tv,tp=tp,tid=tid)
             End Do

          End If
       Else
          Call construction(size,n,tref,tx,tv,tid=tid)
          Do i = n, 2, -1
             Call echanger(1,i,tref,tx,tv,tid=tid)
             size = size - 1
             Call entasser(size,1,tref,tx,tv,tid=tid)
          End Do

       End If
    Else
       If(present(tp)) Then
          If(present(tf)) Then
             Call construction(size,n,tref,tx,tv,tf=tf,tp=tp)
             Do i = n, 2, -1
                Call echanger(1,i,tref,tx,tv,tf=tf,tp=tp)
                size = size - 1
                Call entasser(size,1,tref,tx,tv,tf=tf,tp=tp)
             End Do

          Else
             Call construction(size,n,tref,tx,tv,tp=tp)
             Do i = n, 2, -1
                Call echanger(1,i,tref,tx,tv,tp=tp)
                size = size - 1
                Call entasser(size,1,tref,tx,tv,tp=tp)
             End Do

          End If
       Else
          Call construction(size,n,tref,tx,tv)
          Do i = n, 2, -1
             Call echanger(1,i,tref,tx,tv)
             size = size - 1
             Call entasser(size,1,tref,tx,tv)
          End Do

       End If
    End If
    



  End Subroutine heapsort

  !=======================================================================
  
  Subroutine construction(size,n,tref,tx,tv,tf,tp,tid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=4), Intent(inout),dimension(*)   :: tref 
    Integer(kind=PRI), Intent(inout), dimension(*), optional :: tid
    Real(kind=4), Intent(inout), dimension(*), optional :: tp
    Real(kind=4), Intent(inout), dimension(3,*), optional :: tf

    !Output variable
    Integer(kind=4),Intent(out)                :: size

    ! Input variable
    Integer(kind=4),Intent(in) :: n

    ! Local variables
    Integer(kind=4) :: i

    size = n

    If(present(tid)) Then
       If(present(tp)) Then
          If(present(tf)) Then
             Do i = n/2, 1, -1
                Call entasser(size,i,tref,tx,tv,tf=tf,tp=tp,tid=tid)
             End Do
          Else
             Do i = n/2, 1, -1
                Call entasser(size,i,tref,tx,tv,tp=tp,tid=tid)
             End Do
          End If
       Else
          Do i = n/2, 1, -1
             Call entasser(size,i,tref,tx,tv,tid=tid)
          End Do
       End If
    Else
       If(present(tp)) Then
          If(present(tf)) Then
             Do i = n/2, 1, -1
                Call entasser(size,i,tref,tx,tv,tf=tf,tp=tp)
             End Do
          Else
             Do i = n/2, 1, -1
                Call entasser(size,i,tref,tx,tv,tp=tp)
             End Do
          End If
       Else
          Do i = n/2, 1, -1
             Call entasser(size,i,tref,tx,tv)
          End Do
       End If
    End If
    
    
  End Subroutine construction

  !=======================================================================

  Recursive Subroutine entasser(size,i,tref,tx,tv,tf,tp,tid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=4), Intent(inout),dimension(*)   :: tref
    Integer(kind=PRI), Intent(inout), dimension(*), optional :: tid
    Real(kind=4), Intent(inout), dimension(*), optional :: tp
    Real(kind=4), Intent(inout), dimension(3,*), optional :: tf

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
       If(present(tid)) Then
          If(present(tp)) Then
             If(present(tf)) Then
                Call echanger(i,max,tref,tx,tv,tf=tf, tp=tp, tid=tid)
                Call entasser(size,max,tref,tx,tv,tf=tf,tp=tp,tid=tid)
             Else
                Call echanger(i,max,tref,tx,tv,tp=tp,tid=tid)
                Call entasser(size,max,tref,tx,tv,tp=tp,tid=tid)
             End If
          Else
             Call echanger(i,max,tref,tx,tv,tid=tid)
             Call entasser(size,max,tref,tx,tv,tid=tid)
          End If
       Else
          If(present(tp)) Then
             If(present(tf)) Then
                Call echanger(i,max,tref,tx,tv,tf=tf, tp=tp)
                Call entasser(size,max,tref,tx,tv,tf=tf,tp=tp)
             Else
                Call echanger(i,max,tref,tx,tv,tp=tp)
                Call entasser(size,max,tref,tx,tv,tp=tp)
             End If
          Else
             Call echanger(i,max,tref,tx,tv)
             Call entasser(size,max,tref,tx,tv)
          End If
       End If
    End If

  End Subroutine entasser

  !=======================================================================

  Subroutine echanger(i,j,tref,tx,tv,tf,tp,tid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=4), Intent(inout),dimension(*)   :: tref
    Integer(kind=PRI), Intent(inout), dimension(*), optional :: tid
    Real(kind=4), Intent(inout), dimension(*), optional :: tp
    Real(kind=4), Intent(inout), dimension(3,*), optional :: tf

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
    
    If(present(tf)) Then
       tmpr(:) = tf(:,i)
       tf(:,i) = tf(:,j)
       tf(:,j) = tmpr(:)
    End If

    If(present(tp)) Then
       tmpr(1) = tp(i)
       tp(i) = tp(j)
       tp(j) = tmpr(1)
    End If

  End Subroutine echanger




End Module modsortpart
