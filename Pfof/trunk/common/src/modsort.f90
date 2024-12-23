!> @file
!! This file contains two sort subroutines.

!> This module contains two sort subroutines.
!>
!> Authors: F. Roy, V. Bouillot
Module modsort

  Use modconstant

  Private

  Public :: trirapide, &
       tritas


Contains
  
  !=======================================================================

  !> Quick sort as described in Introduction to algorithms (Cormen, Leiserson, Rivest, Stein)
  Recursive Subroutine trirapide(p,r,tref,tid,tx,tv)
    
    Implicit None

    ! Input parameters
    Integer(kind=4), Intent(in) :: p     !< index of the first element of the array to be sorted
    Integer(kind=4), Intent(in) :: r     !< index of the last element of the array to be sorted

    ! Input parameters modified in subroutine
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx     !< positions array
    Real   (kind=4), Intent(inout),dimension(3,*) :: tv     !< velocities array
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref  !< halo id array
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tid   !< particle id array

    ! Local parameters
    Integer(kind=4) :: q

#ifdef DEBUG2
    Print *,'trirapide begins'
#endif

    If(p<r) Then
       
       Call partition(tref,p,r,q,tid,tx,tv)
       Call trirapide(p,q-1,tref,tid,tx,tv)
       Call trirapide(q+1,r,tref,tid,tx,tv)

    End If

#ifdef DEBUG2
    Print *,'trirapide ends'
#endif

  End Subroutine trirapide


  !=======================================================================

  Subroutine partition(tref,p,r,q,tid,tx,tv)

    Implicit None

    ! Input parameters
    Integer(kind=4),Intent(in)  :: p, r  ! First and last index of the tables to sort
    ! Input parameters modified in subroutine
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,  tv  ! positions and velocities 
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref,tid ! structure id and particle id
    ! Output parameters
    Integer(kind=4), Intent(out) :: q
    ! Local parameters
    Integer(kind=PRI)     :: tmpi,Sref
    Real(kind=4),dimension(3) :: tmpr
    Integer(kind=4) :: i, j

#ifdef DEBUG2
    Print *,'partitions'
#endif

    Sref  = tref(r)

    i = p-1
    Do j = p, r-1
       If(tref(j) <= Sref) Then
          i = i+1
          tmpi = tref(i)
          tref(i) = tref(j)
          tref(j) = tmpi
          tmpi = tid(i)
          tid(i) = tid(j)
          tid(j) = tmpi
          tmpr(:) = tx(:,i)
          tx(:,i) = tx(:,j)
          tx(:,j) = tmpr(:)
          tmpr(:) = tv(:,i)
          tv(:,i) = tv(:,j)
          tv(:,j) = tmpr(:)
       End If
    End Do

    tmpi = tref(i+1)
    tref(i+1) = tref(r)
    tref(r) = tmpi
    tmpi = tid(i+1)
    tid(i+1) = tid(r)
    tid(r) = tmpi
    tmpr(:) = tx(:,i+1)
    tx(:,i+1) = tx(:,r)
    tx(:,r) = tmpr(:)
    tmpr(:) = tv(:,i+1)
    tv(:,i+1) = tv(:,r)
    tv(:,r) = tmpr(:)

    q = i+1

#ifdef DEBUG2
    Print *,'partitions'
#endif

  End Subroutine partition


  !=======================================================================

  !> Heap sort as described in Introduction to algorithms (Cormen, Leiserson, Rivest, Stein)
  Subroutine tritas(tref,tid,tx,tv,n)

    Use modmpicom
    Implicit None

    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid ! structure id and particle id

    ! Input variable
    Integer(kind=4),Intent(in) :: n

    ! Local variables
    Integer(kind=4) :: i,size

#ifdef DEBUG2
    Print *,'tritas ends'
#endif

    Call construction(tref,tid,tx,tv,size,n)

    Do i = n, 2, -1
       Call echanger(tref,tid,tx,tv,1,i)
       size = size - 1
       Call entasser(tref,tid,tx,tv,size,1)
    End Do

#ifdef DEBUG2
    Print *,'tritas ends'
#endif

  End Subroutine tritas


  !=======================================================================

  
  Subroutine construction(tref,tid,tx,tv,size,n)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid ! structure id and particle id

    !Output variable
    Integer(kind=4),Intent(out)                :: size

    ! Input variable
    Integer(kind=4),Intent(in) :: n

    ! Local variables
    Integer(kind=4) :: i

#ifdef DEBUG2
    Print *,'construction ends'
#endif

    size = n
    Do i = n/2, 1, -1
       Call entasser(tref,tid,tx,tv,size,i)
    End Do

#ifdef DEBUG2
    Print *,'construction ends'
#endif

  End Subroutine construction


  !=======================================================================


  Recursive Subroutine entasser(tref,tid,tx,tv,size,i)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid ! structure id and particle id

    ! Input variable
    Integer(kind=4),Intent(in) :: size,i

    ! Local variables
    Integer(kind=4) :: l, r, max

#ifdef DEBUG2
    Print *,'entasser ends'
#endif

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
       Call echanger(tref,tid,tx,tv,i,max)
       Call entasser(tref,tid,tx,tv,size,max)
    End If

#ifdef DEBUG2
    Print *,'entasser ends'
#endif

  End Subroutine entasser


  !=======================================================================


  Subroutine echanger(tref,tid,tx,tv,i,j)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid ! structure id and particle id

    ! Input variable
    Integer(kind=4),Intent(in) :: i,j

    ! Local variables
    Integer(kind=PRI) :: tmpi
    Real(kind=4), dimension(3) :: tmpr

#ifdef DEBUG2
    Print *,'echanger ends'
#endif

    tmpi = tref(i)
    tref(i) = tref(j)
    tref(j) = tmpi

    tmpi = tid(i)
    tid(i) = tid(j)
    tid(j) = tmpi

    tmpr(:) = tx(:,i)
    tx(:,i) = tx(:,j)
    tx(:,j) = tmpr(:)

    tmpr(:) = tv(:,i)
    tv(:,i) = tv(:,j)
    tv(:,j) = tmpr(:)

#ifdef DEBUG2
    Print *,'echanger ends'
#endif

  End Subroutine echanger


End Module modsort
