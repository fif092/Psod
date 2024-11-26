!==============================================================================
! Project: pFoF
! File: common/src/modsort.f90
! Copyright Fabrice Roy (2011)
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

!> @file
!! This file contains sort subroutines.

!> This module contains sort subroutines.
!>
!> Authors: F. Roy

Module modsort

  Use modconstant, only : PRI

  Private

  Public :: quicksort, heapsort

  !> Generic inteface used to sort arrays with quicksort algorithm.
  !> @param[in] p First index of the arrays to sort
  !> @param[in] r Last index of the arrays to sort
  !> @param[in] tref Reference array: this array will be sorted in increasing order, the other arrays will be sorted following this one
  !> @param[in] [tx,tv,tid] Arrays to sort
  !> @param[in] [tp,tf] Optional arrays to sort
  Interface quicksort
     Module procedure quicksort_refI_2x3R_1xI  
     Module procedure quicksort_refI_2x3R_1xR_1xI  
     Module procedure quicksort_refI_3x3R_1xR_1xI  
     Module procedure quicksort_refI_3x3R_1xI  
  End Interface quicksort

  
  !> Generic inteface used to sort arrays with heapsort algorithm.
  !> @param[in] n Lenght of the arrays
  !> @param[in] tref Reference array: this array will be sorted in increasing order, the other arrays will be sorted following this one
  !> @param[in] [tx,tv,tid] Arrays to sort
  !> @param[in] [tp,tf,trid] Optional arrays to sort
  Interface heapsort
     Module procedure heapsort_refI_2x3R_1xI
     Module procedure heapsort_refI_2x3R_1xR_1xI     
     Module procedure heapsort_refI_3x3R_1xR_1xI     
     Module procedure heapsort_refI_3x3R_1xI
     Module procedure heapsort_refI_2x3R_2xI
     Module procedure heapsort_refI_2x3R_1xR_2xI
     Module procedure heapsort_refI_3x3R_1xR_2xI
     Module procedure heapsort_refI_3x3R_2xI
  End Interface heapsort
  

Contains
  !=======================================================================
  !> Quicksort algorith as described in Introduction to algorithms (Cormen, Leiserson, Rivest, Stein)
  Recursive Subroutine quicksort_refI_2x3R_1xI(p,r,tref,tx,tv,tid)

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
    Print *,'quicksort_refI_2x3R_1xI begins between ',p,' and ',r
#endif

    If(p<r) Then
       
       Call partition_refI_2x3R_1xI(p,r,q,tref,tx,tv,tid)
       Call quicksort_refI_2x3R_1xI(p,q-1,tref,tx,tv,tid)
       Call quicksort_refI_2x3R_1xI(q+1,r,tref,tx,tv,tid)

    End If

#ifdef DEBUG2
    Print *,'quicksort ends'
#endif

  End Subroutine quicksort_refI_2x3R_1xI

  !=======================================================================
  !> Quicksort algorith as described in Introduction to algorithms (Cormen, Leiserson, Rivest, Stein)
  Recursive Subroutine quicksort_refI_2x3R_1xR_1xI(p,r,tref,tx,tv,tp,tid)

    Implicit None
    
    ! Input parameters
    Integer(kind=4), Intent(in) :: p     !< index of the first element of the array to be sorted
    Integer(kind=4), Intent(in) :: r     !< index of the last element of the array to be sorted
    
    ! Input parameters modified in subroutine
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx     !< positions array
    Real   (kind=4), Intent(inout),dimension(3,*) :: tv     !< velocities array
    Real   (kind=4), Intent(inout),dimension(*)  :: tp     !< potential array
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref  !< halo id array
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tid   !< particle id array
    
    ! Local parameters
    Integer(kind=4) :: q
    
#ifdef DEBUG2
    Print *,'quicksort_refI_2x3R_1xR_1xI begins'
#endif

    If(p<r) Then
       
       Call partition_refI_2x3R_1xR_1xI(p,r,q,tref,tx,tv,tp,tid)
       Call quicksort_refI_2x3R_1xR_1xI(p,q-1,tref,tx,tv,tp,tid)
       Call quicksort_refI_2x3R_1xR_1xI(q+1,r,tref,tx,tv,tp,tid)

    End If

#ifdef DEBUG2
    Print *,'trirapide ends'
#endif


  End Subroutine quicksort_refI_2x3R_1xR_1xI


  !=======================================================================
  !> Quicksort algorith as described in Introduction to algorithms (Cormen, Leiserson, Rivest, Stein)
  Recursive Subroutine quicksort_refI_3x3R_1xR_1xI(p,r,tref,tx,tv,tf,tp,tid)

    Implicit None
    
    ! Input parameters
    Integer(kind=4), Intent(in) :: p     !< index of the first element of the array to be sorted
    Integer(kind=4), Intent(in) :: r     !< index of the last element of the array to be sorted
    
    ! Input parameters modified in subroutine
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx     !< positions array
    Real   (kind=4), Intent(inout),dimension(3,*) :: tv     !< velocities array
    Real   (kind=4), Intent(inout),dimension(3,*) :: tf     !< gravitational field array
    Real   (kind=4), Intent(inout),dimension(*)   :: tp     !< potential array
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref  !< halo id array
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tid   !< particle id array
    
    ! Local parameters
    Integer(kind=4) :: q
    
#ifdef DEBUG2
    Print *,'quicksort_refI_3x3R_1xR_1xI begins'
#endif

    If(p<r) Then
       
       Call partition_refI_3x3R_1xR_1xI(p,r,q,tref,tx,tv,tf,tp,tid)
       Call quicksort_refI_3x3R_1xR_1xI(p,q-1,tref,tx,tv,tf,tp,tid)
       Call quicksort_refI_3x3R_1xR_1xI(q+1,r,tref,tx,tv,tf,tp,tid)

    End If

#ifdef DEBUG2
    Print *,'quicksort_refI_3x3R_1xR_1xI ends'
#endif


  End Subroutine quicksort_refI_3x3R_1xR_1xI

  !=======================================================================
  !> Quicksort algorith as described in Introduction to algorithms (Cormen, Leiserson, Rivest, Stein)
  Recursive Subroutine quicksort_refI_3x3R_1xI(p,r,tref,tx,tv,tf,tid)

    Implicit None
    
    ! Input parameters
    Integer(kind=4), Intent(in) :: p     !< index of the first element of the array to be sorted
    Integer(kind=4), Intent(in) :: r     !< index of the last element of the array to be sorted
    
    ! Input parameters modified in subroutine
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx     !< positions array
    Real   (kind=4), Intent(inout),dimension(3,*) :: tv     !< velocities array
    Real   (kind=4), Intent(inout),dimension(3,*) :: tf     !< potential array
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref  !< halo id array
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tid   !< particle id array
    
    ! Local parameters
    Integer(kind=4) :: q
    
#ifdef DEBUG2
    Print *,'quicksort_refI_3x3R_1xI begins'
#endif

    If(p<r) Then
       
       Call partition_refI_3x3R_1xI(p,r,q,tref,tx,tv,tf,tid)
       Call quicksort_refI_3x3R_1xI(p,q-1,tref,tx,tv,tf,tid)
       Call quicksort_refI_3x3R_1xI(q+1,r,tref,tx,tv,tf,tid)

    End If

#ifdef DEBUG2
    Print *,'quicksort_refI_3x3R_1xI ends'
#endif

  End Subroutine quicksort_refI_3x3R_1xI
 
  !=======================================================================
  !> Partition function for quicksort 
  Subroutine partition_refI_2x3R_1xI(p,r,q,tref,tx,tv,tid)

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
    Print *,'partition_refI_2x3R_1xI begins'
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
    Print *,'partition_refI_2x3R_1xI ends'
#endif

  End Subroutine partition_refI_2x3R_1xI

  !=======================================================================
  !> Partition function for quicksort 
  Subroutine partition_refI_2x3R_1xR_1xI(p,r,q,tref,tx,tv,tp,tid)

    Implicit None

    ! Input parameters
    Integer(kind=4),Intent(in)  :: p, r  ! First and last index of the tables to sort
    ! Input parameters modified in subroutine
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,  tv  ! positions and velocities 
    Real   (kind=4), Intent(inout), dimension(*) :: tp
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref,tid ! structure id and particle id
    ! Output parameters
    Integer(kind=4), Intent(out) :: q
    ! Local parameters
    Integer(kind=PRI)     :: tmpi,Sref
    Real(kind=4) :: tmpf
    Real(kind=4),dimension(3) :: tmpr
    Integer(kind=4) :: i, j

#ifdef DEBUG2
    Print *,'partition_refI_2x3R_1xR_1xI begins'
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
          tmpf = tp(i)
          tp(i) = tp(j)
          tp(j) = tmpf
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
    tmpf = tp(i+1)
    tp(i+1) = tp(r)
    tp(r) = tmpf

    q = i+1

#ifdef DEBUG2
    Print *,'partition_refI_2x3R_1xR_1xI ends'
#endif

  End Subroutine partition_refI_2x3R_1xR_1xI

  !=======================================================================
  !> Partition function for quicksort 
  Subroutine partition_refI_3x3R_1xR_1xI(p,r,q,tref,tx,tv,tf,tp,tid)

    Implicit None

    ! Input parameters
    Integer(kind=4),Intent(in)  :: p, r  ! First and last index of the tables to sort
    ! Input parameters modified in subroutine
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,  tv  ! positions and velocities 
    Real   (kind=4), Intent(inout), dimension(3,*) :: tf  ! force
    Real   (kind=4), Intent(inout), dimension(*) :: tp    ! potential
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref,tid ! structure id and particle id
    ! Output parameters
    Integer(kind=4), Intent(out) :: q
    ! Local parameters
    Integer(kind=PRI)     :: tmpi,Sref
    Real(kind=4) :: tmpf
    Real(kind=4),dimension(3) :: tmpr
    Integer(kind=4) :: i, j

#ifdef DEBUG2
    Print *,'partition_refI_3x3R_1xR_1xI begins'
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
          tmpr(:) = tf(:,i)
          tf(:,i) = tf(:,j)
          tf(:,j) = tmpr(:)
          tmpf = tp(i)
          tp(i) = tp(j)
          tp(j) = tmpf
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
    tmpr(:) = tf(:,i+1)
    tf(:,i+1) = tf(:,r)
    tf(:,r) = tmpr(:)
    tmpf = tp(i+1)
    tp(i+1) = tp(r)
    tp(r) = tmpf

    q = i+1

#ifdef DEBUG2
    Print *,'partition_refI_3x3R_1xR_1xI ends'
#endif

  End Subroutine partition_refI_3x3R_1xR_1xI

  !=======================================================================
  !> Partition function for quicksort 
  Subroutine partition_refI_3x3R_1xI(p,r,q,tref,tx,tv,tf,tid)

    Implicit None

    ! Input parameters
    Integer(kind=4),Intent(in)  :: p, r  ! First and last index of the tables to sort
    ! Input parameters modified in subroutine
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,  tv  ! positions and velocities 
    Real   (kind=4), Intent(inout), dimension(3,*) :: tf  ! force
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref,tid ! structure id and particle id
    ! Output parameters
    Integer(kind=4), Intent(out) :: q
    ! Local parameters
    Integer(kind=PRI)     :: tmpi,Sref
    Real(kind=4),dimension(3) :: tmpr
    Integer(kind=4) :: i, j

#ifdef DEBUG2
    Print *,'partition_refI_3x3R_1xI begins'
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
          tmpr(:) = tf(:,i)
          tf(:,i) = tf(:,j)
          tf(:,j) = tmpr(:)
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
    tmpr(:) = tf(:,i+1)
    tf(:,i+1) = tf(:,r)
    tf(:,r) = tmpr(:)

    q = i+1

#ifdef DEBUG2
    Print *,'partition_refI_3x3R_1xI ends'
#endif

  End Subroutine partition_refI_3x3R_1xI
  ! End of quicksort
  !=======================================================================


  !=======================================================================
  !> Heap sort as described in Introduction to algorithms (Cormen, Leiserson, Rivest, Stein)
  Subroutine heapsort_refI_2x3R_1xR_1xI(n,tref,tx,tv,tp,tid)

    Implicit None

    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx     !< positions array
    Real   (kind=4), Intent(inout),dimension(3,*) :: tv     !< velocities array
    Real   (kind=4), Intent(inout),dimension(*)   :: tp     !< potential array
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref  !< halo id array
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tid   !< particle id array
    
    ! Input variable
    Integer(kind=4),Intent(in) :: n !< length of the arrays

    ! Local variables
    Integer(kind=4) :: i,size

#ifdef DEBUG2
    Print *,'heapsort_refI_2x3R_1xR_1xI begins'
#endif

    Call heapify_refI_2x3R_1xR_1xI(size,n,tref,tx,tv,tp,tid)

    Do i = n, 2, -1
       Call swap_refI_2x3R_1xR_1xI(1,i,tref,tx,tv,tp,tid)
       size = size - 1
       Call sink_refI_2x3R_1xR_1xI(size,1,tref,tx,tv,tp,tid)
    End Do

#ifdef DEBUG2
    Print *,'heapsort_refI_2x3R_1xR_1xI ends'
#endif

  End Subroutine heapsort_refI_2x3R_1xR_1xI


  !=======================================================================
  !> Heapify routine for heapsort
  Subroutine heapify_refI_2x3R_1xR_1xI(size,n,tref,tx,tv,tp,tid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid ! structure id and particle id
    Real   (kind=4), intent(inout), dimension(*) :: tp

    !Output variable
    Integer(kind=4),Intent(out)                :: size

    ! Input variable
    Integer(kind=4),Intent(in) :: n

    ! Local variables
    Integer(kind=4) :: i

#ifdef DEBUG2
    Print *,'heapify_refI_2x3R_1xR_1xI begins'
#endif

    size = n
    Do i = n/2, 1, -1
       Call sink_refI_2x3R_1xR_1xI(size,i,tref,tx,tv,tp,tid)
    End Do

#ifdef DEBUG2
    Print *,'heapify_refI_2x3R_1xR_1xI ends'
#endif

  End Subroutine heapify_refI_2x3R_1xR_1xI


  !=======================================================================
  !> Sink routine for heapsort
  Recursive Subroutine sink_refI_2x3R_1xR_1xI(size,i,tref,tx,tv,tp, tid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid ! structure id and particle id
    Real   (kind=4), intent(inout), dimension(*) :: tp

    ! Input variable
    Integer(kind=4),Intent(in) :: size,i

    ! Local variables
    Integer(kind=4) :: l, r, max

#ifdef DEBUG2
    Print *,'sink_refI_2x3R_1xR_1xI begins'
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
       Call swap_refI_2x3R_1xR_1xI(i,max,tref,tx,tv,tp,tid)
       Call sink_refI_2x3R_1xR_1xI(size,max,tref,tx,tv,tp,tid)
    End If

#ifdef DEBUG2
    Print *,'sink_refI_2x3R_1xR_1xI ends'
#endif

  End Subroutine sink_refI_2x3R_1xR_1xI


  !=======================================================================
  !> Swap routine for heapsort
  Subroutine swap_refI_2x3R_1xR_1xI(i,j,tref,tx,tv,tp,tid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid ! structure id and particle id
    Real   (kind=4), intent(inout), dimension(*) :: tp

    ! Input variable
    Integer(kind=4),Intent(in) :: i,j

    ! Local variables
    Integer(kind=PRI) :: tmpi
    Real(kind=4), dimension(3) :: tmpr
    Real(kind=4) :: tmps

#ifdef DEBUG2
    Print *,'swap_refI_2x3R_1xR_1xI begins'
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

    tmps = tp(i)
    tp(i) = tp(j)
    tp(j) = tmps

#ifdef DEBUG2
    Print *,'swap_refI_2x3R_1xR_1xI ends'
#endif

  End Subroutine swap_refI_2x3R_1xR_1xI


  !=======================================================================
  !> Heap sort as described in Introduction to algorithms (Cormen, Leiserson, Rivest, Stein)
  Subroutine heapsort_refI_2x3R_1xI(n,tref,tx,tv,tid)

    Implicit None

    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx     !< positions array
    Real   (kind=4), Intent(inout),dimension(3,*) :: tv     !< velocities array
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref  !< halo id array
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tid   !< particle id array
 
    ! Input variable
    Integer(kind=4),Intent(in) :: n  !< length of the arrays

    ! Local variables
    Integer(kind=4) :: i,size

#ifdef DEBUG2
    Print *,'heapsort_refI_2x3R_1xI begins'
#endif

    Call heapify_refI_2x3R_1xI(size,n,tref,tx,tv,tid)

    Do i = n, 2, -1
       Call swap_refI_2x3R_1xI(1,i,tref,tx,tv,tid)
       size = size - 1
       Call sink_refI_2x3R_1xI(size,1,tref,tx,tv,tid)
    End Do

#ifdef DEBUG2
    Print *,'heapsort_refI_2x3R_1xI ends'
#endif

  End Subroutine heapsort_refI_2x3R_1xI


  !=======================================================================
  !> Heapify routine for heapsort
  Subroutine heapify_refI_2x3R_1xI(size,n,tref,tx,tv,tid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid ! structure id and particle id

    !Output variable
    Integer(kind=4),Intent(out) :: size

    ! Input variable
    Integer(kind=4),Intent(in) :: n 

    ! Local variables
    Integer(kind=4) :: i

#ifdef DEBUG2
    Print *,'heapify_refI_2x3R_1xI begins'
#endif

    size = n
    Do i = n/2, 1, -1
       Call sink_refI_2x3R_1xI(size,i,tref,tx,tv,tid)
    End Do

#ifdef DEBUG2
    Print *,'heapify_refI_2x3R_1xI ends'
#endif

  End Subroutine heapify_refI_2x3R_1xI


  !=======================================================================
  !> Sink routine for heapsort
  Recursive Subroutine sink_refI_2x3R_1xI(size,i,tref,tx,tv,tid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid ! structure id and particle id

    ! Input variable
    Integer(kind=4),Intent(in) :: size,i

    ! Local variables
    Integer(kind=4) :: l, r, max

#ifdef DEBUG2
    Print *,'sink_refI_2x3R_1xI begins'
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
       Call swap_refI_2x3R_1xI(i,max,tref,tx,tv,tid)
       Call sink_refI_2x3R_1xI(size,max,tref,tx,tv,tid)
    End If

#ifdef DEBUG2
    Print *,'sink_refI_2x3R_1xI ends'
#endif

  End Subroutine sink_refI_2x3R_1xI


  !=======================================================================
  !> Swap routine for heapsort
  Subroutine swap_refI_2x3R_1xI(i,j,tref,tx,tv,tid)

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
    Print *,'swap_refI_2x3R_1xI begins'
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
    Print *,'swap_refI_2x3R_1xI ends'
#endif

  End Subroutine swap_refI_2x3R_1xI


  !=======================================================================
  !> Heap sort as described in Introduction to algorithms (Cormen, Leiserson, Rivest, Stein)
  Subroutine heapsort_refI_3x3R_1xR_1xI(n,tref,tx,tv,tf,tp,tid)

    Implicit None

    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*)  :: tx     !< positions array
    Real   (kind=4), Intent(inout),dimension(3,*)  :: tv     !< velocities array
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref   !< halo id array
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tid    !< particle id array
    Real   (kind=4), intent(inout), dimension(3,*) :: tf     !< gravitational field array
    Real   (kind=4), intent(inout), dimension(*)   :: tp     !< potential array
    
    ! Input variable
    Integer(kind=4),Intent(in) :: n  !< length of the arrays

    ! Local variables
    Integer(kind=4) :: i,size

#ifdef DEBUG2
    Print *,'heapsort_refI_3x3R_1xR_1xI begins'
#endif

    Call heapify_refI_3x3R_1xR_1xI(size,n,tref,tx,tv,tf,tp,tid)

    Do i = n, 2, -1
       Call swap_refI_3x3R_1xR_1xI(1,i,tref,tx,tv,tf,tp,tid)
       size = size - 1
       Call sink_refI_3x3R_1xR_1xI(size,1,tref,tx,tv,tf,tp,tid)
    End Do

#ifdef DEBUG2
    Print *,'heapsort_refI_3x3R_1xR_1xI ends'
#endif

  End Subroutine heapsort_refI_3x3R_1xR_1xI


  !=======================================================================
  !> Heapify routine for heapsort
  Subroutine heapify_refI_3x3R_1xR_1xI(size,n,tref,tx,tv,tf,tp,tid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Real   (kind=4), intent(inout), dimension(3,*) :: tf  ! force
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid ! structure id and particle id
    Real   (kind=4), intent(inout), dimension(*) :: tp

    !Output variable
    Integer(kind=4),Intent(out)                :: size

    ! Input variable
    Integer(kind=4),Intent(in) :: n 

    ! Local variables
    Integer(kind=4) :: i

#ifdef DEBUG2
    Print *,'heapify_refI_3x3R_1xR_1xI begins'
#endif

    size = n
    Do i = n/2, 1, -1
       Call sink_refI_3x3R_1xR_1xI(size,i,tref,tx,tv,tf,tp,tid)
    End Do

#ifdef DEBUG2
    Print *,'heapify_refI_3x3R_1xR_1xI ends'
#endif

  End Subroutine heapify_refI_3x3R_1xR_1xI


  !=======================================================================
  !> Sink routine for heapsort
  Recursive Subroutine sink_refI_3x3R_1xR_1xI(size,i,tref,tx,tv,tf,tp,tid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Real   (kind=4), intent(inout), dimension(3,*) :: tf  ! force
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid ! structure id and particle id
    Real   (kind=4), intent(inout), dimension(*) :: tp

    ! Input variable
    Integer(kind=4),Intent(in) :: size,i

    ! Local variables
    Integer(kind=4) :: l, r, max

#ifdef DEBUG2
    Print *,'sink_refI_3x3R_1xR_1xI begins'
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
       Call swap_refI_3x3R_1xR_1xI(i,max,tref,tx,tv,tf,tp,tid)
       Call sink_refI_3x3R_1xR_1xI(size,max,tref,tx,tv,tf,tp,tid)
    End If

#ifdef DEBUG2
    Print *,'sink_refI_3x3R_1xR_1xI ends'
#endif

  End Subroutine sink_refI_3x3R_1xR_1xI


  !=======================================================================
  !> Swap routine for heapsort
  Subroutine swap_refI_3x3R_1xR_1xI(i,j,tref,tx,tv,tf,tp,tid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Real   (kind=4), intent(inout), dimension(3,*) :: tf  ! force
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid ! structure id and particle id
    Real   (kind=4), intent(inout), dimension(*) :: tp

    ! Input variable
    Integer(kind=4),Intent(in) :: i,j

    ! Local variables
    Integer(kind=PRI) :: tmpi
    Real(kind=4), dimension(3) :: tmpr
    Real(kind=4) :: tmps

#ifdef DEBUG2
    Print *,'swap_refI_3x3R_1xR_1xI begins'
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

    tmpr(:) = tf(:,i)
    tf(:,i) = tf(:,j)
    tf(:,j) = tmpr(:)

    tmps = tp(i)
    tp(i) = tp(j)
    tp(j) = tmps

#ifdef DEBUG2
    Print *,'swap_refI_3x3R_1xR_1xI ends'
#endif

  End Subroutine swap_refI_3x3R_1xR_1xI


  !=======================================================================
  !> Heap sort as described in Introduction to algorithms (Cormen, Leiserson, Rivest, Stein)
  Subroutine heapsort_refI_3x3R_1xI(n,tref,tx,tv,tf,tid)

    Implicit None

    ! Input/Output variables
    Real   (kind=4), intent(inout),dimension(3,*) :: tx     !< positions array
    Real   (kind=4), intent(inout),dimension(3,*) :: tv     !< velocities array
    Real   (kind=4), intent(inout), dimension(3,*) :: tf    !< gravitational field array
    Integer(kind=PRI),intent(inout),dimension(*)  :: tref   !< halo id array
    Integer(kind=PRI),intent(inout),dimension(*)  :: tid    !< particle id array
    
    ! Input variable
    Integer(kind=4),Intent(in) :: n  !< length of the arrays

    ! Local variables
    Integer(kind=4) :: i,size

#ifdef DEBUG2
    Print *,'heapsort_refI_3x3R_1xI begins'
#endif

    Call heapify_refI_3x3R_1xI(size,n,tref,tx,tv,tf,tid)

    Do i = n, 2, -1
       Call swap_refI_3x3R_1xI(1,i,tref,tx,tv,tf,tid)
       size = size - 1
       Call sink_refI_3x3R_1xI(size,1,tref,tx,tv,tf,tid)
    End Do

#ifdef DEBUG2
    Print *,'heapsort_refI_3x3R_1xI ends'
#endif

  End Subroutine heapsort_refI_3x3R_1xI


  !=======================================================================
  !> Heapify routine for heapsort
  Subroutine heapify_refI_3x3R_1xI(size,n,tref,tx,tv,tf,tid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Real   (kind=4), intent(inout), dimension(3,*) :: tf  ! force
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid ! structure id and particle id

    !Output variable
    Integer(kind=4),Intent(out)                :: size

    ! Input variable
    Integer(kind=4),Intent(in) ::n  !< length of the arrays

    ! Local variables
    Integer(kind=4) :: i

#ifdef DEBUG2
    Print *,'heapify_refI_3x3R_1xI begins'
#endif

    size = n
    Do i = n/2, 1, -1
       Call sink_refI_3x3R_1xI(size,i,tref,tx,tv,tf,tid)
    End Do

#ifdef DEBUG2
    Print *,'heapify_refI_3x3R_1xI ends'
#endif

  End Subroutine heapify_refI_3x3R_1xI


  !=======================================================================
  !> Sink routine for heapsort
  Recursive Subroutine sink_refI_3x3R_1xI(size,i,tref,tx,tv,tf,tid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Real   (kind=4), intent(inout), dimension(3,*) :: tf  ! force
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid ! structure id and particle id

    ! Input variable
    Integer(kind=4),Intent(in) :: size,i

    ! Local variables
    Integer(kind=4) :: l, r, max

#ifdef DEBUG2
    Print *,'sink_refI_3x3R_1xI begins'
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
       Call swap_refI_3x3R_1xI(i,max,tref,tx,tv,tf,tid)
       Call sink_refI_3x3R_1xI(size,max,tref,tx,tv,tf,tid)
    End If

#ifdef DEBUG2
    Print *,'sink_refI_3x3R_1xI ends'
#endif

  End Subroutine sink_refI_3x3R_1xI


  !=======================================================================
  !> Swap routine for heapsort
  Subroutine swap_refI_3x3R_1xI(i,j,tref,tx,tv,tf,tid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Real   (kind=4), intent(inout), dimension(3,*) :: tf  ! force
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid ! structure id and particle id

    ! Input variable
    Integer(kind=4),Intent(in) :: i,j

    ! Local variables
    Integer(kind=PRI) :: tmpi
    Real(kind=4), dimension(3) :: tmpr

#ifdef DEBUG2
    Print *,'swap_refI_3x3R_1xI ends'
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

    tmpr(:) = tf(:,i)
    tf(:,i) = tf(:,j)
    tf(:,j) = tmpr(:)

#ifdef DEBUG2
    Print *,'swap_refI_3x3R_1xI ends'
#endif

  End Subroutine swap_refI_3x3R_1xI


  !=======================================================================  
  !> Heap sort as described in Introduction to algorithms (Cormen, Leiserson, Rivest, Stein)
  Subroutine heapsort_refI_2x3R_1xR_2xI(n,tref,tx,tv,tp,tid,trid)

    Implicit None

    ! Input/Output variables
    Real   (kind=4), intent(inout),dimension(3,*) :: tx     !< positions array   
    Real   (kind=4), intent(inout),dimension(3,*) :: tv     !< velocities array
    Integer(kind=PRI),intent(inout),dimension(*)  :: tref   !< halo id array
    Integer(kind=PRI),intent(inout),dimension(*)  :: tid    !< particle id array
    Integer(kind=PRI),intent(inout),dimension(*)  :: trid   !< ramses particle id array (for cone)
    Real   (kind=4), intent(inout), dimension(*)  :: tp     !< potential array
    
    ! Input variable
    Integer(kind=4),Intent(in) :: n  !< length of the arrays

    ! Local variables
    Integer(kind=4) :: i,size

#ifdef DEBUG2
    Print *,'heapsort_refI_2x3R_1xR_2xI begins'
#endif

    Call heapify_refI_2x3R_1xR_2xI(size,n,tref,tx,tv,tp,tid,trid)

    Do i = n, 2, -1
       Call swap_refI_2x3R_1xR_2xI(1,i,tref,tx,tv,tp,tid,trid)
       size = size - 1
       Call sink_refI_2x3R_1xR_2xI(size,1,tref,tx,tv,tp,tid,trid)
    End Do

#ifdef DEBUG2
    Print *,'heapsort_refI_2x3R_1xR_2xI ends'
#endif

  End Subroutine heapsort_refI_2x3R_1xR_2xI


  !=======================================================================  
  !> Heapify routine for heapsort
  Subroutine heapify_refI_2x3R_1xR_2xI(size,n,tref,tx,tv,tp,tid,trid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid,trid ! structure id and particle id
    Real   (kind=4), intent(inout), dimension(*) :: tp

    !Output variable
    Integer(kind=4),Intent(out)                :: size

    ! Input variable
    Integer(kind=4),Intent(in) :: n 

    ! Local variables
    Integer(kind=4) :: i

#ifdef DEBUG2
    Print *,'heapify_refI_2x3R_1xR_2xI begins'
#endif

    size = n
    Do i = n/2, 1, -1
       Call sink_refI_2x3R_1xR_2xI(size,i,tref,tx,tv,tp,tid,trid)
    End Do

#ifdef DEBUG2
    Print *,'heapify_refI_2x3R_1xR_2xI ends'
#endif

  End Subroutine heapify_refI_2x3R_1xR_2xI


  !=======================================================================
  !> Sink routine for heapsort
  Recursive Subroutine sink_refI_2x3R_1xR_2xI(size,i,tref,tx,tv,tp,tid,trid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid,trid ! structure id and particle id
    Real   (kind=4), intent(inout), dimension(*) :: tp

    ! Input variable
    Integer(kind=4),Intent(in) :: size,i

    ! Local variables
    Integer(kind=4) :: l, r, max

#ifdef DEBUG2
    Print *,'sink_refI_2x3R_1xR_2xI begins'
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
       Call swap_refI_2x3R_1xR_2xI(i,max,tref,tx,tv,tp,tid,trid)
       Call sink_refI_2x3R_1xR_2xI(size,max,tref,tx,tv,tp,tid,trid)
    End If

#ifdef DEBUG2
    Print *,'sink_refI_2x3R_1xR_2xI ends'
#endif

  End Subroutine sink_refI_2x3R_1xR_2xI


  !=======================================================================
  !> Swap routine for heapsort
  Subroutine swap_refI_2x3R_1xR_2xI(i,j,tref,tx,tv,tp,tid,trid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid,trid ! structure id and particle id
    Real   (kind=4), intent(inout), dimension(*) :: tp

    ! Input variable
    Integer(kind=4),Intent(in) :: i,j

    ! Local variables
    Integer(kind=PRI) :: tmpi
    Real(kind=4), dimension(3) :: tmpr
    Real(kind=4) :: tmps

#ifdef DEBUG2
    Print *,'swap_refI_2x3R_1xR_2xI begins'
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

    tmps = tp(i)
    tp(i) = tp(j)
    tp(j) = tmps

    tmpi = trid(i)
    trid(i) = trid(j)
    trid(j) = tmpi

#ifdef DEBUG2
    Print *,'swap_refI_2x3R_1xR_2xI ends'
#endif

  End Subroutine swap_refI_2x3R_1xR_2xI


  !=======================================================================
  !> Heap sort as described in Introduction to algorithms (Cormen, Leiserson, Rivest, Stein)
  Subroutine heapsort_refI_2x3R_2xI(n,tref,tx,tv,tid,trid)

    Implicit None

    ! Input/Output variables
    Real   (kind=4), intent(inout),dimension(3,*) :: tx     !< positions array   
    Real   (kind=4), intent(inout),dimension(3,*) :: tv     !< velocities array
    Integer(kind=PRI),intent(inout),dimension(*)  :: tref   !< halo id array
    Integer(kind=PRI),intent(inout),dimension(*)  :: tid    !< particle id array
    Integer(kind=PRI),intent(inout),dimension(*)  :: trid   !< ramses particle id array (for cone)
    
    ! Input variable
    Integer(kind=4),Intent(in) :: n  !< length of the arrays

    ! Local variables
    Integer(kind=4) :: i,size

#ifdef DEBUG2
    Print *,'heapsort_refI_2x3R_2xI begins'
#endif

    Call heapify_refI_2x3R_2xI(size,n,tref,tx,tv,tid,trid)

    Do i = n, 2, -1
       Call swap_refI_2x3R_2xI(1,i,tref,tx,tv,tid,trid)
       size = size - 1
       Call sink_refI_2x3R_2xI(size,1,tref,tx,tv,tid,trid)
    End Do

#ifdef DEBUG2
    Print *,'heapsort_refI_2x3R_2xI ends'
#endif

  End Subroutine heapsort_refI_2x3R_2xI


  !=======================================================================
  !> Heapify routine for heapsort
  Subroutine heapify_refI_2x3R_2xI(size,n,tref,tx,tv,tid,trid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid,trid ! structure id and particle id

    !Output variable
    Integer(kind=4),Intent(out)                :: size

    ! Input variable
    Integer(kind=4),Intent(in) :: n 

    ! Local variables
    Integer(kind=4) :: i

#ifdef DEBUG2
    Print *,'heapify_refI_2x3R_2xI begins'
#endif

    size = n
    Do i = n/2, 1, -1
       Call sink_refI_2x3R_2xI(size,i,tref,tx,tv,tid,trid)
    End Do

#ifdef DEBUG2
    Print *,'heapify_refI_2x3R_2xI ends'
#endif

  End Subroutine heapify_refI_2x3R_2xI


  !=======================================================================
  !> Sink routine for heapsort
  Recursive Subroutine sink_refI_2x3R_2xI(size,i,tref,tx,tv,tid,trid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid,trid ! structure id and particle id

    ! Input variable
    Integer(kind=4),Intent(in) :: size,i

    ! Local variables
    Integer(kind=4) :: l, r, max

#ifdef DEBUG2
    Print *,'sink_refI_2x3R_2xI begins'
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
       Call swap_refI_2x3R_2xI(i,max,tref,tx,tv,tid,trid)
       Call sink_refI_2x3R_2xI(size,max,tref,tx,tv,tid,trid)
    End If

#ifdef DEBUG2
    Print *,'sink_refI_2x3R_2xI ends'
#endif

  End Subroutine sink_refI_2x3R_2xI


  !=======================================================================
  !> Swap routine for heapsort
  Subroutine swap_refI_2x3R_2xI(i,j,tref,tx,tv,tid,trid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid,trid ! structure id and particle id

    ! Input variable
    Integer(kind=4),Intent(in) :: i,j

    ! Local variables
    Integer(kind=PRI) :: tmpi
    Real(kind=4), dimension(3) :: tmpr

#ifdef DEBUG2
    Print *,'swap_refI_2x3R_2xI begins'
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

    tmpi = trid(i)
    trid(i) = trid(j)
    trid(j) = tmpi

#ifdef DEBUG2
    Print *,'swap_refI_2x3R_2xI ends'
#endif

  End Subroutine swap_refI_2x3R_2xI


  !=======================================================================
  !> Heap sort as described in Introduction to algorithms (Cormen, Leiserson, Rivest, Stein)
  Subroutine heapsort_refI_3x3R_1xR_2xI(n,tref,tx,tv,tf,tp,tid,trid)

    Implicit None

    ! Input/Output variables
    Real   (kind=4), intent(inout),dimension(3,*) :: tx     !< positions array   
    Real   (kind=4), intent(inout),dimension(3,*) :: tv     !< velocities array
    Integer(kind=PRI),intent(inout),dimension(*)  :: tref   !< halo id array
    Integer(kind=PRI),intent(inout),dimension(*)  :: tid    !< particle id array
    Integer(kind=PRI),intent(inout),dimension(*)  :: trid   !< ramses particle id array (for cone)
    Real   (kind=4), intent(inout), dimension(*)  :: tp     !< potential array
    Real   (kind=4), intent(inout), dimension(3,*) :: tf    !< gravitational field array

    ! Input variable
    Integer(kind=4),Intent(in) :: n  !< length of the arrays

    ! Local variables
    Integer(kind=4) :: i,size

#ifdef DEBUG2
    Print *,'heapsort_refI_3x3R_1xR_2xI begins'
#endif

    Call heapify_refI_3x3R_1xR_2xI(size,n,tref,tx,tv,tf,tp,tid,trid)

    Do i = n, 2, -1
       Call swap_refI_3x3R_1xR_2xI(1,i,tref,tx,tv,tf,tp,tid,trid)
       size = size - 1
       Call sink_refI_3x3R_1xR_2xI(size,1,tref,tx,tv,tf,tp,tid,trid)
    End Do

#ifdef DEBUG2
    Print *,'heapsort_refI_3x3R_1xR_2xI ends'
#endif

  End Subroutine heapsort_refI_3x3R_1xR_2xI


  !=======================================================================
  !> Heapify routine for heapsort
  Subroutine heapify_refI_3x3R_1xR_2xI(size,n,tref,tx,tv,tf,tp,tid,trid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Real   (kind=4), intent(inout), dimension(3,*) :: tf  ! force
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid, trid ! structure id and particle id
    Real   (kind=4), intent(inout), dimension(*) :: tp

    !Output variable
    Integer(kind=4),Intent(out)                :: size

    ! Input variable
    Integer(kind=4),Intent(in) :: n

    ! Local variables
    Integer(kind=4) :: i

#ifdef DEBUG2
    Print *,'heapify_refI_3x3R_1xR_2xI begins'
#endif

    size = n
    Do i = n/2, 1, -1
       Call sink_refI_3x3R_1xR_2xI(size,i,tref,tx,tv,tf,tp,tid,trid)
    End Do

#ifdef DEBUG2
    Print *,'heapify_refI_3x3R_1xR_2xI ends'
#endif

  End Subroutine heapify_refI_3x3R_1xR_2xI


  !=======================================================================
  !> Sink routine for heapsort
  Recursive Subroutine sink_refI_3x3R_1xR_2xI(size,i,tref,tx,tv,tf,tp,tid,trid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Real   (kind=4), intent(inout), dimension(3,*) :: tf  ! force
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid, trid ! structure id and particle id
    Real   (kind=4), intent(inout), dimension(*) :: tp

    ! Input variable
    Integer(kind=4),Intent(in) :: size,i

    ! Local variables
    Integer(kind=4) :: l, r, max

#ifdef DEBUG2
    Print *,'sink_refI_3x3R_1xR_2xI begins'
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
       Call swap_refI_3x3R_1xR_2xI(i,max,tref,tx,tv,tf,tp,tid,trid)
       Call sink_refI_3x3R_1xR_2xI(size,max,tref,tx,tv,tf,tp,tid,trid)
    End If

#ifdef DEBUG2
    Print *,'sink_refI_3x3R_1xR_2xI ends'
#endif

  End Subroutine sink_refI_3x3R_1xR_2xI


  !=======================================================================
  !> Swap routine for heapsort
  Subroutine swap_refI_3x3R_1xR_2xI(i,j,tref,tx,tv,tf,tp,tid,trid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Real   (kind=4), intent(inout), dimension(3,*) :: tf  ! force
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid, trid ! structure id and particle id
    Real   (kind=4), intent(inout), dimension(*) :: tp

    ! Input variable
    Integer(kind=4),Intent(in) :: i,j

    ! Local variables
    Integer(kind=PRI) :: tmpi
    Real(kind=4), dimension(3) :: tmpr
    Real(kind=4) :: tmps

#ifdef DEBUG2
    Print *,'swap_refI_3x3R_1xR_2xI begins'
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

    tmpr(:) = tf(:,i)
    tf(:,i) = tf(:,j)
    tf(:,j) = tmpr(:)

    tmps = tp(i)
    tp(i) = tp(j)
    tp(j) = tmps

    tmpi = trid(i)
    trid(i) = trid(j)
    trid(j) = tmpi

#ifdef DEBUG2
    Print *,'swap_refI_3x3R_1xR_2xI ends'
#endif

  End Subroutine swap_refI_3x3R_1xR_2xI


  !=======================================================================
  !> Heap sort as described in Introduction to algorithms (Cormen, Leiserson, Rivest, Stein)
  Subroutine heapsort_refI_3x3R_2xI(n,tref,tx,tv,tf,tid,trid)

    Implicit None

    ! Input/Output variables
    Real   (kind=4), intent(inout),dimension(3,*) :: tx     !< positions array   
    Real   (kind=4), intent(inout),dimension(3,*) :: tv     !< velocities array
    Integer(kind=PRI),intent(inout),dimension(*)  :: tref   !< halo id array
    Integer(kind=PRI),intent(inout),dimension(*)  :: tid    !< particle id array
    Integer(kind=PRI),intent(inout),dimension(*)  :: trid   !< ramses particle id array (for cone)
    Real   (kind=4), intent(inout), dimension(3,*) :: tf    !< gravitational field array

    ! Input variable
    Integer(kind=4),Intent(in) :: n  !< length of the arrays

    ! Local variables
    Integer(kind=4) :: i,size

#ifdef DEBUG2
    Print *,'heapsort_refI_3x3R_2xI begins'
#endif

    Call heapify_refI_3x3R_2xI(size,n,tref,tx,tv,tf,tid,trid)

    Do i = n, 2, -1
       Call swap_refI_3x3R_2xI(1,i,tref,tx,tv,tf,tid,trid)
       size = size - 1
       Call sink_refI_3x3R_2xI(size,1,tref,tx,tv,tf,tid,trid)
    End Do

#ifdef DEBUG2
    Print *,'heapsort_refI_3x3R_2xI ends'
#endif

  End Subroutine heapsort_refI_3x3R_2xI


  !=======================================================================
  !> Heapify routine for heapsort
  Subroutine heapify_refI_3x3R_2xI(size,n,tref,tx,tv,tf,tid,trid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Real   (kind=4), intent(inout), dimension(3,*) :: tf  ! force
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid, trid ! structure id and particle id

    !Output variable
    Integer(kind=4),Intent(out)                :: size

    ! Input variable
    Integer(kind=4),Intent(in) :: n

    ! Local variables
    Integer(kind=4) :: i

#ifdef DEBUG2
    Print *,'heapify_refI_3x3R_2xI begins'
#endif

    size = n
    Do i = n/2, 1, -1
       Call sink_refI_3x3R_2xI(size,i,tref,tx,tv,tf,tid,trid)
    End Do

#ifdef DEBUG2
    Print *,'heapify_refI_3x3R_2xI ends'
#endif

  End Subroutine heapify_refI_3x3R_2xI


  !=======================================================================
  !> Sink routine for heapsort
  Recursive Subroutine sink_refI_3x3R_2xI(size,i,tref,tx,tv,tf,tid,trid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Real   (kind=4), intent(inout), dimension(3,*) :: tf  ! force
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid, trid ! structure id and particle id

    ! Input variable
    Integer(kind=4),Intent(in) :: size,i

    ! Local variables
    Integer(kind=4) :: l, r, max

#ifdef DEBUG2
    Print *,'sink_refI_3x3R_2xI begins'
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
       Call swap_refI_3x3R_2xI(i,max,tref,tx,tv,tf,tid,trid)
       Call sink_refI_3x3R_2xI(size,max,tref,tx,tv,tf,tid,trid)
    End If

#ifdef DEBUG2
    Print *,'sink_refI_3x3R_2xI ends'
#endif

  End Subroutine sink_refI_3x3R_2xI


  !=======================================================================
  !> Swap routine for heapsort
  Subroutine swap_refI_3x3R_2xI(i,j,tref,tx,tv,tf,tid,trid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Real   (kind=4), intent(inout), dimension(3,*) :: tf  ! force
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref, tid, trid ! structure id and particle id

    ! Input variable
    Integer(kind=4),Intent(in) :: i,j

    ! Local variables
    Integer(kind=PRI) :: tmpi
    Real(kind=4), dimension(3) :: tmpr

#ifdef DEBUG2
    Print *,'swap_refI_3x3R_2xI begins'
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

    tmpr(:) = tf(:,i)
    tf(:,i) = tf(:,j)
    tf(:,j) = tmpr(:)

    tmpi = trid(i)
    trid(i) = trid(j)
    trid(j) = tmpi

#ifdef DEBUG2
    Print *,'swap_refI_3x3R_2xI ends'
#endif

  End Subroutine swap_refI_3x3R_2xI

End Module modsort
