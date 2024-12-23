!=======================================================================
! Author: Fabrice Roy & Vincent Bouillot (LUTH/CNRS/Observatoire de Paris)
! Fabrice.Roy@obspm.fr 
!=======================================================================
Module modsortinterf

  Use modconstant, only : PRI

  Private

  Public :: quicksort, heapsort

  Interface quicksort
     Module procedure quicksort_refI_2x3R_1xI  
     Module procedure quicksort_refI_2x3R_1xR_1xI  
     Module procedure quicksort_refI_3x3R_1xR_1xI  
     Module procedure quicksort_refI_3x3R_1xI  
  End Interface quicksort

  
  Interface heapsort
     Module procedure heapsort_refI_2x3R_1xI
     Module procedure heapsort_refI_2x3R_1xR_1xI     
     Module procedure heapsort_refI_3x3R_1xR_1xI     
     Module procedure heapsort_refI_3x3R_1xI     
  End Interface heapsort
  

Contains

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

  Recursive Subroutine quicksort_refI_2x3R_1xR_1xI(p,r,tref,tx,tv,tp,tid)

    Implicit None
    
    ! Input parameters
    Integer(kind=4), Intent(in) :: p     !< index of the first element of the array to be sorted
    Integer(kind=4), Intent(in) :: r     !< index of the last element of the array to be sorted
    
    ! Input parameters modified in subroutine
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx     !< positions array
    Real   (kind=4), Intent(inout),dimension(3,*) :: tv     !< velocities array
    Real   (kind=4), Intent(inout),dimension(*) :: tp     !< potential array
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref  !< halo id array
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tid   !< particle id array
    
    ! Local parameters
    Integer(kind=4) :: q
    
#ifdef DEBUG2
    Print *,'quicksort_refI_2x3R_1xI begins'
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

  Recursive Subroutine quicksort_refI_3x3R_1xR_1xI(p,r,tref,tx,tv,tf,tp,tid)

    Implicit None
    
    ! Input parameters
    Integer(kind=4), Intent(in) :: p     !< index of the first element of the array to be sorted
    Integer(kind=4), Intent(in) :: r     !< index of the last element of the array to be sorted
    
    ! Input parameters modified in subroutine
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx     !< positions array
    Real   (kind=4), Intent(inout),dimension(3,*) :: tv     !< velocities array
    Real   (kind=4), Intent(inout),dimension(3,*) :: tf     !< potential array
    Real   (kind=4), Intent(inout),dimension(*)   :: tp     !< potential array
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tref  !< halo id array
    Integer(kind=PRI),Intent(inout),dimension(*)   :: tid   !< particle id array
    
    ! Local parameters
    Integer(kind=4) :: q
    
#ifdef DEBUG2
    Print *,'quicksort_refI_2x3R_1xI begins'
#endif

    If(p<r) Then
       
       Call partition_refI_3x3R_1xR_1xI(p,r,q,tref,tx,tv,tf,tp,tid)
       Call quicksort_refI_3x3R_1xR_1xI(p,q-1,tref,tx,tv,tf,tp,tid)
       Call quicksort_refI_3x3R_1xR_1xI(q+1,r,tref,tx,tv,tf,tp,tid)

    End If

#ifdef DEBUG2
    Print *,'trirapide ends'
#endif


  End Subroutine quicksort_refI_3x3R_1xR_1xI

  !=======================================================================

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

  End Subroutine partition_refI_2x3R_1xI

  !=======================================================================

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
    Print *,'partitions'
#endif

  End Subroutine partition_refI_2x3R_1xR_1xI

  !=======================================================================
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
    Print *,'partitions'
#endif

  End Subroutine partition_refI_3x3R_1xR_1xI

  !=======================================================================
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
    Print *,'partitions'
#endif

  End Subroutine partition_refI_3x3R_1xI

  !=======================================================================
  !> Heap sort as described in Introduction to algorithms (Cormen, Leiserson, Rivest, Stein)
  Subroutine heapsort_refI_2x3R_1xR_1xI(n,tref,tx,tv,tp,tid)

    Use modmpicom
    Implicit None

    ! Input/Output variables
    Real   (kind=4), intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=PRI),intent(inout),dimension(*)  :: tref, tid ! structure id and particle id
    Real   (kind=4), intent(inout), dimension(*) :: tp
    
    ! Input variable
    Integer(kind=4),Intent(in) :: n

    ! Local variables
    Integer(kind=4) :: i,size

#ifdef DEBUG2
    Print *,'tritas ends'
#endif

    Call construction_refI_2x3R_1xR_1xI(size,n,tref,tx,tv,tp,tid)

    Do i = n, 2, -1
       Call echanger_refI_2x3R_1xR_1xI(1,i,tref,tx,tv,tp,tid)
       size = size - 1
       Call entasser_refI_2x3R_1xR_1xI(size,1,tref,tx,tv,tp,tid)
    End Do

#ifdef DEBUG2
    Print *,'tritas ends'
#endif

  End Subroutine heapsort_refI_2x3R_1xR_1xI


  !=======================================================================

  
  Subroutine construction_refI_2x3R_1xR_1xI(size,n,tref,tx,tv,tp,tid)

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
    Print *,'construction ends'
#endif

    size = n
    Do i = n/2, 1, -1
       Call entasser_refI_2x3R_1xR_1xI(size,i,tref,tx,tv,tp,tid)
    End Do

#ifdef DEBUG2
    Print *,'construction ends'
#endif

  End Subroutine construction_refI_2x3R_1xR_1xI


  !=======================================================================


  Recursive Subroutine entasser_refI_2x3R_1xR_1xI(size,i,tref,tx,tv,tp, tid)

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
       Call echanger_refI_2x3R_1xR_1xI(i,max,tref,tx,tv,tp,tid)
       Call entasser_refI_2x3R_1xR_1xI(size,max,tref,tx,tv,tp,tid)
    End If

#ifdef DEBUG2
    Print *,'entasser ends'
#endif

  End Subroutine entasser_refI_2x3R_1xR_1xI


  !=======================================================================


  Subroutine echanger_refI_2x3R_1xR_1xI(i,j,tref,tx,tv,tp,tid)

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

    tmps = tp(i)
    tp(i) = tp(j)
    tp(j) = tmps

#ifdef DEBUG2
    Print *,'echanger ends'
#endif

  End Subroutine echanger_refI_2x3R_1xR_1xI


  !=======================================================================

  Subroutine heapsort_refI_2x3R_1xI(n,tref,tx,tv,tid)

    Use modmpicom
    Implicit None

    ! Input/Output variables
    Real   (kind=4), intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Integer(kind=PRI),intent(inout),dimension(*)  :: tref, tid ! structure id and particle id
    
    ! Input variable
    Integer(kind=4),Intent(in) :: n

    ! Local variables
    Integer(kind=4) :: i,size

#ifdef DEBUG2
    Print *,'tritas ends'
#endif

    Call construction_refI_2x3R_1xI(size,n,tref,tx,tv,tid)

    Do i = n, 2, -1
       Call echanger_refI_2x3R_1xI(1,i,tref,tx,tv,tid)
       size = size - 1
       Call entasser_refI_2x3R_1xI(size,1,tref,tx,tv,tid)
    End Do

#ifdef DEBUG2
    Print *,'tritas ends'
#endif

  End Subroutine heapsort_refI_2x3R_1xI


  !=======================================================================

  Subroutine construction_refI_2x3R_1xI(size,n,tref,tx,tv,tid)

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
       Call entasser_refI_2x3R_1xI(size,i,tref,tx,tv,tid)
    End Do

#ifdef DEBUG2
    Print *,'construction ends'
#endif

  End Subroutine construction_refI_2x3R_1xI


  !=======================================================================

  Recursive Subroutine entasser_refI_2x3R_1xI(size,i,tref,tx,tv,tid)

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
       Call echanger_refI_2x3R_1xI(i,max,tref,tx,tv,tid)
       Call entasser_refI_2x3R_1xI(size,max,tref,tx,tv,tid)
    End If

#ifdef DEBUG2
    Print *,'entasser ends'
#endif

  End Subroutine entasser_refI_2x3R_1xI


  !=======================================================================

  Subroutine echanger_refI_2x3R_1xI(i,j,tref,tx,tv,tid)

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

  End Subroutine echanger_refI_2x3R_1xI


  !=======================================================================
  Subroutine heapsort_refI_3x3R_1xR_1xI(n,tref,tx,tv,tf,tp,tid)

    Use modmpicom
    Implicit None

    ! Input/Output variables
    Real   (kind=4), intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Real   (kind=4), intent(inout), dimension(3,*) :: tf  ! force
    Integer(kind=PRI),intent(inout),dimension(*)  :: tref, tid ! structure id and particle id
    Real   (kind=4), intent(inout), dimension(*) :: tp
    
    ! Input variable
    Integer(kind=4),Intent(in) :: n

    ! Local variables
    Integer(kind=4) :: i,size

#ifdef DEBUG2
    Print *,'tritas ends'
#endif

    Call construction_refI_3x3R_1xR_1xI(size,n,tref,tx,tv,tf,tp,tid)

    Do i = n, 2, -1
       Call echanger_refI_3x3R_1xR_1xI(1,i,tref,tx,tv,tf,tp,tid)
       size = size - 1
       Call entasser_refI_3x3R_1xR_1xI(size,1,tref,tx,tv,tf,tp,tid)
    End Do

#ifdef DEBUG2
    Print *,'tritas ends'
#endif

  End Subroutine heapsort_refI_3x3R_1xR_1xI


  !=======================================================================

  Subroutine construction_refI_3x3R_1xR_1xI(size,n,tref,tx,tv,tf,tp,tid)

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
    Print *,'construction ends'
#endif

    size = n
    Do i = n/2, 1, -1
       Call entasser_refI_3x3R_1xR_1xI(size,i,tref,tx,tv,tf,tp,tid)
    End Do

#ifdef DEBUG2
    Print *,'construction ends'
#endif

  End Subroutine construction_refI_3x3R_1xR_1xI


  !=======================================================================


  Recursive Subroutine entasser_refI_3x3R_1xR_1xI(size,i,tref,tx,tv,tf,tp, tid)

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
       Call echanger_refI_3x3R_1xR_1xI(i,max,tref,tx,tv,tf,tp,tid)
       Call entasser_refI_3x3R_1xR_1xI(size,max,tref,tx,tv,tf,tp,tid)
    End If

#ifdef DEBUG2
    Print *,'entasser ends'
#endif

  End Subroutine entasser_refI_3x3R_1xR_1xI


  !=======================================================================


  Subroutine echanger_refI_3x3R_1xR_1xI(i,j,tref,tx,tv,tf,tp,tid)

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

    tmpr(:) = tf(:,i)
    tf(:,i) = tf(:,j)
    tf(:,j) = tmpr(:)

    tmps = tp(i)
    tp(i) = tp(j)
    tp(j) = tmps

#ifdef DEBUG2
    Print *,'echanger ends'
#endif

  End Subroutine echanger_refI_3x3R_1xR_1xI


  !=======================================================================
  Subroutine heapsort_refI_3x3R_1xI(n,tref,tx,tv,tf,tid)

    Use modmpicom
    Implicit None

    ! Input/Output variables
    Real   (kind=4), intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Real   (kind=4), intent(inout), dimension(3,*) :: tf  ! force
    Integer(kind=PRI),intent(inout),dimension(*)  :: tref, tid ! structure id and particle id
    
    ! Input variable
    Integer(kind=4),Intent(in) :: n

    ! Local variables
    Integer(kind=4) :: i,size

#ifdef DEBUG2
    Print *,'tritas ends'
#endif

    Call construction_refI_3x3R_1xI(size,n,tref,tx,tv,tf,tid)

    Do i = n, 2, -1
       Call echanger_refI_3x3R_1xI(1,i,tref,tx,tv,tf,tid)
       size = size - 1
       Call entasser_refI_3x3R_1xI(size,1,tref,tx,tv,tf,tid)
    End Do

#ifdef DEBUG2
    Print *,'tritas ends'
#endif

  End Subroutine heapsort_refI_3x3R_1xI


  !=======================================================================

  Subroutine construction_refI_3x3R_1xI(size,n,tref,tx,tv,tf,tid)

    Implicit None
    
    ! Input/Output variables
    Real   (kind=4), Intent(inout),dimension(3,*) :: tx,tv     ! positions and velocities
    Real   (kind=4), intent(inout), dimension(3,*) :: tf  ! force
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
       Call entasser_refI_3x3R_1xI(size,i,tref,tx,tv,tf,tid)
    End Do

#ifdef DEBUG2
    Print *,'construction ends'
#endif

  End Subroutine construction_refI_3x3R_1xI


  !=======================================================================


  Recursive Subroutine entasser_refI_3x3R_1xI(size,i,tref,tx,tv,tf,tid)

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
       Call echanger_refI_3x3R_1xI(i,max,tref,tx,tv,tf,tid)
       Call entasser_refI_3x3R_1xI(size,max,tref,tx,tv,tf,tid)
    End If

#ifdef DEBUG2
    Print *,'entasser ends'
#endif

  End Subroutine entasser_refI_3x3R_1xI


  !=======================================================================


  Subroutine echanger_refI_3x3R_1xI(i,j,tref,tx,tv,tf,tid)

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

    tmpr(:) = tf(:,i)
    tf(:,i) = tf(:,j)
    tf(:,j) = tmpr(:)

#ifdef DEBUG2
    Print *,'echanger ends'
#endif

  End Subroutine echanger_refI_3x3R_1xI


End Module modsortinterf
