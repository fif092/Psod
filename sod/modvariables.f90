Module modvariables
  Use modconstant,only:Type_inforamses
  Implicit None

#ifdef LONGINT
  Integer(kind=4), parameter :: PRI=8
#else
  Integer(kind=4), parameter :: PRI=4
#endif

  Real(kind=4), dimension(:,:), allocatable :: pos
  Real(kind=4), dimension(:,:), allocatable :: vel
  Real(kind=4), dimension(:), allocatable :: pot
  Integer(kind=PRI), dimension(:), allocatable :: id
  Real(kind=4), dimension(:,:), allocatable :: for
  Real(kind=4), dimension(:), allocatable :: mass ! FT contains mass of particles

  Type(Type_inforamses) :: inforamses

End Module modvariables
