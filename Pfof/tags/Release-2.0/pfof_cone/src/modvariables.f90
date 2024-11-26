!> @file
!!This file contains the common variables used in pFoF_cone and a subroutine used to deallocate all allocatable arrays.

!> This module contains the common variables used in pFoF_cone and a subroutine used to deallocate all allocatable arrays.
!>
!> Authors: F. Roy, V. Bouillot


Module modvariables

  Use modconstant
  Implicit None

  Integer(kind=4) :: shell_nb
  Integer(kind=4), dimension(:,:), allocatable :: indexcube
  Integer(kind=4), dimension(:), allocatable :: ictable
  Integer(kind=4) :: shellcubes_size
  Real(kind=8), dimension(6) :: boundaries 
  Integer(kind=4), dimension(:), allocatable :: cubepershell
  Integer(kind=4) :: mynpart
  Integer(kind=PRI) :: npart

  Real(kind=PR), dimension(:,:), allocatable :: pos
  Real(kind=PR), dimension(:,:), allocatable :: vel
  Integer(kind=PRI), dimension(:), allocatable :: id
  Integer(kind=PRI), dimension(:), allocatable :: structure
  Real(kind=PR), dimension(:), allocatable :: pot !(here for compatibility reason with pfof)
  Real(kind=PR), dimension(:,:), allocatable :: for !(here for compatibility reason with pfof)

  Real(kind=4), dimension(:,:), allocatable :: posf
  Real(kind=4), dimension(:,:), allocatable :: velf
  Integer(kind=PRI), dimension(:), allocatable :: idf
  Integer(kind=PRI), dimension(:), allocatable :: stf
  Real(kind=PR), dimension(:), allocatable :: potf !(here for compatibility reason with pfof)
  Real(kind=PR), dimension(:,:), allocatable :: forf !(here for compatibility reason with pfof)

  Real(kind=PR), dimension(3) :: partmin, partmax

Contains

  Subroutine deallocall()

    Implicit None

    If(Allocated(indexcube)) Deallocate(indexcube)
    If(Allocated(ictable)) Deallocate(ictable)
    If(Allocated(cubepershell)) Deallocate(cubepershell)
    If(Allocated(pos)) Deallocate(pos)
    If(Allocated(vel)) Deallocate(vel)
    If(Allocated(id)) Deallocate(id)

  End Subroutine deallocall

End Module modvariables
