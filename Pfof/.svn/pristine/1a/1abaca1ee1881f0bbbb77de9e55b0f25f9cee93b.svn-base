!==============================================================================
! Project: pFoF
! File: pfof_cone/src/modvariables.f90
! Copyright Fabrice Roy and Vincent Bouillot (2011)
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
!!This file contains the common variables used in pFoF_cone and a subroutine used to deallocate all allocatable arrays.

!> This module contains the common variables used in pFoF_cone and a subroutine used to deallocate all allocatable arrays.
!>
!> Authors: F. Roy, V. Bouillot


Module modvariables

  Use modvarcommons
  Use modconstant
  Implicit None

  Integer(kind=4) :: shell_nb
  Integer(kind=4), dimension(:,:), allocatable :: indexcube
  Integer(kind=4), dimension(:), allocatable :: ictable
  Integer(kind=4) :: shellcubes_size
  Real(kind=8), dimension(6) :: boundaries 
  Integer(kind=4), dimension(:), allocatable :: cubepershell

  Real(kind=PR), dimension(3) :: partmin, partmax

  Integer(kind=4) :: ncube
  Integer(kind=4), dimension(:), allocatable :: npart_tab

  Type(Type_info_ramses) :: inforamses, inforamseslast
  Type(Type_info_cone_part) :: infocone, infoconelast
  Type(Type_parameter_pfof_cone) :: param
  Type(Type_common_metadata) :: common_meta

  Private
  Public :: shell_nb, &
       indexcube, &
       ictable, &
       shellcubes_size, &
       boundaries, &
       cubepershell, &
       ncube, &
       npart_tab, &
       partmin, &
       partmax, &
       inforamses, inforamseslast, &
       infocone, infoconelast, &
       param, &
       common_meta, &
       deallocall

Contains

  Subroutine deallocall()

    Implicit None

    If(Allocated(indexcube)) Deallocate(indexcube)
    If(Allocated(ictable)) Deallocate(ictable)
    If(Allocated(cubepershell)) Deallocate(cubepershell)
    If(Allocated(position)) Deallocate(position)
    If(Allocated(velocity)) Deallocate(velocity)
    If(Allocated(pfof_id)) Deallocate(pfof_id)
    If(Allocated(ramses_id)) Deallocate(ramses_id)
    If(Allocated(potential)) Deallocate(potential)
    If(Allocated(field)) Deallocate(field)

  End Subroutine deallocall

End Module modvariables
