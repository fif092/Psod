!==============================================================================
! Project: pFoF
! File: tools/conemapper/src/modvariables.f90
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

!> @file
!!This file contains some global variables declarations.

!> This module contains some global variables declarations.
Module modvariables
  Use modconstant

  Implicit None

  Type(Type_parameter_pfof_cone) :: parameters !< Input parameters read from the pfof_cone.nml file.

  Integer(kind=4) :: shell_nb  !< Number of shell files read.
  Integer(kind=4), dimension(:), allocatable :: npartcube  !< Number of particles in each shell cube.
  Integer(kind=8), dimension(:), allocatable :: npartcubefof  ! nb of particles in each pfof cube/process
  Integer(kind=4), dimension(:,:), allocatable :: pictable  ! list of shell cube indices in each pfof cube/process
  Integer(kind=4), dimension(:,:), allocatable :: neigh  ! process ID of each neighbour for each pfof cube/process
  Integer(kind=8), dimension(:), allocatable :: my_npart_tab  ! array of nb of particles in each pfof cube/process
  Integer(kind=4), dimension(3) :: nctab  ! nb of shell cubes in each dimension in the global mapping
  Integer(kind=4) :: ncube  ! nb of shell cubes
  Integer(kind=8) :: npart  ! total number of particles in the cone
  Real(kind=8) :: cubesize  ! length of the shell cube edge
  Real(kind=8), dimension(:,:), allocatable :: boundaries ! (x,y,z)min, (x,y,z)max of each pfof cube
  Logical :: fullsky  ! is it a fullsky cone?

End Module modvariables
