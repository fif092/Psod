!==============================================================================
! Project: pFoF
! File: tools/conepartcreator/src/modvariables.f90
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

Module modvariables

  Use modconstant
  Use mpi

  Type(Type_parameter_conecreator_part) :: param

  Integer(kind=PRI), dimension(:), allocatable :: id        ! id of the particles in the halo just read
  Real(kind=4), dimension(:,:), allocatable :: pos  ! position of the particles in the halo just read
  Real(kind=4), dimension(:,:), allocatable :: vel  ! velocities of the particles in the halo just read
  Integer(kind=PRI), dimension(:), allocatable :: ramsespartid
  Real(kind=4), dimension(:), allocatable :: pot
  Real(kind=4), dimension(:,:), allocatable :: field


  Integer(kind=4), dimension(:), allocatable :: npartcubeloc, npartcube
  Integer(kind=4), dimension(:), allocatable :: idcube

  Integer(kind=4) :: ncx
  Integer(kind=4) :: ncy
  Integer(kind=4) :: ncz
  Integer(kind=4) :: ncube

  Real(kind=8) :: hy
  Real(kind=8) :: hz

  !! MPI Variables
  Integer :: req_sumnpc
  Integer(kind=4), dimension(Mpi_Status_Size) :: mpistat

  Type(Type_info_cone_part) :: infocone
  Type(Type_info_ramses) :: inforamses

End Module modvariables
