! Copyright (c) 2007-2016 CNRS, Fabrice Roy, Vincent Bouillot
! Author: Fabrice Roy (LUTH/CNRS/PSL), fabrice.roy@obspm.fr
! Vincent Bouillot (LUTH/CNRS/PSL)
! Distributed under the CeCILL 2.0 license. For full terms see the file LICENSE.

!> @file 
!! This file contains the definition of the 3 arrays containing the position, 
!! the velocity and the id of the particles.
!! @author Fabrice Roy
!! @author Vincent Bouillot

Module modvariables

  Use modconstant
  Implicit None

  Integer(kind=4) :: mynpart       !< local particle number
  Integer(kind=4) :: nres          !< 1-D resolution: number of grid points in each dimension ; nres = 2 ** lmin
  Integer(kind=PRI) :: npart       !< total particle number: npart = nres ** 3
  Integer(kind=PRI) :: ngrid       !< total number of grid points: ngrid = nres ** 3
  Real   (kind=4)   :: xmin, xmax, ymin, ymax, zmin, zmax  !< min and max (x,y,z) for each process

  Integer(kind=4) :: neighbours(6)  !< array containing the id of the neighbours of the local process in the MpiCube communicator

  Real(kind=4), dimension(:,:), allocatable :: pos                !< position of the particles
  Real(kind=4), dimension(:,:), allocatable :: vel                !< velocity of the particles
  Real(kind=4), dimension(:,:), allocatable :: for                !< force on the particles
  Real(kind=4), dimension(:), allocatable :: pot                  !< potential at the position of the particles (optional)
  Integer(kind=PRI), dimension(:), allocatable :: id            !< RAMSES id of the particles; 
                                                                !! Integer8 if particle number exceeds 2^31-1
  Integer(kind=PRI), dimension(:), allocatable :: structure     !< Halo ID of the partickes;
                                                                !! Integer8 if particle number exceeds 2^31-1

  Real(kind=4), dimension(:,:), allocatable :: posf               !< position of the particles
  Real(kind=4), dimension(:,:), allocatable :: velf               !< velocity of the particles
  Real(kind=4), dimension(:,:), allocatable :: forf               !< force on the particles
  Real(kind=4), dimension(:), allocatable :: potf                 !< potential at the position of the particles (optional)
  Integer(kind=PRI), dimension(:), allocatable :: idf           !< RAMSES id of the particles; 
                                                                !! Integer8 if particle number exceeds 2^31-1
  Integer(kind=PRI), dimension(:), allocatable :: stf           !< Halo ID of the partickes;
                                                                !! Integer8 if particle number exceeds 2^31-1

  Integer(kind=PRI), dimension(:), allocatable :: ramsesid, ramsesidf ! just for compatibility
  Character(len=50) :: origin 

  ! information read from Ramses info file:
  Type(Type_inforamses) :: inforamses


End Module modvariables
