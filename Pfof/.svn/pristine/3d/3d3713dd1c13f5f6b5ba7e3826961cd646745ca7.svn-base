!==============================================================================
! Project: pFoF
! File: pfof_snap/src/modvariables.f90
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
!! This file contains some global variables declarations.

!> This module contains some global variables declarations.
Module modvariables

  Use modconstant
  Implicit None

  Integer(kind=PRI) :: ngrid       !< total number of grid points: ngrid = nres ** 3
  Real   (kind=4)   :: xmin, xmax, ymin, ymax, zmin, zmax  !< min and max (x,y,z) for each process
  Type(Type_parameter_pfof_snap) :: param !< Input parameters for pfof_snap

  ! information read from Ramses info file:
  Type(Type_info_ramses) :: inforamses !< Information read from RAMSES info file.

  Type(Type_common_metadata) :: common_meta !< Metadata common to every HDF5 file created by pfof_snap.

  Private

  Public :: ngrid, &
       xmin, xmax, ymin, ymax, zmin, zmax, &
       param, &
       inforamses,&
       common_meta

End Module modvariables
