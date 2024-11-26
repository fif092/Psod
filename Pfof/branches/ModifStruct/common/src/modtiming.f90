!==============================================================================
! Project: pFoF
! File: common/src/modtiming.f90
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
!! This file contains variables used for timings.

!> This module contains variables used for timings.
!>
!> Authors: F. Roy, V. Bouillot

Module modtiming

  Implicit None
  Real(kind=8) :: time0 = 0.d0                  !< beginning of the execution
  Real(kind=8) :: timeInt = 0.d0                !< intermediate time
  Real(kind=8) :: tReadFile = 0.d0              !< time used to read particles files (cubes or Ramses files)
  Real(kind=8) :: tTailPart = 0.d0              !< time used to distribute the particles if read from Ramses files
  Real(kind=8) :: tInitRead = 0.d0              !< time used to initialize reading from Ramses files
  Real(kind=8) :: tRead = 0.d0                  !< global reading time
  Real(kind=8) :: tObs = 0.d0                   !< time used to compute observables for haloes
  Real(kind=8) :: tOut = 0.d0                   !< time used to write output files
  Real(kind=8) :: tFoF = 0.d0                   !< global time for friends-of-friends halo detection
  Real(kind=8) :: tFoFloc = 0.d0                !< time for local fof algorithm
  Real(kind=8) :: tFoFinit = 0.d0               !< time for fof initialization (construction of the tree)
  Real(kind=8) :: tRaccord = 0.d0               !< time for halo merging procedure 
  Real(kind=8) :: tGatherhalo = 0.d0            !< time for gathering particles following their haloID
  Real(kind=8) :: tOuthalopart = 0.d0           !< time for writing halo part files
  Real(kind=8) :: tOutmass = 0.d0               !< time for writing halo mass files
  Real(kind=8) :: tSelectHalo = 0.d0            !< time for selecting halo with mass > Mmin
  Real(kind=8) :: tSort = 0.d0                  !< time for sorting particles following their haloID

End Module modtiming
