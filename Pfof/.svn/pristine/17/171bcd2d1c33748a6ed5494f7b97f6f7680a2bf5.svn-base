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
