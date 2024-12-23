!> @file 
!! This file contains the definition of the 3 arrays containing the position, the velocity and the id of the particles.

Module modvariable

  Use modconst
  Real(   kind=SP ), dimension(:,:), allocatable, target :: x   !< position of the particles
  Real(   kind=SP ), dimension(:,:), allocatable, target :: v   !< velocity of the particles
  Integer(kind=PRI), dimension(:),   allocatable, target :: id  !< RAMSES id of the particles; 
                                                                !! Integer8 if particle number exceeds 2^31-1
End Module modvariable
