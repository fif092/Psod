!=======================================================================
! Author: Fabrice Roy (LUTH/CNRS/Observatoire de Paris)
! Fabrice.Roy@obspm.fr
!=======================================================================

Module modvariables

  Implicit None
  
  Integer(kind=4) :: shell_nb  ! nb of shell files
  Integer(kind=4), dimension(:), allocatable :: npartcube  ! nb of particles in each shell cube
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
