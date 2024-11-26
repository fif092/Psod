!> @file
!! This file contains the definition of the input parameters.

!> This module contains the definition of the input parameters.
!> Authors: F. Roy, V. Bouillot
Module modparameters

  Use modconstant
  Implicit None

!!$  ! Input parameters
!!$  Character(len=3)   :: code_index           !< version of RAMSES used to produce the RAMSES files to analyze
!!$  Character(len=200) :: input_path           !< path to the directory containing the RAMSES files
!!$  Character(len=200) :: part_inputfile       !< base name of the files containing the particles information
!!$  Character(len=200) :: info_inputfile       !< name of the RAMSES info file
!!$  Integer(kind=4)    :: grpsize              !< size of the I/O group used for the RAMSES simulation
!!$  Logical(kind=4)    :: do_skip_star         !< do the RAMSES files contain stars?
!!$  Logical(kind=4)    :: do_skip_metal        !< do the RAMSES files contain metalicities?
!!$  Logical(kind=4)    :: do_read_potential    !< do the RAMSES files contain potential?
!!$  Logical(kind=4)    :: do_read_force        !< do the RAMSES files contain force?
!!$  Logical(kind=4)    :: do_read_from_cube    !< should pFOF read particles information from previously created cube files?
!!$  Integer(kind=4)    :: gatherread_factor    !< gather parameter for parallel hdf5 input
!!$
!!$  ! Friend of friend parameters
!!$  Real(kind=4)       :: percolation_length          !< FOF percolation length
!!$  Integer(kind=4)    :: mmin            !< minimum mass of the haloes
!!$  Integer(kind=4)    :: mmax            !< maximum mass of the haloes
!!$  Logical(kind=4)    :: do_fof          !< should pFOF perform FOF halo finding? (if not, it can still write cube files)
!!$  Logical(kind=4)    :: do_unbinding    !< Unbinding process? (not implemented yet)
!!$  Logical(kind=4)    :: do_subHalo      !< SubHalo detection? (not implemented yet)
!!$
!!$  ! Output parameters
!!$  Character(len=200) :: simulation_name      !< Simulation name to be written in output files
!!$  Integer(kind=4)    :: snapshot             !< Ramses output number to be written in output files
!!$  Logical(kind=4)    :: do_write_cube        !< should pFOF write cube files?
!!$  Integer(kind=4)    :: gatherwrite_factor   !< gather parameter for parallel hdf5 output
!!$  Logical(kind=4)    :: do_sort_cube         !< sort the particles and write a 'sorted cube'
!!$  Logical(kind=4)    :: do_timings      !< should pFOF perform timings (imply extra synchronizations)

  Type(Type_parameter_psod_snap) :: param

!!$  ! Namelist for input file
!!$  Namelist / input_parameters  / code_index, input_path, part_inputfile, info_inputfile,  grpsize, &
!!$       do_skip_star, do_skip_metal, do_read_potential, do_read_force, do_read_from_cube, &
!!$       gatherread_factor
!!$  Namelist / fof_parameters    / percolation_length, mmin, mmax, do_fof, do_unbinding, do_subhalo
!!$  Namelist / output_parameters / simulation_name, snapshot, do_write_cube,&
!!$       gatherwrite_factor,  do_sort_cube, do_timings


End Module modparameters
