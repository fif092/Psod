!==============================================================================
! Project: pFoF
! File: common/src/modmpicommons.f90
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
!! This file contains MPI related variables, type declarations and initialization procedures useful for different codes.

!> Common MPI related variables, type declarations and initialization procedures.
!>
!> Authors: F. Roy

Module modmpicommons

  Use modconstant
  Use mpi

  Implicit none

  Integer(kind=4) :: procID !< process id in the global communicator
  Integer(kind=4) :: procNB !< process number in the global communicator
!  Integer(kind=4) :: neighbours(6)  !< array containing the id of the neighbours of the local process in the MpiCube communicator

  Integer(kind=4) :: Mpi_Type_info_ramses !< MPI type corresponding to Type_info_ramses
  Integer(kind=4) :: Mpi_Type_info_cone_part !< MPI type corresponding to Type_info_cone_part
  Integer(kind=4) :: Mpi_Type_info_cone_grav !< MPI type corresponding to Type_info_cone_grav
  Integer(kind=4) :: Mpi_Type_parameter_pfof_snap !< MPI type corresponding to Type_parameter_pfof_snap
  Integer(kind=4) :: Mpi_Type_parameter_pfof_cone !< MPI type corresponding to Type_parameter_pfof_cone
  Integer(kind=4) :: Mpi_Type_parameter_cone_part !< MPI type corresponding to Type_parameter_conecreator_part
  Integer(kind=4) :: Mpi_Type_parameter_cone_grav !< MPI type corresponding to Type_parameter_conecreator_grav

  Type :: Type_info_communicator
     Integer(kind=4) :: name
     Integer(kind=4) :: pid
     Integer(kind=4) :: size
     Integer(kind=4) :: color
     Integer(kind=4), dimension(3) :: dims
     Logical(kind=4), dimension(3) :: periods
     Integer(kind=4), dimension(3) :: coords
     Integer(kind=4), dimension(6) :: neighbours
  End type Type_info_communicator

  Type :: Type_info_process
     Type(Type_info_communicator) :: global_comm
     Type(Type_info_communicator) :: read_comm
     Type(Type_info_communicator) :: write_comm
  End type Type_info_process
  
  Private

  Public :: procID, &
       procNB, &
 !      neighbours, &
       Type_info_communicator, &
       Type_info_process, &
       Mpi_Type_info_ramses, &
       Mpi_Type_info_cone_part, &
       Mpi_Type_info_cone_grav, &
       Mpi_Type_parameter_pfof_snap, &
       Mpi_Type_parameter_pfof_cone, &
       Mpi_Type_parameter_cone_part, &
       Mpi_Type_parameter_cone_grav, &
       create_mpi_type_info_cone_part, &
       create_mpi_type_info_cone_grav, &
       create_mpi_type_info_ramses, &
       create_mpi_type_param_pfof_snap, &
       create_mpi_type_param_pfof_cone, &
       create_mpi_type_param_cone_part, &
       create_mpi_type_param_cone_grav, &
       EmergencyStop

Contains

  !=======================================================================
  !< Create MPI type corresponding to Type_info_cone_part
  Subroutine create_mpi_type_info_cone_part

    Implicit none

    Type(Type_info_cone_part) :: fooconepart
    Integer(kind=4), allocatable, dimension(:) :: blocklen, type 
    Integer(kind=MPI_ADDRESS_KIND), allocatable, dimension(:) :: disp
    Integer(kind=MPI_ADDRESS_KIND) :: base
    Integer(kind=4) :: mpierr

    Allocate(disp(39), type(39), blocklen(39))

    Call Mpi_Get_Address(fooconepart%ncpu, disp(1), mpierr)
    Call Mpi_Get_Address(fooconepart%nstride, disp(2), mpierr)
    Call Mpi_Get_Address(fooconepart%nstep_coarse, disp(3), mpierr)
    Call Mpi_Get_Address(fooconepart%cone_id, disp(4), mpierr)
    Call Mpi_Get_Address(fooconepart%nglobalfile, disp(5), mpierr)
    Call Mpi_Get_Address(fooconepart%isfullsky, disp(6), mpierr)
    Call Mpi_Get_Address(fooconepart%future, disp(7), mpierr)
    Call Mpi_Get_Address(fooconepart%npart, disp(8), mpierr)
    Call Mpi_Get_Address(fooconepart%aexp, disp(9), mpierr)
    Call Mpi_Get_Address(fooconepart%observer_x, disp(10), mpierr)
    Call Mpi_Get_Address(fooconepart%observer_y, disp(11), mpierr)
    Call Mpi_Get_Address(fooconepart%observer_z, disp(12), mpierr)
    Call Mpi_Get_Address(fooconepart%observer_rds, disp(13), mpierr)
    Call Mpi_Get_Address(fooconepart%cone_zlim, disp(14), mpierr)
    Call Mpi_Get_Address(fooconepart%amax, disp(15), mpierr)
    Call Mpi_Get_Address(fooconepart%amin, disp(16), mpierr)
    Call Mpi_Get_Address(fooconepart%zmax, disp(17), mpierr)
    Call Mpi_Get_Address(fooconepart%zmin, disp(18), mpierr)
    Call Mpi_Get_Address(fooconepart%dmax, disp(19), mpierr)
    Call Mpi_Get_Address(fooconepart%dmin, disp(20), mpierr)
    Call Mpi_Get_Address(fooconepart%dtol, disp(21), mpierr)
    Call Mpi_Get_Address(fooconepart%thetay, disp(22), mpierr)
    Call Mpi_Get_Address(fooconepart%thetaz, disp(23), mpierr)
    Call Mpi_Get_Address(fooconepart%theta, disp(24), mpierr)
    Call Mpi_Get_Address(fooconepart%phi, disp(25), mpierr)
    Call Mpi_Get_Address(fooconepart%aendconem2, disp(26), mpierr)
    Call Mpi_Get_Address(fooconepart%aendconem1, disp(27), mpierr)
    Call Mpi_Get_Address(fooconepart%aendcone, disp(28), mpierr)
    Call Mpi_Get_Address(fooconepart%aexpold, disp(29), mpierr)
    Call Mpi_Get_Address(fooconepart%zendconem2, disp(30), mpierr)
    Call Mpi_Get_Address(fooconepart%zendconem1, disp(31), mpierr)
    Call Mpi_Get_Address(fooconepart%zendcone, disp(32), mpierr)
    Call Mpi_Get_Address(fooconepart%zexpold, disp(33), mpierr)
    Call Mpi_Get_Address(fooconepart%zexp, disp(34), mpierr)
    Call Mpi_Get_Address(fooconepart%dendconem2, disp(35), mpierr)
    Call Mpi_Get_Address(fooconepart%dendconem1, disp(36), mpierr)
    Call Mpi_Get_Address(fooconepart%dendcone, disp(37), mpierr)
    Call Mpi_Get_Address(fooconepart%dexpold, disp(38), mpierr)
    Call Mpi_Get_Address(fooconepart%dexp, disp(39), mpierr)

    base = disp(1) 
    disp(:) = disp(:) - base 
    blocklen(:) = 1 
    type(1:7) = MPI_INTEGER 
    type(8:8) = MPI_INTEGER8
    type(9:39) = MPI_DOUBLE_PRECISION 
 
    Call Mpi_Type_Create_Struct(39, blocklen, disp, type, Mpi_Type_info_cone_part, mpierr) 
    Call Mpi_Type_Commit(Mpi_Type_info_cone_part, mpierr) 

    Deallocate(disp, type, blocklen)

  End Subroutine create_mpi_type_info_cone_part


  !=======================================================================
  !< Create MPI type corresponding to Type_info_cone_grav
  Subroutine create_mpi_type_info_cone_grav

    Implicit none

    Type(Type_info_cone_grav) :: fooconegrav
    Integer(kind=4), allocatable, dimension(:) :: blocklen, type 
    Integer(kind=MPI_ADDRESS_KIND), allocatable, dimension(:) :: disp
    Integer(kind=MPI_ADDRESS_KIND) :: base
    Integer(kind=4) :: mpierr

    Allocate(disp(37), type(37), blocklen(37))

    Call Mpi_Get_Address(fooconegrav%ncpu, disp(1), mpierr)
    Call Mpi_Get_Address(fooconegrav%nstride, disp(2), mpierr)
    Call Mpi_Get_Address(fooconegrav%nstep_coarse, disp(3), mpierr)
    Call Mpi_Get_Address(fooconegrav%cone_id, disp(4), mpierr)
    Call Mpi_Get_Address(fooconegrav%nglobalfile, disp(5), mpierr)
    Call Mpi_Get_Address(fooconegrav%isfullsky, disp(6), mpierr)
    Call Mpi_Get_Address(fooconegrav%future, disp(7), mpierr)
    Call Mpi_Get_Address(fooconegrav%nlevel, disp(8), mpierr)
    Call Mpi_Get_Address(fooconegrav%levelmin, disp(9), mpierr)
    Call Mpi_Get_Address(fooconegrav%levelmax, disp(10), mpierr)    
    Call Mpi_Get_Address(fooconegrav%nglobalcell, disp(11), mpierr)
    Call Mpi_Get_Address(fooconegrav%aexp, disp(12), mpierr)
    Call Mpi_Get_Address(fooconegrav%observer_x, disp(13), mpierr)
    Call Mpi_Get_Address(fooconegrav%observer_y, disp(14), mpierr)
    Call Mpi_Get_Address(fooconegrav%observer_z, disp(15), mpierr)
    Call Mpi_Get_Address(fooconegrav%observer_rds, disp(16), mpierr)
    Call Mpi_Get_Address(fooconegrav%cone_zlim, disp(17), mpierr)
    Call Mpi_Get_Address(fooconegrav%amax, disp(18), mpierr)
    Call Mpi_Get_Address(fooconegrav%amin, disp(19), mpierr)
    Call Mpi_Get_Address(fooconegrav%zmax, disp(20), mpierr)
    Call Mpi_Get_Address(fooconegrav%zmin, disp(21), mpierr)
    Call Mpi_Get_Address(fooconegrav%dmax, disp(22), mpierr)
    Call Mpi_Get_Address(fooconegrav%dmin, disp(23), mpierr)
    Call Mpi_Get_Address(fooconegrav%dtol, disp(24), mpierr)
    Call Mpi_Get_Address(fooconegrav%thetay, disp(25), mpierr)
    Call Mpi_Get_Address(fooconegrav%thetaz, disp(26), mpierr)
    Call Mpi_Get_Address(fooconegrav%theta, disp(27), mpierr)
    Call Mpi_Get_Address(fooconegrav%phi, disp(28), mpierr)
    Call Mpi_Get_Address(fooconegrav%aendconem2, disp(29), mpierr)
    Call Mpi_Get_Address(fooconegrav%aendconem1, disp(30), mpierr)
    Call Mpi_Get_Address(fooconegrav%aendcone, disp(31), mpierr)
    Call Mpi_Get_Address(fooconegrav%zendconem2, disp(32), mpierr)
    Call Mpi_Get_Address(fooconegrav%zendconem1, disp(33), mpierr)
    Call Mpi_Get_Address(fooconegrav%zendcone, disp(34), mpierr)
    Call Mpi_Get_Address(fooconegrav%dendconem2, disp(35), mpierr)
    Call Mpi_Get_Address(fooconegrav%dendconem1, disp(36), mpierr)
    Call Mpi_Get_Address(fooconegrav%dendcone, disp(37), mpierr)

    base = disp(1) 
    disp(:) = disp(:) - base 
    blocklen(:) = 1 
    type(1:10) = MPI_INTEGER 
    type(11:11) = MPI_INTEGER8
    type(12:37) = MPI_DOUBLE_PRECISION 
 
    Call Mpi_Type_Create_Struct(37, blocklen, disp, type, Mpi_Type_info_cone_grav, mpierr) 
    Call Mpi_Type_Commit(Mpi_Type_info_cone_grav, mpierr) 

    Deallocate(disp, type, blocklen)

  End Subroutine create_mpi_type_info_cone_grav


  !=======================================================================
  !< Create MPI type corresponding to Type_info_ramses
  Subroutine create_mpi_type_info_ramses

    Implicit none

    Type(Type_info_ramses) :: fooramses
    Integer(kind=4), allocatable, dimension(:) :: blocklen, type 
    Integer(kind=MPI_ADDRESS_KIND), allocatable, dimension(:) :: disp
    Integer(kind=MPI_ADDRESS_KIND) :: base
    Integer(kind=4) :: mpierr


    Allocate(disp(17), type(17), blocklen(17))

    Call Mpi_Get_Address(fooramses%ncpu, disp(1), mpierr) 
    Call Mpi_Get_Address(fooramses%ndim, disp(2), mpierr) 
    Call Mpi_Get_Address(fooramses%levelmin, disp(3), mpierr) 
    Call Mpi_Get_Address(fooramses%levelmax, disp(4), mpierr) 
    Call Mpi_Get_Address(fooramses%ngridmax, disp(5), mpierr) 
    Call Mpi_Get_Address(fooramses%nstep_coarse, disp(6), mpierr) 
    Call Mpi_Get_Address(fooramses%boxlen, disp(7), mpierr) 
    Call Mpi_Get_Address(fooramses%time, disp(8), mpierr) 
    Call Mpi_Get_Address(fooramses%aexp, disp(9), mpierr) 
    Call Mpi_Get_Address(fooramses%h0, disp(10), mpierr) 
    Call Mpi_Get_Address(fooramses%omega_m, disp(11), mpierr) 
    Call Mpi_Get_Address(fooramses%omega_l, disp(12), mpierr) 
    Call Mpi_Get_Address(fooramses%omega_k, disp(13), mpierr) 
    Call Mpi_Get_Address(fooramses%omega_b, disp(14), mpierr) 
    Call Mpi_Get_Address(fooramses%unit_l, disp(15), mpierr) 
    Call Mpi_Get_Address(fooramses%unit_d, disp(16), mpierr) 
    Call Mpi_Get_Address(fooramses%unit_t, disp(17), mpierr) 

    base = disp(1) 
    disp(:) = disp(:) - base 
    blocklen(:) = 1 
    type(1:6) = MPI_INTEGER 
    type(7:17) = MPI_DOUBLE_PRECISION 
 
    Call Mpi_Type_Create_Struct(17, blocklen, disp, type, Mpi_Type_info_ramses, mpierr) 
    Call Mpi_Type_Commit(Mpi_Type_info_ramses, mpierr) 

    Deallocate(disp, type, blocklen)

  End Subroutine create_mpi_type_info_ramses


  !=======================================================================
  !< Create MPI type corresponding to Type_parameter_pfof_snap
  Subroutine create_mpi_type_param_pfof_snap

    Implicit none

    Type(Type_parameter_pfof_snap) :: foopfofsnap
    Integer(kind=4), allocatable, dimension(:) :: blocklen, type 
    Integer(kind=MPI_ADDRESS_KIND), allocatable, dimension(:) :: disp
    Integer(kind=MPI_ADDRESS_KIND) :: base
    Integer(kind=4) :: mpierr

    Allocate(disp(23), type(23), blocklen(23))

    Call Mpi_Get_Address(foopfofsnap%input_path, disp(1), mpierr)
    Call Mpi_Get_Address(foopfofsnap%part_input_file, disp(2), mpierr)
    Call Mpi_Get_Address(foopfofsnap%info_input_file, disp(3), mpierr)
    Call Mpi_Get_Address(foopfofsnap%simulation_name, disp(4), mpierr)
    Call Mpi_Get_Address(foopfofsnap%mmin, disp(5), mpierr)
    Call Mpi_Get_Address(foopfofsnap%mmax, disp(6), mpierr)
    Call Mpi_Get_Address(foopfofsnap%do_read_potential, disp(7), mpierr)
    Call Mpi_Get_Address(foopfofsnap%do_read_gravitational_field, disp(8), mpierr)
    Call Mpi_Get_Address(foopfofsnap%do_unbinding, disp(9), mpierr)
    Call Mpi_Get_Address(foopfofsnap%do_subHalo, disp(10), mpierr)
    Call Mpi_Get_Address(foopfofsnap%do_timings, disp(11), mpierr)
    Call Mpi_Get_Address(foopfofsnap%percolation_length, disp(12), mpierr)

    Call Mpi_Get_Address(foopfofsnap%grpsize, disp(13), mpierr)
    Call Mpi_Get_Address(foopfofsnap%gatherread_factor, disp(14), mpierr)
    Call Mpi_Get_Address(foopfofsnap%snapshot, disp(15), mpierr)
    Call Mpi_Get_Address(foopfofsnap%gatherwrite_factor, disp(16), mpierr)
    Call Mpi_Get_Address(foopfofsnap%do_skip_star, disp(17), mpierr)
    Call Mpi_Get_Address(foopfofsnap%do_skip_metal, disp(18), mpierr)
    Call Mpi_Get_Address(foopfofsnap%do_read_from_cube, disp(19), mpierr)
    Call Mpi_Get_Address(foopfofsnap%do_fof, disp(20), mpierr)
    Call Mpi_Get_Address(foopfofsnap%do_write_cube, disp(21), mpierr)
    Call Mpi_Get_Address(foopfofsnap%do_sort_cube, disp(22), mpierr)
    Call Mpi_Get_Address(foopfofsnap%code_index, disp(23), mpierr)

    base = disp(1) 
    disp(:) = disp(:) - base 
    blocklen(1:4) = 200 
    blocklen(5:22) = 1
    blocklen(23:23) = 3
    type(1:4) = MPI_CHARACTER 
    type(5:6) = MPI_INTEGER
    type(7:11) = MPI_LOGICAL
    type(12:12) = Mpi_Real
    type(13:16) = Mpi_Integer
    type(17:22) = Mpi_Logical
    type(23) = MPI_CHARACTER

    Call Mpi_Type_Create_Struct(23, blocklen, disp, type, Mpi_Type_parameter_pfof_snap, mpierr) 
    Call Mpi_Type_Commit(Mpi_Type_parameter_pfof_snap, mpierr) 

    Deallocate(disp, type, blocklen)

  End Subroutine create_mpi_type_param_pfof_snap


  !=======================================================================
  !< Create MPI type corresponding to Type_parameter_pfof_cone
  Subroutine create_mpi_type_param_pfof_cone

    Implicit none

    Type(Type_parameter_pfof_cone) :: foopfofcone
    Integer(kind=4), allocatable, dimension(:) :: blocklen, type 
    Integer(kind=MPI_ADDRESS_KIND), allocatable, dimension(:) :: disp
    Integer(kind=MPI_ADDRESS_KIND) :: base
    Integer(kind=4) :: mpierr

  
    Allocate(disp(15), type(15), blocklen(15))

    Call Mpi_Get_Address(foopfofcone%input_path, disp(1), mpierr)
    Call Mpi_Get_Address(foopfofcone%part_input_file, disp(2), mpierr)
    Call Mpi_Get_Address(foopfofcone%simulation_name, disp(3), mpierr)

    Call Mpi_Get_Address(foopfofcone%mmin, disp(4), mpierr)
    Call Mpi_Get_Address(foopfofcone%mmax, disp(5), mpierr)

    Call Mpi_Get_Address(foopfofcone%do_read_ramses_part_id, disp(6), mpierr)
    Call Mpi_Get_Address(foopfofcone%do_read_potential, disp(7), mpierr)
    Call Mpi_Get_Address(foopfofcone%do_read_gravitational_field, disp(8), mpierr)
    Call Mpi_Get_Address(foopfofcone%do_unbinding, disp(9), mpierr)
    Call Mpi_Get_Address(foopfofcone%do_subHalo, disp(10), mpierr)
    Call Mpi_Get_Address(foopfofcone%do_timings, disp(11), mpierr)
    Call Mpi_Get_Address(foopfofcone%do_gather_halo, disp(12), mpierr) 

    Call Mpi_Get_Address(foopfofcone%percolation_length, disp(13), mpierr)

    Call Mpi_Get_Address(foopfofcone%shell_first_id, disp(14), mpierr)
    Call Mpi_Get_Address(foopfofcone%shell_last_id, disp(15), mpierr)

    base = disp(1) 
    disp(:) = disp(:) - base 
    blocklen(1:3) = 200 
    blocklen(4:15) = 1
    type(1:3) = MPI_CHARACTER 
    type(4:5) = MPI_INTEGER
    type(6:12) = MPI_LOGICAL
    type(13:13) = Mpi_Real
    type(14:15) = Mpi_Integer

    Call Mpi_Type_Create_Struct(15, blocklen, disp, type, Mpi_Type_parameter_pfof_cone, mpierr) 
    Call Mpi_Type_Commit(Mpi_Type_parameter_pfof_cone, mpierr) 

    Deallocate(disp, type, blocklen)

  End Subroutine create_mpi_type_param_pfof_cone


  !=======================================================================
  !< Create MPI type corresponding to Type_parameter_conecreator_part
  Subroutine create_mpi_type_param_cone_part

    Implicit none

    Type(Type_parameter_conecreator_part) :: foo
    Integer(kind=4), allocatable, dimension(:) :: blocklen, type 
    Integer(kind=MPI_ADDRESS_KIND), allocatable, dimension(:) :: disp
    Integer(kind=MPI_ADDRESS_KIND) :: base
    Integer(kind=4) :: mpierr

    Allocate(disp(13), type(13), blocklen(13))

    Call Mpi_Get_Address(foo%input_path, disp(1), mpierr)
    Call Mpi_Get_Address(foo%simulation_name, disp(2), mpierr)
    Call Mpi_Get_Address(foo%cone_input_file, disp(3), mpierr)
    Call Mpi_Get_Address(foo%info_cone_input_file, disp(4), mpierr)
    Call Mpi_Get_Address(foo%info_ramses_input_file, disp(5), mpierr)
    Call Mpi_Get_Address(foo%nfile, disp(6), mpierr)
    Call Mpi_Get_Address(foo%first_file, disp(7), mpierr)
    Call Mpi_Get_Address(foo%last_file, disp(8), mpierr)
    Call Mpi_Get_Address(foo%cone_max_radius, disp(9), mpierr)
    Call Mpi_Get_Address(foo%cube_size, disp(10), mpierr)
    Call Mpi_Get_Address(foo%do_read_ramses_part_id, disp(11), mpierr)
    Call Mpi_Get_Address(foo%do_read_potential, disp(12), mpierr)
    Call Mpi_Get_Address(foo%do_read_gravitational_field, disp(13), mpierr)

    base = disp(1) 
    disp(:) = disp(:) - base 
    blocklen(1:5) = 200
    blocklen(6:13) = 1
    type(1:5) = MPI_CHARACTER 
    type(6:8) = MPI_INTEGER
    type(9:10) = Mpi_Double_Precision
    type(11:13) = Mpi_Logical

    Call Mpi_Type_Create_Struct(13, blocklen, disp, type, Mpi_Type_parameter_cone_part, mpierr) 
    Call Mpi_Type_Commit(Mpi_Type_parameter_cone_part, mpierr) 
    
    Deallocate(disp, type, blocklen)

  End Subroutine create_mpi_type_param_cone_part


  !=======================================================================
  !< Create MPI type corresponding to Type_parameter_conecreator_grav
  Subroutine create_mpi_type_param_cone_grav

    Implicit none

    Type(Type_parameter_conecreator_grav) :: foo
    Integer(kind=4), allocatable, dimension(:) :: blocklen, type 
    Integer(kind=MPI_ADDRESS_KIND), allocatable, dimension(:) :: disp
    Integer(kind=MPI_ADDRESS_KIND) :: base
    Integer(kind=4) :: mpierr

    Allocate(disp(12), type(12), blocklen(12))

    Call Mpi_Get_Address(foo%input_path, disp(1), mpierr)
    Call Mpi_Get_Address(foo%simulation_name, disp(2), mpierr)
    Call Mpi_Get_Address(foo%cone_input_file, disp(3), mpierr)
    Call Mpi_Get_Address(foo%info_cone_input_file, disp(4), mpierr)
    Call Mpi_Get_Address(foo%info_ramses_input_file, disp(5), mpierr)
    Call Mpi_Get_Address(foo%nfile, disp(6), mpierr)
    Call Mpi_Get_Address(foo%first_file, disp(7), mpierr)
    Call Mpi_Get_Address(foo%last_file, disp(8), mpierr)
    Call Mpi_Get_Address(foo%nlevel, disp(9), mpierr)
    Call Mpi_Get_Address(foo%levelmin, disp(10), mpierr)
    Call Mpi_Get_Address(foo%cone_max_radius, disp(11), mpierr)
    Call Mpi_Get_Address(foo%cube_size, disp(12), mpierr)

    base = disp(1) 
    disp(:) = disp(:) - base 
    blocklen(1:5) = 200
    blocklen(6:12) = 1
    type(1:5) = MPI_CHARACTER 
    type(6:10) = MPI_INTEGER
    type(11:12) = Mpi_Double_Precision

    Call Mpi_Type_Create_Struct(12, blocklen, disp, type, Mpi_Type_parameter_cone_grav, mpierr) 
    Call Mpi_Type_Commit(Mpi_Type_parameter_cone_grav, mpierr) 
    
    Deallocate(disp, type, blocklen)

  End Subroutine create_mpi_type_param_cone_grav


  !=======================================================================
  !> The subroutine is used to abort the execution if something wrong is happening.
  Subroutine EmergencyStop(message,errcode)

    Implicit none

    Character(len=*), Intent(in) :: message
    Integer, Intent(in) :: errcode

    Integer(kind=4) :: mpierr

    Write(*,1000) trim(message),procID

1000 Format('*** ',A,' on process ',I5.5,'. ***')

    If(procID==0) Close(50)
    Close(54)

    Call Mpi_Abort(Mpi_Comm_World,errcode,mpierr)
    Stop

  End Subroutine EmergencyStop
  
End Module modmpicommons
