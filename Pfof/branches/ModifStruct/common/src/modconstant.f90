!==============================================================================
! Project: pFoF
! File: common/src/modconstant.f90
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
!! This file contains some constant parameters definition

!> Constant parameters definition
!>
!> @author  F. Roy, V. Bouillot

Module modconstant

  Use iso_c_binding
  Use mpi, only : Mpi_Integer, Mpi_Integer8

#include "svnrev.h"

#ifdef LONGREAL
  Integer, parameter :: PR=8   !< Precision for real arrays read from Ramses simulations (position/velocities)
#else
  Integer, parameter :: PR=4   !< Precision for real arrays read from Ramses simulations (position/velocities)
#endif

#ifdef LONGINT
  Integer, parameter :: PRI = 8   !< Precision for integer arrays (id)
  Integer, parameter :: MPI_PRI = Mpi_Integer8  !< MPI precision for integer arrays
#else
  Integer, parameter :: PRI = 4   !< Precision for integer arrays (id)
  Integer, parameter :: MPI_PRI = Mpi_Integer   !< MPI precision for integer arrays
#endif

  ! Name of the codes to write as metadata in HDF5 files
  Character(len=16), parameter :: NAME_CONECREATOR_PART='conecreator_part'
  Character(len=16), parameter :: NAME_CONECREATOR_GRAV='conecreator_grav'
  Character(len=16), parameter :: NAME_CONEMAPPER='conemapper'
  Character(len=16), parameter :: NAME_PFOF_SNAP='pfof_snap'
  Character(len=16), parameter :: NAME_PFOF_CONE='pfof_cone'
  Character(len=16), parameter :: NAME_PSOD_SNAP='psod_snap'
  

  ! Output Units 
  Integer, parameter :: Ulog=50   !< I/O unit for text log file
  Integer, parameter :: Ucub=51   !< I/O unit for binary cube file
  Integer, parameter :: Umas=52   !< I/O unit for binary mass file
  Integer, parameter :: Ustr=53   !< I/O unit for binary haloes file
  Integer, parameter :: Uopa=54   !< I/O unit for text input parameter log file

  ! SVN Revision 
  Character(len=32) :: svn_version=SVNREV


  ! Type declaration for information read from cone info file
  Type :: Type_info_cone !< Type for information read from cone info file
     Integer(kind=4) :: ncpu
     Integer(kind=4) :: nstride
     Integer(kind=4) :: nstep_coarse
     Integer(kind=4) :: cone_id
     Integer(kind=4) :: nglobalfile
     Integer(kind=4) :: isfullsky
     Integer(kind=4) :: future
     Real(kind=8) :: aexp
     Real(kind=8) :: observer_x
     Real(kind=8) :: observer_y
     Real(kind=8) :: observer_z
     Real(kind=8) :: observer_rds
     Real(kind=8) :: cone_zlim
     Real(kind=8) :: amax
     Real(kind=8) :: amin
     Real(kind=8) :: zmax
     Real(kind=8) :: zmin
     Real(kind=8) :: dmax
     Real(kind=8) :: dmin
     Real(kind=8) :: dtol
     Real(kind=8) :: thetay
     Real(kind=8) :: thetaz
     Real(kind=8) :: theta
     Real(kind=8) :: phi
     Real(kind=8) :: aendconem2
     Real(kind=8) :: aendconem1
     Real(kind=8) :: aendcone
     Real(kind=8) :: zendconem2
     Real(kind=8) :: zendconem1
     Real(kind=8) :: zendcone
     Real(kind=8) :: dendconem2
     Real(kind=8) :: dendconem1
     Real(kind=8) :: dendcone
  End Type Type_info_cone

  Type, extends(Type_info_cone) :: Type_info_cone_part !< Extension of Type_info_cone for particles cone
     Integer(kind=8) :: npart          !! nglobalcell in the file
     Real(kind=8) :: aexpold
     Real(kind=8) :: zexpold
     Real(kind=8) :: zexp
     Real(kind=8) :: dexpold
     Real(kind=8) :: dexp
  End Type Type_info_cone_part

  Type, extends(Type_info_cone) :: Type_info_cone_grav !< Extension of Type_info_cone for cells cone
     Integer(kind=4) :: nlevel
     Integer(kind=4) :: levelmin
     Integer(kind=4) :: levelmax
     Integer(kind=8) :: nglobalcell
  End type Type_info_cone_grav

  ! Type declaration for information read from Ramses info file
  Type :: Type_info_ramses !< Type for information read from Ramses info file
     Integer(kind=4) :: ncpu
     Integer(kind=4) :: ndim          !< number of dimensions
     Integer(kind=4) :: levelmin          !< minimum mesh refinement level
     Integer(kind=4) :: levelmax
     Integer(kind=4) :: ngridmax
     Integer(kind=4) :: nstep_coarse
     Real(kind=8) :: boxlen
     Real(kind=8) :: time
     Real(kind=8) :: aexp
     Real(kind=8) :: h0
     Real(kind=8) :: omega_m
     Real(kind=8) :: omega_l
     Real(kind=8) :: omega_k
     Real(kind=8) :: omega_b
     Real(kind=8) :: unit_l
     Real(kind=8) :: unit_d
     Real(kind=8) :: unit_t
  End Type Type_info_ramses

  Type :: Type_parameter_halofinder
     Character(len=200) :: input_path      !< path to the directory containing the RAMSES files
     Character(len=200) :: part_input_file  !< name of the files containing the particles data
     Character(len=200) :: simulation_name  !< Simulation name to be written in output files
     Integer(kind=4)    :: mmin            !< minimum mass of the haloes
     Integer(kind=4)    :: mmax            !< maximum mass of the haloes
     Logical(kind=4)    :: do_read_potential    !< do the RAMSES files contain potential?
     Logical(kind=4)    :: do_read_gravitational_field  !< do the RAMSES files contain force?
     Logical(kind=4)    :: do_unbinding    !< Unbinding process? (not implemented yet)
     Logical(kind=4)    :: do_subHalo      !< SubHalo detection? (not implemented yet)
     Logical(kind=4)    :: do_timings  !< should pFOF perform timings (imply extra synchronizations)
  End type Type_parameter_halofinder

  Type, extends(Type_parameter_halofinder) :: Type_parameter_psod_snap
     Logical(kind=4) :: do_sod
  End Type Type_parameter_psod_snap
  
  Type, extends(Type_parameter_halofinder) :: Type_parameter_pfof !< Type for pfof input parameters common to snapshot and cone versions
     Real(kind=4)       :: percolation_length  !< FOF percolation length
  End Type Type_parameter_pfof

  Type, extends(Type_parameter_pfof)  :: Type_parameter_pfof_snap
     Integer(kind=4)    :: grpsize     !< size of the I/O group used for the RAMSES simulation
     Integer(kind=4)    :: gatherread_factor    !< gather parameter for parallel hdf5 input
     Integer(kind=4)    :: snapshot             !< Ramses output number to be written in output files
     Integer(kind=4)    :: gatherwrite_factor   !< gather parameter for parallel hdf5 output
     Logical(kind=4)    :: do_skip_star         !< do the RAMSES files contain stars?
     Logical(kind=4)    :: do_skip_metal        !< do the RAMSES files contain metalicities?
     Logical(kind=4)    :: do_read_from_cube    !< should pFOF read particles from cube files?
     Logical(kind=4)    :: do_fof   !< should pFOF perform FOF halo finding?
     Logical(kind=4)    :: do_write_cube        !< should pFOF write cube files?
     Logical(kind=4)    :: do_sort_cube         !< sort the particles and write a 'sorted cube'
     Character(len=3)   :: code_index  !< version of RAMSES used to produce the RAMSES files
     Character(len=200) :: info_input_file  !< name of the RAMSES info file
  End Type Type_parameter_pfof_snap

  Type, extends(Type_parameter_pfof) :: Type_parameter_pfof_cone !< Extension of Type_parameter_pfof for the cone version
     Logical(kind=4) :: do_read_ramses_part_id
     Integer(kind=4) :: shell_first_id
     Integer(kind=4) :: shell_last_id
     Logical(kind=4) :: do_gather_halo
  End Type Type_parameter_pfof_cone

  Type :: Type_parameter_conecreator !< Type for conecreator input parameters common to particles and cells versions
     Character(len=200) :: input_path
     Character(len=200) :: simulation_name
     Character(len=200) :: cone_input_file
     Character(len=200) :: info_cone_input_file
     Character(len=200) :: info_ramses_input_file
     Integer(kind=4) :: nfile
     Integer(kind=4) :: first_file
     Integer(kind=4) :: last_file
     Real(kind=8) :: cone_max_radius
     Real(kind=8) :: cube_size    
  End Type Type_parameter_conecreator
  
  Type, extends(Type_parameter_conecreator) :: Type_parameter_conecreator_part !< Extension of Type_parameter_conecreator for the particles version
     Logical(kind=4) :: do_read_ramses_part_id
     Logical(kind=4) :: do_read_potential
     Logical(kind=4) :: do_read_gravitational_field
  End type Type_parameter_conecreator_part

  Type, extends(Type_parameter_conecreator) :: Type_parameter_conecreator_grav !< Extension of Type_parameter_conecreator for the cells version
     Integer(kind=4) :: nlevel
     Integer(kind=4) :: levelmin
  End type Type_parameter_conecreator_grav


  Type :: Type_common_metadata !< Type for metadata common to every kind of hdf5 output file
     Character(len=16) :: created_by
     Integer(kind=4) :: svn_version
     Character(len=16) :: simulation_code
     Character(len=16) :: particle_type
     Integer(kind=4) :: constant_mass
     Character(len=16) :: units
     Integer(kind=8) :: npart_simulation
  End type Type_common_metadata


  Private
  Public :: Type_parameter_conecreator, &
       Type_parameter_conecreator_part, &
       Type_parameter_conecreator_grav, &
       Type_parameter_halofinder, &
       Type_parameter_pfof, &
       Type_parameter_pfof_cone, &
       Type_parameter_pfof_snap, &
       Type_parameter_psod_snap, &
       Type_info_ramses, &
       Type_info_cone, &
       Type_info_cone_grav, &
       Type_info_cone_part, &
       Type_common_metadata, &
       PR, &
       PRI, &
       MPI_PRI, &
       NAME_CONECREATOR_PART, &
       NAME_CONECREATOR_GRAV, &
       NAME_CONEMAPPER, &
       NAME_PFOF_SNAP, &
       NAME_PFOF_CONE, &
       NAME_PSOD_SNAP, &
       Ulog, &
       Ucub, &
       Umas, &
       Ustr, &
       Uopa, &
       svn_version

End Module modconstant