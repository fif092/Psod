!==============================================================================
! Project: pFoF
! File: common/src/modiocommons.f90
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
!! This file contains input/output procedures useful for different codes

!> Common input/output procedures
!>
!> Authors: F. Roy

Module modiocommons

  Use modvariables
  Use modhdf5

  Private
  Public :: h5write_meta_common, &
       h5write_meta_pfof_parameter, &
       h5write_meta_conecreator_parameter, &
       h5write_meta_info_ramses, &
       h5write_meta_info_cone, &
       h5write_info_haloanalyzer, &
       h5read_meta_common, & 
       h5read_meta_pfof_parameter, &
       h5read_meta_info_cone, &
       h5read_meta_info_ramses, &
       h5read_meta_conecreator_parameter!, &
!       h5read_meta_common_metadata


Contains

  !=======================================================================
  !> Write metadata common to every hdf5 ouput files of the pfof toolbox
  Subroutine h5write_meta_common(file_id, npart)

    Implicit none

    Integer(kind=hid_t), intent(in) :: file_id
    Integer(kind=8), intent(in) :: npart

    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: aname
    Character(len=H5STRLEN) :: adata
    Character(len=H5STRLEN) :: dsetname

    Integer(kind=hid_t) :: gr_id
    Integer(kind=4) :: tmpint4

    groupname = 'metadata'
    Call hdf5_open_group(file_id, groupname, gr_id)

    ! Common metadata
    aname = 'created_by'
    Call hdf5_write_attr(gr_id, aname, common_metadata%created_by)

    aname = 'svn_version'
    Call hdf5_write_attr(gr_id, aname, common_metadata%svn_version)

    aname = 'simulation_code'
    adata = 'ramses'
    Call hdf5_write_attr(gr_id, aname, adata)

    aname = 'particle_type'
    adata = 'dark_matter'
    Call hdf5_write_attr(gr_id, aname, adata)

    aname = 'constant_mass'
    tmpint4 = 1
    Call hdf5_write_attr(gr_id, aname, tmpint4)

    aname = 'units'
    adata = 'ramses'
    Call hdf5_write_attr(gr_id, aname, adata) 

    ! Write the number of particles in this simulation as an integer(kind=8) dataset
    dsetname = 'npart_simulation'
    Call hdf5_write_data(gr_id, dsetname, npart)

    Call hdf5_close_group(gr_id)

  End Subroutine h5write_meta_common

  !=======================================================================
  !> Write pfof input parameters as metadata 
  Subroutine h5write_meta_pfof_parameter(file_id, param_pfof)

    Implicit none

    Integer(kind=hid_t), intent(in) :: file_id
    Class(Type_parameter_pfof), intent(in) :: param_pfof

    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: aname

    Integer(kind=hid_t) :: gr_id
    Integer(kind=hid_t) :: gr_pfof_id
    Integer(kind=hid_t) :: gr_input_id
    Integer(kind=hid_t) :: gr_fof_id
    Integer(kind=hid_t) :: gr_output_id

    Integer(kind=4) :: tmpint4
    Character(len=16) :: c_b_name


    Select Type (param_pfof)
    Type is (Type_parameter_pfof_snap)
       c_b_name = NAME_PFOF_SNAP
    Type is (Type_parameter_pfof_cone)
       c_b_name = NAME_PFOF_CONE
    End Select

    groupname = 'metadata'
    Call hdf5_open_group(file_id,groupname, gr_id)
    ! pfof parameters:
    groupname = trim(c_b_name)//'_parameters'
    Call hdf5_create_group(gr_id, groupname, gr_pfof_id)
    groupname = 'input_parameters'
    Call hdf5_create_group(gr_pfof_id, groupname, gr_input_id)
    aname = 'input_path'
    Call hdf5_write_attr(gr_input_id, aname, param_pfof%input_path)
    aname = 'part_input_file'
    Call hdf5_write_attr(gr_input_id, aname, param_pfof%part_input_file)

    ! Write potential logical as an integer attribute (1=true, 0=false)
    aname = 'do_read_potential'
    If(param_pfof%do_read_potential) Then
       tmpint4=1
    Else
       tmpint4=0
    End If
    Call hdf5_write_attr(gr_input_id, aname, tmpint4)

    ! Write force logical as an integer attribute (1=true, 0=false)
    aname = 'do_read_gravitational_field'
    If(param_pfof%do_read_gravitational_field) Then
       tmpint4=1
    Else
       tmpint4=0
    End If
    Call hdf5_write_attr(gr_input_id, aname, tmpint4)

    Select Type (param_pfof)
    Type Is (Type_parameter_pfof_snap)
       aname = 'info_input_file'
       Call hdf5_write_attr(gr_input_id, aname, param_pfof%info_input_file)
       aname = 'grpsize'
       Call hdf5_write_attr(gr_input_id, aname, param_pfof%grpsize)
       aname = 'code_index'
       Call hdf5_write_attr(gr_input_id, aname, param_pfof%code_index)
       aname = 'do_skip_star'
       If(param_pfof%do_skip_star) Then
          tmpint4=1
       Else
          tmpint4=0
       End If
       Call hdf5_write_attr(gr_input_id, aname, tmpint4)
       aname = 'do_skip_metal'
       If(param_pfof%do_skip_metal) Then
          tmpint4=1
       Else
          tmpint4=0
       End If
       Call hdf5_write_attr(gr_input_id, aname, tmpint4)       
       aname='do_read_from_cube'
       tmpint4 = 0
       If(param_pfof%do_read_from_cube) tmpint4 = 1
       Call hdf5_write_attr(gr_input_id, aname, tmpint4)
       aname = 'gatherread_factor'
       Call hdf5_write_attr(gr_input_id, aname, param_pfof%gatherread_factor)

    Type Is (Type_parameter_pfof_cone)
       aname = 'shell_first_id'
       Call hdf5_write_attr(gr_input_id, aname, param_pfof%shell_first_id)
       aname = 'shell_last_id'
       Call hdf5_write_attr(gr_input_id, aname, param_pfof%shell_last_id)
    End Select
    Call hdf5_close_group(gr_input_id)
#ifdef SOD
    groupname = 'sod_parameters'
#else
    groupname = 'fof_parameters'
#endif
    Call hdf5_create_group(gr_pfof_id, groupname, gr_fof_id)
    ! Write percolation parameter b as attribute
    aname = 'percolation_length'
    Call hdf5_write_attr(gr_fof_id, aname, param_pfof%percolation_length)
    ! Write minimum halo mass Mmin as attribute
    aname = 'npart_halo_min'
    Call hdf5_write_attr(gr_fof_id, aname, param_pfof%mmin)    
    aname = 'npart_halo_max'
    Call hdf5_write_attr(gr_fof_id, aname, param_pfof%mmax)

    ! Write doUnbinding as attribute (not implemented yet) 
    aname = 'do_unbinding'
    If(param_pfof%do_unbinding) Then
       tmpint4=1
    Else
       tmpint4=0
    End If
    Call hdf5_write_attr(gr_fof_id, aname, tmpint4)

    ! Write doSubHalo as attribute (not implemented yet) 
    aname = 'do_subhalo'
    If(param_pfof%do_subhalo) Then
       tmpint4=1
    Else
       tmpint4=0
    End If
    Call hdf5_write_attr(gr_fof_id, aname, tmpint4)

    Select Type (param_pfof)
    Type Is (Type_parameter_pfof_snap)
#ifdef SOD
       aname = 'do_sod'
#else
       aname = 'do_fof'
#endif
       If(param_pfof%do_fof) Then
          tmpint4=1
       Else
          tmpint4=0
       End If
       Call hdf5_write_attr(gr_fof_id, aname, tmpint4)
    End Select
    Call hdf5_close_group(gr_fof_id)       

    groupname = 'output_parameters'
    Call hdf5_create_group(gr_pfof_id, groupname, gr_output_id)
    ! Write the simulation name as attribute
    aname = 'simulation_name'
    Call hdf5_write_attr(gr_output_id, aname, param_pfof%simulation_name)

    ! Write do_timings as attribute
    aname = 'do_timings'
    If(param_pfof%do_timings) Then
       tmpint4=1
    Else
       tmpint4=0
    End If
    Call hdf5_write_attr(gr_fof_id, aname, tmpint4)

    Select Type (param_pfof)
    Type Is (Type_parameter_pfof_snap)
       ! write the snapshot number: if halo computed from a lightcone then outputNB = -1
       aname = 'snapshot'
       Call hdf5_write_attr(gr_output_id, aname, param_pfof%snapshot)
       tmpint4 = 0
       If(param_pfof%do_write_cube) tmpint4 = 1
       aname = 'do_write_cube'
       Call hdf5_write_attr(gr_output_id, aname, tmpint4)       
       tmpint4 = 0
       If(param_pfof%do_sort_cube) tmpint4 = 1
       aname = 'do_sort_cube'
       Call hdf5_write_attr(gr_output_id, aname, tmpint4)
       aname = 'gatherwrite_factor'
       Call hdf5_write_attr(gr_output_id, aname, param_pfof%gatherwrite_factor)
    End Select
    Call hdf5_close_group(gr_output_id)
    Call hdf5_close_group(gr_pfof_id)
    Call hdf5_close_group(gr_id)
  End Subroutine h5write_meta_pfof_parameter

  !=======================================================================
  !> Write conecreator input parameters as metadata
  Subroutine h5write_meta_conecreator_parameter(file_id, param_cone)

    Implicit none

    Integer(kind=hid_t), intent(in) :: file_id
    Class(Type_parameter_conecreator), intent(in) :: param_cone

    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: aname

    Integer(kind=hid_t) :: gr_id
    Integer(kind=hid_t) :: gr_param_id
    Integer(kind=hid_t) :: gr_input_id
    Integer(kind=hid_t) :: gr_output_id

    Integer(kind=4) :: tmpint4
    Character(len=16) :: c_b_name

    Select Type(param_cone)
    Type Is(Type_parameter_conecreator_part)
       c_b_name = NAME_CONECREATOR_PART
    Type Is(Type_parameter_conecreator_grav)
       c_b_name = NAME_CONECREATOR_GRAV
    End Select

    groupname = 'metadata'
    Call hdf5_open_group(file_id, groupname, gr_id)

    groupname = c_b_name//'_parameters'
    Call hdf5_create_group(gr_id, groupname, gr_param_id)

    groupname = 'input_parameters'
    Call hdf5_create_group(gr_param_id, groupname, gr_input_id)

    aname = 'input_path'
    Call hdf5_write_attr(gr_input_id, aname, param_cone%input_path)

    aname = 'cone_input_file'
    Call hdf5_write_attr(gr_input_id, aname, param_cone%cone_input_file)

    aname = 'info_cone_input_file'
    Call hdf5_write_attr(gr_input_id, aname, param_cone%info_cone_input_file)

    aname = 'info_ramses_input_file'
    Call hdf5_write_attr(gr_input_id, aname, param_cone%info_ramses_input_file)

    aname = 'nfile'
    Call hdf5_write_attr(gr_input_id, aname, param_cone%nfile)

    aname = 'first_file'
    Call hdf5_write_attr(gr_input_id, aname, param_cone%first_file)

    aname = 'last_file'
    Call hdf5_write_attr(gr_input_id, aname, param_cone%last_file)

    aname = 'cone_max_radius'
    Call hdf5_write_attr(gr_input_id, aname, param_cone%cone_max_radius)

    Select Type(param_cone)
    Type Is (Type_parameter_conecreator_part)
       aname = 'do_read_ramses_part_id'
       tmpint4 = 0 
       If(param_cone%do_read_ramses_part_id) tmpint4 = 1
       Call hdf5_write_attr(gr_input_id, aname, tmpint4)
       aname = 'do_read_potential'
       tmpint4 = 0 
       If(param_cone%do_read_potential) tmpint4 = 1
       Call hdf5_write_attr(gr_input_id, aname, tmpint4)
       aname = 'do_read_gravitational_field'
       tmpint4 = 0 
       If(param_cone%do_read_gravitational_field) tmpint4 = 1
       Call hdf5_write_attr(gr_input_id, aname, tmpint4)
    Type Is (Type_parameter_conecreator_grav)
       aname = 'nlevel'
       Call hdf5_write_attr(gr_input_id, aname, param_cone%nlevel)
       aname = 'level_min'
       Call hdf5_write_attr(gr_input_id, aname, param_cone%level_min)
    End Select

    Call hdf5_close_group(gr_input_id)

    groupname = 'output_parameters'
    Call hdf5_create_group(gr_param_id, groupname, gr_output_id)

    aname = 'simulation_name'
    Call hdf5_write_attr(gr_output_id, aname, param_cone%simulation_name)

    aname = 'cube_size'
    Call hdf5_write_attr(gr_input_id, aname, param_cone%cube_size)

    Call hdf5_close_group(gr_output_id)
    Call hdf5_close_group(gr_param_id)
    Call hdf5_close_group(gr_id)

  End Subroutine h5write_meta_conecreator_parameter


  !=======================================================================
  !> Write ramses info as metadata
  Subroutine h5write_meta_info_ramses(file_id, simu_info)

    Implicit none

    Integer(kind=hid_t), intent(in) :: file_id
    Type(Type_info_ramses), intent(in) :: simu_info

    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: aname

    Integer(kind=hid_t) :: gr_id
    Integer(kind=hid_t) :: gr_ramses_id

    groupname = 'metadata'
    Call hdf5_open_group(file_id, groupname, gr_id)

    ! Ramses Info Metadata
    groupname = 'ramses_info'
    Call hdf5_create_group(gr_id,groupname,gr_ramses_id)

    aname = 'ncpu'
    Call hdf5_write_attr(gr_ramses_id,aname,simu_info%ncpu)

    aname = 'ndim'
    Call hdf5_write_attr(gr_ramses_id,aname,simu_info%ndim)

    aname = 'levelmin'
    Call hdf5_write_attr(gr_ramses_id,aname,simu_info%levelmin)

    aname = 'levelmax'
    Call hdf5_write_attr(gr_ramses_id,aname,simu_info%levelmax)

    aname = 'ngridmax'
    Call hdf5_write_attr(gr_ramses_id,aname,simu_info%ngridmax)

    aname = 'nstep_coarse'
    Call hdf5_write_attr(gr_ramses_id,aname,simu_info%nstep_coarse)

    aname = 'boxlen'
    Call hdf5_write_attr(gr_ramses_id,aname,simu_info%boxlen)

    aname = 'time'
    Call hdf5_write_attr(gr_ramses_id,aname,simu_info%time)

    aname = 'aexp'
    Call hdf5_write_attr(gr_ramses_id,aname,simu_info%aexp)

    aname = 'h0'
    Call hdf5_write_attr(gr_ramses_id,aname,simu_info%h0)

    aname = 'omega_m'
    Call hdf5_write_attr(gr_ramses_id,aname,simu_info%omega_m)

    aname = 'omega_l'
    Call hdf5_write_attr(gr_ramses_id,aname,simu_info%omega_l)

    aname = 'omega_k'
    Call hdf5_write_attr(gr_ramses_id,aname,simu_info%omega_k)

    aname = 'omega_b'
    Call hdf5_write_attr(gr_ramses_id,aname,simu_info%omega_b)

    aname = 'unit_l'
    Call hdf5_write_attr(gr_ramses_id,aname,simu_info%unit_l)

    aname = 'unit_d'
    Call hdf5_write_attr(gr_ramses_id,aname,simu_info%unit_d)

    aname = 'unit_t'
    Call hdf5_write_attr(gr_ramses_id,aname,simu_info%unit_t)

    Call hdf5_close_group(gr_ramses_id)

    Call hdf5_close_group(gr_id)

  End Subroutine h5write_meta_info_ramses


  !=======================================================================
  !> Write light cone info as metadata
  Subroutine h5write_meta_info_cone(file_id, infocone)

    Implicit none

    Integer(kind=hid_t), intent(in) :: file_id
    Class(Type_info_cone), intent(in) :: infocone

    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: aname

    Integer(kind=hid_t) :: gr_id
    Integer(kind=hid_t) :: gr_cone_id

    groupname = 'metadata'
    Call hdf5_open_group(file_id, groupname, gr_id)

    groupname = 'cone_info'
    Call hdf5_create_group(gr_id, groupname, gr_cone_id)

    aname = 'ncpu'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%ncpu)

    aname = 'nstride'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%ncpu)

    aname = 'nstep_coarse'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%nstep_coarse)

    aname = 'cone_id'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%cone_id)

    aname = 'nglobalfile'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%nglobalfile)

    aname = 'isfullsky'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%isfullsky)

    aname = 'aexp'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%aexp)

    aname = 'observer_x'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%observer_x)

    aname = 'observer_y'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%observer_y)

    aname = 'observer_z'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%observer_z)

    aname = 'observer_rds'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%observer_rds)

    aname = 'cone_zlim'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%cone_zlim)

    aname = 'amax'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%amax)

    aname = 'amin'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%amin)

    aname = 'zmax'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%zmax)

    aname = 'zmin'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%zmin)

    aname = 'dmax'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%dmax)

    aname = 'dmin'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%dmin)

    aname = 'dtol'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%dtol)

    aname = 'thetay'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%thetay)

    aname = 'thetaz'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%thetaz)

    aname = 'theta'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%theta)

    aname = 'phi'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%phi)

    aname = 'aendconem2'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%aendconem2)

    aname = 'aendconem1'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%aendconem1)

    aname = 'aendcone'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%aendcone)

    aname = 'zendconem2'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%zendconem2)

    aname = 'zendconem1'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%zendconem1)

    aname = 'zendcone'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%zendcone)

    aname = 'dendconem2'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%dendconem2)

    aname = 'dendconem1'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%dendconem1)

    aname = 'dendcone'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%dendcone)

    aname = 'future'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%future)

    Select Type(infocone)
    Type Is(Type_info_cone_part)
       aname = 'npart'
       Call hdf5_write_data(gr_cone_id, aname, infocone%npart)
       aname = 'aexpold'
       Call hdf5_write_attr(gr_cone_id, aname, infocone%aexpold)
       aname = 'zexpold'
       Call hdf5_write_attr(gr_cone_id, aname, infocone%zexpold)
       aname = 'zexp'
       Call hdf5_write_attr(gr_cone_id, aname, infocone%zexp)
       aname = 'dexpold'
       Call hdf5_write_attr(gr_cone_id, aname, infocone%dexpold)
       aname = 'dexp'
       Call hdf5_write_attr(gr_cone_id, aname, infocone%dexp)
    Type Is(Type_info_cone_grav)
       aname = 'nglobalcell'
       Call hdf5_write_data(gr_cone_id, aname, infocone%nglobalcell)
       aname = 'nlevel'
       Call hdf5_write_attr(gr_cone_id, aname, infocone%nlevel)
       aname = 'level_min'
       Call hdf5_write_attr(gr_cone_id, aname, infocone%level_min)
       aname = 'level_max'
       Call hdf5_write_attr(gr_cone_id, aname, infocone%level_max)
    End Select


    Call hdf5_close_group(gr_cone_id)

    Call hdf5_close_group(gr_id)

  End Subroutine h5write_meta_info_cone

  !=======================================================================
  !> Write units and constants used in the haloanalyzer code
  Subroutine h5write_info_haloanalyzer(file_id)

    Implicit none

    Integer(kind=hid_t), intent(in) :: file_id
    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: aname
    Character(len=H5STRLEN) :: adata

    Integer(kind=hid_t) :: gr_id, def_id, const_id,param_id,un_id,propdef_id, profdef_id,propun_id, profun_id,aparam_id,inp_id,out_id
    Integer(kind=4) :: tmpint4
    
    groupname = 'metadata'
    Call hdf5_open_group(file_id, groupname, gr_id)

    groupname = 'haloanalyzer_definitions'
    Call hdf5_create_group(gr_id, groupname, def_id)

    
    groupname = 'properties'
    Call hdf5_create_group(def_id, groupname, propdef_id)

    aname = 'identity_halo'
    adata = 'Halo identity '
    Call hdf5_write_attr(propdef_id, aname, adata)

    aname = 'number_particles_halo'
    adata = 'Number of particles in the halo '
    Call hdf5_write_attr(propdef_id, aname, adata)

    aname = 'mass_halo'
    adata = 'Mass of the halo '
    Call hdf5_write_attr(propdef_id, aname, adata)


    aname = 'normalisation_mass_halo'
    adata = 'Mass normalisation for other quantities'
    Call hdf5_write_attr(propdef_id, aname, adata)

    aname = 'normalisation_radius_physical_halo'
    adata = 'Radius normalisation for other quantities (given in physical). Usually deduced from mass normalisation and delta_m'
    Call hdf5_write_attr(propdef_id, aname, adata)

    aname = 'normalisation_radius_comoving_halo'
    adata = 'Radius normalisation for other quantities (given in comoving). Usually deduced from mass normalisation and delta_m'
    Call hdf5_write_attr(propdef_id, aname, adata)

    aname = 'normalisation_velocity_halo'
    adata = 'Velocity normalisation for other quantities. Usually deduced from mass normalisation and delta_m '
    Call hdf5_write_attr(propdef_id, aname, adata)

    
    aname = 'position_center_of_mass_halo'
    adata = 'Center of mass position '
    Call hdf5_write_attr(propdef_id, aname, adata)

    aname = 'velocity_center_of_mass_halo'
    adata = 'Velocity of center of mass '
    Call hdf5_write_attr(propdef_id, aname, adata)

    aname = 'radius_maximum_halo'
    adata = 'Maximum distance of particle to the center of mass'
    Call hdf5_write_attr(propdef_id, aname, adata)

    aname = 'velocity_maximum_halo'
    adata = 'Maximum relative velocity from center of mass velocity '
    Call hdf5_write_attr(propdef_id, aname, adata)

    aname = 'dispersion_position_halo'
    adata = 'Dispersion of position '
    Call hdf5_write_attr(propdef_id, aname, adata)

    aname = 'dispersion_velocity_halo'
    adata = 'Dispersion of velocity '
    Call hdf5_write_attr(propdef_id, aname, adata)

    

    aname = 'energy_self_binding_halo'
    adata = 'Self binding energy '
    Call hdf5_write_attr(propdef_id, aname, adata)

    aname = 'energy_kinetic_halo'
    adata = 'Kinetic energy '
    Call hdf5_write_attr(propdef_id, aname, adata)

    

    aname = 'angular_momentum_lambda_prime_halo'
    adata = 'Normalized angular momentum '
    Call hdf5_write_attr(propdef_id, aname, adata)

    aname = 'angular_momentum_lambda_halo'
    adata = 'Another normalized angular momentum '
    Call hdf5_write_attr(propdef_id, aname, adata)

    aname = 'inertia_minor_axis_vector_halo'
    adata = 'Eigen Vector of the minor axis '
    Call hdf5_write_attr(propdef_id, aname, adata)

    aname = 'inertia_eigen_values_halo'
    adata = 'Eigen value of the inertia tensor '
    Call hdf5_write_attr(propdef_id, aname, adata)

    aname = 'cosine_angle_angular_momentum_minor_axis_halo'
    adata = 'Cosine of the angle between the angular momentum direction and minor axis '
    Call hdf5_write_attr(propdef_id, aname, adata)

    aname = 'position_most_bounded_halo'
    adata = 'Position of the most bounded particles x,y,z '
    Call hdf5_write_attr(propdef_id, aname, adata)

    aname = 'circular_values_halos'
    adata = 'Mass, radius, velocity and enclosed density at position of maximum circular velocity: M_circ_max, R_circ_max, V_circ_max, Rho_circ_max'
    Call hdf5_write_attr(propdef_id, aname, adata)

    aname = 'delta_values_halo'
    adata = 'Mass and radius where overdensity reaches Delta_matter. Only halo particles are considered : Mdeltam, Rdeltam'
    Call hdf5_write_attr(propdef_id, aname, adata)

    aname = 'volume_halo'
    adata = 'Rough estimate of the halo volume. Method: pave halo volume with small cubes and sum volume of pencil excluding outer empty cells.'
    Call hdf5_write_attr(propdef_id, aname, adata)


    Call hdf5_close_group(propdef_id)

    groupname = 'profiles'
    Call hdf5_create_group(def_id, groupname, profdef_id)
    
    aname = 'profile_radial_bins_halo'
    adata = 'Radial bins for the profile. Center of logarithmic bins (for local profile). Outer shell (for integrated or circular profile) is a factor exp(delta_lr/2) away.'
    Call hdf5_write_attr(profdef_id, aname, adata)

    aname = 'profile_radial_bins_rdelta_halo'
    adata = 'Radial bins relative to Rdeltam for the profile. Center of logarithmic bins (for local profile). Outer shell (for integrated or circular profile) is a factor exp(delta_lr178/2) away.'
    Call hdf5_write_attr(profdef_id, aname, adata)

    aname = 'profile_circular_velocity_halo'
    adata = 'Circular velocity profile deduced from enclosed mass '
    Call hdf5_write_attr(profdef_id, aname, adata)

    aname = 'profile_density_halo'
    adata = 'Density profile computed using halo particles. Bins are proportional to ln R'
    Call hdf5_write_attr(profdef_id, aname, adata)

    aname = 'profile_density_integrated_halo'
    adata = 'Enclosed overdensity computed using halo particles. Bins are proportional to ln R '
    Call hdf5_write_attr(profdef_id, aname, adata)

    aname = 'profile_density_bins_rdelta_halo'
    adata = 'Density profile computed using halo particles. Bins are proportional to ln R/Rdeltam '
    Call hdf5_write_attr(profdef_id, aname, adata)

    aname = 'profile_density_integrated_bins_rdelta_halo'
    adata = 'Enclosed overdensity computed using halo particles. Bins are proportional to ln R/Rdeltam '
    Call hdf5_write_attr(profdef_id, aname, adata)

    
    Call hdf5_close_group(profdef_id)


    Call hdf5_close_group(def_id)
    Call hdf5_close_group(gr_id)


    !write units used in the haloanalyzer code
    
    groupname = 'metadata'
    Call hdf5_open_group(file_id, groupname, gr_id)

    groupname = 'haloanalyzer_units'
    Call hdf5_create_group(gr_id, groupname, un_id)

    
    groupname = 'properties'
    Call hdf5_create_group(un_id, groupname, propun_id)

    aname = 'identity_halo'
    adata = 'dimensionless'
    Call hdf5_write_attr(propun_id, aname, adata)

    aname = 'number_particles_halo'
    adata = 'dimensionless'
    Call hdf5_write_attr(propun_id, aname, adata)

    aname = 'mass_halo'
    adata = 'Msun/h'
    Call hdf5_write_attr(propun_id, aname, adata)


    aname = 'normalisation_mass_halo'
    adata = 'Msun/h'
    Call hdf5_write_attr(propun_id, aname, adata)


    aname = 'normalisation_radius_physical_halo'
    adata = 'kpc/h physical'
    Call hdf5_write_attr(propun_id, aname, adata)

    aname = 'normalisation_radius_comoving_halo'
    adata = 'kpc/h comoving'
    Call hdf5_write_attr(propun_id, aname, adata)

    aname = 'normalisation_velocity_halo'
    adata = 'km/s'
    Call hdf5_write_attr(propun_id, aname, adata)

    
    aname = 'position_center_of_mass_halo'
    adata = 'Mpc/h comoving'
    Call hdf5_write_attr(propun_id, aname, adata)

    aname = 'velocity_center_of_mass_halo'
    adata = 'km/s'
    Call hdf5_write_attr(propun_id, aname, adata)

    aname = 'radius_maximum_halo'
    adata = 'normalisation_radius_halo'
    Call hdf5_write_attr(propun_id, aname, adata)

    aname = 'velocity_maximum_halo'
    adata = 'normalisation_velocity_halo'
    Call hdf5_write_attr(propun_id, aname, adata)

    aname = 'dispersion_position_halo'
    adata = 'normalisation_radius_halo'
    Call hdf5_write_attr(propun_id, aname, adata)

    aname = 'dispersion_velocity_halo'
    adata = 'normalisation_velocity_halo'
    Call hdf5_write_attr(propun_id, aname, adata)

    

    aname = 'energy_self_binding_halo'
    adata = '0.5*normalisation_mass_halo*normalisation_velocity_halo^2'
    Call hdf5_write_attr(propun_id, aname, adata)

    aname = 'energy_kinetic_halo'
    adata = '0.5*normalisation_mass_halo*normalisation_velocity_halo^2'
    Call hdf5_write_attr(propun_id, aname, adata)

    

    aname = 'angular_momentum_lambda_prime_halo'
    adata = 'sqrt(2)*normalisation_mass_halo*normalisation_velocity_halo*normalisation_radius_halo'
    Call hdf5_write_attr(propun_id, aname, adata)

    aname = 'angular_momentum_lambda_halo'
    adata = 'G*normalisation_mass_halo^2.5)/sqrt(abs(energy_self_binding_halo))'
    Call hdf5_write_attr(propun_id, aname, adata)

    aname = 'inertia_minor_axis_vector_halo'
    adata = 'dimensionless'
    Call hdf5_write_attr(propun_id, aname, adata)

    aname = 'inertia_eigen_values_halo'
    adata = 'normalized by normalisation_radius_halo'
    Call hdf5_write_attr(propun_id, aname, adata)

    aname = 'cosine_angle_angular_momentum_minor_axis_halo'
    adata = 'dimensionless'
    Call hdf5_write_attr(propun_id, aname, adata)

    aname = 'position_most_bounded_halo'
    adata = 'Mpc/h comoving'
    Call hdf5_write_attr(propun_id, aname, adata)

    aname = 'circular_values_halos'
    adata = 'M_circ_max(Msun/h), R_circ_max(kpc/h physical), V_circ_max(km/s), Rho_circ_max (background matter)'
    Call hdf5_write_attr(propun_id, aname, adata)

    aname = 'delta_values_halo'
    adata = 'Mdeltam(Msun/h), Rdeltam(kpc/h physical)'
    Call hdf5_write_attr(propun_id, aname, adata)

    aname = 'volume_halo'
    adata = ',normalized by 4/3*pi*normalisation_radius^3'
    Call hdf5_write_attr(propun_id, aname, adata)

    Call hdf5_close_group(propun_id)

    groupname = 'profiles'
    Call hdf5_create_group(un_id, groupname, profun_id)
    
    aname = 'profile_radial_bins_halo'
    adata = 'kpc/h physical'
    Call hdf5_write_attr(profun_id, aname, adata)

    aname = 'profile_radial_bins_rdelta_halo'
    adata = 'normalized by Rdeltam'
    Call hdf5_write_attr(profun_id, aname, adata)

    aname = 'profile_circular_velocity_halo'
    adata = 'km/s'
    Call hdf5_write_attr(profun_id, aname, adata)

    aname = 'profile_density_halo'
    adata = 'background matter'
    Call hdf5_write_attr(profun_id, aname, adata)

    aname = 'profile_density_integrated_halo'
    adata = 'background matter'
    Call hdf5_write_attr(profun_id, aname, adata)

    aname = 'profile_density_bins_rdelta_halo'
    adata = 'background matter'
    Call hdf5_write_attr(profun_id, aname, adata)

    aname = 'profile_density_integrated_bins_rdelta_halo'
    adata = 'background matter'
    Call hdf5_write_attr(profun_id, aname, adata)

    
    Call hdf5_close_group(profun_id)


    Call hdf5_close_group(un_id)
    Call hdf5_close_group(gr_id)

    !write constants used in the haloanalyzer code
    groupname = 'metadata'
    Call hdf5_open_group(file_id, groupname, gr_id)
    
    groupname = 'haloanalyzer_constants'
    Call hdf5_create_group(gr_id, groupname, const_id)


    aname = 'msun'
    Call hdf5_write_attr(const_id, aname, msun)

    aname = 'mpc'
    Call hdf5_write_attr(const_id, aname, mpc) 

    aname = 'rho_crit0'
    Call hdf5_write_attr(const_id, aname, rho_crit0) 

    aname = 'pi'
    Call hdf5_write_attr(const_id, aname, pi)

    aname = 'g'
    Call hdf5_write_attr(const_id, aname, g)

    aname = 'clight'
    Call hdf5_write_attr(const_id, aname, clight)
    

    Call hdf5_close_group(const_id)
    Call hdf5_close_group(gr_id)




!write parameters used in the haloanalyzer code
    groupname = 'metadata'
    Call hdf5_open_group(file_id, groupname, gr_id)
    
    groupname = 'haloanalyzer_parameters'
    Call hdf5_create_group(gr_id, groupname, param_id)

    groupname = 'input_parameters'
    Call hdf5_create_group(param_id, groupname, inp_id)
    
    aname = 'file_cosmo'
    Call hdf5_write_attr(inp_id, aname, file_cosmo)
    
    aname = 'filelist(1)'
    Call hdf5_write_attr(inp_id, aname, filelist(1))
    
    aname = 'filelistsize'
    Call hdf5_write_attr(inp_id, aname, filelistsize)
    
    Call hdf5_close_group(inp_id)

    groupname = 'output_parameters'
    Call hdf5_create_group(param_id, groupname, out_id)
    
    aname = 'filenameh5'
    Call hdf5_write_attr(out_id, aname, filenameh5)

    Call hdf5_close_group(out_id)    

    groupname = 'analyzer_parameters'
    Call hdf5_create_group(param_id, groupname, aparam_id)

    aname = 'do_lightcone'
    If(do_lightcone) Then
       tmpint4=1
    Else
       tmpint4=0
    End If
    Call hdf5_write_attr(aparam_id, aname, tmpint4)

    aname = 'do_periodic_boundary'
    If(do_periodic_boundary) Then
       tmpint4=1
    Else
       tmpint4=0
    End If
    Call hdf5_write_attr(aparam_id, aname, tmpint4)

    aname = 'do_center_mass'
    If(do_center_mass) Then
       tmpint4=1
    Else
       tmpint4=0
    End If
    Call hdf5_write_attr(aparam_id, aname, tmpint4)

    aname = 'delta_matter'
    Call hdf5_write_attr(aparam_id, aname, delta_matter)
    
    aname = 'delta_critz'
    Call hdf5_write_attr(aparam_id, aname, delta_critz)
    
    aname = 'min_r'
    Call hdf5_write_attr(aparam_id, aname, min_r)
    
    aname = 'min_r178'
    Call hdf5_write_attr(aparam_id, aname, min_r178)

    aname = 'delta_lr'
    Call hdf5_write_attr(aparam_id, aname, delta_lr)

    aname = 'delta_lr178'
    Call hdf5_write_attr(aparam_id, aname, delta_lr178)

    aname = 'nbin1'
    Call hdf5_write_attr(aparam_id, aname, nbin1)

    aname = 'nbin2'
    Call hdf5_write_attr(aparam_id, aname, nbin2)

    aname = 'box_len_mpc'
    Call hdf5_write_attr(aparam_id, aname, box_len_mpc)

    aname = 'simu_info_delta_vir'
    Call hdf5_write_attr(aparam_id, aname, simu_info%delta_vir)
    
    
    Call hdf5_close_group(aparam_id)
    Call hdf5_close_group(param_id)
    Call hdf5_close_group(gr_id)



    

    
  End Subroutine h5write_info_haloanalyzer
  !=======================================================================
  !> Read ramses info as metadata
  Subroutine h5read_meta_info_ramses(file_id, simu_info)

    Implicit none
    Integer(kind=hid_t), intent(in) :: file_id
    Type(Type_info_ramses), intent(out) :: simu_info
    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: aname
    Integer(kind=hid_t) :: gr_id
    Integer(kind=hid_t) :: gr_ramses_id

    groupname = 'metadata'
    Call hdf5_open_group(file_id, groupname, gr_id)
    ! Ramses Info Metadata
    groupname = 'ramses_info'
    Call hdf5_open_group(gr_id,groupname,gr_ramses_id)

    aname = 'ncpu'
    Call hdf5_read_attr(gr_ramses_id,aname,simu_info%ncpu)
    aname = 'ndim'
    Call hdf5_read_attr(gr_ramses_id,aname,simu_info%ndim)
    aname = 'levelmin'
    Call hdf5_read_attr(gr_ramses_id,aname,simu_info%levelmin)
    aname = 'levelmax'
    Call hdf5_read_attr(gr_ramses_id,aname,simu_info%levelmax)
    aname = 'ngridmax'
    Call hdf5_read_attr(gr_ramses_id,aname,simu_info%ngridmax)
    aname = 'nstep_coarse'
    Call hdf5_read_attr(gr_ramses_id,aname,simu_info%nstep_coarse)
    aname = 'boxlen'
    Call hdf5_read_attr(gr_ramses_id,aname,simu_info%boxlen)
    aname = 'time'
    Call hdf5_read_attr(gr_ramses_id,aname,simu_info%time)
    aname = 'aexp'
    Call hdf5_read_attr(gr_ramses_id,aname,simu_info%aexp)
    aname = 'h0'
    Call hdf5_read_attr(gr_ramses_id,aname,simu_info%h0)
    aname = 'omega_m'
    Call hdf5_read_attr(gr_ramses_id,aname,simu_info%omega_m)
    aname = 'omega_l'
    Call hdf5_read_attr(gr_ramses_id,aname,simu_info%omega_l)
    aname = 'omega_k'
    Call hdf5_read_attr(gr_ramses_id,aname,simu_info%omega_k)
    aname = 'omega_b'
    Call hdf5_read_attr(gr_ramses_id,aname,simu_info%omega_b)
    aname = 'unit_l'
    Call hdf5_read_attr(gr_ramses_id,aname,simu_info%unit_l)
    aname = 'unit_d'
    Call hdf5_read_attr(gr_ramses_id,aname,simu_info%unit_d)
    aname = 'unit_t'
    Call hdf5_read_attr(gr_ramses_id,aname,simu_info%unit_t)
    Call hdf5_close_group(gr_ramses_id)
    Call hdf5_close_group(gr_id)
  End Subroutine h5read_meta_info_ramses

  !=======================================================================
  !> Read light cone info as metadata 
  Subroutine h5read_meta_info_cone(file_id, infocone)

    Implicit none

    Integer(kind=hid_t), intent(in) :: file_id
    Class(Type_info_cone), intent(out) :: infocone

    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: aname

    Integer(kind=hid_t) :: gr_id
    Integer(kind=hid_t) :: gr_cone_id

    groupname = 'metadata'
    Call hdf5_open_group(file_id, groupname, gr_id)

    groupname = 'cone_info'
    Call hdf5_open_group(gr_id, groupname, gr_cone_id)

    groupname = 'cone_info'
    Call hdf5_open_group(gr_id,groupname, gr_cone_id)

    aname= "ncpu"
    Call hdf5_read_attr(gr_cone_id, aname, infocone%ncpu)

    aname = 'nstride'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%ncpu)

    aname = 'nstep_coarse'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%nstep_coarse)

    aname = 'cone_id'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%cone_id)

    aname = 'nglobalfile'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%nglobalfile)

    aname = 'isfullsky'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%isfullsky)

    aname = 'aexp'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%aexp)

    aname = 'observer_x'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%observer_x)

    aname = 'observer_y'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%observer_y)

    aname = 'observer_z'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%observer_z)

    aname = 'observer_rds'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%observer_rds)

    aname = 'cone_zlim'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%cone_zlim)

    aname = 'amax'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%amax)

    aname = 'amin'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%amin)

    aname = 'zmax'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%zmax)

    aname = 'zmin'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%zmin)

    aname = 'dmax'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%dmax)

    aname = 'dmin'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%dmin)

    aname = 'dtol'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%dtol)

    aname = 'thetay'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%thetay)

    aname = 'thetaz'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%thetaz)

    aname = 'theta'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%theta)

    aname = 'phi'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%phi)

    aname = 'aendconem2'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%aendconem2)

    aname = 'aendconem1'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%aendconem1)

    aname = 'aendcone'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%aendcone)

    aname = 'zendconem2'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%zendconem2)

    aname = 'zendconem1'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%zendconem1)

    aname = 'zendcone'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%zendcone)

    aname = 'dendconem2'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%dendconem2)

    aname = 'dendconem1'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%dendconem1)

    aname = 'dendcone'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%dendcone)

    aname = 'future'
    Call hdf5_read_attr(gr_cone_id, aname, infocone%future)


    Select Type(infocone)
    Type Is(Type_info_cone_part)
       aname = 'npart'
       Call hdf5_read_data(gr_cone_id, aname, infocone%npart)
       aname = 'aexpold'
       Call hdf5_read_attr(gr_cone_id, aname, infocone%aexpold)
       aname = 'zexpold'
       Call hdf5_read_attr(gr_cone_id, aname, infocone%zexpold)
       aname = 'zexp'
       Call hdf5_read_attr(gr_cone_id, aname, infocone%zexp)
       aname = 'dexpold'
       Call hdf5_read_attr(gr_cone_id, aname, infocone%dexpold)
       aname = 'dexp'
       Call hdf5_read_attr(gr_cone_id, aname, infocone%dexp)
    Type Is(Type_info_cone_grav)
       aname = 'nglobalcell'
       Call hdf5_read_data(gr_cone_id, aname, infocone%nglobalcell)
       aname = 'nlevel'
       Call hdf5_read_attr(gr_cone_id, aname, infocone%nlevel)
       aname = 'level_min'
       Call hdf5_read_attr(gr_cone_id, aname, infocone%level_min)
       aname = 'level_max'
       Call hdf5_read_attr(gr_cone_id, aname, infocone%level_max)
    End Select

    Call hdf5_close_group(gr_cone_id)

    Call hdf5_close_group(gr_id)

  End Subroutine h5read_meta_info_cone

  !=======================================================================
  !> Read pfof input parameters as metadata 
  Subroutine h5read_meta_pfof_parameter(file_id, param_pfof)

    Implicit none

    Integer(kind=hid_t), intent(in) :: file_id
    Class(Type_parameter_pfof), intent(out) :: param_pfof

    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: aname

    Integer(kind=hid_t) :: gr_id
    Integer(kind=hid_t) :: gr_pfof_id
    Integer(kind=hid_t) :: gr_input_id
    Integer(kind=hid_t) :: gr_fof_id
    Integer(kind=hid_t) :: gr_output_id

    Integer(kind=4) :: tmpint4
    Character(len=16) :: c_b_name

    Select Type (param_pfof)
    Type is (Type_parameter_pfof_snap)
       c_b_name = NAME_PFOF_SNAP
    Type is (Type_parameter_pfof_cone)
       c_b_name = NAME_PFOF_CONE
    End Select

    groupname = 'metadata'
    Call hdf5_open_group(file_id,groupname, gr_id)
    ! pfof parameters:
    groupname = trim(c_b_name)//'_parameters'
    Call hdf5_open_group(gr_id, groupname, gr_pfof_id)
    groupname = 'input_parameters'
    Call hdf5_open_group(gr_pfof_id, groupname, gr_input_id)
    aname = 'input_path'
    Call hdf5_read_attr(gr_input_id, aname, len(param_pfof%input_path), &
         param_pfof%input_path)
    aname = 'part_input_file'
    Call hdf5_read_attr(gr_input_id, aname, len(param_pfof%part_input_file),&
         param_pfof%part_input_file)
    aname = 'do_read_potential'
    Call hdf5_read_attr(gr_input_id, aname, tmpint4)
    param_pfof%do_read_potential = .false.
    If(tmpint4==1) param_pfof%do_read_potential = .true.
    aname = 'do_read_gravitational_field'
    Call hdf5_read_attr(gr_input_id, aname, tmpint4)
    param_pfof%do_read_gravitational_field = .false.
    If(tmpint4==1) param_pfof%do_read_gravitational_field = .true.

    Select Type (param_pfof)
    Type Is (Type_parameter_pfof_snap)
       aname = 'info_input_file'
       Call hdf5_read_attr(gr_input_id, aname, len(param_pfof%info_input_file), &
            param_pfof%info_input_file)

       aname = 'grpsize'
       Call hdf5_read_attr(gr_input_id, aname, param_pfof%grpsize)
       aname = 'code_index'
       Call hdf5_read_attr(gr_input_id, aname, len(param_pfof%code_index), param_pfof%code_index)
       aname = 'do_skip_star'
       Call hdf5_read_attr(gr_input_id, aname, tmpint4)
       param_pfof%do_skip_star = .false.
       If(tmpint4==1) param_pfof%do_skip_star = .true.
       aname = 'do_skip_metal'
       Call hdf5_read_attr(gr_input_id, aname, tmpint4)       
       param_pfof%do_skip_metal = .false.
       If(tmpint4==1) param_pfof%do_skip_metal = .true.
       aname='do_read_from_cube'
       Call hdf5_read_attr(gr_input_id, aname, tmpint4)
       param_pfof%do_read_from_cube = .false.
       If(tmpint4==1) param_pfof%do_read_from_cube = .true.
       aname = 'gatherread_factor'
       Call hdf5_read_attr(gr_input_id, aname, param_pfof%gatherread_factor)

    Type Is (Type_parameter_pfof_cone)
       aname = 'shell_first_id'
       Call hdf5_read_attr(gr_input_id, aname, param_pfof%shell_first_id)
       aname = 'shell_last_id'
       Call hdf5_read_attr(gr_input_id, aname, param_pfof%shell_last_id)
    End Select
    Call hdf5_close_group(gr_input_id)
#ifdef SOD
    groupname = 'sod_parameters'
#else
    groupname = 'fof_parameters'
#endif
    Call hdf5_open_group(gr_pfof_id, groupname, gr_fof_id)
    ! Write percolation parameter b as attribute
    aname = 'percolation_length'
    Call hdf5_read_attr(gr_fof_id, aname, param_pfof%percolation_length)
    ! Write minimum halo mass Mmin as attribute
    aname = 'npart_halo_min'
    Call hdf5_read_attr(gr_fof_id, aname, param_pfof%mmin)    
    aname = 'npart_halo_max'
    Call hdf5_read_attr(gr_fof_id, aname, param_pfof%mmax)
    ! Write doUnbinding as attribute (not implemented yet) 
    aname = 'do_unbinding'
    Call hdf5_read_attr(gr_fof_id, aname, tmpint4)
    param_pfof%do_unbinding = .false.
    If(tmpint4==1) param_pfof%do_unbinding = .true.
    ! Write doSubHalo as attribute (not implemented yet) 
    aname = 'do_subhalo'
    Call hdf5_read_attr(gr_fof_id, aname, tmpint4)
    param_pfof%do_subhalo = .false.
    If(tmpint4==1) param_pfof%do_subhalo = .true.

    Select Type (param_pfof)
    Type Is (Type_parameter_pfof_snap)
#ifdef SOD
       aname = 'do_sod'
#else
       aname = 'do_fof'
#endif
       Call hdf5_read_attr(gr_fof_id, aname, tmpint4)
       param_pfof%do_fof = .false.
       If(tmpint4==1) param_pfof%do_fof = .true.
    End Select
    Call hdf5_close_group(gr_fof_id)       

    groupname = 'output_parameters'
    Call hdf5_open_group(gr_pfof_id, groupname, gr_output_id)
    ! Write the simulation name as attribute
    aname = 'simulation_name'
    Call hdf5_read_attr(gr_output_id, aname, len(param_pfof%simulation_name), &
         param_pfof%simulation_name)
    ! Write do_timings as attribute
#ifndef SOD
    aname = 'do_timings'
    Call hdf5_read_attr(gr_fof_id, aname, tmpint4)
    param_pfof%do_timings = .false.
    If(tmpint4==1) param_pfof%do_timings = .true.
    !RY: not sure why there is an error message on irene skylake 21/11/23 so I assume do_timings=.false. for sod
#endif
    
    Select Type (param_pfof)
    Type Is (Type_parameter_pfof_snap)
       ! write the snapshot number: if halo computed from a lightcone then outputNB = -1
       aname = 'snapshot'
       Call hdf5_read_attr(gr_output_id, aname, param_pfof%snapshot)
       aname = 'do_write_cube'
       Call hdf5_read_attr(gr_output_id, aname, tmpint4)       
       param_pfof%do_write_cube = .false.
       If(tmpint4==1) param_pfof%do_write_cube = .true.
       aname = 'do_sort_cube'
       Call hdf5_read_attr(gr_output_id, aname, tmpint4)
       param_pfof%do_sort_cube = .false.
       If(tmpint4==1) param_pfof%do_sort_cube = .true.
       aname = 'gatherwrite_factor'
       Call hdf5_read_attr(gr_output_id, aname, param_pfof%gatherwrite_factor)
    End Select
    Call hdf5_close_group(gr_output_id)
    Call hdf5_close_group(gr_pfof_id)
    Call hdf5_close_group(gr_id)
  End Subroutine h5read_meta_pfof_parameter

  !=======================================================================
  !> Read conecreator input parameters as metadata
  Subroutine h5read_meta_conecreator_parameter(file_id, param_cone)

    Implicit none

    Integer(kind=hid_t), intent(in) :: file_id
    Class(Type_parameter_conecreator), intent(out) :: param_cone

    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: aname

    Integer(kind=hid_t) :: gr_id
    Integer(kind=hid_t) :: gr_param_id
    Integer(kind=hid_t) :: gr_input_id
    Integer(kind=hid_t) :: gr_output_id

    Integer(kind=4) :: tmpint4
    Character(len=16) :: c_b_name

    Select Type(param_cone)
    Type Is(Type_parameter_conecreator_part)
       c_b_name = NAME_CONECREATOR_PART
    Type Is(Type_parameter_conecreator_grav)
       c_b_name = NAME_CONECREATOR_GRAV
    End Select

    groupname = 'metadata'
    Call hdf5_open_group(file_id, groupname, gr_id)

    groupname = c_b_name//'_parameters'
    Call hdf5_open_group(gr_id, groupname, gr_param_id)

    groupname = 'input_parameters'
    Call hdf5_open_group(gr_param_id, groupname, gr_input_id)

    aname = 'input_path'
    Call hdf5_read_attr(gr_input_id, aname, len(param_cone%input_path), &
         param_cone%input_path)

    aname = 'cone_input_file'
    Call hdf5_read_attr(gr_input_id, aname, len(param_cone%cone_input_file), &
         param_cone%cone_input_file)

    aname = 'info_cone_input_file'
    Call hdf5_read_attr(gr_input_id, aname, len(param_cone%info_cone_input_file), &
         param_cone%info_cone_input_file)

    aname = 'info_ramses_input_file'
    Call hdf5_read_attr(gr_input_id, aname, len(param_cone%info_ramses_input_file), &
         param_cone%info_ramses_input_file)

    aname = 'nfile'
    Call hdf5_read_attr(gr_input_id, aname, param_cone%nfile)

    aname = 'first_file'
    Call hdf5_read_attr(gr_input_id, aname, param_cone%first_file)

    aname = 'last_file'
    Call hdf5_read_attr(gr_input_id, aname, param_cone%last_file)

    aname = 'cone_max_radius'
    Call hdf5_read_attr(gr_input_id, aname, param_cone%cone_max_radius)

    Select Type(param_cone)
    Type Is (Type_parameter_conecreator_part)
       aname = 'do_read_ramses_part_id'
       Call hdf5_read_attr(gr_input_id, aname, tmpint4)
       param_cone%do_read_ramses_part_id = .false.
       If(tmpint4==1) param_cone%do_read_ramses_part_id = .true.

    Type Is (Type_parameter_conecreator_grav)
       aname = 'nlevel'
       Call hdf5_read_attr(gr_input_id, aname, param_cone%cone_max_radius)
       aname = 'level_min'
       Call hdf5_read_attr(gr_input_id, aname, param_cone%cone_max_radius)
    End Select

    Call hdf5_close_group(gr_input_id)

    groupname = 'output_parameters'
    Call hdf5_open_group(gr_param_id, groupname, gr_output_id)

    aname = 'simulation_name'
    Call hdf5_read_attr(gr_output_id, aname, len(param_cone%simulation_name), &
         param_cone%simulation_name)

    aname = 'cube_size'
    Call hdf5_read_attr(gr_input_id, aname, param_cone%cube_size)

    Call hdf5_close_group(gr_output_id)
    Call hdf5_close_group(gr_param_id)
    Call hdf5_close_group(gr_id)


  End Subroutine h5read_meta_conecreator_parameter


  !=======================================================================
  !> Read metadata common to every hdf5 ouput files of the pfof toolbox
  Subroutine h5read_meta_common(file_id, common)
        
    Implicit none
    Integer(kind=hid_t), intent(in) :: file_id
    Type(Type_common_metadata), intent(out) :: common
    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: name
    Integer(kind=hid_t) :: meta_id
    Integer(kind=4) :: length

    groupname = 'metadata'
    Call hdf5_open_group(file_id, groupname, meta_id)
    length = len(common%created_by)
    name = 'created_by'
    Call hdf5_read_attr(meta_id, name, length, common%created_by)
    name = 'svn_version'
    Call hdf5_read_attr(meta_id, name, common%svn_version)
    length = len(common%simulation_code)
    name = 'simulation_code'
    Call hdf5_read_attr(meta_id, name, length, common%simulation_code)
    length = len(common%particle_type)
    name = 'particle_type'
    Call hdf5_read_attr(meta_id, name, length, common%particle_type)
    name = 'constant_mass'
    Call hdf5_read_attr(meta_id, name, common%constant_mass)
    length = len(common%units)
    name = 'units'
    Call hdf5_read_attr(meta_id, name, length, common%units) 
    name = 'npart_simulation'
    Call hdf5_read_data(meta_id, name, common%npart_simulation)
    Call hdf5_close_group(meta_id)
  End Subroutine h5read_meta_common

End Module modiocommons
