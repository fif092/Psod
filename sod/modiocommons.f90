Module modiocommons

  Use modconstant
  Use modhdf5
  Implicit none

  ! Some metadata constants whose value may change for some software or options 
  Character(len=H5STRLEN) :: PARTICLE_TYPE = 'dark_matter'
  Character(len=H5STRLEN) :: PFOF_CELL_ORDER = 'none'

Contains

  !=======================================================================
  Subroutine h5write_meta_common(file_id, codename, npart, process_id)

    Implicit none

    Integer(kind=hid_t), intent(in) :: file_id
    Character(len=H5STRLEN), intent(in) :: codename
    Integer(kind=8), intent(in) :: npart
    Integer(kind=4), intent(in), optional :: process_id

    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: aname
    Character(len=H5STRLEN) :: adata
    Character(len=H5STRLEN) :: dsetname

    Integer(kind=hid_t) :: gr_id
    Integer(kind=4) :: tmpint4

    groupname = 'metadata'
    Call hdf5_open_group(file_id,groupname, gr_id)

    ! Common metadata
    aname = 'created_by'
    Call hdf5_write_attr(gr_id, aname, codename)

    aname = 'svn_version'
    Call hdf5_write_attr(gr_id, aname, svn_version)

    aname = 'simulation_code'
    adata = 'ramses'
    Call hdf5_write_attr(gr_id, aname,adata)

    aname = 'particle_type'
    adata = PARTICLE_TYPE
    Call hdf5_write_attr(gr_id, aname, adata)

    aname = 'process_id'
    If(present(process_id)) Then
       tmpint4=process_id
    Else
       tmpint4=0
    End If
    Call hdf5_write_attr(gr_id, aname, tmpint4)

    aname = 'constant_mass'
    tmpint4 = 1
    Call hdf5_write_attr(gr_id, aname, tmpint4)

    aname = 'units'
    adata = 'ramses'
    Call hdf5_write_attr(gr_id, aname, adata) 

    aname = 'pfof_cells_order'
    adata = PFOF_CELL_ORDER
    Call hdf5_write_attr(gr_id, aname, adata)

    ! Write the number of particles in this simulation as an integer(kind=8) dataset
    dsetname = 'npart_simulation'
    Call hdf5_write_data(gr_id, dsetname, npart)

    Call hdf5_close_group(gr_id)

  End Subroutine h5write_meta_common

  !=======================================================================
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
    Character(len=16) :: c_b_name,hfname

    
    Select Type (param_pfof)
    Type is (Type_parameter_pfof_snap)
       c_b_name = NAME_PFOF_SNAP
       hfname='fof'
    Type is (Type_parameter_pfof_cone)
       c_b_name = NAME_PFOF_CONE
       hfname='fof'
    Type is (Type_parameter_psod_snap)
       c_b_name = NAME_PSOD_SNAP
       hfname='sod'
    End Select

    groupname = 'metadata'
    Call hdf5_open_group(file_id,groupname, gr_id)

    ! pfof parameters:
    groupname = trim(c_b_name)//'_parameters'
    Call hdf5_create_group(gr_id, groupname, gr_pfof_id)

    ! input parameter
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


    Type Is (Type_parameter_psod_snap)

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
    ! End of input parameters

    ! Friend of friend parameters
    groupname = trim(hfname)//'_parameters'
    Call hdf5_create_group(gr_pfof_id, groupname, gr_fof_id)

    ! Write percolation parameter b as attribute
    aname = 'percolation_length'
    Call hdf5_write_attr(gr_fof_id, aname, param_pfof%percolation_length)
    

    ! Write sod detection threshold as attribute
    Select Type (param_pfof)
    Type Is (Type_parameter_psod_snap)
       aname = 'matter_overdensity_threshold'
       Call hdf5_write_attr(gr_fof_id, aname, param_pfof%matter_overdensity_threshold)
    End Select

    

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

       aname = 'do_fof'
       If(param_pfof%do_fof) Then
          tmpint4=1
       Else
          tmpint4=0
       End If
       Call hdf5_write_attr(gr_fof_id, aname, tmpint4)

    Type Is (Type_parameter_psod_snap)

       aname = 'do_sod'
       If(param_pfof%do_sod) Then
          tmpint4=1
       Else
          tmpint4=0
       End If
       Call hdf5_write_attr(gr_fof_id, aname, tmpint4)

    End Select

    Call hdf5_close_group(gr_fof_id)
    ! end of friend of friend parameters

    ! output parameters
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
    Call hdf5_write_attr(gr_output_id, aname, tmpint4)

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


    Type Is (Type_parameter_psod_snap)

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
    ! end of output parameters
    
    Call hdf5_close_group(gr_pfof_id)

    Call hdf5_close_group(gr_id)


  End Subroutine h5write_meta_pfof_parameter


  !=======================================================================
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
       aname = 'levelmin'
       Call hdf5_write_attr(gr_input_id, aname, param_cone%levelmin)
    End Select

    Call hdf5_close_group(gr_input_id)

    groupname = 'output_parameters'
    Call hdf5_create_group(gr_param_id, groupname, gr_output_id)

    aname = 'simulation_name'
    Call hdf5_write_attr(gr_output_id, aname, param_cone%simulation_name)

    aname = 'cube_size'
    Call hdf5_write_attr(gr_output_id, aname, param_cone%cube_size)

    Call hdf5_close_group(gr_output_id)
    Call hdf5_close_group(gr_param_id)
    Call hdf5_close_group(gr_id)

  End Subroutine h5write_meta_conecreator_parameter


  !=======================================================================
  Subroutine h5write_meta_info_ramses(file_id, inforamses, islast)

    Implicit none

    Integer(kind=hid_t), intent(in) :: file_id
    Type(Type_inforamses), intent(in) :: inforamses
    Logical(kind=4), intent(in) :: islast

    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: aname

    Integer(kind=hid_t) :: gr_id
    Integer(kind=hid_t) :: gr_ramses_id

    groupname = 'metadata'
    Call hdf5_open_group(file_id, groupname, gr_id)

    ! Ramses Info Metadata
    groupname = 'ramses_info'
    If(islast) Then
       groupname = 'ramses_info_last'
    End If
    Call hdf5_create_group(gr_id,groupname,gr_ramses_id)

    aname = 'ncpu'
    Call hdf5_write_attr(gr_ramses_id,aname,inforamses%ncpu)

    aname = 'ndim'
    Call hdf5_write_attr(gr_ramses_id,aname,inforamses%ndim)

    aname = 'levelmin'
    Call hdf5_write_attr(gr_ramses_id,aname,inforamses%levelmin)

    aname = 'levelmax'
    Call hdf5_write_attr(gr_ramses_id,aname,inforamses%levelmax)

    aname = 'ngridmax'
    Call hdf5_write_attr(gr_ramses_id,aname,inforamses%ngridmax)

    aname = 'nstep_coarse'
    Call hdf5_write_attr(gr_ramses_id,aname,inforamses%nstep_coarse)

    aname = 'boxlen'
    Call hdf5_write_attr(gr_ramses_id,aname,inforamses%boxlen)

    aname = 'time'
    Call hdf5_write_attr(gr_ramses_id,aname,inforamses%time)

    aname = 'aexp'
    Call hdf5_write_attr(gr_ramses_id,aname,inforamses%aexp)

    aname = 'h0'
    Call hdf5_write_attr(gr_ramses_id,aname,inforamses%h0)

    aname = 'omega_m'
    Call hdf5_write_attr(gr_ramses_id,aname,inforamses%omega_m)

    aname = 'omega_l'
    Call hdf5_write_attr(gr_ramses_id,aname,inforamses%omega_l)

    aname = 'omega_k'
    Call hdf5_write_attr(gr_ramses_id,aname,inforamses%omega_k)

    aname = 'omega_b'
    Call hdf5_write_attr(gr_ramses_id,aname,inforamses%omega_b)

    aname = 'unit_l'
    Call hdf5_write_attr(gr_ramses_id,aname,inforamses%unit_l)

    aname = 'unit_d'
    Call hdf5_write_attr(gr_ramses_id,aname,inforamses%unit_d)

    aname = 'unit_t'
    Call hdf5_write_attr(gr_ramses_id,aname,inforamses%unit_t)

    Call hdf5_close_group(gr_ramses_id)

    Call hdf5_close_group(gr_id)

  End Subroutine h5write_meta_info_ramses


  !=======================================================================
  Subroutine h5write_meta_info_cone(file_id, infocone, islast)

    Implicit none

    Integer(kind=hid_t), intent(in) :: file_id
    Class(Type_infocone), intent(in) :: infocone
    Logical(kind=4), intent(in) :: islast

    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: aname

    Integer(kind=hid_t) :: gr_id
    Integer(kind=hid_t) :: gr_cone_id

    groupname = 'metadata'
    Call hdf5_open_group(file_id, groupname, gr_id)

    groupname = 'cone_info'
    If(islast) Then
       groupname = 'cone_info_last'
    End If
    Call hdf5_create_group(gr_id, groupname, gr_cone_id)

    aname = 'ncpu'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%ncpu)

    aname = 'nstride'
    Call hdf5_write_attr(gr_cone_id, aname, infocone%nstride)

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
    Type Is(Type_infocone_part)
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
    Type Is(Type_infocone_grav)
       aname = 'nglobalcell'
       Call hdf5_write_data(gr_cone_id, aname, infocone%nglobalcell)
       aname = 'nlevel'
       Call hdf5_write_attr(gr_cone_id, aname, infocone%nlevel)
       aname = 'levelmin'
       Call hdf5_write_attr(gr_cone_id, aname, infocone%levelmin)
       aname = 'levelmax'
       Call hdf5_write_attr(gr_cone_id, aname, infocone%levelmax)
    End Select


    Call hdf5_close_group(gr_cone_id)

    Call hdf5_close_group(gr_id)

  End Subroutine h5write_meta_info_cone


  !=======================================================================
  Subroutine h5readcommonmetadata(file_id, islast, inforamses, infocone)

    Implicit none

    Integer(kind=hid_t), intent(in) :: file_id
    Type(Type_inforamses), intent(out), optional :: inforamses
    Class(Type_infocone), intent(out), optional :: infocone
    Logical(kind=4), intent(in) :: islast

    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: aname
    Character(len=H5STRLEN) :: dsetname

    Integer(kind=hid_t) :: gr_id
    Integer(kind=hid_t) :: gr_ramses_id
    Integer(kind=hid_t) :: gr_cone_id
    Integer(kind=8) :: npart8
    Integer(kind=4) :: tmpint4

    groupname = 'metadata'
    Call hdf5_open_group(file_id,groupname, gr_id)

    aname = 'constant_mass'
    tmpint4 = 1
    Call hdf5_read_attr(gr_id, aname, tmpint4)

    aname = 'pfof_cells_order'
    Call hdf5_read_attr(gr_id, aname, H5STRLEN, pfof_cell_order)
    
    ! Read the number of particles in this simulation as an integer(kind=8) dataset
    dsetname = 'npart_simulation'
    Call hdf5_read_data(gr_id, dsetname, npart8)


    If(present(inforamses)) Then
       ! Ramses Info Metadata
#ifdef DEBUG
       Print *,'Read inforamses'
#endif       
       groupname = 'ramses_info'
       If(islast) Then
          groupname = 'ramses_info_last'
       End If
       Call hdf5_open_group(gr_id,groupname,gr_ramses_id)

       aname = 'ncpu'
       Call hdf5_read_attr(gr_ramses_id,aname,inforamses%ncpu)

       aname = 'ndim'
       Call hdf5_read_attr(gr_ramses_id,aname,inforamses%ndim)

       aname = 'levelmin'
       Call hdf5_read_attr(gr_ramses_id,aname,inforamses%levelmin)
#ifdef DEBUG
       Print *,'levelmin=',inforamses%levelmin
#endif       

       aname = 'levelmax'
       Call hdf5_read_attr(gr_ramses_id,aname,inforamses%levelmax)

       aname = 'ngridmax'
       Call hdf5_read_attr(gr_ramses_id,aname,inforamses%ngridmax)

       aname = 'nstep_coarse'
       Call hdf5_read_attr(gr_ramses_id,aname,inforamses%nstep_coarse)

       aname = 'boxlen'
       Call hdf5_read_attr(gr_ramses_id,aname,inforamses%boxlen)

       aname = 'time'
       Call hdf5_read_attr(gr_ramses_id,aname,inforamses%time)

       aname = 'aexp'
       Call hdf5_read_attr(gr_ramses_id,aname,inforamses%aexp)

       aname = 'h0'
       Call hdf5_read_attr(gr_ramses_id,aname,inforamses%h0)

       aname = 'omega_m'
       Call hdf5_read_attr(gr_ramses_id,aname,inforamses%omega_m)

       aname = 'omega_l'
       Call hdf5_read_attr(gr_ramses_id,aname,inforamses%omega_l)

       aname = 'omega_k'
       Call hdf5_read_attr(gr_ramses_id,aname,inforamses%omega_k)

       aname = 'omega_b'
       Call hdf5_read_attr(gr_ramses_id,aname,inforamses%omega_b)

       aname = 'unit_l'
       Call hdf5_read_attr(gr_ramses_id,aname,inforamses%unit_l)

       aname = 'unit_d'
       Call hdf5_read_attr(gr_ramses_id,aname,inforamses%unit_d)

       aname = 'unit_t'
       Call hdf5_read_attr(gr_ramses_id,aname,inforamses%unit_t)

       
       Call hdf5_close_group(gr_ramses_id)
    End If

    If(present(infocone)) Then
       groupname = 'cone_info'
       If(islast) Then
          groupname = 'cone_info_last'
       End If
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
       Type Is(Type_infocone_part)
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
       Type Is(Type_infocone_grav)
          aname = 'nglobalcell'
          Call hdf5_read_data(gr_cone_id, aname, infocone%nglobalcell)
          aname = 'nlevel'
          Call hdf5_read_attr(gr_cone_id, aname, infocone%nlevel)
          aname = 'levelmin'
          Call hdf5_read_attr(gr_cone_id, aname, infocone%levelmin)
          aname = 'levelmax'
          Call hdf5_read_attr(gr_cone_id, aname, infocone%levelmax)
       End Select


       Call hdf5_close_group(gr_cone_id)
    End If

    Call hdf5_close_group(gr_id)

  End Subroutine h5readcommonmetadata


  !=======================================================================
  Subroutine h5read_meta_info_ramses(file_id, inforamses, islast)

    Implicit none

    Integer(kind=hid_t), intent(in) :: file_id
    Type(Type_inforamses), intent(out) :: inforamses
    Logical(kind=4), intent(in) :: islast

    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: aname

    Integer(kind=hid_t) :: gr_id
    Integer(kind=hid_t) :: gr_ramses_id

    groupname = 'metadata'
    Call hdf5_open_group(file_id, groupname, gr_id)

!!$    ! Ramses Info Metadata
!!$    groupname = 'ramses_info'
!!$    Call hdf5_open_group(gr_id,groupname,gr_ramses_id)

    ! Ramses Info Metadata
    groupname = 'ramses_info'
    If(islast) Then
       groupname = 'ramses_info_last'
    End If
    Call hdf5_open_group(gr_id,groupname,gr_ramses_id)

    aname = 'ncpu'
    Call hdf5_read_attr(gr_ramses_id,aname,inforamses%ncpu)

    aname = 'ndim'
    Call hdf5_read_attr(gr_ramses_id,aname,inforamses%ndim)

    aname = 'levelmin'
    Call hdf5_read_attr(gr_ramses_id,aname,inforamses%levelmin)

    aname = 'levelmax'
    Call hdf5_read_attr(gr_ramses_id,aname,inforamses%levelmax)

    aname = 'ngridmax'
    Call hdf5_read_attr(gr_ramses_id,aname,inforamses%ngridmax)

    aname = 'nstep_coarse'
    Call hdf5_read_attr(gr_ramses_id,aname,inforamses%nstep_coarse)

    aname = 'boxlen'
    Call hdf5_read_attr(gr_ramses_id,aname,inforamses%boxlen)

    aname = 'time'
    Call hdf5_read_attr(gr_ramses_id,aname,inforamses%time)

    aname = 'aexp'
    Call hdf5_read_attr(gr_ramses_id,aname,inforamses%aexp)

    aname = 'h0'
    Call hdf5_read_attr(gr_ramses_id,aname,inforamses%h0)

    aname = 'omega_m'
    Call hdf5_read_attr(gr_ramses_id,aname,inforamses%omega_m)

    aname = 'omega_l'
    Call hdf5_read_attr(gr_ramses_id,aname,inforamses%omega_l)

    aname = 'omega_k'
    Call hdf5_read_attr(gr_ramses_id,aname,inforamses%omega_k)

    aname = 'omega_b'
    Call hdf5_read_attr(gr_ramses_id,aname,inforamses%omega_b)

    aname = 'unit_l'
    Call hdf5_read_attr(gr_ramses_id,aname,inforamses%unit_l)

    aname = 'unit_d'
    Call hdf5_read_attr(gr_ramses_id,aname,inforamses%unit_d)

    aname = 'unit_t'
    Call hdf5_read_attr(gr_ramses_id,aname,inforamses%unit_t)

    Call hdf5_close_group(gr_ramses_id)

    Call hdf5_close_group(gr_id)

  End Subroutine h5read_meta_info_ramses

  !=======================================================================
  Subroutine h5read_meta_info_cone(file_id, infocone, islast)

    Implicit none

    Integer(kind=hid_t), intent(in) :: file_id
    Class(Type_infocone), intent(out) :: infocone
    Logical(kind=4), intent(in) :: islast

    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: aname

    Integer(kind=hid_t) :: gr_id
    Integer(kind=hid_t) :: gr_cone_id

    groupname = 'metadata'
    Call hdf5_open_group(file_id, groupname, gr_id)

    groupname = 'cone_info'
    If(islast) Then
       groupname = 'cone_info_last'
    End If
    Call hdf5_open_group(gr_id, groupname, gr_cone_id)

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
    Type Is(Type_infocone_part)
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
    Type Is(Type_infocone_grav)
       aname = 'nglobalcell'
       Call hdf5_read_data(gr_cone_id, aname, infocone%nglobalcell)
       aname = 'nlevel'
       Call hdf5_read_attr(gr_cone_id, aname, infocone%nlevel)
       aname = 'levelmin'
       Call hdf5_read_attr(gr_cone_id, aname, infocone%levelmin)
       aname = 'levelmax'
       Call hdf5_read_attr(gr_cone_id, aname, infocone%levelmax)
    End Select

    Call hdf5_close_group(gr_cone_id)

    Call hdf5_close_group(gr_id)

  End Subroutine h5read_meta_info_cone

  !=======================================================================
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

    groupname = 'fof_parameters'
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

       aname = 'do_fof'
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
    aname = 'do_timings'
    Call hdf5_read_attr(gr_fof_id, aname, tmpint4)
    param_pfof%do_timings = .false.
    If(tmpint4==1) param_pfof%do_timings = .true.

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
       Call hdf5_read_attr(gr_input_id, aname, param_cone%nlevel)
       aname = 'levelmin'
       Call hdf5_read_attr(gr_input_id, aname, param_cone%levelmin)
    End Select

    Call hdf5_close_group(gr_input_id)

    groupname = 'output_parameters'
    Call hdf5_open_group(gr_param_id, groupname, gr_output_id)

    aname = 'simulation_name'
    Call hdf5_read_attr(gr_output_id, aname, len(param_cone%simulation_name), &
         param_cone%simulation_name)

    aname = 'cube_size'
    Call hdf5_read_attr(gr_output_id, aname, param_cone%cube_size)

    Call hdf5_close_group(gr_output_id)
    Call hdf5_close_group(gr_param_id)
    Call hdf5_close_group(gr_id)


  End Subroutine h5read_meta_conecreator_parameter

End Module modiocommons
