Module modreadparameters
  
Contains
  
  Subroutine read_pfof_cone_parameters(filename, param)
    
    Use modconstant
    Implicit none
    
    Character(len=400), intent(in) :: filename
    Type(Type_parameter_pfof_cone), intent(out) :: param
    Integer(kind=4) :: ioerr

    Character(len=200) :: input_path      !< path to the directory containing the RAMSES files
    Character(len=200) :: part_input_file  !< name of the files containing the particles data
!    Character(len=200) :: info_input_file  !< name of the RAMSES info file
    Character(len=200) :: simulation_name  !< Simulation name to be written in output files
    Real(kind=4)       :: percolation_length  !< FOF percolation length
    Integer(kind=4)    :: mmin            !< minimum mass of the haloes
    Integer(kind=4)    :: mmax            !< maximum mass of the haloes
    Logical(kind=4)    :: do_read_ramses_part_id
    Logical(kind=4)    :: do_read_potential    !< do the RAMSES files contain potential?
    Logical(kind=4)    :: do_read_gravitational_field  !< do the RAMSES files contain force?
    Logical(kind=4)    :: do_unbinding    !< Unbinding process? (not implemented yet)
    Logical(kind=4)    :: do_subHalo      !< SubHalo detection? (not implemented yet)
    Logical(kind=4)    :: do_timings  !< should pFOF perform timings (imply extra synchronizations)
    Logical(kind=4)    :: do_gather_halo !< should pFoF gather halos while writing halopart files
    Integer(kind=4) :: shell_first_id
    Integer(kind=4) :: shell_last_id

    Namelist / input_parameters / input_path, part_input_file, do_read_ramses_part_id, &
         do_read_potential, do_read_gravitational_field, shell_first_id, shell_last_id
    Namelist / fof_parameters / percolation_length, mmin, mmax, do_unbinding, do_subHalo
    Namelist / output_parameters / simulation_name, do_gather_halo, do_timings
    

    Open(Unit=10, file=filename, iostat=ioerr) 
    If(ioerr>0) Then
!       Print *,'** Error opening input file pfof_cone.nml. Please check this file. **'
       Stop '** Error opening input file pfof_cone.nml. Please check this file. **'
    End If
    Read(10, nml=input_parameters)
    Read(10, nml=fof_parameters)
    Read(10, nml=output_parameters)
    Close(10)

    param%input_path = input_path
    param%part_input_file = part_input_file
    param%simulation_name = simulation_name
    param%mmin = mmin
    param%mmax = mmax
    param%do_read_ramses_part_id = do_read_ramses_part_id
    param%do_read_potential = do_read_potential
    param%do_read_gravitational_field = do_read_gravitational_field
    param%do_unbinding = do_unbinding
    param%do_subHalo = do_subHalo
    param%do_gather_halo = do_gather_halo
    param%do_timings = do_timings
    param%shell_first_id = shell_first_id
    param%shell_last_id = shell_last_id

  End Subroutine read_pfof_cone_parameters

End Module modreadparameters
