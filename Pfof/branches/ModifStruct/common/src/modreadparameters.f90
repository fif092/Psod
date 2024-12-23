!==============================================================================
! Project: pFoF
! File: common/src/modreadparameters.f90
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
!! This file contains subroutines to read input parameters for snapshot and light cone versions of pfof.

!> This module contains subroutines to read input parameters for snapshot and light cone versions of pfof.
!>
!> Authors: F. Roy

Module modreadparameters

  Private
  Public :: read_pfof_cone_parameters, &
       read_pfof_snap_parameters
Contains
  
  !=======================================================================
  !> Read pfof cone input parameters from pfof_cone.nml and save them if requested.
  Subroutine read_pfof_cone_parameters(filename, param, ioerr, errormessage, do_save)
    
    Use modconstant
    Implicit none

    ! Input/output variables
    Character(len=400), intent(in) :: filename
    Type(Type_parameter_pfof_cone), intent(out) :: param
    Integer(kind=4), intent(out)    :: ioerr
    Character(len=500), intent(out) :: errormessage
    Logical(kind=4), intent(in), optional :: do_save

    ! Local variable 
    Character(len=230) :: filesave

    !! pfof parameters
    ! Input parameters
    Character(len=200) :: input_path      !< path to the directory containing the RAMSES files
    Character(len=200) :: part_input_file  !< name of the files containing the particles data
    Logical(kind=4)    :: do_read_ramses_part_id
    Logical(kind=4)    :: do_read_potential    !< do the RAMSES files contain potential?
    Logical(kind=4)    :: do_read_gravitational_field  !< do the RAMSES files contain force?
    Integer(kind=4)    :: shell_first_id
    Integer(kind=4)    :: shell_last_id

    ! Friend of friend parameters
    Real(kind=4)       :: percolation_length  !< FOF percolation length
    Integer(kind=4)    :: mmin            !< minimum mass of the haloes
    Integer(kind=4)    :: mmax            !< maximum mass of the haloes
    Logical(kind=4)    :: do_unbinding    !< Unbinding process? (not implemented yet)
    Logical(kind=4)    :: do_subHalo      !< SubHalo detection? (not implemented yet)

    ! Output parameters
    Character(len=200) :: simulation_name  !< Simulation name to be written in output files
    Logical(kind=4)    :: do_gather_halo   !< should pFOF gather halos in output files
    Logical(kind=4)    :: do_timings  !< should pFOF perform timings (imply extra synchronizations)

    Namelist / input_parameters / input_path, part_input_file, do_read_ramses_part_id, &
         do_read_potential, do_read_gravitational_field, shell_first_id, shell_last_id
    Namelist / fof_parameters / percolation_length, mmin, mmax, do_unbinding, do_subHalo
    Namelist / output_parameters / simulation_name, do_gather_halo, do_timings
    

    Open(Unit=10, file=filename, iostat=ioerr) 

    If(ioerr>0) Then
       errormessage='** Error opening input file '//trim(filename)//'. Please check this file. **'
       Return
    End If
    
    Read(10, nml=input_parameters, iostat=ioerr)
    If(ioerr>0) Then
       errormessage='** Error reading input parameters in '//trim(filename)//'. Please check this file. **'
       Return
    End If

    Read(10, nml=fof_parameters, iostat=ioerr)
    If(ioerr>0) Then
       errormessage='** Error reading fof parameters in '//trim(filename)//'. Please check this file. **'
       Return
    End If

    Read(10, nml=output_parameters, iostat=ioerr)
    If(ioerr>0) Then
       errormessage='** Error reading output parameters in '//trim(filename)//'. Please check this file. **'
       Return
    End If

    Close(10)

    param%input_path = input_path
    param%part_input_file = part_input_file
    param%simulation_name = simulation_name
    param%mmin = mmin
    param%mmax = mmax
    param%percolation_length = percolation_length 
    param%do_read_ramses_part_id = do_read_ramses_part_id
    param%do_read_potential = do_read_potential
    param%do_read_gravitational_field = do_read_gravitational_field
    param%do_unbinding = do_unbinding
    param%do_subHalo = do_subHalo
    param%do_gather_halo = do_gather_halo
    param%do_timings = do_timings
    param%shell_first_id = shell_first_id
    param%shell_last_id = shell_last_id

    If(present(do_save)) Then
       If(do_save) Then
          filesave = 'pfof_cone_parameters_'//trim(simulation_name)//'.nml'
          Open(Unit=10,file=filesave)
          Write(10, nml=input_parameters)
          Write(10, nml=fof_parameters)
          Write(10, nml=output_parameters)
          Close(10)
       End If
    End If

  End Subroutine read_pfof_cone_parameters


  !=======================================================================
  !> Read pfof snap input parameters from pfof_snap.nml and save them if requested.
  Subroutine read_pfof_snap_parameters(filename, param, ioerr, errormessage, do_save)
  
    Use modconstant

    Implicit none

    ! I/O variables
    Character(len=400), intent(in) :: filename
    Type(Type_parameter_pfof_snap), intent(out) :: param
    Integer(kind=4), intent(out)    :: ioerr
    Character(len=500), intent(out) :: errormessage
    Logical(kind=4), intent(in), optional :: do_save

    ! Local variable 
    Character(len=230) :: filesave

    ! Input parameters
    Character(len=3)   :: code_index
    Character(len=200) :: input_path
    Character(len=200) :: part_input_file
    Character(len=200) :: info_input_file
    Integer(kind=4)    :: grpsize
    Logical(kind=4)    :: do_skip_star
    Logical(kind=4)    :: do_skip_metal
    Logical(kind=4)    :: do_read_potential
    Logical(kind=4)    :: do_read_gravitational_field
    Logical(kind=4)    :: do_read_from_cube
    Integer(kind=4)    :: gatherread_factor

    ! Friend of friend parameters
    Real(kind=4)       :: percolation_length
    Integer(kind=4)    :: mmin
    Integer(kind=4)    :: mmax
    Logical(kind=4)    :: do_fof
    Logical(kind=4)    :: do_unbinding
    Logical(kind=4)    :: do_subHalo
    
    ! Output parameters
    Character(len=200) :: simulation_name
    Integer(kind=4)    :: snapshot
    Logical(kind=4)    :: do_write_cube
    Integer(kind=4)    :: gatherwrite_factor
    Logical(kind=4)    :: do_sort_cube
    Logical(kind=4)    :: do_timings
    
    
    ! Namelist for input file
    Namelist / input_parameters  / code_index, input_path, part_input_file, info_input_file, &
         grpsize, do_skip_star, do_skip_metal, do_read_potential, do_read_gravitational_field, &
         do_read_from_cube, gatherread_factor
    Namelist / fof_parameters    / percolation_length, mmin, mmax, do_fof, do_unbinding, do_subhalo
    Namelist / output_parameters / simulation_name, snapshot, do_write_cube,&
         gatherwrite_factor,  do_sort_cube, do_timings
    

    ! Read input parameters'
    Open(10, file=filename, iostat=ioerr, status='old') 
    
    If(ioerr /= 0) Then
       errormessage='** Error opening input file '//trim(filename)//'. Please check this file. **'
       Return
    End If
    
    Read(10, nml=input_parameters, iostat=ioerr)
    If(ioerr /= 0) Then
       errormessage='** Error reading input parameters in '//trim(filename)//'. Please check this file. **'
       Return
    End If

    Read(10, nml=fof_parameters, iostat=ioerr)
    If(ioerr /= 0) Then
       errormessage='** Error reading fof parameters in '//trim(filename)//'. Please check this file. **'
       Return
    End If

    Read(10, nml=output_parameters, iostat=ioerr)
    If(ioerr /= 0) Then
       errormessage='** Error reading output parameters in '//trim(filename)//'. Please check this file. **'
       Return
    End If

    Close(10)
        
    param%code_index = code_index 
    param%input_path = input_path 
    param%part_input_file = part_input_file 
    param%info_input_file = info_input_file 
    param%grpsize = grpsize 
    param%do_skip_star = do_skip_star 
    param%do_skip_metal = do_skip_metal 
    param%do_read_potential = do_read_potential 
    param%do_read_gravitational_field = do_read_gravitational_field 
    param%do_read_from_cube = do_read_from_cube 
    param%gatherread_factor = gatherread_factor
    param%percolation_length = percolation_length 
    param%mmin = mmin 
    param%mmax = mmax 
    param%do_fof = do_fof 
    param%do_unbinding = do_unbinding 
    param%do_subhalo = do_subhalo
    param%simulation_name = simulation_name 
    param%snapshot = snapshot 
    param%do_write_cube = do_write_cube
    param%gatherwrite_factor = gatherwrite_factor 
    param%do_sort_cube = do_sort_cube 
    param%do_timings = do_timings
    
    If(present(do_save)) Then
       If(do_save) Then
          filesave = 'pfof_snap_parameters_'//trim(simulation_name)//'.nml'
          Print *,'FILE=',filesave
          Open(Unit=10,file=filesave)
          Write(10, nml=input_parameters)
          Write(10, nml=fof_parameters)
          Write(10, nml=output_parameters)
          Close(10)
       End If
    End If


  End Subroutine read_pfof_snap_parameters
  
End Module modreadparameters
