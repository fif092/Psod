!==============================================================================
! Project: pFoF
! File: tools/conemapper/src/modio.f90
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
!! This file contains some input/output subroutines for the conemapper tool.

!> This module contains some input/output subroutines for the conemapper tool.
Module modio
  
  Use modhdf5
  Use modvariables
  Use modreadparameters
  
  Integer :: ioerr

Contains

  !=======================================================================

  Subroutine theend()

    Implicit none
    
    Print *,' '
    Print *,'Run Completed!'
    Print *,' '

  End Subroutine theend

  !=======================================================================
  !> This subroutine reads the parameters from the file pfof_cone.nml
  Subroutine readparameters()

    Implicit None

    Character(len=400) :: parameters_file
    Integer(kind=4) :: errcode
    Character(len=500) :: errmessage

    parameters_file = 'pfof_cone.nml'

    Call read_pfof_cone_parameters(parameters_file, parameters, errcode, errmessage)

    If(errcode > 0) Then
       Print *,trim(errmessage)
       Stop 
       !Call exit(errcode)
    End If
 
    ! nb of shell files 
    shell_nb = parameters%shell_last_id - parameters%shell_first_id + 1
    
  End Subroutine readparameters


  !=======================================================================
  !> This subroutine reads the cone size and the number of particles in each
  !> cubic group from each shell file, then perform a reduction to get the 
  !> total number of particles for each cubic group.
  Subroutine h5readshellinfo()
    
    Implicit None

    Character(len=H5STRLEN) :: dataname
    Character(len=H5STRLEN) :: groupname
    Character(len=5) :: charish
    Character(len=400) :: shellname
    Integer :: ish, ic
    Integer(kind=4) :: isfullsky
    Integer, dimension(:), allocatable :: nptemp
    
    Integer(kind=hid_t) :: file_id
    Integer(kind=hid_t) :: gr_m_id
    Integer(kind=hid_t) :: gr_ci_id
    Integer(kind=hid_t) :: gr_pp_id
    Integer(kind=hid_t) :: gr_op_id
    
    Do ish = parameters%shell_first_id, parameters%shell_last_id
       
       Write(charish(1:5),'(I5.5)') ish
       shellname = trim(parameters%input_path)//trim(parameters%part_input_file)//'_'//charish//'.h5'
       Call hdf5_open_file(shellname, file_id)
       
       groupname='metadata'
       Call hdf5_open_group(file_id, groupname, gr_m_id)

       groupname='cone_info'
       Call hdf5_open_group(gr_m_id, groupname, gr_ci_id)

       fullsky=.false.
       dataname='isfullsky'
       Call hdf5_read_attr(gr_ci_id, dataname, isfullsky)
       If(isfullsky==1) fullsky=.true.

       Call hdf5_close_group(gr_ci_id)
    
       groupname='conecreator_part_parameters'
       Call hdf5_open_group(gr_m_id, groupname, gr_pp_id)
       groupname='output_parameters'
       Call hdf5_open_group(gr_pp_id, groupname, gr_op_id)
       
       dataname = 'cube_size'
       Call hdf5_read_attr(gr_op_id, dataname, cubesize)

       Call hdf5_close_group(gr_op_id)
       Call hdf5_close_group(gr_pp_id)

       dataname = 'ncube_array'
       Call hdf5_read_attr(gr_m_id, dataname, 3, nctab)
       
       dataname = 'npart_cube_array'
       ncube = nctab(1)*nctab(2)*nctab(3)
       
       If(.not.Allocated(npartcube)) Then
          Allocate(npartcube(ncube))
          npartcube = 0
       End If
       
       If(.not.Allocated(nptemp)) Allocate(nptemp(ncube))
       Call hdf5_read_data(gr_m_id, dataname, ncube, nptemp)

       npartcube = npartcube + nptemp

       Call hdf5_close_group(gr_m_id)
       Call hdf5_close_file(file_id)
       
    End Do
    
    npart = sum(npartcube)
    
    npart = 0
    Do ic = 1, ncube
       npart = npart + npartcube(ic)
    End Do


    Deallocate(nptemp)

  End Subroutine h5readshellinfo


  !=======================================================================
  !> This subroutine writes a map for pfof_cone, this map contains the 
  !> id of the shell cubes and the procID of the neighbours for each 
  !> pfof_cone process.
  Subroutine h5writemap(factor, procNB, dims, coords)

    Implicit None

    Integer(kind=4), intent(in) :: factor
    Integer(kind=4), intent(in) :: procNB
    Integer(kind=4), dimension(3), intent(in) :: dims
    Integer(kind=4), dimension(3,procNB), intent(in) :: coords
    Character(len=400) :: filename
    Character(len=H5STRLEN) :: groupname
    Character(len=H5STRLEN) :: dataname
    Character(len=6) :: charpnb
    Character(len=4) :: charfac
    Character(len=8) :: charpid
    Integer(kind=hid_t) :: file_id, gr_id
    Integer(kind=4) :: factorcube
    Integer(kind=4) :: ip

    factorcube = size(pictable,1)

    Write(charfac(1:4),'(I4.4)') factor
    Write(charpnb(1:6),'(I6.6)') procNB
    filename = 'procmap_'//charfac//'_'//charpnb//'.h5'
    Print '(A)','Name of the map file: ', trim(filename)
    Print *,' '

    Call hdf5_create_file(filename, file_id)

    dataname = 'process_number'
    Call hdf5_write_attr(file_id, dataname, procNB)

    dataname='ncube_shell'
    Call hdf5_write_attr(file_id, dataname, factorcube)

    dataname='commdims'
    Call hdf5_write_attr(file_id, dataname, 3, dims)
    

    Do ip = 1, procNB
       Write(charpid(1:8),'(I8.8)') ip-1
       groupname = 'process_'//charpid
       Call hdf5_create_group(file_id, groupname, gr_id)

       dataname = 'boundaries'
       Call hdf5_write_attr(gr_id, dataname, 6, boundaries(:,ip))

       dataname='npart_process'
       Call hdf5_write_attr(gr_id, dataname, int(my_npart_tab(ip),kind=4))
       
       dataname = 'idc_array'
       Call hdf5_write_data(gr_id, dataname, factorcube, pictable(:,ip))
       dataname = 'neighbours'
       Call hdf5_write_data(gr_id, dataname, 6, neigh(:,ip))

       dataname = 'coords'
       Call hdf5_write_data(gr_id, dataname, 3, coords(:,ip))

       Call hdf5_close_group(gr_id)
    End Do

    Call hdf5_close_file(file_id)

  End Subroutine h5writemap

End Module modio